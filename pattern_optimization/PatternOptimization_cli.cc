////////////////////////////////////////////////////////////////////////////////
// PatternOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure parameters to bring the structure's
//      homogenized elasticity tensor closer to a target tensor.
//
//      Expects the following per-element scalar fields:
//      lookup table fields:
//          fitted_poisson_yx
//          fitted_young_x
//          fitted_young_y
//          shear_xy
//      target fields:
//          young
//          poisson
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/12/2014 01:15:28
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <MSHFieldWriter.hh>
#include <MSHFieldParser.hh>
#include <JSFieldWriter.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>

#include <WireInflator2D.h>

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
using namespace std;
using namespace PeriodicHomogenization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: PatternOptimization_cli [options] input.msh" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("msh",       po::value<string>(),                     "input msh")
        ;
    po::positional_options_description p;
    p.add("msh",                1);

    po::options_description visible_opts;
    visible_opts.add_options()("help", "Produce this help message")
        ("pattern,p",   po::value<string>(), "Pattern wire mesh (.obj)")
        ("material,m",  po::value<string>(), "base material")
        ("output,o",    po::value<string>(), "output .js mesh + fields")
        ("max_area,a",  po::value<double>()->default_value(0.0001), "max_area parameter for wire inflator")
        ("step,s",      po::value<double>()->default_value(0.0001), "gradient step size")
        ("nIters,n",    po::value<size_t>()->default_value(20), "number of iterations")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;
    if (vm.count("msh") == 0) {
        cout << "Error: must specify input .msh file" << endl;
        fail = true;
    }

    if (vm.count("pattern") == 0) {
        cout << "Error: must specify pattern mesh" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using ETensor = typename LinearElasticityND<_N>::ETensor;
template<size_t _N>
using VField = typename LinearElasticityND<_N>::VField;
typedef ScalarField<Real> SField;

template<size_t _N, class Simulator>
void writeDescentVectors(JSFieldWriter<_N> &writer,
        const string &name, const SField &v_n, const Simulator &sim) {
    size_t numBE = sim.mesh().numBoundaryElements();
    VField<_N> direction(numBE);
    for (size_t be = 0; be < numBE; ++be)
        direction(be) = v_n[be] * sim.mesh().boundaryElement(be)->normal();
    writer.addField(name + " direction", direction,
            JSFieldWriter<_N>::PER_BDRY_ELEM);
}

////////////////////////////////////////////////////////////////////////////////
/*! Computes grad(1/2 sum_ijkl (curr_ijkl - target_ijlk|)^2) =
//      (curr_ijkl - target_ijlk) * grad(curr_ikjl)) =
//  @param[in]  target, curr target and current rank 4 tensors
//  @param[in]  gradCurr     Shape derivative of the current tensor. This is a
//                           per-boundary-edge tensor field. Each component
//                           gives the average normal velocity a boundary edge
//                           should be advected with to most rapidly increase
//                           the corresponding component of the current tensor.
//  @return     SField       Per-boundary-edge scalar field giving the steepest
//                           ascent normal velocity perturbation for the fitting
//                           objective.
*///////////////////////////////////////////////////////////////////////////////
template<class _ETensor>
SField targetTensorShapeDerivative(const _ETensor &target,
                const _ETensor &curr, const vector<_ETensor> &gradCurr) {
    auto diff = curr - target;
    SField v_n(gradCurr.size());

    for (size_t be = 0; be < gradCurr.size(); ++be)
        v_n[be] = diff.quadrupleContract(gradCurr[be]);

    return v_n;
}

template<size_t _N>
void step(WireInflator2D &inflator, TessellationParameters &t_params,
        CellParameters &p_params, const ETensor<_N> &targetS,
        const string outName, Real alpha) {
    size_t nParams = p_params.numberOfParameters();
    cout << "p = [" << p_params.parameter(0);
    for (size_t p = 1; p < nParams; ++p)
        cout << ", " << p_params.parameter(p);
    cout << "]" << endl;

    if (!inflator.patternGenerator().parametersValid(p_params))
        throw runtime_error("Invalid parameters specified.");

    WireInflator2D::OutMeshType inflatedMesh;
    inflator.generatePattern(p_params, t_params, inflatedMesh);


    typename LinearElasticityND<_N>:: template
        HomogenousSimulator<Materials::Constant> sim(inflatedMesh.elements,
                                                     inflatedMesh.nodes);

    // Output geometry .msh
    MeshIO::save(outName + ".msh", sim.mesh());

    JSFieldWriter<_N> *writer = NULL;
    if (outName != "")
        writer = new JSFieldWriter<_N>(outName, sim.mesh());

    std::vector<VField<_N>> w_ij;
    solveCellProblems(w_ij, sim);
    ETensor<_N> C = homogenizedElasticityTensor(w_ij, sim);

    ETensor<_N> S = C.inverse();
    vector<Real> moduli(flatLen(_N));

    // Shear moduli are multiplied by 4 in flattened compliance tensor...
    for (size_t i = 0; i < flatLen(_N); ++i)
        moduli[i] = ((i < _N) ? 1.0 : 0.25) / S.D(i, i);

    vector<Real> poisson = { -S.D(0, 1) / S.D(1, 1),   // v_yx
                             -S.D(1, 0) / S.D(0, 0) }; // v_xy

    cout << setprecision(16) << endl;
    cout << "Current moduli:\t"  << moduli[0] << "\t" << moduli[1] << "\t"
         << poisson[0] << "\t" << moduli[2] << endl;

    std::vector<ETensor<_N>> gradEh = homogenizedTensorGradient(w_ij, sim);
    std::vector<ETensor<_N>> gradS; gradS.reserve(gradEh.size());
    for (const auto &G : gradEh)
        gradS.push_back(-S.doubleDoubleContract(G));

    SField complianceFitGrad = targetTensorShapeDerivative(targetS,
            S, gradS);

    size_t numBE = sim.mesh().numBoundaryElements();
    std::vector<SField> vn_p(nParams, SField(numBE));
    std::vector<Real> gradp_JS(nParams, 0.0);
    for (size_t bei = 0; bei < numBE; ++bei) {
        auto be = sim.mesh().boundaryElement(bei);
        auto edge = make_pair(be.tip(). volumeVertex().index(),
                              be.tail().volumeVertex().index());
        const auto &field = inflatedMesh.edge_fields.at(edge);
        assert(field.size() == nParams);
        for (size_t p = 0; p < nParams; ++p) {
            vn_p[p][bei] = field[p];
            gradp_JS[p] += vn_p[p][bei] * complianceFitGrad[bei] * be->area();
        }
    }

    auto diffS = S - targetS;
    cout << "J_S = " << diffS.quadrupleContract(diffS) << endl;
    cout << "grad_p(J_S) = [" << gradp_JS[0];
    for (size_t p = 1; p < nParams; ++p)
        cout << ", " << gradp_JS[p];
    cout << "]" << endl;

    Real normSq = 0;
    for (size_t p = 0; p < nParams; ++p) normSq += gradp_JS[p] * gradp_JS[p];
    cout << "||grad_p||:\t" << sqrt(normSq) << endl;

    // Step parameters
    for (size_t p = 0; p < nParams; ++p)
        p_params.parameter(p) -= alpha * gradp_JS[p];

    // Compute effective boundary velocity for visualization
    SField projectedNormalVelocity(numBE);
    for (size_t bei = 0; bei < numBE; ++bei) {
        projectedNormalVelocity[bei] = 0;
        for (size_t p = 0; p < nParams; ++p)
            projectedNormalVelocity[bei] += gradp_JS[p] * vn_p[p][bei];
    }

    if (writer)  {
        writer->addField("gradFit", complianceFitGrad, JSFieldWriter<_N>::PER_BDRY_ELEM);
        writeDescentVectors(*writer, "gradFit", complianceFitGrad, sim);

        writer->addField("projectedVn", projectedNormalVelocity,
                         JSFieldWriter<_N>::PER_BDRY_ELEM);
        writeDescentVectors(*writer, "projectedVn", projectedNormalVelocity, sim);

        for (size_t p = 0; p < nParams; ++p) {
            string name("vn_" + to_string(p));
            writer->addField(name, vn_p[p], JSFieldWriter<_N>::PER_BDRY_ELEM);
            writeDescentVectors(*writer, name, vn_p[p], sim);
        }

        ////////////////////////////////////////////////////////////////////////
        // Moduli Derivatives - for orthotropic only!
        ////////////////////////////////////////////////////////////////////////
        // Gradient of the compliance tensor--useful for moduli derivatives.
        std::vector<ETensor<_N>> gradEhinv;
        for (const auto &G : gradEh)
            gradEhinv.push_back(-S.doubleDoubleContract(G));

        SField gradEx(numBE), gradEy(numBE), gradVyx(numBE), gradVxy(numBE),
               gradmu(numBE);
        for (size_t i = 0; i < numBE; ++i) {
             gradEx[i] = -moduli[0] * moduli[0] * gradEhinv[i].D(0, 0);
             gradEy[i] = -moduli[1] * moduli[1] * gradEhinv[i].D(1, 1);
            gradVyx[i] = -gradEy[i] * S.D(0, 1) - moduli[1] * gradEhinv[i].D(0, 1);
            gradVxy[i] = -gradEx[i] * S.D(1, 0) - moduli[0] * gradEhinv[i].D(1, 0);
             gradmu[i] =  gradEh[i].D(2, 2);
        }
        writer->addField( "gradEx",  gradEx, JSFieldWriter<_N>::PER_BDRY_ELEM);
        writer->addField( "gradEy",  gradEy, JSFieldWriter<_N>::PER_BDRY_ELEM);
        writer->addField("gradVyx", gradVyx, JSFieldWriter<_N>::PER_BDRY_ELEM);
        writer->addField("gradVxy", gradVxy, JSFieldWriter<_N>::PER_BDRY_ELEM);
        writer->addField( "gradmu",  gradmu, JSFieldWriter<_N>::PER_BDRY_ELEM);

        writeDescentVectors(*writer, "gradEx",  gradEx, sim);
        writeDescentVectors(*writer, "gradEy",  gradEy, sim);
        writeDescentVectors(*writer,"gradVyx", gradVyx, sim);
        writeDescentVectors(*writer,"gradVxy", gradVxy, sim);
        writeDescentVectors(*writer, "gradmu",  gradmu, sim);
    }

    if (writer) delete writer;
}

template<size_t _N>
void execute(const po::variables_map &args,
             const vector<MeshIO::IOVertex> &inVertices, 
             const vector<MeshIO::IOElement> &inElements);
template<>
void execute<2>(const po::variables_map &args,
                const vector<MeshIO::IOVertex> &inVertices, 
                const vector<MeshIO::IOElement> &inElements)
{
    constexpr size_t _N = 2;
    MSHFieldParser<2> parser(args["msh"].as<string>());

    SField lut_Ex    = parser.scalarField("fitted_young_x"),
           lut_Ey    = parser.scalarField("fitted_young_y"),
           lut_nu_yx = parser.scalarField("fitted_poisson_yx"),
           lut_mu    = parser.scalarField("shear_xy"),
           tgt_E     = parser.scalarField("young"),
           tgt_nu    = parser.scalarField("poisson");

    // Focus on the worst fitted tensor in the object
    size_t optElement = inElements.size();
    Real worstDist = 0;
    ETensor<2> currentC, targetC;
    for (size_t i = 0; i < inElements.size(); ++i) {
        ETensor<2> lutC, tgtC;
        lutC.setOrthotropic2D(lut_Ex[i], lut_Ey[i], lut_nu_yx[i], lut_mu[i]);
        tgtC.setIsotropic(tgt_E[i], tgt_nu[i]);
        auto diffS = lutC.inverse() - tgtC.inverse();
        Real dist = diffS.quadrupleContract(diffS);
        if (dist > worstDist) {
            worstDist = dist;
            currentC = lutC;
             targetC = tgtC;
            optElement = i;
        }
    }

    cout << "Optimizing worst-fit on element " << optElement << endl
         << "initial distance:\t" << worstDist << endl;

    cout << "LUT tensors:" << endl << currentC << endl << endl << currentC.inverse() << endl << endl;
    cout << "Target tensor:" << endl << targetC << endl << endl << targetC.inverse() << endl << endl;
    cout << "LUT moduli: "
         << lut_Ex[optElement] << ", " << lut_Ey[optElement] << ", "
         << lut_nu_yx[optElement] << ", " << lut_mu[optElement] << endl;
    cout << "Target moduli (table): " << tgt_E[optElement] << ", "
         << tgt_nu[optElement] << endl;

    Real tgt_Ex, tgt_Ey, tgt_nuyx, tgt_mu;
    targetC.getOrthotropic2D(tgt_Ex, tgt_Ey, tgt_nuyx, tgt_mu);
    cout << "Target moduli (tensor):"
         << "\t" << tgt_Ex << "\t" << tgt_Ey << "\t" << tgt_nuyx
         << "\t" << tgt_mu << endl;

	WireInflator2D inflator(args["pattern"].as<string>());
    TessellationParameters t_params;
    t_params.max_area = args["max_area"].as<double>();
    CellParameters         p_params = inflator.createParameters();

    // Current mapping of parameters...
    // vertex_orbit_0_thickness: 1
    // vertex_orbit_1_thickness: 3
    // vertex_orbit_2_thickness: 0
    // vertex_orbit_3_thickness: 2
    // vertex_orbit_4_thickness: 4
    // average vertex_orbit_0_offset_[01]: 6
    // average vertex_orbit_2_offset_[01]: 5
    // average vertex_orbit_4_offset_[01]: 7 and 8
    Real scale = (1 / 5.0) / 2.0;
    p_params.parameter(1) = scale * parser.scalarField("vertex_orbit_0_thickness")[optElement];
    p_params.parameter(3) = scale * parser.scalarField("vertex_orbit_1_thickness")[optElement];
    p_params.parameter(0) = scale * parser.scalarField("vertex_orbit_2_thickness")[optElement];
    p_params.parameter(2) = scale * parser.scalarField("vertex_orbit_3_thickness")[optElement];
    p_params.parameter(4) = scale * parser.scalarField("vertex_orbit_4_thickness")[optElement];
    p_params.parameter(5) = 0.5 * (parser.scalarField("vertex_orbit_2_offset_0")[optElement] + parser.scalarField("vertex_orbit_2_offset_1")[optElement]);
    p_params.parameter(6) = 0.5 * (parser.scalarField("vertex_orbit_0_offset_0")[optElement] + parser.scalarField("vertex_orbit_0_offset_1")[optElement]);
    p_params.parameter(7) = 0.5 * (parser.scalarField("vertex_orbit_4_offset_0")[optElement] + parser.scalarField("vertex_orbit_4_offset_1")[optElement]) / sqrt(2);
    p_params.parameter(8) = p_params.parameter(7);

    // Set up material
    auto &mat = LinearElasticityND<_N>::
        template homogenousMaterial<Materials::Constant>();
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    for (size_t i = 0; i < args["nIters"].as<size_t>(); ++i) {
        string outName = "";
        if (args.count("output"))
            outName = to_string(i) + "_" + args["output"].as<string>();
        step<_N>(inflator, t_params, p_params, targetC.inverse(),
                 outName, args["step"].as<double>());
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    po::variables_map args = parseCmdLine(argc, argv);

    vector<MeshIO::IOVertex>  inVertices;
    vector<MeshIO::IOElement> inElements;
    MeshIO::MeshIO_MSH io;
    string mshPath = args["msh"].as<string>();
    ifstream infile(mshPath);
    if (!infile.is_open()) throw runtime_error("Couldn't open " + mshPath);
    auto type = io.load(infile, inVertices, inElements, MeshIO::MESH_GUESS);

    // Infer dimension from mesh type.
    size_t dim;
    if (type == MeshIO::MESH_QUAD) dim = 2;
    else    throw std::runtime_error("Only 2D quad input is supported.");

    // Look up and run appropriate optimization instantiation.

    execute<2>(args, inVertices, inElements);

    return 0;
}
