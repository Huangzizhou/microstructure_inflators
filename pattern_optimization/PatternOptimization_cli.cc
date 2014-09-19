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
        ("pattern,p",       po::value<string>(), "Pattern wire mesh (.obj)")
        ("material,m",      po::value<string>(), "base material")
        ("output,o",        po::value<string>(), "output .js mesh + fields")
        ("max_area,a",      po::value<double>()->default_value(0.0001), "max_area parameter for wire inflator")
        ("step,s", po::value<double>(), "gradient step size")
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

// Writes steepest descent direction for
//      1/2 sum_ijkl (target_ijkl - Eh_ijlk|)^2
// That is, -grad(1/2 sum_ijkl (target_ijkl - Eh_ijlk|)^2) =
//  -(target_ijkl - Eh_ijlk) * -grad(Eh_ikjl)) =
//   (target_ijkl - Eh_ijlk) *  grad(Eh_ikjl))
template<size_t _N, class Simulator>
void writeTargetTensorShapeDerivative(JSFieldWriter<_N> &writer,
                const string &name, const ETensor<_N> &target,
                const ETensor<_N> &current, const vector<ETensor<_N>> &gradEh,
                const std::vector<VField<_N>> &w_ij, const Simulator &sim) {
    size_t numBE = sim.mesh().numBoundaryElements();
    assert(gradEh.size() == numBE);
    ETensor<_N> diff = target - current;
    SField v_n(numBE);

    for (size_t be = 0; be < numBE; ++be)
        v_n[be] = diff.quadrupleContract(gradEh[be]);

    writer.addField(name + " v_n", v_n, JSFieldWriter<_N>::PER_BDRY_ELEM);
    writeDescentVectors(writer, name, v_n, sim);
}

// Writes steepest descent direction for
//      1/2 sum_ijkl (target_ijkl - Einv_ijlk|)^2
// That is, -grad(1/2 sum_ijkl (target_ijkl - Einv_ijlk|)^2) =
//  -(target_ijkl - Einv_ijlk) * -grad(Einv_ikjl)) =
//   (target_ijkl - Einv_ijlk) *  grad(Einv_ikjl))
template<size_t _N, class Simulator>
void writeTargetTensorShapeDerivative(JSFieldWriter<_N> &writer,
                const string &name, const ETensor<_N> &target,
                const ETensor<_N> &current, const vector<ETensor<_N>> &gradEh,
                const std::vector<VField<_N>> &w_ij, const Simulator &sim) {
    size_t numBE = sim.mesh().numBoundaryElements();
    assert(gradEh.size() == numBE);
    ETensor<_N> diff = target - current;
    SField v_n(numBE);

    for (size_t be = 0; be < numBE; ++be)
        v_n[be] = diff.quadrupleContract(gradEh[be]);

    writer.addField(name + " v_n", v_n, JSFieldWriter<_N>::PER_BDRY_ELEM);
    writeDescentVectors(writer, name, v_n, sim);
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

    // Focus on the best fitted tensor in the object
    size_t optElement = inElements.size();
    Real bestDist = 1e10;
    ETensor<2> bestLUTTensor, bestTargetTensor;
    ETensor<2> lutTensor, targetTensor;
    for (size_t i = 0; i < inElements.size(); ++i) {
        lutTensor.setOrthotropic2D(lut_Ex[i], lut_Ey[i], lut_nu_yx[i], lut_mu[i]);
        targetTensor.setIsotropic(tgt_E[i], tgt_nu[i]);
        auto diff = lutTensor.inverse() - targetTensor.inverse();
        Real dist = diff.quadrupleContract(diff);
        if (dist < bestDist) {
            bestDist = dist;
            bestLUTTensor = lutTensor;
            bestTargetTensor = targetTensor;
            optElement = i;
        }
    }

    cout << "Optimizing best-fit on element " << optElement << endl
         << "initial distance:\t" << bestDist << endl;

    cout << "LUT tensor:" << bestLUTTensor << endl;
    cout << "LUT tensor:" << endl << bestLUTTensor << endl <<endl;
    cout << "Target tensor:" << endl << bestTargetTensor << endl <<endl;
    cout << "LUT moduli: "
         << lut_Ex[optElement] << ", " << lut_Ey[optElement] << ", "
         << lut_nu_yx[optElement] << ", " << lut_mu[optElement] << endl;
    cout << "Target moduli: " << tgt_E[optElement] << ", "
         << tgt_nu[optElement] << endl;

    // Inflate pattern.
    // WireInflator2D inflator(args["pattern"].as<string>());
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

    if (!inflator.patternGenerator().parametersValid(p_params))
        throw runtime_error("Invalid parameters specified.");

    WireInflator2D::OutMeshType mesh;
    inflator.generatePattern(p_params, t_params, mesh);

    auto &mat = LinearElasticityND<_N>::
        template homogenousMaterial<Materials::Constant>();
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());
    typename LinearElasticityND<_N>:: template
        HomogenousSimulator<Materials::Constant> sim(mesh.elements, mesh.nodes);

    JSFieldWriter<_N> *writer = NULL;
    if (args.count("output"))
        writer = new JSFieldWriter<_N>(args["output"].as<string>(), sim.mesh());

    std::vector<VField<_N>> w_ij;
    // solveCellProblems(w_ij, sim, writer);
    solveCellProblems(w_ij, sim, (JSFieldWriter<_N> *) NULL);
    ETensor<_N> Eh = homogenizedElasticityTensor(w_ij, sim);

    cout << setprecision(16) << endl;
    cout << "Homogenized elasticity tensor:" << endl;
    cout << Eh << endl << endl;

    ETensor<_N> Einv = Eh.inverse();
    auto moduli((1.0 / Einv.diag().array()).eval());

    vector<Real> poisson;
    if (_N == 2) poisson = { -Einv.D(0, 1) / Einv.D(1, 1),   // v_yx
                             -Einv.D(1, 0) / Einv.D(0, 0) }; // v_xy
    else         poisson = { -Einv.D(0, 1) / Einv.D(1, 1),   // v_yx
                             -Einv.D(0, 2) / Einv.D(2, 2),   // v_zx
                             -Einv.D(1, 2) / Einv.D(2, 2),   // v_zy
                             -Einv.D(1, 0) / Einv.D(0, 0),   // v_xy
                             -Einv.D(2, 0) / Einv.D(0, 0),   // v_xz
                             -Einv.D(2, 1) / Einv.D(1, 1) }; // v_zy

    if (_N == 2)  {
        cout << "Approximate Young moduli:\t"  << moduli[0] << "\t" << moduli[1] << endl;
        cout << "Approximate shear modulus:\t" << moduli[2] << endl;

        cout << "v_yx, v_xy:\t" << poisson[0] << "\t" << poisson[1] << endl;
    }
    else {
        cout << "Approximate Young moduli:\t" << moduli[0] << "\t" << moduli[1] << "\t"
             << moduli[2] << endl;
        cout << "Approximate shear moduli:\t" << moduli[3] << "\t" << moduli[4] << "\t"
             << moduli[5] << endl;

        cout << "v_yx, v_zx, v_zy:\t" << poisson[0] << "\t" << poisson[1] << "\t" << poisson[2] << endl;
        cout << "v_xy, v_xz, v_yz:\t" << poisson[3] << "\t" << poisson[4] << "\t" << poisson[5] << endl;
    }

    if (writer)  {
        std::vector<ETensor<_N>> gradEh = homogenizedTensorGradient(w_ij, sim);
        ////////////////////////////////////////////////////////////////////////
        // Moduli Derivatives - for orthotropic only!
        ////////////////////////////////////////////////////////////////////////
        // Gradient of the compliance tensor--useful for moduli derivatives.
        std::vector<ETensor<_N>> gradEhinv;
        for (const auto &G : gradEh)
            gradEhinv.push_back(-Einv.doubleDoubleContract(G));

        size_t numBE = sim.mesh().numBoundaryElements();
        SField gradEx(numBE), gradEy(numBE), gradVyx(numBE), gradVxy(numBE),
               gradmu(numBE);
        for (size_t i = 0; i < numBE; ++i) {
             gradEx[i] = -moduli[0] * moduli[0] * gradEhinv[i].D(0, 0);
             gradEy[i] = -moduli[1] * moduli[1] * gradEhinv[i].D(1, 1);
            gradVyx[i] = -gradEy[i] * Einv.D(0, 1) - moduli[1] * gradEhinv[i].D(0, 1);
            gradVxy[i] = -gradEx[i] * Einv.D(1, 0) - moduli[0] * gradEhinv[i].D(1, 0);
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
    // if (writer && args.count("parameterStep")) {
    //     Real parameterStep = args["parameterStep"].as<double>();
    //     cout << "Parameter step: " << parameterStep << endl;

    //     ETensor<_N> ETargetinv(Einv);
    //     Real currentPoisson = -ETargetinv.D(0, 1) / ETargetinv.D(1, 1);
    //     Real targetPoisson = currentPoisson + parameterStep * (-0.5 - currentPoisson);
    //     ETargetinv.D(0, 1) = -targetPoisson * ETargetinv.D(1, 1);
    //     writeTargetTensorShapeDerivative(*writer, "v_yx decreasing",
    //             ETargetinv.inverse(), Eh, gradEh, w_ij, sim);

    //     ETargetinv = Einv;
    //     currentPoisson = -ETargetinv.D(0, 1) / ETargetinv.D(0, 0);
    //     targetPoisson = currentPoisson + parameterStep * (-0.5 - currentPoisson);
    //     ETargetinv.D(0, 1) = -targetPoisson * ETargetinv.D(1, 1);
    //     writeTargetTensorShapeDerivative(*writer, "v_xy decreasing",
    //             ETargetinv.inverse(), Eh, gradEh, w_ij, sim);

    //     // Step a bit toward the tensor with Ex halved.
    //     ETargetinv = Einv;
    //     ETargetinv.D(0, 0) *= 2;
    //     ETargetinv.D(1, 0) *= ETargetinv.D(0, 0) / Einv.D(0, 0);
    //     // ETargetinv.D(0, 0) += parameterStep * ETargetinv.D(0, 0);
    //     // ETargetinv.D(1, 0) *= ETargetinv.D(0, 0) / Einv.D(0, 0);
    //     writeTargetTensorShapeDerivative(*writer, "Ex decreasing, v_xy fixed",
    //             ETargetinv.inverse(), Eh, gradEh, w_ij, sim);

    //     ETargetinv = Einv;
    //     ETargetinv.D(1, 1) *= 2;
    //     ETargetinv.D(0, 1) *= ETargetinv.D(1, 1) / Einv.D(1, 1);
    //     // ETargetinv.D(1, 1) += parameterStep * ETargetinv.D(1, 1);
    //     // ETargetinv.D(0, 1) *= ETargetinv.D(1, 1) / Einv.D(1, 1);
    //     writeTargetTensorShapeDerivative(*writer, "Ey decreasing, v_yx fixed",
    //             ETargetinv.inverse(), Eh, gradEh, w_ij, sim);

    //     ETargetinv = Einv;
    //     ETargetinv.D(2, 2) += parameterStep * ETargetinv.D(2, 2);
    //     writeTargetTensorShapeDerivative(*writer, "mu decreasing",
    //             ETargetinv.inverse(), Eh, gradEh, w_ij, sim);
    // 
    // }

    if (writer) delete writer;
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
