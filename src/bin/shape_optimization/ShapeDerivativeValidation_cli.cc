////////////////////////////////////////////////////////////////////////////////
// ShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validate the discrete shape derivative of stress quantities.
//
//      Boundary vertices are offset in the normal direction to create a
//      perturbed mesh, and post- and pre-perturbation quantities are
//      subtracted to compute forward/centered difference (material)
//      derivatives.
//
//      The BoundaryPerturbationInflator is used to ease the creation of a
//      perturbed mesh.
*/
//  Based on DiscreteShapeDerivativeValidation_cli for worst_case_stress (which
//   was implemented by Julian Panetta)
//  Author:  Davi Colli Tozoni (dctozoni), davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  02/05/2018
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/Materials.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <iomanip>
#include <MeshFEM/Laplacian.hh>
#include <MeshFEM/BoundaryConditions.hh>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <inflators/wrappers/ConstrainedInflator.hh>
#include <inflators/wrappers/BoundaryPerturbationInflator.hh>
#include <isosurface_inflator/ShapeVelocityInterpolator.hh>
#include "MicroscopicStress.hh"
#include "ParametersMask.hh"

namespace po = boost::program_options;
using namespace std;

template<size_t _N> using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;
template<size_t _N> using ETensor = ElasticityTensor<Real, _N>;

[[ noreturn ]] void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: ShapeDerivativeValidation_cli [options] mesh.msh" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
            ("mesh", po::value<string>(), "input mesh")
            ;
    po::positional_options_description p;
    p.add("mesh", 1);

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
            ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp global worst case stress measure")
            ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
            ("StressWeight",    po::value<double>()->default_value(1.0),         "Weight for the WCS term of the objective")
            ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
            ("material,m",   po::value<string>(),                    "Base material")
            ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
            ;

    po::options_description simulationOptions;
    simulationOptions.add_options()
            ("boundaryConditions,b", po::value<string>(),                    "boundary conditions")
            ("zeroPerturbationAreas,z", po::value<string>(),                 "areas where there is no perturbation")
            ;

    po::options_description generalOptions;
    generalOptions.add_options()
            ("help,h",                                               "Produce this help message")
            ("output,o",     po::value<string>(),                    "Output the Lagrangian derivatives computed by forward difference and the discrete shape derivative.")
            ("fullDegreeFieldOutput,D",                              "Output full-degree nodal fields (don't do piecewise linear subsample)")
            ("perturbationAmplitude,a", po::value<double>()->default_value(0.001), "Amplitude of boundary perturbation")
            ("perturbationFrequency,f", po::value<double>()->default_value(1.0),  "Frequency of boundary perturbation")
            ;

    po::options_description visibleOptions;
    visibleOptions.add(objectiveOptions).add(elasticityOptions).add(generalOptions).add(simulationOptions);

    po::options_description cli_opts;
    cli_opts.add(visibleOptions).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visibleOptions);
    }

    bool fail = false;

    if (vm.count("mesh") == 0) {
        cout << "Error: must specify input mesh" << endl;
        fail = true;
    }

    if (vm.count("output") == 0) {
        cout << "Error: must specify output mesh" << endl;
        fail = true;
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

template<size_t _N, class Simulator, class VField>
PthRootObjective<IntegratedMicroscopicStressObjective<_N, MicroscopicStressIntegrandLp<Simulator>, Simulator>> buildStressObjective(const po::variables_map &args, const Simulator &sim, VField u, const ETensor<_N> CBase) {
    PthRootObjective<IntegratedMicroscopicStressObjective<_N, MicroscopicStressIntegrandLp<Simulator>, Simulator>> objective;

    objective.integrand.p = args["pnorm"].as<double>();
    objective.p = args.count("usePthRoot") ? 2.0 * args["pnorm"].as<double>() : 1.0;

    const bool major_symmetry = CBase.MajorSymmetry;

    MicroscopicStress<_N, Simulator> microStress = MicroscopicFrobeniusStress<CBase.Dim, major_symmetry, Simulator>(CBase, sim.averageStressField(u));
    objective.setPointwiseStress(sim.mesh(), microStress);
    return objective;
}

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args,
             const std::vector<MeshIO::IOVertex>  &vertices,
             const std::vector<MeshIO::IOElement> &elements)
{
    using Mesh      = typename LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using VField    = typename Simulator::VField;

    string bcondsPath = args["boundaryConditions"].as<string>();
    string zeroPerturbationPath = args["zeroPerturbationAreas"].as<string>();
    bool no_rigid_motion;
    Simulator tmp(elements, vertices);
    vector<CondPtr<_N> > zeroPerturbationConds = readBoundaryConditions<_N>(zeroPerturbationPath, tmp.mesh().boundingBox(), no_rigid_motion);

    // Original mesh and simulator
    std::vector<bool> paramsMask = ParametersMask::generateParametersMask<_N>(vertices, elements, zeroPerturbationPath);
    std::unique_ptr<BoundaryPerturbationInflator<_N>> originalBpi = Future::make_unique<BoundaryPerturbationInflator<_N>>(vertices, elements, false);
    ConstrainedInflator<_N> bpi(std::move(originalBpi), paramsMask);
    vector<Real> perturbParams(bpi.numParameters(), 0.0);
    bpi.inflate(perturbParams);
    Simulator sim(bpi.elements(), bpi.vertices());

    BoundaryPerturbationInflator<_N> * baseInflator = dynamic_cast<BoundaryPerturbationInflator<_N> *>(bpi.m_infl.get());

    vector<CondPtr<_N> > bconds = readBoundaryConditions<_N>(bcondsPath, sim.mesh().boundingBox(), no_rigid_motion);

    // Perturb mesh
    auto &currentMesh = sim.mesh();
    auto normals = baseInflator->boundaryVertexNormals();
    VField perturbation(sim.mesh().numBoundaryVertices());
    Real A = args["perturbationAmplitude"].as<Real>();
    Real f = args["perturbationFrequency"].as<Real>();
    for (auto bv : currentMesh.boundaryVertices()) {
        auto pt = bv.node().volumeNode()->p;
        Real a = A * cos(M_PI * f * pt[0]) * cos(M_PI * f * pt[1]);
        perturbation(bv.index()) = a * normals(bv.index());
    }
    bpi.paramsFromBoundaryVField(perturbation).getFlattened(perturbParams);

    // Perturbed meshes and simulators
    bpi.inflate(perturbParams);
    Simulator perturbed_sim(bpi.elements(), bpi.vertices());
    for (Real &p : perturbParams) p *= -1.0;
    bpi.inflate(perturbParams);
    Simulator neg_perturbed_sim(bpi.elements(), bpi.vertices());

    // Determine change in each vertex's position.
    VField delta_p(sim.mesh().numVertices());
    for (auto v : sim.mesh().vertices()) {
        delta_p(v.index()) = perturbed_sim.mesh().vertex(v.index()).node()->p;
        delta_p(v.index()) -= v.node()->p;
    }

    // Set up simulators' (base) material
    Materials::Constant<_N> &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    NonPeriodicCellOperations<Simulator> cell_operations(sim, bconds);
    NonPeriodicCellOperations<Simulator> perturbed_cell_operations(perturbed_sim, bconds);
    NonPeriodicCellOperations<Simulator> neg_perturbed_cell_operations(neg_perturbed_sim, bconds);

    cell_operations.m_solveCellProblems(sim, bconds);
    VField u = cell_operations.displacement();
    //auto delta_u = PeriodicHomogenization::deltaDisplacements(sim, u, delta_p); //TODO: implement forward version

    VField perturbed_u, neg_perturbed_u;
    perturbed_cell_operations.m_solveCellProblems(perturbed_sim, bconds);
    perturbed_u = perturbed_cell_operations.displacement();
    neg_perturbed_cell_operations.m_solveCellProblems(neg_perturbed_sim, bconds);
    neg_perturbed_u = neg_perturbed_cell_operations.displacement();

    auto delta_u = cell_operations.deltaDisplacements(u, delta_p);
    VField delta_u_forward_diff = perturbed_u;
    VField delta_u_centered_diff = perturbed_u;
    delta_u_forward_diff  -= u;
    delta_u_centered_diff -= neg_perturbed_u;
    delta_u_centered_diff *= 0.5;

    string output = args["output"].as<string>();
    bool linearSubsampleFields = args.count("fullDegreeFieldOutput") == 0;
    MSHFieldWriter writer(output, sim.mesh(), linearSubsampleFields);

    writer.addField("u", u);
    writer.addField("forward u", perturbed_u);
    writer.addField("backward u", neg_perturbed_u);
    writer.addField("delta u ", delta_u);
    writer.addField("forward difference delta u", delta_u_forward_diff);
    writer.addField("centered difference delta u", delta_u_centered_diff);

    auto origStrain = sim.averageStrainField(u);
    auto forwardDiffStrain   = perturbed_sim.averageStrainField(perturbed_u);
    auto centeredDiffStrain  = forwardDiffStrain;
    forwardDiffStrain  -= origStrain;
    centeredDiffStrain -= neg_perturbed_sim.averageStrainField(neg_perturbed_u);
    centeredDiffStrain *= 0.5;

    auto origNeumannLoad = sim.neumannLoad();
    auto deltaNeumannLoad = sim.deltaNeumannLoad(delta_p);
    auto perturbedNeumannLoad = perturbed_sim.neumannLoad();
    auto negPerturbedNeumannLoad = neg_perturbed_sim.neumannLoad();
    auto centeredDiffLoad = perturbedNeumannLoad - negPerturbedNeumannLoad;
    centeredDiffLoad *= 0.5;

    writer.addField("delta neumann load", deltaNeumannLoad);
    writer.addField("centered difference delta neumann load", centeredDiffLoad);

    auto origArea = sim.neumannBoundaryArea();
    auto perturbedArea = perturbed_sim.neumannBoundaryArea();
    auto negPerturbedArea = neg_perturbed_sim.neumannBoundaryArea();
    std::cout << "area: " << origArea << std::endl;
    std::cout << "centered difference delta area: " << 0.5*(perturbedArea - negPerturbedArea) << std::endl;

    writer.addField("strain u", origStrain);
    writer.addField("delta strain u", sim.deltaAverageStrainField(u, delta_u, delta_p));
    writer.addField("forward difference delta strain u",  forwardDiffStrain);
    writer.addField("centered difference delta strain u", centeredDiffStrain);

    auto origStress = sim.averageStressField(u);
    auto forwardDiffStress   = perturbed_sim.averageStressField(perturbed_u);
    auto centeredDiffStress  = forwardDiffStress;
    forwardDiffStress  -= origStress;
    centeredDiffStress -= neg_perturbed_sim.averageStressField(neg_perturbed_u);
    centeredDiffStress *= 0.5;

    writer.addField("stress u", origStress);
    writer.addField("forward difference delta stress u",  forwardDiffStress);
    writer.addField("centered difference delta stress u", centeredDiffStress);

    MSHFieldWriter perturbed_writer(output + ".perturbed.msh", perturbed_sim.mesh(), linearSubsampleFields);
    perturbed_writer.addField("u", perturbed_u);
    perturbed_writer.addField("strain u", perturbed_sim.averageStrainField(perturbed_u));

    MSHFieldWriter neg_perturbed_writer(output + ".neg_perturbed.msh", neg_perturbed_sim.mesh(), linearSubsampleFields);
    neg_perturbed_writer.addField("u", neg_perturbed_u);
    neg_perturbed_writer.addField("strain u", perturbed_sim.averageStrainField(neg_perturbed_u));

    auto origStressObjective          = buildStressObjective(args, sim, u, mat.getTensor());
    auto perturbedStressObjective     = buildStressObjective(args, perturbed_sim, perturbed_u, mat.getTensor());
    auto neg_perturbedStressObjective = buildStressObjective(args, neg_perturbed_sim, neg_perturbed_u, mat.getTensor());

    cout << "Stress:\t" << origStressObjective.evaluate() << endl;
    cout << "Perturbed Stress:\t" << perturbedStressObjective.evaluate() << endl;
    cout << "Neg Perturbed Stress:\t" << neg_perturbedStressObjective.evaluate() << endl;

    cout << "Forward  difference Stress:\t" << perturbedStressObjective.evaluate() - origStressObjective.evaluate() << endl;
    cout << "Centered difference Stress:\t" << 0.5 * (perturbedStressObjective.evaluate() - neg_perturbedStressObjective.evaluate()) << endl;

    const auto &mesh = sim.mesh();
    ShapeVelocityInterpolator interpolator(sim);
    VField bdry_svel(sim.mesh().numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices())
        bdry_svel(bv.index()) = delta_p(bv.volumeVertex().index());

    cout << "Forward shape derivative Stress:\t" << origStressObjective.deltaJ(sim, u, delta_p, cell_operations) << endl;
    OneForm<Real, _N> dJ = origStressObjective.adjointDeltaJ(cell_operations);
    cout << "Adjoint discrete shape derivative Stress (volume):\t" << dJ[interpolator.interpolate(sim, bdry_svel)] << endl;
    OneForm<Real, _N> dJbdry = interpolator.adjoint(sim, dJ);
    cout << "Adjoint discrete shape derivative Stress (boundary):\t" << dJbdry[bdry_svel] << endl;
    cout << "Adjoint discrete shape derivative Stress (delta p):\t" << dJ[delta_p] << endl;

    VField dJ_field = SDConversions::descent_from_diff_vol(dJ, sim);
    writer.addField("dJ field", dJ_field);
    VField dJbdry_field = SDConversions::descent_from_diff_bdry(dJbdry, sim);
    writer.addField("dJ bdry field", dJ_field);
    writer.addField("delta p", delta_p);
}

int main(int argc, const char *argv[])
{
    po::variables_map args = parseCmdLine(argc, argv);

    vector<MeshIO::IOVertex>  inVertices;
    vector<MeshIO::IOElement> inElements;
    string meshPath = args["mesh"].as<string>();
    auto type = load(meshPath, inVertices, inElements, MeshIO::FMT_GUESS,
                     MeshIO::MESH_GUESS);

    // Infer dimension from mesh type.
    size_t dim;
    if      (type == MeshIO::MESH_TET) dim = 3;
    else if (type == MeshIO::MESH_TRI) dim = 2;
    else    throw std::runtime_error("Mesh must be triangle or tet.");

    // Look up and run appropriate homogenizer instantiation.
    int deg = args["degree"].as<size_t>();
    auto exec = (dim == 3) ? ((deg == 2) ? execute<3, 2> : execute<3, 1>)
                           : ((deg == 2) ? execute<2, 2> : execute<2, 1>);

    cout << setprecision(19) << endl;

    exec(args, inVertices, inElements);

    return 0;
}
