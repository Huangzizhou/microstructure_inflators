////////////////////////////////////////////////////////////////////////////////
// DiscreteShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validate the discrete shape derivative of periodic homogenization and
//      worst case stress quantities.
//
//      Boundary vertices are offset in the normal direction to create a
//      perturbed mesh, and post- and pre-perturbation quantities are
//      subtracted to compute finite difference (material) derivatives.
//
//      The BoundaryPerturbationInflator is used to ease the creation of a
//      perturbed mesh.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  03/16/2016 18:14:07
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include <GlobalBenchmark.hh>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "../pattern_optimization/BoundaryPerturbationInflator.hh"
#include "WCSOptimization.hh"
#include "WCStressOptimizationConfig.hh"

namespace po = boost::program_options;
using namespace std;
using namespace WCStressOptimization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: DiscreteShapeDerivativeValidation_cli [options] mesh.msh" << endl;
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
        ("WCSWeight",    po::value<double>()->default_value(1.0),         "Weight for the WCS term of the objective")
        ("JSWeight",     po::value<double>()->default_value(0.0),         "Weight for the JS term of the objective")
        ("JVolWeight",   po::value<double>()->default_value(0.0),         "Weight for the JVol term of the objective")
        ("LaplacianRegWeight,r", po::value<double>()->default_value(0.0), "Weight for the boundary Laplacian regularization term")
        ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
        ("material,m",   po::value<string>(),                    "Base material")
        ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                               "Produce this help message")
        ("output,o",     po::value<string>(),                    "Output the Lagrangian derivatives computed by forward difference and the discrete shape derivative.")
        ("fullDegreeFieldOutput,D",                              "Output full-degree nodal fields (don't do piecewise linear subsample)")
        ("perturbationAmplitude,a", po::value<double>()->default_value(0.01), "Amplitude of boundary perturbation")
        ("perturbationFrequency,f", po::value<double>()->default_value(1.0),  "Frequency of boundary perturbation")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(objectiveOptions).add(elasticityOptions).add(generalOptions);

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

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args,
             const std::vector<MeshIO::IOVertex>  &vertices,
             const std::vector<MeshIO::IOElement> &elements)
{
    using Mesh      = typename LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using ETensor   = typename Simulator::ETensor;
    using VField    = typename Simulator::VField;
    using Vector    = VectorND<_N>;

    // Original mesh and simulator
    BoundaryPerturbationInflator<_N> bpi(vertices, elements);
    vector<Real> perturbParams(bpi.numParameters());
    bpi.inflate(perturbParams);
    Simulator sim(bpi.elements(), bpi.vertices());

    // Perturb mesh
    {
        auto &mesh = sim.mesh();
        auto normals = bpi.normals();
        vector<Vector> perturbation(sim.mesh().numBoundaryVertices());
        Real A = args["perturbationAmplitude"].as<Real>();
        Real f = args["perturbationFrequency"].as<Real>();
        for (auto bv : mesh.boundaryVertices()) {
            auto pt = bv.node().volumeNode()->p;
            Real a = A * sin(M_PI * f * pt[0]) * sin(M_PI * f * pt[1]);
            perturbation[bv.index()] = a * normals(bv.index());
        }
        perturbParams = bpi.extractParamsFromBoundaryValues(perturbation);
    }

    // Perturbed mesh and simulator
    bpi.inflate(perturbParams);
    Simulator perturbed_sim(bpi.elements(), bpi.vertices());

    // Determine change in each vertex's position.
    VField delta_p(sim.mesh().numVertices());
    for (auto v : sim.mesh().vertices()) {
        delta_p(v.index()) = perturbed_sim.mesh().vertex(v.index()).node()->p;
        delta_p(v.index()) -= v.node()->p;
    }

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    // Configure WCS Objective
    auto &wcsConfig = WCStressOptimization::Config::get();
    wcsConfig.globalObjectivePNorm = args["pnorm"].as<double>();
    if (args.count("usePthRoot"))
        wcsConfig.globalObjectiveRoot = 2.0 * wcsConfig.globalObjectivePNorm;
    wcsConfig.useVtxNormalPerturbationGradientVersion = args.count("vtxNormalPerturbationGradient");

    // Create scalarized multi-objective with weights specified by the
    // arguments.
    WCStressOptimization::Objective<_N> fullObjective(ETensor(),
                                args[  "JSWeight"].as<double>(),
                                args[ "WCSWeight"].as<double>(),
                                args["JVolWeight"].as<double>(),
                                args["LaplacianRegWeight"].as<double>());

    std::vector<VField> w, delta_w;
    PeriodicHomogenization::solveCellProblems(w, sim);
    
    using SMatrix = typename Simulator::SMatrix;

    delta_w.reserve(w.size());
    for (size_t ij = 0; ij < w.size(); ++ij) {
        auto rhs = sim.deltaConstantStrainLoad(-SMatrix::CanonicalBasis(ij), delta_p);
        rhs     -= sim.applyDeltaStiffnessMatrix(w[ij], delta_p);
        delta_w.push_back(sim.solve(rhs));
    }

    std::vector<VField> perturbed_w;
    PeriodicHomogenization::solveCellProblems(perturbed_w, perturbed_sim);
    std::vector<VField> delta_w_finited_diff = perturbed_w;
    for (size_t ij = 0; ij < w.size(); ++ij)
        delta_w_finited_diff[ij] -= w[ij];

    string output = args["output"].as<string>();
    bool linearSubsampleFields = args.count("fullDegreeFieldOutput") == 0;
    MSHFieldWriter writer(output, sim.mesh(), linearSubsampleFields);

    for (size_t ij = 0; ij < w.size(); ++ij) {
        writer.addField("w " + std::to_string(ij), w[ij]);
        writer.addField("delta w " + std::to_string(ij), delta_w[ij]);
        writer.addField("finite difference delta w " + std::to_string(ij), delta_w_finited_diff[ij]);
    }

    for (size_t ij = 0; ij < w.size(); ++ij) {
        auto origStrain = sim.averageStrainField(w[ij]);
        auto finiteDiffStrain = perturbed_sim.averageStrainField(perturbed_w[ij]);
        for (size_t i = 0; i < finiteDiffStrain.domainSize(); ++i)
            finiteDiffStrain(i) -= origStrain(i);

        writer.addField("strain w " + std::to_string(ij), origStrain);
        writer.addField("delta strain w " + std::to_string(ij), sim.deltaAverageStrainField(w[ij], delta_w[ij], delta_p));
        writer.addField("finite difference delta strain w " + std::to_string(ij), finiteDiffStrain);
    }

    MSHFieldWriter perturbed_writer(output + ".perturbed.msh", perturbed_sim.mesh(), linearSubsampleFields);
    for (size_t ij = 0; ij < w.size(); ++ij) {
        perturbed_writer.addField("w "+ std::to_string(ij), perturbed_w[ij]);
        perturbed_writer.addField("strain w " + std::to_string(ij), perturbed_sim.averageStrainField(perturbed_w[ij]));
    }

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

    exec(args, inVertices, inElements);

    return 0;
}
