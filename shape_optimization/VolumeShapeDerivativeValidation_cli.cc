////////////////////////////////////////////////////////////////////////////////
// VolumeShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validate the discrete shape derivative of volume quantities.
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
//  Created:  03/15/2018
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <objective_terms/TargetVolume.hh>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "../pattern_optimization/inflators/BoundaryPerturbationInflator.hh"

namespace po = boost::program_options;
using namespace std;

template<size_t _N> using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;
template<size_t _N> using ETensor = ElasticityTensor<Real, _N>;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: VolumeShapeDerivativeValidation_cli [options] mesh.msh" << endl;
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
            ("volWeight",    po::value<double>()->default_value(1.0),         "Weight for the volume target term of the objective function")
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
    visibleOptions.add(objectiveOptions).add(generalOptions);

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

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args,
             const std::vector<MeshIO::IOVertex>  &vertices,
             const std::vector<MeshIO::IOElement> &elements)
{
    using Mesh      = typename LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using VField    = typename Simulator::VField;

    // Original mesh and simulator
    BoundaryPerturbationInflator<_N> bpi(vertices, elements);
    vector<Real> perturbParams(bpi.numParameters(), 0.0);
    bpi.inflate(perturbParams);
    Simulator sim(bpi.elements(), bpi.vertices());

    // Perturb mesh
    auto &mesh = sim.mesh();
    auto normals = bpi.boundaryVertexNormals();
    VField perturbation(sim.mesh().numBoundaryVertices());
    Real A = args["perturbationAmplitude"].as<Real>();
    Real f = args["perturbationFrequency"].as<Real>();
    for (auto bv : mesh.boundaryVertices()) {
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

    string output = args["output"].as<string>();
    bool linearSubsampleFields = args.count("fullDegreeFieldOutput") == 0;
    MSHFieldWriter writer(output, sim.mesh(), linearSubsampleFields);

    Real volume = sim.mesh().volume();
    Real perturbed_volume = perturbed_sim.mesh().volume();
    Real neg_perturbed_volume = neg_perturbed_sim.mesh().volume();
    cout << "Volume:\t" << volume << endl;
    cout << "Perturbed Volume:\t" << perturbed_volume << endl;
    cout << "Neg Perturbed Volume:\t" << neg_perturbed_volume << endl;

    cout << "Forward  difference Volume:\t" << perturbed_volume - volume << endl;
    cout << "Centered difference Volume:\t" << 0.5 * (perturbed_volume - neg_perturbed_volume) << endl;

    ShapeVelocityInterpolator interpolator(sim);
    VField bdry_svel(sim.mesh().numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices())
        bdry_svel(bv.index()) = delta_p(bv.volumeVertex().index());

    //cout << "Shape derivative Volume:\t" << YYY << endl;
    OneForm<Real, _N> dV = sim.deltaVolumeForm();
    cout << "Shape derivative Volume (using one form):\t" << dV[delta_p] << endl;

    VField dV_field = SDConversions::descent_from_diff_vol(dV, sim);
    writer.addField("dV field", dV_field);
    writer.addField("delta p", delta_p);

    Real targetVol = 0.5;
    Real perturbed_cost = PatternOptimization::ObjectiveTerms::TargetVolumeTerm<Simulator>::compute_cost(targetVol, perturbed_volume);
    Real neg_perturbed_cost = PatternOptimization::ObjectiveTerms::TargetVolumeTerm<Simulator>::compute_cost(targetVol, neg_perturbed_volume);
    OneForm<Real, _N> differential = PatternOptimization::ObjectiveTerms::TargetVolumeTerm<Simulator>::compute_differential(targetVol, volume, sim, dV);
    cout << "Centered difference (Objective Term):\t" << 0.5 * (perturbed_cost - neg_perturbed_cost) << endl;
    cout << "Applied Differential (Objective Term):\t" << differential[bdry_svel] << endl;
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
    auto exec = (dim == 3) ? execute<3, 2> : execute<2, 2>;

    cout << setprecision(19) << endl;

    exec(args, inVertices, inElements);

    return 0;
}
