////////////////////////////////////////////////////////////////////////////////
// ScaleInvariantSmoothingShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validate the discrete shape derivative of smoothing objective term.
//
//      Boundary vertices are offset in the normal direction to create a
//      perturbed mesh, and post- and pre-perturbation quantities are
//      subtracted to compute forward/centered difference (material)
//      derivatives.
//
//      The BoundaryPerturbationInflator is used to ease the creation of a
//      perturbed mesh.
*/
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <pattern_optimization/objective_terms/ScaleInvariantSmoothingRegularization.hh>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <inflators/wrappers/BoundaryPerturbationInflator.hh>

#include <CLI/CLI.hpp>

namespace po = boost::program_options;
using namespace std;

template<size_t _N> using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;
template<size_t _N> using ETensor = ElasticityTensor<Real, _N>;

struct Args {
    std::string mesh;
    std::string output;
    double amplitude;
    double frequency;
};

template<size_t _N, size_t _FEMDegree>
void execute(Args args)
{
    using Mesh      = typename LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using VField    = typename Simulator::VField;
    using SField    = typename Simulator::SField;
    using ScaleInvariantSmoothingTerm = typename PatternOptimization::ObjectiveTerms::ScaleInvariantSmoothingRegularizationTerm<Simulator>;

    std::vector<MeshIO::IOVertex>  vertices;
    std::vector<MeshIO::IOElement> elements;

    string meshPath = args.mesh;
    auto type = load(args.mesh, vertices, elements, MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

    // Original mesh and simulator
    BoundaryPerturbationInflator<_N> bpi(vertices, elements, false);
    vector<Real> perturbParams(bpi.numParameters(), 0.0);
    bpi.inflate(perturbParams);
    Simulator sim(bpi.elements(), bpi.vertices());

    // Perturb mesh
    auto &mesh = sim.mesh();
    auto normals = bpi.boundaryVertexNormals();
    VField perturbation(sim.mesh().numBoundaryVertices());
    Real A = args.amplitude;
    Real f = args.frequency;
    for (auto bv : mesh.boundaryVertices()) {
        auto pt = bv.node().volumeNode()->p;
        Real a = A * cos(M_PI * f * pt[0]) * cos(M_PI * f * pt[1]);
        perturbation(bv.index()) = a * normals(bv.index());

        /*PointND<2> target;
        target << 0.0, 2.0;
        if ((target - pt).norm() < 1e-5) {
            //std::cout << "Point being moved: " << bv.node().volumeNode()->p;
            perturbation(bv.index()) = A * normals(bv.index());
        }
        else {
            perturbation(bv.index()) = 0 * normals(bv.index());
        }*/

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

    string output = args.output;
    MSHFieldWriter writer(output, sim.mesh());

    ShapeVelocityInterpolator interpolator(sim);
    VField bdry_svel(sim.mesh().numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices())
        bdry_svel(bv.index()) = delta_p(bv.volumeVertex().index());

    OneForm<Real, _N> differential = ScaleInvariantSmoothingTerm::computeDifferential(sim);
    OneForm<Real, _N> volumeDifferential = ScaleInvariantSmoothingTerm::computeVolumeDifferential(sim);

    /*std::cout << "Volume differential: " << volumeDifferential.domainSize() << std::endl;
    for (size_t i=0; i < volumeDifferential.domainSize(); i++) {
        std::cout << volumeDifferential(i)[0] << ", " << volumeDifferential(i)[1] << std::endl;
    }*/

    VField dV_field = SDConversions::descent_from_diff_vol(volumeDifferential, sim);
    writer.addField("dV field", dV_field, DomainType::PER_NODE);

    /*for (size_t i=0; i < dV_field.domainSize(); i++) {
        std::cout << dV_field(i)[0] << ", " << dV_field(i)[1] << std::endl;
    }
    for (auto v : sim.mesh().vertices()) {
        std::cout << "Point being moved from " << v.node()->p[0] << ", " << v.node()->p[1] << " by " << delta_p(v.index())[0] << ", " << delta_p(v.index())[1] << std::endl;
    }*/

    writer.addField("delta p", delta_p, DomainType::PER_NODE);

    SField costField = ScaleInvariantSmoothingTerm::computeSmoothingField(sim);
    writer.addField("smoothing cost", costField, DomainType::PER_NODE);

    Real cost = ScaleInvariantSmoothingTerm::computeCost(sim);
    Real perturbed_cost = ScaleInvariantSmoothingTerm::computeCost(perturbed_sim);
    Real neg_perturbed_cost = ScaleInvariantSmoothingTerm::computeCost(neg_perturbed_sim);

    std::cout << "Cost: " << cost << std::endl;
    std::cout << "Perturbed Cost: " << perturbed_cost << std::endl;
    std::cout << "Negative Perturbed Cost: " << neg_perturbed_cost << std::endl;

    cout << "Centered difference (Objective Term):\t" << 0.5 * (perturbed_cost - neg_perturbed_cost) << endl;
    cout << "Applied Differential (Objective Term):\t" << differential[bdry_svel] << endl;
    cout << "Volume Differential (Objective Term):\t" << volumeDifferential[delta_p] << endl;
}

int main(int argc, const char *argv[])
{
    Args args;
    args.amplitude = 0.0001;
    args.frequency = 1.0;

    // Parse arguments
    CLI::App app{"FancySmoothingShapeDerivativeValidation"};

    app.add_option("mesh",           args.mesh,         "mesh path")->required()->check(CLI::ExistingFile);
    app.add_option("--output",       args.output,       "output");
    app.add_option("--amplitude",    args.amplitude,    "amplitude used in finite differences check");
    app.add_option("--frequency",    args.frequency,    "frequency used in finite differences check");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    execute<2,2>(args);

    return 0;
}
