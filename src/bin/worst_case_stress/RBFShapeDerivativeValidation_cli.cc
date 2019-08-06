////////////////////////////////////////////////////////////////////////////////
// RBFShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//   Checks shape derivative and velocity by using changes of parameters of the RBF Inflator.
//
*/
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>

#include <inflators/wrappers/RBFInflator.hh>
#include <isosurface_inflator/ShapeVelocityInterpolator.hh>
#include <pattern_optimization/SDConversions.hh>

#include <CLI/CLI.hpp>

using namespace std;

template<size_t _N> using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;
template<size_t _N> using ETensor = ElasticityTensor<Real, _N>;

struct Args {
    std::string pngPath;
    std::string outputTable;
};

vector<Real> perturbOriginalParams(vector<Real> originalParams, int p, double perturbation = 1e-3) {
    vector<Real> result(originalParams);

    result[p] += perturbation;

    return result;
}

template<size_t _N, size_t _FEMDegree>
void execute(Args args)
{
    using Mesh      = typename LinearElasticity::Mesh<2, 2, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using ETensor   = typename Simulator::ETensor;
    using VField    = typename Simulator::VField;


    // Create inflator
    size_t d = 5;
    Real epsilon = d/2;
    RBFInflator inflator(args.pngPath, epsilon, d, d);

    // Create simulator
    std::vector<Real> params = inflator.defaultParameters();
    inflator.inflate(params);
    Simulator sim(inflator.elements(), inflator.vertices());


    // Compute volume and shape derivative form
    Real originalVolume = sim.mesh().volume();
    OneForm<Real, _N> diff_vol = sim.deltaVolumeForm();

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;

    // Compute volume velocities
    vector<VectorField<Real, _N>> vvels = inflator.volumeShapeVelocities();

    // gets material and original objective function
    const ETensor CBase = mat.getTensor();
    MSHFieldWriter comp_writer("rbfValidation.msh", sim.mesh(), true);

    auto dJ_field = SDConversions::descent_from_diff_vol(diff_vol, sim);
    comp_writer.addField("diff_vol", dJ_field, DomainType::PER_NODE);

    // Reads perturbations
    vector<Real> perturbations = {1e-6};

    ofstream * ofs;
    if (!args.outputTable.empty()) {
        ofs = new ofstream(args.outputTable, std::ofstream::out);
        (*ofs) << setw(7) << "Param";
        (*ofs) << setw(20) << "dJ[v]";
        for (auto perturbation : perturbations) {
            (*ofs) << setw(20) << perturbation;
        }
        (*ofs) << endl;
    }

    /* For each parameter, compute the velocity vector v generated from inflating with slighly different value
     * and use the dJ[v]. For comparisons, also compute the finite difference, simply by computing the J for initial
     * parameters and then J for the perturbed parameters. Compare the difference!! */
    for (size_t p = 0; p < params.size(); p++) {
        cout << "Running test for parameter: " << p << std::endl;

        if (!args.outputTable.empty())
            (*ofs) << setw(7) << p;

        comp_writer.addField("vvels-" + std::to_string(p), vvels[p], DomainType::PER_NODE);

        // compute velocity fields
        VField bdry_svel(sim.mesh().numBoundaryVertices());
        for (auto bv : sim.mesh().boundaryVertices())
            bdry_svel(bv.index()) = vvels[p](bv.volumeVertex().index());

        ShapeVelocityInterpolator interpolator(sim);
        OneForm<Real, _N> dJ = diff_vol;
        cout << "Adjoint discrete shape derivative Stress (volume):\t" << dJ[interpolator.interpolate(sim, bdry_svel)] << endl;
        OneForm<Real, _N> dJbdry = interpolator.adjoint(sim, dJ);
        Real dJv = dJbdry[bdry_svel];
        cout << "Adjoint discrete shape derivative Stress (boundary):\t" << dJv << endl;

        if (!args.outputTable.empty())
            (*ofs) << setw(20) << dJv;

        for (auto perturbation : perturbations) {
            cout << "Running test for perturbation: " << perturbation << std::endl;

            try {
                // Compute perturbed parameters and corresponding simulator
                vector<Real> perturbedParams = perturbOriginalParams(params, p, perturbation);
                inflator.inflate(perturbedParams);
                Simulator perturbed_sim(inflator.elements(), inflator.vertices());
                Real perturbedVolume = perturbed_sim.mesh().volume();

                // Compute negative perturbed parameters and corresponding simulator
                vector<Real> negPerturbedParams = perturbOriginalParams(params, p, -perturbation);
                inflator.inflate(negPerturbedParams);
                Simulator neg_perturbed_sim(inflator.elements(), inflator.vertices());
                Real negPerturbedVolume = neg_perturbed_sim.mesh().volume();

                Real finitDiff = (perturbedVolume - originalVolume) / (perturbedParams[p] - params[p]);
                Real centeredDiff = (perturbedVolume - negPerturbedVolume) / (perturbedParams[p] - negPerturbedParams[p]);
                cout << "Forward  difference Stress:\t" << finitDiff << endl;
                cout << "Centered  difference Stress:\t" << centeredDiff << endl;

                if (!args.outputTable.empty())
                    (*ofs) << setw(20) << centeredDiff;
            }
            catch (...) {
                if (!args.outputTable.empty())
                    (*ofs) << setw(20) << "fail";
            }
        }

        if (!args.outputTable.empty())
            (*ofs) << endl;
    }
}


int main(int argc, const char *argv[]) {
    Args args;

    // Parse arguments
    CLI::App app{"RBFShapeDerivativeValidation"};

    app.add_option("pngPath",  args.pngPath,  "png path")->required()->check(CLI::ExistingFile);
    app.add_option("--outputTable", args.outputTable, "output table");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    execute<2,2>(args);

    return 0;
}


