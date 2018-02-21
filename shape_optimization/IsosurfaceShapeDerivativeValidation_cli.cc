////////////////////////////////////////////////////////////////////////////////
// IsosurfaceShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//   Checks shape derivative using changes of parameters of the Isosurface Inflator.
//
*/
//  Author:  Davi Colli Tozoni (dctozoni) davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  1/13/18
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>

#include <set>

#include <json.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <Inflator.hh>
#include <MakeInflator.hh>
#include <PatternOptimizationJob.hh>

#include <IterateFactory.hh>
#include <IterateManager.hh>

#include "ShapeOptimizationIterate.hh"

#include "StressObjectiveTerm.hh"
#include "ParametersMask.h"
#include "inflators/IsoinflatorWrapper.hh"
#include "../pattern_optimization/ShapeVelocityInterpolator.hh"

namespace po = boost::program_options;
using namespace std;

template<size_t _N> using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;
template<size_t _N> using ETensor = ElasticityTensor<Real, _N>;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: IsosurfaceShapeDerivativeValidation_cli [options]" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description patternOptions;
    patternOptions.add_options()
            ("pattern,p",    po::value<string>(),                              "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
            ("params",       po::value<string>(),                              "Initial params (overrides those specified in job file).")
            ;

    po::options_description simulationOptions;
    simulationOptions.add_options()
            ("boundaryConditions,b", po::value<string>(),                    "boundary conditions")
            ("perturbations", po::value<string>()->default_value("1e-3"),    "perturbations to be considered in tests")
            ;

    po::options_description meshingOptions;
    meshingOptions.add_options()
            ("meshingOptions,M",      po::value<string>(), "Meshing options configuration file")
            ;

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
            ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp stress measure")
            ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
            ("stressWeight",    po::value<double>()->default_value(1.0),         "Weight for the Microscopic stress term of the objective")
            ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
            ("material,m",   po::value<string>(),                    "Base material")
            ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
            ;

    po::options_description generalOptions;
    generalOptions.add_options()
            ("help,h",                                    "Produce this help message")
            ("output,o",     po::value<string>(),         "Output")
            ("outputTable,o",     po::value<string>(),         "Output table")
            ;

    po::options_description visibleOptions;
    visibleOptions.add(patternOptions).add(simulationOptions).add(meshingOptions)
            .add(objectiveOptions).add(elasticityOptions).add(generalOptions);

    po::options_description cli_opts;
    cli_opts.add(visibleOptions);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(cli_opts).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visibleOptions);
    }

    bool fail = false;
    if (vm.count("pattern") == 0) {
        cout << "Error: must specify pattern mesh" << endl;
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

// Parse parameters
vector<Real> parseParams(string pstring) {
    boost::trim(pstring);
    vector<string> tokens;
    boost::split(tokens, pstring, boost::is_any_of("\t "), boost::token_compress_on);
    vector<Real> pvals;
    for (string &s : tokens) pvals.push_back(std::stod(s));
    return  pvals;
}

// Filter parameters
vector<Real> originalToOptimizedParameters(std::vector<bool> parametersMask, vector<Real> originalParams) {
    int originalIdx = 0;
    vector<Real> result;

    for (; originalIdx < originalParams.size(); originalIdx++) {
        if (!parametersMask[originalIdx]) {
            result.push_back(originalParams[originalIdx]);
        }
    }

    return result;
}

string paramsToString(vector<Real> params) {
    std::string result;

    for (unsigned i=0; i<(params.size()-1); i++) {
        result += std::to_string(params[i]) + string(", ");
    }
    result += std::to_string(params[params.size()-1]);

    return result;
}

// Deal with parameters mask
string paramsMaskToString(vector<bool> mask) {
    std::string result;

    for (unsigned i=0; i<(mask.size()-1); i++) {
        if (mask[i])
            result += "1, ";
        else
            result += "0, ";
    }
    if (mask[mask.size()-1])
        result += "1";
    else
        result += "0";

    return result;
}


template<size_t _N, class Simulator, class VField>
PthRootObjective<IntegratedMicroscopicStressObjective<_N, MicroscopicStressIntegrandLp<Simulator>, Simulator>> buildStressObjective(po::variables_map &args, const Simulator &sim, VField u, const ETensor<_N> CBase) {
    PthRootObjective<IntegratedMicroscopicStressObjective<_N, MicroscopicStressIntegrandLp<Simulator>, Simulator>> objective;

    objective.integrand.p = args["pnorm"].as<double>();
    objective.p = args.count("usePthRoot") ? 2.0 * args["pnorm"].as<double>() : 1.0;

    const bool major_symmetry = CBase.MajorSymmetry;

    MicroscopicStress<_N, Simulator> microStress = MicroscopicFrobeniusStress<CBase.Dim, major_symmetry, Simulator>(CBase, sim.averageStressField(u));
    objective.setPointwiseStress(sim.mesh(), microStress);
    return objective;
};

vector<Real> perturbOriginalParams(vector<Real> originalParams, int p, double perturbation = 1e-3) {
    vector<Real> result(originalParams);

    result[p] += perturbation;

    return result;
}

template<size_t _N, size_t _FEMDegree>
void execute(po::variables_map &args)
{
    using Mesh      = typename LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using ETensor   = typename Simulator::ETensor;
    using VField    = typename Simulator::VField;
    using Vector    = VectorND<_N>;

    string bcondsPath = args["boundaryConditions"].as<string>();

    // If requested, override the initial parameters set in the job file
    vector<double> originalParams = parseParams(args["params"].as<string>());
    vector<bool> paramsMask = ParametersMask::generateParametersMask(args["pattern"].as<string>(), args["params"].as<string>(), args["boundaryConditions"].as<string>());

    assert(originalParams.size() == paramsMask.size());

    /*for (unsigned i = 0; i < paramsMask.size(); i++) {
        if (paramsMask[i])
            cout << i << " " << endl;
    }*/

    vector<Real> params = originalToOptimizedParameters(paramsMask, originalParams);

    // Create inflator
    std::vector<std::string> parameterConstraints;
    string symmetry_name = "non_periodic";
    //IsoinflatorWrapper<_N> inflator(args["pattern"].as<string>(), "non_periodic", true, paramsMask, originalParams, 2);
    args.insert(std::make_pair("symmetry", po::variable_value(symmetry_name, true)));
    args.insert(std::make_pair("vertexThickness", po::variable_value(true, true)));
    args.insert(std::make_pair("params", po::variable_value(paramsToString(originalParams), true)));
    args.insert(std::make_pair("paramsMask", po::variable_value(paramsMaskToString(paramsMask), true)));
    po::notify(args);
    std::unique_ptr<InflatorBase> inflatorBase =make_inflator("Isosurface2D", filterInflatorOptions(args), parameterConstraints);
    InflatorBase* temp = inflatorBase.get();
    IsoinflatorWrapper<_N>* inflator = (IsoinflatorWrapper<_N> *) temp;

    // Create simulator
    inflator->inflate(params);
    Simulator sim(inflator->elements(), inflator->vertices());

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    // Read boundary conditions
    bool no_rigid_motion;
    vector<CondPtr<_N> > bconds = readBoundaryConditions<_N>(bcondsPath, sim.mesh().boundingBox(), no_rigid_motion);

    bool linearSubsampleFields = true; // TODO: what does it mean?

    // Compute volume velocities
    vector<VectorField<Real, _N>> vvels = inflator->volumeShapeVelocities();
    //vector<VectorField<Real, _N>> bvels = inflator->shapeVelocities(vvels);

    // Find displacement for original mesh
    NonPeriodicCellOperations<Simulator> cell_operations(sim, bconds);
    cell_operations.m_solveCellProblems(sim, bconds);
    VField u = cell_operations.displacement();

    // gets material and original objective function
    const ETensor CBase = mat.getTensor();
    auto origStressObjective = buildStressObjective(args, sim, u, CBase);

    MSHFieldWriter comp_writer("isosurfaceValidation.msh", sim.mesh(), true);
    comp_writer.addField("Pointwise Stress", origStressObjective.microStress.sqrtStressMeasure());

    auto diff_vol = origStressObjective.adjointDeltaJ(cell_operations);
    auto dJ_field = cell_operations.descent_from_diff_vol(diff_vol);
    comp_writer.addField("diff_vol", dJ_field, DomainType::PER_NODE);

    // Reads perturbations
    vector<Real> perturbations = parseParams(args["perturbations"].as<string>());

    ofstream * ofs;
    if (args.count("outputTable") > 0) {
        ofs = new ofstream(args["outputTable"].as<string>(), std::ofstream::out);
        (*ofs) << setw(7) << "Param";
        (*ofs) << setw(20) << "dJ[v]";
        for (auto perturbation : perturbations) {
            (*ofs) << setw(20) << perturbation;
        }
        (*ofs) << endl;

        //(*ofs) << setw(7) << "";
        //for (auto perturbation : perturbations) {
        //    (*ofs) << setw(10) << "FinDiff";
        //}
        //(*ofs) << endl;
        //cout << setw(5) << 1 << setw(15) << "January" << setw(15) << "Abhilash" << endl;
    }

    // For next computations, use cheap mode:
    inflator->enableCheapPostprocess();

    /* For each remaining parameter, compute the velocity vector v generated from inflating with slighly different value
     * and use the dJ[v]. For comparisons, also compute the finite difference, simply by computing the J for initial
     * parameters and then J for the perturbed parameters. Compare the difference!! */
    for (unsigned p = 0; p < params.size(); p++) {
        cout << "Running test for parameter: " << p << std::endl;

        if (args.count("outputTable") > 0)
            (*ofs) << setw(7) << p;

        // compute velocity fields
        VField bdry_svel(sim.mesh().numBoundaryVertices());
        for (auto bv : sim.mesh().boundaryVertices())
            bdry_svel(bv.index()) = vvels[p](bv.volumeVertex().index());

        ShapeVelocityInterpolator interpolator(sim);
        OneForm<Real, _N> dJ = origStressObjective.adjointDeltaJ(cell_operations);
        cout << "Adjoint discrete shape derivative Stress (volume):\t"
             << dJ[interpolator.interpolate(sim, bdry_svel)] << endl;
        OneForm<Real, _N> dJbdry = interpolator.adjoint(sim, dJ);
        Real dJv = dJbdry[bdry_svel];
        cout << "Adjoint discrete shape derivative Stress (boundary):\t" << dJv << endl;

        //VField dJ_field = SDConversions::descent_from_diff_vol(dJ, sim);
        //MSHFieldWriter dJ_writer("isosurface.dJ.msh", sim.mesh(), true);
        //dJ_writer.addField("dJ field", dJ_field);

        if (args.count("outputTable") > 0)
            (*ofs) << setw(20) << dJv;

        for (auto perturbation : perturbations) {
            cout << "Running test for perturbation: " << perturbation << std::endl;

            try {
                // Compute perturbed parameters and corresponding simulator
                vector<Real> perturbedParams = perturbOriginalParams(params, p, perturbation);
                inflator->inflate(perturbedParams);
                Simulator perturbed_sim(inflator->elements(), inflator->vertices());

                // Compute negative perturbed parameters and corresponding simulator
                vector<Real> negPerturbedParams = perturbOriginalParams(params, p, -perturbation);
                inflator->inflate(negPerturbedParams);
                Simulator neg_perturbed_sim(inflator->elements(), inflator->vertices());

                // Create helper for non periodic cell simulation computations
                NonPeriodicCellOperations<Simulator> perturbed_cell_operations(perturbed_sim, bconds);
                NonPeriodicCellOperations<Simulator> neg_perturbed_cell_operations(neg_perturbed_sim, bconds);

                // Find displacement for perturbed mesh
                perturbed_cell_operations.m_solveCellProblems(perturbed_sim, bconds);
                VField perturbed_u = perturbed_cell_operations.displacement();

                neg_perturbed_cell_operations.m_solveCellProblems(neg_perturbed_sim, bconds);
                VField neg_perturbed_u = neg_perturbed_cell_operations.displacement();


                auto origStrain = sim.averageStrainField(u);
                auto perturbedStrain = perturbed_sim.averageStrainField(perturbed_u);
                auto negPerturbedStrain = neg_perturbed_sim.averageStrainField(neg_perturbed_u);

                // set output name to include parameter number
                if (args.count("output") > 0) {
                    string output = args["output"].as<string>() + string(" ") + to_string(p);
                    MSHFieldWriter writer(output, sim.mesh(), linearSubsampleFields);

                    // add displacement field and perturbed field
                    writer.addField("u", u);
                    writer.addField("perturbed u", perturbed_u);
                    writer.addField("neg perturbed u", neg_perturbed_u);

                    // add also the strain values
                    writer.addField("strain u", origStrain);
                    writer.addField("neg perturbed u", negPerturbedStrain);
                }

                auto perturbedStressObjective = buildStressObjective(args, perturbed_sim, perturbed_u, CBase);
                auto negPerturbedStressObjective = buildStressObjective(args, neg_perturbed_sim, neg_perturbed_u, CBase);

                //cout << "Stress:\t" << origStressObjective.evaluate() << endl;
                //cout << "Perturbed Stress:\t" << perturbedStressObjective.evaluate() << endl;
                //cout << "Neg Perturbed Stress:\t" << negPerturbedStressObjective.evaluate() << endl;

                Real finitDiff = (perturbedStressObjective.evaluate() - origStressObjective.evaluate()) /
                                 (perturbedParams[p] - params[p]);
                Real centeredDiff = (perturbedStressObjective.evaluate() - negPerturbedStressObjective.evaluate()) /
                                 (perturbedParams[p] - negPerturbedParams[p]);
                cout << "Forward  difference Stress:\t" << finitDiff << endl;
                cout << "Centered  difference Stress:\t" << centeredDiff << endl;

                if (args.count("outputTable") > 0)
                    (*ofs) << setw(20) << centeredDiff;
            }
            catch (...) {
                if (args.count("outputTable") > 0)
                    (*ofs) << setw(20) << "fail";
            }
        }

        if (args.count("outputTable") > 0)
            (*ofs) << endl;
    }
}

int main(int argc, const char *argv[]) {
    po::variables_map args = parseCmdLine(argc, argv);

#if HAS_TBB
    size_t np = tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(np);
#endif

    cout << setprecision(16);

    size_t deg = args["degree"].as<size_t>();
    bool is2D = true; //TODO: fix for 3D cases
    if (is2D) {
        if (deg == 1) execute<2, 1>(args);
        if (deg == 2) execute<2, 2>(args);
    }
    else {
        if (deg == 1) execute<3, 1>(args);
        if (deg == 2) execute<3, 2>(args);
    }

    return 0;
}


