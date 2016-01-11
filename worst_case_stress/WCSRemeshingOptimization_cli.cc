////////////////////////////////////////////////////////////////////////////////
// WCSRemeshingOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure boundary to bring the structure's
//      homogenized elasticity tensor closer to a target while also minimizing
//      worst-case stress.
//
//      The microstructure is remeshed after each step to prevent collapsed
//      triangles from blowing up the computation. Remeshing means that the
//      number of optimization variables (number of boundary vertices) may
//      change during the optimization, so fancier black-box optimizers cannot
//      be used.
//
//      Only 2D is supported for now.
//
//      TODO: subdivision of the inflated mesh for WCS/gradient computation.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/28/2015 18:11:43
////////////////////////////////////////////////////////////////////////////////
#include <GlobalBenchmark.hh>
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include <Future.hh>

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <Inflator.hh>
#include <PatternOptimizationJob.hh>
#include <PatternOptimizationConfig.hh>

#include "WCStressOptimizationIterate.hh"
#include "BoundaryPerturbationIterate.hh"

#include <LpHoleInflator.hh>
#include <BoundaryPerturbationInflator.hh>
#include <BoundaryPerturbationInflatorRemesher.hh>

namespace po = boost::program_options;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: PatternOptimization_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file")
        ;
    po::positional_options_description p;
    p.add("job", 1);

    po::options_description meshingOptions;
    meshingOptions.add_options()
        ("subdivide,S",  po::value<size_t>()->default_value(0),      "number of subdivisions to run before creating simulator")
        ("max_volume,v", po::value<double>()->default_value(0.01),   "maximum element area for remeshing")
        ("noRemesh",                                                 "disable remeshing")
        ("preRemeshOutput",                                          "Output the pre-remesh stats/msh file.")
        ("remeshMergeThreshold", po::value<double>(),                "edge merging threshold for remeshing")
        ("remeshSplitThreshold", po::value<double>(),                "edge splitting threshold for remeshing")
        ;

    po::options_description optimizerOptions;
    optimizerOptions.add_options()
        ("nIters,n",     po::value<size_t>()->default_value(20),     "number of iterations")
        ("step,s",       po::value<double>()->default_value(0.0001), "gradient step size")
        ("maxNormStep",  po::value<double>(),                        "Choose step size so step has specified Linv norm")
        ("vtxNormalPerturbationGradient,N",                          "use the vertex-normal-based version of the boundary perturbation gradient")
        ;

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
        ("ignoreShear",                                              "Ignore the shear components in the isotropic tensor fitting")
        ("alpha",        po::value<double>()->default_value(1.0),    "Trade-off between fitting and WCS minimization. 1.0 = tensor fit, 0.0 = WCS minimization.")
        ("pnorm,P",      po::value<double>()->default_value(1.0),    "the pnorm used in the Lp global worst case stress measure")
        ("JVolWeight,V", po::value<double>()->default_value(0.0),    "Weight of the volume constraint objective term.")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                                   "Produce this help message")
        ("pattern,p",    po::value<string>(),                        "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
        ("material,m",   po::value<string>(),                        "base material")
        ("degree,d",     po::value<size_t>()->default_value(2),      "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""),     "output .js mesh + fields at each iteration")
        ("outputEvery,O",po::value<size_t>()->default_value(1),      "output every O iterations (default n=1)")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(meshingOptions).add(optimizerOptions)
                  .add(objectiveOptions).add(generalOptions);

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
    if (vm.count("job") == 0) {
        cout << "Error: must specify input job.opt file" << endl;
        fail = true;
    }

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

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _FEMDegree>
void execute(const po::variables_map &args, const PatternOptimization::Job<2> *job) {
    constexpr size_t _N = 2;

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();
    Real targetSNormSq = targetS.quadrupleContract(targetS);
    std::cout << "targetSNormSq:\t" << targetSNormSq << std::endl;

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    cout << "target tensor: " << targetC << endl;

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    auto &wcsConfig = WCStressOptimization::Config::get();
    wcsConfig.globalObjectivePNorm = args["pnorm"].as<double>();
    wcsConfig.useVtxNormalPerturbationGradientVersion = args.count("vtxNormalPerturbationGradient");
    
    auto &poConfig = PatternOptimization::Config::get();
    poConfig.ignoreShear = args.count("ignoreShear");
    if (poConfig.ignoreShear) cout << "Ignoring shear components" << endl;
    if (args.count("remeshMergeThreshold")) poConfig.remeshMergeThreshold = args["remeshMergeThreshold"].as<Real>();
    if (args.count("remeshSplitThreshold")) poConfig.remeshSplitThreshold = args["remeshSplitThreshold"].as<Real>();

    Real   stepSize = args[      "step"].as<double>();
    string  outName = args[    "output"].as<string>();
    Real      alpha = args[     "alpha"].as<double>();
    size_t   niters = args[    "nIters"].as<size_t>();
    Real JVolWeight = args["JVolWeight"].as<double>();
    poConfig.fem2DSubdivRounds = args["subdivide"].as<size_t>();

    std::vector<MeshIO::IOVertex>  vertices;
    std::vector<MeshIO::IOElement> elements;
    MeshIO::load(args["pattern"].as<string>(), vertices, elements);

    using _Iterate = WCStressOptimization::BoundaryPerturbationIterate<Simulator>;
    using BPI = BoundaryPerturbationInflator<_N>;
    auto bpi = Future::make_unique<BPI>(vertices, elements);

    SField params(bpi->numParameters());
    params.clear();

    bool hasVolumeTarget = bool(job->targetVolume);
    Real targetVolSq = 0;
    if (hasVolumeTarget) {
        targetVolSq = *job->targetVolume;
        targetVolSq *= targetVolSq;
    }

    auto processIterate = [=](const _Iterate &it) ->SField {
        cout << "moduli:\t";
        it.elasticityTensor().printOrthotropic(cout);
        cout << "anisotropy:\t" << it.elasticityTensor().anisotropy() << endl;

        BENCHMARK_START_TIMER("Evaluate JS/WCS");
        Real JS = it.evaluateJS(), WCS = it.evaluateWCS();
        cout << "JS:\t"     << JS  << endl;
        cout << "WCS:\t"    << WCS << endl;
        if (hasVolumeTarget)
            cout << "JVol:\t"   << it.evaluateJVol()  << endl;
        cout << "Volume:\t" << it.mesh().volume() << endl;
        BENCHMARK_STOP_TIMER("Evaluate JS/WCS");

        SField   JS_p = it.gradp_JS();
        SField  WCS_p = it.gradientWCS_adjoint();
        SField JVol_p(WCS_p.domainSize());
        JVol_p.clear();

        cout << "||grad_p JS||:\t" << JS_p.norm() << endl;
        cout << "||grad_p WCS||:\t" << WCS_p.norm() << endl;;

        if (hasVolumeTarget) {
            JVol_p = it.gradp_JVol();
            cout << "||grad_p Jvol||:\t" << JVol_p.norm() << endl;
        }

        cout << "Normalized JS:\t" << it.evaluateJS() / targetSNormSq << endl;
        if (hasVolumeTarget)
            cout << "Normalized JVol:\t" << it.evaluateJVol() / targetVolSq << endl;

        // Compute composite objective and gradient
        Real J = JS * alpha + WCS * (1 - alpha);
        if (hasVolumeTarget)
            J += it.evaluateJVol() * JVolWeight;

        SField gradp = JS_p * alpha + WCS_p * (1 - alpha);
        if (hasVolumeTarget)
            gradp += JVol_p * JVolWeight;

        // Output composite iterate stats.
        cout << "J_full:\t" << J << endl;
        cout << "||grad_p J_full||:\t" << gradp.norm() << endl;
        return gradp;
    };

    // Custom gradient descent: number of parameters changes at each iteration.
    for (size_t i = 0; i < niters; ++i) {
        bool out = ((outName != "") &&
                    (i % args["outputEvery"].as<size_t>() == 0));
        SField gradp;

        std::unique_ptr<_Iterate> it;
        if (args.count("preRemeshOutput")) {
            it = Future::make_unique<_Iterate>(*bpi, params.size(), params.data(), targetS);
            if (hasVolumeTarget) it->setTargetVolume(*job->targetVolume);
            if (out) it->writeMeshAndFields(outName + "_" + std::to_string(i) + ".msh");
            std::cout << "Pre-remesh iterate:" << std::endl;
            gradp = processIterate(*it);
            std::cout << std::endl;
        }
        else {
            // For efficiency, don't homogenize the pre-remesh iterate.
            // (But we still need to inflate it so the remesher can read the
            //  perturbed boundary).
            std::vector<Real> pvector(params.domainSize());
            for (size_t i = 0; i < pvector.size(); ++i)
                pvector[i] = params[i];
            bpi->inflate(pvector);
        }

        if (!args.count("noRemesh")) {
            // post-remesh inflator and iterate
            remeshPerturbedShape(*bpi, args["max_volume"].as<Real>(), vertices, elements);
            bpi = Future::make_unique<BPI>(vertices, elements);
            params.resizeDomain(bpi->numParameters());
            params.clear();
            bpi->setNoPerturb(true); // Use mesh from Triangle directly.
            it = Future::make_unique<_Iterate>(*bpi, params.size(), params.data(), targetS);
            bpi->setNoPerturb(false);
            if (out) it->writeMeshAndFields(outName + "_" + std::to_string(i) + ".remeshed.msh");
            if (hasVolumeTarget) it->setTargetVolume(*job->targetVolume);

            std::cout << "Post-remesh iterate:" << std::endl;
            gradp = processIterate(*it);
            std::cout << std::endl;
        }

        if (args.count("maxNormStep"))
            stepSize = args["maxNormStep"].as<Real>() / abs(gradp.maxMag());

        params -= gradp * stepSize;

        // Apply bound constraints
        params.minRelax(job->translationBounds[1]);
        params.maxRelax(job->translationBounds[0]);
    }

    BENCHMARK_REPORT();
}


int main(int argc, const char *argv[])
{
    cout << setprecision(16);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = PatternOptimization::parseJobFile(args["job"].as<string>());

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<PatternOptimization::Job<2> *>(job)) {
        if (deg == 1) execute<1>(args, job2D);
        if (deg == 2) execute<2>(args, job2D);
    }
    else if (/* auto job3D = */ dynamic_cast<PatternOptimization::Job<3> *>(job)) {
        throw std::runtime_error("3D unsupported.");
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
