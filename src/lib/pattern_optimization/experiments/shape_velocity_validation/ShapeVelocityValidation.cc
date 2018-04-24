////////////////////////////////////////////////////////////////////////////////
// ShapeVelocityValidation.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validates the shape velocity by computing the shape derivative of the
//      trivial object volume functional:
//          dVol[v] = d[v] int_omega 1 dV = int_domega v.n dA
//      In other words, we consider the integral of the normal shape velocity
//      scalar field over the whole surface. By printing both this quantity and
//      the volume for a fine 1D sweep over pattern parameters, we can see how
//      well it approximates the true derivative of volume.
//  Output:
//      iterate    paramValue    volume    dVol
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  11/21/2015 14:04:43
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <MSHFieldWriter.hh>
#include <LinearElasticity.hh>
#include <GlobalBenchmark.hh>

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

#include "../../Inflator.hh"
#include "../../PatternOptimizationJob.hh"
#include "../../PatternOptimizationConfig.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: ShapeVelocityValidation [options] job.opt component_idx -l'lower_bd' -u'upper_bd' nsamples" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file")
        ("component_idx", po::value<size_t>(), "index of component to sweep")
        ("nsamples", po::value<size_t>(), "number of steps to sweep")
        ;
    po::positional_options_description p;
    p.add("job", 1);
    p.add("component_idx", 1);
    p.add("nsamples", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",      "Produce this help message")
        ("lower_bd,l", po::value<double>(), "sweep lower bound (must be positional to support negative values)")
        ("upper_bd,u", po::value<double>(), "sweep upper bound (must be positional to support negative values)")
        ("pattern,p",       po::value<string>(), "Pattern wire mesh (.obj|wire)")
        ("volumeMeshOut",   po::value<string>()->default_value(""),       "output volume mesh at each iteration")
        ("velocityOut,o",   po::value<string>(),                          "output shape velocity scalar field at each iteration")
        ("subdivide,S",  po::value<size_t>()->default_value(0),           "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"), "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                             "maximum element volume parameter for wire inflator")
        ("cell_size,c",  po::value<double>()->default_value(5.0),         "Inflation cell size (3D only)")
        ("directSDNormalVelocity,n",                                      "For isosurface inflator: use signed distance normal velocity directly (instead of multiplying by vertex normal . face normal)")
        ("isotropicParameters,I",                                         "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                             "Use vertex thickness instead of edge thickness (3D only)")
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
    if (vm.count("job") + vm.count("component_idx")  + vm.count("lower_bd") + vm.count("upper_bd") + vm.count("nsamples") != 5) {
        cout << "Error: must specify input job.opt file, sweep component index, sweep bounds, and sweep samples" << endl;
        fail = true;
    }

    if (vm.count("pattern") == 0) {
        cout << "Error: must specify pattern mesh" << endl;
        fail = true;
    }

    set<string> subdivisionAlgorithms = {"simple", "loop"};
    if (subdivisionAlgorithms.count(vm["sub_algorithm"].as<string>()) == 0) {
        cout << "Illegal subdivision algorithm specified" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
void execute(const po::variables_map &args, const Job<_N> *job)
{
    ConstrainedInflator<_N> inflator(job->parameterConstraints,
                args["pattern"].as<string>(),
                args["cell_size"].as<double>(),
                0.5 * sqrt(2),
                args.count("isotropicParameters"),
                args.count("vertexThickness"));
    if (args.count("max_volume"))
        inflator.setMaxElementVolume(args["max_volume"].as<double>());
    if (_N == 3) {
        inflator.configureSubdivision(args["sub_algorithm"].as<string>(),
                                      args["subdivide"].as<size_t>());
    }

    // Only relevant to isosurface inflator...
    Config::get().useSDNormalShapeVelocityDirectly = args.count("directSDNormalVelocity");

    if (job->numParams() != inflator.numParameters()) {
        for (size_t i = 0; i < inflator.numParameters(); ++i) {
            cout << "param " << i << " role: " <<
                (inflator.parameterType(i) == ParameterType::Offset ? "Offset" : "Thickness")
                << endl;
        }
        throw runtime_error("Invalid number of parameters.");
    }

    // Need to use LinearElasticity for the periodic cell boundary condtions--in
    // the future that code should probably be moved into the general LinearFEM
    // class.
    using Mesh      = LinearElasticity::Mesh<_N, 1, HMG>;
    using Simulator = LinearElasticity::Simulator<Mesh>;

    vector<Real> params(job->initialParams);
    string volumeMeshOut = args["volumeMeshOut"].as<string>();
    size_t compIdx = args["component_idx"].as<size_t>();
    assert(compIdx < params.size());
    double lowerBound = args["lower_bd"].as<double>();
    double upperBound = args["upper_bd"].as<double>();
    size_t nSamples = args["nsamples"].as<size_t>();

    for (size_t i = 0; i < nSamples; ++i) {
        params[compIdx] = lowerBound + ((nSamples == 1) ? 0.0
                        : (upperBound - lowerBound) * (double(i) / (nSamples - 1)));
        inflator.inflate(params);
        Simulator sim(inflator.elements(), inflator.vertices());
        sim.applyPeriodicConditions();
        const auto &mesh = sim.mesh();
        auto vn_p = inflator.computeShapeNormalVelocities(mesh);
        assert(vn_p.size() == params.size());
        const auto &vn = vn_p[compIdx];
        assert(vn.size() == mesh.numBoundaryElements());

        Real normalVelocityIntegral = 0;
        for (auto be : mesh.boundaryElements())
            normalVelocityIntegral += vn.at(be.index()).integrate(be->volume());

        cout << i << "\t"
             << params[compIdx] << "\t"
             << mesh.volume() << "\t"
             << normalVelocityIntegral << "\t"
             << endl;

        if (volumeMeshOut != "") MeshIO::save(volumeMeshOut + "_" + std::to_string(i) + ".msh", mesh);
        if (args.count("velocityOut")) {
            ScalarField<Real> pboundaryIndicator(mesh.numBoundaryElements());
            for (auto be : mesh.boundaryElements()) {
                pboundaryIndicator[be.index()] = be->isInternal ? 1.0 : 0.0;
            }
            MSHBoundaryFieldWriter writer(args["velocityOut"].as<string>() + "_" + std::to_string(i) + ".msh", mesh);
            writer.addField("vn", vn);
            writer.addField("periodic boundary", pboundaryIndicator);
        }
    }

    BENCHMARK_REPORT();
}

int main(int argc, const char *argv[])
{
    cout << setprecision(20);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = parseJobFile(args["job"].as<string>());

    if (auto job2D = dynamic_cast<Job<2> *>(job))
        execute<2>(args, job2D);
    else if (auto job3D = dynamic_cast<Job<3> *>(job))
        execute<3>(args, job3D);
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
