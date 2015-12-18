////////////////////////////////////////////////////////////////////////////////
// LpHoleGradientValidation.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validates a single component of the Lp hole pattern optimization
//      gradient. Does this by doing a fine 1D sweep of the objective around the
//      sample point, writing the objective and gradient component for each
//      sample.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/11/2015 02:25:06
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>

#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include <GlobalBenchmark.hh>

#include <LpHoleInflator.hh>
#include <PatternOptimizationJob.hh>
#include <PatternOptimizationConfig.hh>
#include "WCStressOptimizationIterate.hh"

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

namespace po = boost::program_options;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: LpHoleInflator_cli [options] radius p" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("radius", po::value<Real>(), "Radius parameter")
        ("p",      po::value<Real>(), "Lp hole p parameter")
        ("output",    po::value<string>(), "output path")
        ;
    po::positional_options_description p;
    p.add("radius", 1);
    p.add("p", 1);
    p.add("output", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",        "Produce this help message")
        ("material,m",      po::value<string>(), "base material")
        ("degree,d",        po::value<size_t>()->default_value(2),        "FEM Degree")
        ("nsubdiv,n",  po::value<size_t>()->default_value(64),            "number of subdivisions of Lp hole boundary")
        ("max_area,a", po::value<double>()->default_value(0.001),         "maximum triangle area for meshing")
        ("pnorm,P",    po::value<double>()->default_value(1.0),           "the pnorm used in the Lp global worst case stress measure")
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
    if (vm.count("radius") + vm.count("p")  + vm.count("output") != 3) {
        cout << "Error: must specify radius, p, and output path" << endl;
        fail = true;
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _FEMDegree>
void execute(const po::variables_map &args)
{
    constexpr size_t _N = 2;

    LpHoleInflator inflator;
    inflator.setMaxElementVolume(args["max_area"].as<double>());
    inflator.setNumSubdiv(args["nsubdiv"].as<size_t>());

    WCStressOptimization::Config::get().globalObjectivePNorm = args["pnorm"].as<double>();

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    typename Simulator::ETensor targetS(1.0, 0.0);

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    SField params(2);
    params[0] = args["radius"].as<Real>();
    params[1] = args["p"].as<Real>();

    WCStressOptimization::Iterate<Simulator> it(inflator, params.domainSize(), &params[0], targetS);
    it.writeMeshAndFields(args["output"].as<string>());

    BENCHMARK_REPORT();
}

int main(int argc, const char *argv[])
{
    cout << setprecision(20);

    po::variables_map args = parseCmdLine(argc, argv);
    size_t deg = args["degree"].as<size_t>();

    if (deg == 1) execute<1>(args);
    if (deg == 2) execute<2>(args);

    return 0;
}
