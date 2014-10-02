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
#include <EdgeFields.hh>
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
#include <set>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "PatternOptimization.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

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
        ("pattern,p",   po::value<string>(), "Pattern wire mesh (.obj)")
        ("material,m",  po::value<string>(), "base material")
        ("output,o",    po::value<string>(), "output .js mesh + fields")
        ("max_area,a",  po::value<double>()->default_value(0.0001), "max_area parameter for wire inflator")
        ("solver",      po::value<string>()->default_value("gradient_descent"), "solver to use: gradient_descent, bfgs, lbfgs, levenberg_marquardt")
        ("step,s",      po::value<double>()->default_value(0.0001), "gradient step size")
        ("nIters,n",    po::value<size_t>()->default_value(20), "number of iterations")
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

    set<string> solvers = {"gradient_descent", "bfgs", "lbfgs", "levenberg_marquardt"};
    if (solvers.count(vm["solver"].as<string>()) == 0) {
        cout << "Illegal solver specified" << endl;
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

    // Focus on the worst fitted tensor in the object
    size_t optElement = inElements.size();
    Real worstDist = 0;
    ETensor<2> currentC, targetC;
    for (size_t i = 0; i < inElements.size(); ++i) {
        ETensor<2> lutC, tgtC;
        lutC.setOrthotropic2D(lut_Ex[i], lut_Ey[i], lut_nu_yx[i], lut_mu[i]);
        tgtC.setIsotropic(tgt_E[i], tgt_nu[i]);
        auto diffS = lutC.inverse() - tgtC.inverse();
        Real dist = diffS.quadrupleContract(diffS);
        if (dist > worstDist) {
            worstDist = dist;
            currentC = lutC;
             targetC = tgtC;
            optElement = i;
        }
    }
    ETensor<2> targetS = targetC.inverse();

    cout << setprecision(16) << endl;

    cout << "Optimizing worst-fit on element " << optElement << endl
         << "initial distance:\t" << worstDist << endl;

    cout << "LUT tensors:" << endl << currentC << endl << endl << currentC.inverse() << endl << endl;
    cout << "Target tensor:" << endl << targetC << endl << endl << targetS << endl << endl;
    cout << "LUT moduli: "
         << lut_Ex[optElement] << ", " << lut_Ey[optElement] << ", "
         << lut_nu_yx[optElement] << ", " << lut_mu[optElement] << endl;
    cout << "Target moduli (table): " << tgt_E[optElement] << ", "
         << tgt_nu[optElement] << endl;

    Real tgt_Ex, tgt_Ey, tgt_nuyx, tgt_mu;
    targetC.getOrthotropic2D(tgt_Ex, tgt_Ey, tgt_nuyx, tgt_mu);
    cout << "Target moduli (tensor):"
         << "\t" << tgt_Ex << "\t" << tgt_Ey << "\t" << tgt_nuyx
         << "\t" << tgt_mu << endl;

	WireInflator2D inflator(args["pattern"].as<string>());
    TessellationParameters t_params;
    t_params.max_area = args["max_area"].as<double>();

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
    size_t nParams = inflator.patternGenerator().numberOfParameters();
    SField params(nParams);
    params[1] = scale * parser.scalarField("vertex_orbit_0_thickness")[optElement];
    params[3] = scale * parser.scalarField("vertex_orbit_1_thickness")[optElement];
    params[0] = scale * parser.scalarField("vertex_orbit_2_thickness")[optElement];
    params[2] = scale * parser.scalarField("vertex_orbit_3_thickness")[optElement];
    params[4] = scale * parser.scalarField("vertex_orbit_4_thickness")[optElement];
    params[5] = 0.5 * (parser.scalarField("vertex_orbit_2_offset_0")[optElement] + parser.scalarField("vertex_orbit_2_offset_1")[optElement]);
    params[6] = 0.5 * (parser.scalarField("vertex_orbit_0_offset_0")[optElement] + parser.scalarField("vertex_orbit_0_offset_1")[optElement]);
    params[7] = 0.5 * (parser.scalarField("vertex_orbit_4_offset_0")[optElement] + parser.scalarField("vertex_orbit_4_offset_1")[optElement]) / sqrt(2);
    params[8] = params[7];

    // Set up material
    auto &mat = LinearElasticityND<_N>::
        template homogenousMaterial<Materials::Constant>();
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    Optimizer<_N> optimizer(inflator, t_params);
    string solver = args["solver"].as<string>();
    if (solver == "levenberg_marquardt")
        optimizer.optimize_lm(params, targetS, args["output"].as<string>());
    else if (solver == "gradient_descent")
        optimizer.optimize_gd(params, targetS, args["nIters"].as<size_t>(),
                          args["step"].as<double>(), args["output"].as<string>());
    else if (solver == "bfgs")
        optimizer.optimize_bfgs(params, targetS, args["output"].as<string>());
    else if (solver == "lbfgs")
        optimizer.optimize_bfgs(params, targetS, args["output"].as<string>(), 10);
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
