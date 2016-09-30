////////////////////////////////////////////////////////////////////////////////
// WCS_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Compute the worst case stress of a periodic microstructure mesh.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/12/2016 16:42:04
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include <OrthotropicHomogenization.hh>
#include <MSHFieldWriter.hh>
#include <iomanip>

#include "WorstCaseStress.hh"

#include <vector>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace PH = PeriodicHomogenization;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: WCS_cli [options] mesh" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("mesh",       po::value<string>(),                     "input mesh")
        ;
    po::positional_options_description p;
    p.add("mesh",                1);

    po::options_description visible_opts;
    visible_opts.add_options()("help", "Produce this help message")
        ("material,m",               po::value<string>(),                 "base material")
        ("degree,d",                 po::value<int>()->default_value(2),  "degree of finite elements")
        ("fieldOutput,o",            po::value<string>(),                 "Dump fluctation stress and strain fields to specified msh file")
        ("orthotropicCell,O",                                             "Analyze the orthotropic symmetry base cell only")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;
    if (vm.count("mesh") == 0) {
        cout << "Error: must specify input mesh" << endl;
        fail = true;
    }

    int d = vm["degree"].as<int>();
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

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args,
             const vector<MeshIO::IOVertex> &inVertices, 
             const vector<MeshIO::IOElement> &inElements) {
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    using Mesh = LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Sim  = LinearElasticity::Simulator<Mesh>;
    Sim sim(inElements, inVertices);
    using ETensor = typename Sim::ETensor;
    using VField  = typename Sim::VField;
    // using WCSObjective = PthRootObjective<IntegratedWorstCaseObjective<Sim::N, WCStressIntegrandLp>>;
    const bool orthoCell = args.count("orthotropicCell");

    // Periodic Homogenization
    vector<VField> w;
    ETensor Eh;
    if (!orthoCell) {
        PH::solveCellProblems(w, sim, 1e-7);
        Eh = PH::homogenizedElasticityTensorDisplacementForm(w, sim);
    }
    else {
        PH::Orthotropic::solveCellProblems(w, sim, 1e-7);
        Eh = PH::Orthotropic::homogenizedElasticityTensorDisplacementForm(w, sim);
    }

    ETensor Sh = Eh.inverse();

    // Pointwise worst-case stress quantities are all computed the same way for
    // orthotropic and triply periodic cells
    auto m2m = PH::macroStrainToMicroStrainTensors(w, sim);
    auto wcs = worstCaseFrobeniusStress(mat.getTensor(), Sh, m2m);
    // WCSObjective wcsObjective(mesh, wcs);
    auto wcsValues = wcs.sqrtStressMeasure();

    Real maxMag = 0;
    size_t argmaxMag = 0;
    for (size_t i = 0; i < wcsValues.domainSize(); ++i) {
        Real val = std::abs(wcsValues[i]);
        if (val > maxMag) {
            maxMag = val;
            argmaxMag = i;
        }
    }

    typename Sim::SMatrix peakMacroStress = wcs.wcMacroStress(argmaxMag);

    cout << "Peak WCS Stress:\t" << maxMag << endl;
    cout << "Corresponding Macro Stress:" << endl;
    cout << peakMacroStress  << endl;
    cout << "Corresponding Macro strain:" << endl;
    cout << Sh.doubleContract(peakMacroStress) << endl;

    if (args.count("fieldOutput")) {
        MSHFieldWriter writer(args["fieldOutput"].as<string>(), sim.mesh());
        
        writer.addField("WC Macro Stress", wcs.wcMacroStress);
        writer.addField("WC Micro Stress", wcs.wcMicroStress());
        writer.addField("Pointwise WCS",   wcsValues);
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
    else    throw runtime_error("Mesh must be triangle or tet.");

    int deg = args["degree"].as<int>();
    auto exec = (dim == 3) ? ((deg == 2) ? execute<3, 2> : execute<3, 1>)
                           : ((deg == 2) ? execute<2, 2> : execute<2, 1>);

    cout << setprecision(19);

    exec(args, inVertices, inElements);

    return 0;
}
