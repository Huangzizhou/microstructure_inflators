////////////////////////////////////////////////////////////////////////////////
// CreateJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
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
//  Created:  10/05/2014 12:59:00
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/MSHFieldParser.hh>
#include "PatternOptimizationJob.hh"
#include <MeshFEM/MSHFieldParser.hh>
#include <CSGFEM/utils.hh>

#include <boost/program_options.hpp>

#include <string>

using namespace std;

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "usage: ./CreateJob fields.msh outJob.opt [options]" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("fields", po::value<string>(), "msh with fields")
        ("job",    po::value<string>(), "out job file")
        ;
    po::positional_options_description p;
    p.add("fields", 1);
    p.add("job", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help", "Produce this help message")
        ("percentile,p", po::value<Real>(), "Fit percentile to select (default: 1.0)")
        ("element,e",    po::value<size_t>(), "Element to select")
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
    if (vm.count("job") == 0) {
        cout << "Error: must specify input fields.msh and output job.opt file" << endl;
        fail = true;
    }

    if (vm.count("percentile") && vm.count("element")) {
        cout << "Error: must specify only one of --percentile and --element" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
template<size_t _N>
using VField = VectorField<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _N>
void execute(const po::variables_map &args);
template<>
void execute<2>(const po::variables_map &args) {
    MSHFieldParser<2> parser(args["fields"].as<string>());

    SField lut_Ex    = parser.scalarField("fitted_young_x"),
           lut_Ey    = parser.scalarField("fitted_young_y"),
           lut_nu_yx = parser.scalarField("fitted_poisson_yx"),
           lut_mu    = parser.scalarField("shear_xy"),
           tgt_E     = parser.scalarField("young"),
           tgt_nu    = parser.scalarField("poisson");

    size_t optElement;
    std::vector<Real> dist;
    size_t nElems = parser.numElements();
    for (size_t i = 0; i < nElems; ++i) {
        ETensor<2> lutC, tgtC;
        lutC.setOrthotropic2D(lut_Ex[i], lut_Ey[i], lut_nu_yx[i], lut_mu[i]);
        tgtC.setIsotropic(tgt_E[i], tgt_nu[i]);
        auto diffS = lutC.inverse() - tgtC.inverse();
        dist.push_back(diffS.quadrupleContract(diffS));
    }

    if (args.count("element"))
        optElement = args["element"].as<size_t>();
    else {
        Real percentile = 1.0;
        if (args.count("percentile"))
            percentile = args["percentile"].as<Real>();
        vector<size_t> perm;
        sortPermutation(dist, perm, true); // Descending order sort permutation
        size_t idx = ceil(percentile * (dist.size() - 1));
        optElement = perm.at(idx);
    }

    if (optElement >= nElems) throw runtime_error("Illegal element index");

    ETensor<2> currentC, targetC;
    currentC.setOrthotropic2D(lut_Ex[optElement], lut_Ey[optElement],
                          lut_nu_yx[optElement], lut_mu[optElement]);
    targetC.setIsotropic(tgt_E[optElement], tgt_nu[optElement]);
    ETensor<2> targetS = targetC.inverse();

    cout << setprecision(16);

    cout << "Creating job to optimize fit on element " << optElement << endl
         << "initial distance:\t" << dist[optElement] << endl;

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
    SField params(9);
    params[1] = scale * parser.scalarField("vertex_orbit_0_thickness")[optElement];
    params[3] = scale * parser.scalarField("vertex_orbit_1_thickness")[optElement];
    params[0] = scale * parser.scalarField("vertex_orbit_2_thickness")[optElement];
    params[2] = scale * parser.scalarField("vertex_orbit_3_thickness")[optElement];
    params[4] = scale * parser.scalarField("vertex_orbit_4_thickness")[optElement];
    params[5] = 0.5 * (parser.scalarField("vertex_orbit_2_offset_0")[optElement] + parser.scalarField("vertex_orbit_2_offset_1")[optElement]);
    params[6] = 0.5 * (parser.scalarField("vertex_orbit_0_offset_0")[optElement] + parser.scalarField("vertex_orbit_0_offset_1")[optElement]);
    params[7] = 0.5 * (parser.scalarField("vertex_orbit_4_offset_0")[optElement] + parser.scalarField("vertex_orbit_4_offset_1")[optElement]) / sqrt(2);
    params[8] = params[7];

    string outPath(args["job"].as<string>());
    ofstream of(outPath);
    of << setprecision(16);
    if (!of.is_open()) throw runtime_error("Couldn't open file " + outPath);
    of << "{" << endl;
    of << "\t\"dim\": 2," << endl;
    of << "\t\"target\": {" << endl;
    of << "\t\t\"type\": \"isotropic\"," << endl;
    of << "\t\t\"young\": " << tgt_E[optElement]  << "," << endl;
    of << "\t\t\"poisson\": " << tgt_nu[optElement] << endl;
    of << "\t}," << endl;
    of << "\t\"initial_params\":"
       << " [" << params[1] 
       << ", " << params[3] 
       << ", " << params[0] 
       << ", " << params[2] 
       << ", " << params[4] 
       << ", " << params[5] 
       << ", " << params[6] 
       << ", " << params[7] 
       << ", " << params[8] << "]," << endl;
    of << "\t\"radiusBounds\": [0.05, 0.1]," << endl;
    of << "\t\"translationBounds\": [-0.15, 0.15]" << endl;
    of << "}" << endl;
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

    string mshPath(args["fields"].as<string>());

    vector<MeshIO::IOVertex>  inVertices;
    vector<MeshIO::IOElement> inElements;
    MeshIO::MeshIO_MSH io;
    ifstream infile(mshPath);
    if (!infile.is_open()) throw runtime_error("Couldn't open " + mshPath);
    auto type = io.load(infile, inVertices, inElements, MeshIO::MESH_GUESS);

    // Infer dimension from mesh type.
    size_t dim;
    if (type == MeshIO::MESH_QUAD) dim = 2;
    else    throw std::runtime_error("Only 2D quad input is supported.");

    if (dim == 2) execute<2>(args);
    return 0;
}
