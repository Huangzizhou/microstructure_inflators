#include <iostream>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>

#include <CLI/CLI.hpp>

using namespace std;
using namespace Eigen;

// Parser of parameters
int parseCmdLine(int argc, char *argv[], string &input, string &output) {

    CLI::App app{"Glue and clean mesh"};

    app.add_option("input", input, "input triangle mesh file (.obj)")->required()->check(CLI::ExistingFile);
    app.add_option("output", output, "output glued and clean mesh (.obj)")->required();

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    return 0;
}


int main(int argc, char *argv[])
{
    string input;
    string output;
    parseCmdLine(argc, argv, input, output);

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    // Load mesh
    igl::readOBJ(input, V, T);

    Eigen::MatrixXd SV;
    Eigen::MatrixXi SVI;
    Eigen::MatrixXi SVJ;
    Eigen::MatrixXi ST;
    igl::remove_duplicate_vertices(V, T, 1e-5, SV, SVI, SVJ, ST);
    
    Eigen::MatrixXd NV;
    Eigen::MatrixXi NT;
    Eigen::MatrixXi IM;
    igl::remove_unreferenced(SV, ST, NV, NT, IM);

    igl::writeOBJ(output, NV, NT);

    return 0;
}
