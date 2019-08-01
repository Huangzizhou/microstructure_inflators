#include <bitset>
#include <string>

#include <CLI/CLI.hpp>
#include <inflators/wrappers/RBFInflator.hh>

using namespace std;

struct Args {
    std::string inputPath;
    std::string outputPath;
    std::string moptsPath;
};

void parseRBFFile(string path, double &epsilon, size_t &d1, size_t &d2, std::vector<double> &coeffs) {
    ifstream infile(path);

    string line;
    size_t idx = 0;
    while (getline(infile, line))
    {
        std::istringstream l(line);
        if (idx == 0) {
            if (line.compare("rbf") == 0) {
                cout << "RBF" << endl;
            }
            else {
                throw runtime_error("Method " + line + "is not implemented");
            }
        }
        else if (idx == 1) {
            epsilon = stof(line);
            cout << "epsilon: " << epsilon << endl;
        }
        else if (idx == 2) {
            l >> d1 >> d2;

            cout << "d1: " << d1 << endl;
            cout << "d2: " << d2 << endl;
        }
        else {
            size_t total = d1 * d2;
            double value;
            size_t so_far = 0;
            coeffs.clear();

            std::cout << "coeffs: ";
            while (l >> value && so_far < total) {
                coeffs.push_back(value);
                cout << value << " ";
                so_far++;
            }
            cout << endl;

            if (so_far != total) {
                throw runtime_error("Not the expected number of coefficients");
            }
        }

        idx++;
    }
}

int main(int argc, char ** argv) {
    Args args;

    // Parse arguments
    CLI::App app{"rbf2mesh"};

    app.add_option("inputPath",  args.inputPath,  "input  path")->required()->check(CLI::ExistingFile);
    app.add_option("outputPath", args.outputPath, "output path")->required();
    app.add_option("--mopts",    args.moptsPath,  "mesh options");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    // read input coefficient file
    double epsilon;
    size_t d1, d2;
    vector<double> coeffs;
    parseRBFFile(args.inputPath, epsilon, d1, d2, coeffs);

    RBFInflator inflator(epsilon, d1, d2);
    if (!args.moptsPath.empty())  inflator.meshingOptions().load(args.moptsPath);

    inflator.inflate(coeffs);

    MeshIO::save(args.outputPath, inflator.vertices(), inflator.elements());
}