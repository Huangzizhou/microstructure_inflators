#include <bitset>
#include <string>

#include <CLI/CLI.hpp>
#include <inflators/wrappers/RBFInflator.hh>
#include <inflators/wrappers/RBFOrthoInflator.hh>

using namespace std;

struct Args {
    std::string frame1;
    std::string frame2;
    std::string output;
    std::string moptsPath;
    size_t dim;
    bool png = false;
    bool ortho = false;
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
    args.dim = 10;

    // Parse arguments
    CLI::App app{"rbf_interpolation"};

    app.add_option("frame1",  args.frame1,  "input1 path")->required()->check(CLI::ExistingFile);
    app.add_option("frame2",  args.frame2,  "input2 path")->required()->check(CLI::ExistingFile);
    app.add_option("output",  args.output,  "output path")->required();
    app.add_option("--dim",   args.dim,     "dim of rbf representation");
    app.add_flag(  "--ortho", args.ortho,   "orthotropic symmetry (uses less parameters)");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    if (args.ortho) {
        size_t d = args.dim;
        double epsilon = (d + d - 1.0) / 2.0;

        RBFOrthoInflator inflator1(args.frame1, epsilon, d);
        vector<double> coeffs1 = inflator1.defaultParameters();

        RBFOrthoInflator inflator2(args.frame2, epsilon, d);
        vector<double> coeffs2 = inflator2.defaultParameters();

        vector<double> final_coeffs(d*d, 0.0);;
        for (size_t i = 0; i < final_coeffs.size(); i++) {
            final_coeffs[i] = (coeffs1[i] + coeffs2[i]) / 2.0;
            std::cout << "Coeff 1: " << coeffs1[i] << std::endl;
            std::cout << "Coeff 2: " << coeffs2[i] << std::endl;
            std::cout << "Final Coeffs: " << final_coeffs[i] << std::endl;
        }

        inflator1.savePng(final_coeffs, args.output);
        inflator1.inflate(final_coeffs);
        //MeshIO::save(args.output, inflator1.vertices(), inflator1.elements());
    }
    else {
        size_t d = args.dim;
        double epsilon = d / 2.0;

        RBFInflator inflator1(args.frame1, epsilon, d);
        vector<double> coeffs1 = inflator1.defaultParameters();

        RBFInflator inflator2(args.frame2, epsilon, d);
        vector<double> coeffs2 = inflator2.defaultParameters();

        vector<double> final_coeffs(d*d, 0.0);;
        for (size_t i = 0; i < final_coeffs.size(); i++) {
            final_coeffs[i] = (coeffs1[i] + coeffs2[i]) / 2.0;
        }

        inflator1.savePng(final_coeffs, args.output);
        //inflator1.inflate(final_coeffs);
        //MeshIO::save(args.output, inflator1.vertices(), inflator1.elements());
    }


}