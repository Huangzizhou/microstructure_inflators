#include <CLI/CLI.hpp>
#define STB_IMAGE_IMPLEMENTATION
#include "read_png.hh"
#include <bitset>
#include <inflators/wrappers/VoxelsInflator.hh>

struct Args {
    std::string inputPath;
    std::string outputPath;
    std::string moptsPath;
};

std::vector<Real>
matToVec(const std::vector<std::vector<Real>> &mat) {
    std::vector<Real> result;

    for (size_t i = 0; i < mat.size(); i++) {
        result.insert(result.end(), mat[i].begin(), mat[i].end() );
    }

    return result;
}

int main(int argc, char ** argv) {
    Args args;

    // Parse arguments
    CLI::App app{"voxel2mesh"};

    app.add_option("inputPath",  args.inputPath,  "input  path")->required()->check(CLI::ExistingFile);
    app.add_option("outputPath", args.outputPath, "output path")->required();
    app.add_option("--mopts", args.moptsPath, "mesh options");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    // read input png file, transforming it into a png_bytep pointer
    int width, height;
    int channels;
    unsigned char * png_matrix = read_png_file(args.inputPath, width, height, channels);

    // Construct density matrix
    std::vector<std::vector<double>> density_matrix(height);

    if (channels == 1) {
        for (size_t y = 0; y < height; y++) {
            unsigned char *row = &png_matrix[y * width];
            density_matrix[y] = std::vector<double>(width);
            //std::cout << "width: " << width << std::endl;

            for (size_t x = 0; x < width; x++) {
                unsigned char byte = row[x];
                density_matrix[y][x] = 1.0 - byte / 255.0;
                // std::cout << density_matrix[y][x];
            }
            // std::cout << std::endl;
        }
    }

    // Print densities
    /*for (size_t i = 0; i < density_matrix.size(); i++) {
        for (size_t j = 0; j < density_matrix[i].size(); j++) {
            std::cout << density_matrix[i][j];
        }
        std::cout << std::endl;
    }*/

    VoxelsInflator inflator(density_matrix);
    if (!args.moptsPath.empty())  inflator.meshingOptions().load(args.moptsPath);

    std::vector<double> params = matToVec(density_matrix);
    inflator.inflate(params);

    MeshIO::save(args.outputPath, inflator.vertices(), inflator.elements());
}
