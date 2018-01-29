//
// Created by Davi Colli Tozoni on 1/28/18.
//

#ifndef EXAMPLE_SIMULATION_TOOLS_H
#define EXAMPLE_SIMULATION_TOOLS_H

#include <boost/filesystem.hpp>
#include <igl/write_triangle_mesh.h>
#include "mesh_tools.h"

// Json
#include "json.hpp"
using json = nlohmann::json;

void read_packed_mesh(const std::string &filename, Mesh &M) {
    std::ifstream file(filename, std::ios::binary);
    std::vector<int32_t> size(4);
    file.read(reinterpret_cast<char*>(size.data()), 4 * sizeof(int32_t));

    int num_vert = size[0];
    int num_tris = size[1];
    int num_vattr = size[2];
    int num_fattr = size[3];
    M.clear();
    M.vertices.resize(num_vert, 3);
    M.facets.resize(num_tris, 3);

    // Check size, otherwise we need to use an tmp buffer
    assert(sizeof(int) == sizeof(int32_t));
    assert(sizeof(double) == sizeof(int64_t));
    file.read(reinterpret_cast<char*>(M.vertices.data()), M.vertices.size() * sizeof(double));
    file.read(reinterpret_cast<char*>(M.facets.data()), M.facets.size() * sizeof(int32_t));

    // Attributes
    assert(num_vattr == 3);
    for (auto va : {"u", "load", "Ku"}) {
        auto &V = M.vertexAttrs[va] = Eigen::MatrixXd(num_vert, 3);
        file.read(reinterpret_cast<char*>(V.data()), V.size() * sizeof(double));
    }
    assert(num_fattr == 2);
    for (auto fa : {"strain", "stress"}) {
        auto &F = M.facetAttrs[fa] = Eigen::MatrixXd(num_tris, 9);
        file.read(reinterpret_cast<char*>(F.data()), F.size() * sizeof(double));
    }
}

void simulate(Mesh &M, double young, double poisson) {
    using namespace boost;
    auto tmp_dir = filesystem::temp_directory_path();
    auto f_tri = tmp_dir / filesystem::unique_path("so_%%%%-%%%%-%%%%-%%%%.obj");
    auto f_msh = tmp_dir / filesystem::unique_path("so_%%%%-%%%%-%%%%-%%%%.msh");
    auto f_bin = tmp_dir / filesystem::unique_path("so_%%%%-%%%%-%%%%-%%%%.bin");
    auto f_mat = tmp_dir / filesystem::unique_path("so_%%%%-%%%%-%%%%-%%%%.json");
    std::string app(PYTHON_DIR "simulate_load.py");

    std::cout << "Young:" << young << std::endl;
    std::cout << "poisson:" << poisson << std::endl;

    std::string cmd = app + " " + f_tri.string() + " " + f_msh.string() + " -e " + f_bin.string();
    {
        json material = {
                {"type", "isotropic_material"},
                {"dim", 2},
                {"young", young},
                {"poisson", poisson},
        };
        std::ofstream out(f_mat.string());
        out << material << std::endl;
        cmd += " -m " + f_mat.string();
    }
    igl::write_triangle_mesh(f_tri.c_str(), M.vertices, M.facets);
    std::cout << "Running command:\n" + cmd << std::endl;
    if (::system(cmd.c_str()) == 0) {
        read_packed_mesh(f_bin.string(), M);
        //filesystem::remove(f_msh);
        filesystem::remove(f_bin);
    }
    filesystem::remove(f_tri);
    filesystem::remove(f_mat);

    // Display stress attribute
    M._selected_fa = "stress";

    // Set max displacement
    M.set_max_displacement();
}

void von_mises_stress(const Eigen::MatrixXd &S, Eigen::VectorXd &M) {
    assert(S.cols() == 9);
    M.resize(S.rows());
    for (int i = 0; i < S.rows(); ++i) {
        // stress: {σx, σy, σz, τxy, τyz, τzx}
        const Eigen::Matrix<double, 1, 9> buffer = S.row(i);
        const Eigen::Matrix3d stress(buffer.data());
        const Eigen::Vector3d sigma(stress(0,0), stress(1,1), stress(2,2));
        const Eigen::Vector3d tau(stress(0,1), stress(1,2), stress(2,0));
        const Eigen::Vector3d shifted(sigma(1), sigma(2), sigma(0));
        M(i) = std::sqrt(0.5 * (
                (sigma - shifted).squaredNorm() + 6.0 * tau.squaredNorm()
        ));
    }
}

#endif //EXAMPLE_SIMULATION_TOOLS_H
