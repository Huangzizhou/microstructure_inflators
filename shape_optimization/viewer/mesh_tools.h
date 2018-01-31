//
// Created by Davi Colli Tozoni on 1/28/18.
// Based on code from quadfoam project
//

#ifndef EXAMPLE_MESH_TOOLS_H
#define EXAMPLE_MESH_TOOLS_H

#include <boost/filesystem.hpp>

// Json
#include "json.hpp"
using json = nlohmann::json;

struct Mesh {
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi facets;
    Eigen::MatrixXd colors;

    std::map<std::string, Eigen::MatrixXd> vertexAttrs;
    std::map<std::string, Eigen::MatrixXd> facetAttrs;

    json params;

#ifndef NO_VIEWER
    // Viewer stuff
    Eigen::Affine3d _model_matrix;
    std::string _selected_va;
    std::string _selected_fa;
    int _selected_quad;
    bool _show_bbox;
    bool _show_lines;
    bool _show_fitted_quads;
    bool _show_jacobian;
    bool _show_orientation;
    bool _color_normalize;
    float _color_range[2];
    mutable float _attr_range[2];
    float _displacements;
    float _displacements_max;

    // Memorize viewer core state
    mutable std::vector<char> _viewer_core;
#endif

    Mesh(
        const Eigen::MatrixXd &V = Eigen::MatrixXd(0, 3),
        const Eigen::MatrixXi &F = Eigen::MatrixXi(0, 3)
    ) : vertices(V) , facets(F)
    {
#ifndef NO_VIEWER
        _model_matrix.setIdentity();
        _selected_quad = -1;
        _show_bbox = false;
        _show_lines = true;
        _show_fitted_quads = false;
        _show_jacobian = true;
        _show_orientation = false;
        _color_normalize = false;
        _color_range[0] = 0;
        _color_range[1] = 1;
        _attr_range[0] = 0;
        _attr_range[1] = 1;
        _displacements = 0;
        _displacements_max = 1;
#endif
    }

    // Clear mesh data (but doesn't reset the viewer variables)
    void clear() {
        vertices.resize(0, 3);
        facets.resize(0, 3);
        vertexAttrs.clear();
        facetAttrs.clear();
    }

    void normalize_mesh();

    bool is_triangle_mesh() const { return facets.cols() == 3; }
    bool is_quad_mesh() const { return facets.cols() == 4; }

#ifndef NO_VIEWER
    void set_max_displacement();

    void display(igl::viewer::Viewer &viewer, bool align_camera = false) const;

    void save_state(const igl::viewer::Viewer &viewer);
#endif
};

void Mesh::normalize_mesh() {
    const Eigen::RowVector3d minV = vertices.colwise().minCoeff();
    const Eigen::RowVector3d maxV = vertices.colwise().maxCoeff();
    double eps = std::numeric_limits<double>::epsilon();
    vertices = vertices.rowwise() - minV;
    vertices = vertices.array().rowwise() / (maxV - minV).cwiseMax(eps).array();
}

#ifndef NO_VIEWER
void Mesh::set_max_displacement() {
    // Set max displacement
    if (vertexAttrs.count("u")) {
        double max_disp = vertexAttrs.at("u").cwiseAbs().maxCoeff();
        const auto & minV = vertices.colwise().minCoeff().array();
        const auto & maxV = vertices.colwise().maxCoeff().array();
        double max_bbox = (maxV - minV).maxCoeff();
        _displacements_max = std::max((float) (max_bbox / max_disp), 1e-6f);
    }
}
#endif

void append_meshes(Mesh &mesh1, Mesh &mesh2, Mesh &scene) {
    // Create one huge mesh containing both meshes
    igl::cat(1, mesh1.vertices, mesh2.vertices, scene.vertices);
    igl::cat(1, mesh1.facets, Eigen::MatrixXi(mesh2.facets.array() + mesh1.vertices.rows()), scene.facets);
    igl::cat(1, mesh1.colors, mesh2.colors, scene.colors);

    //std::cout << "After merge... " << std::endl;
    //std::cout << "Vertices: " << scene.vertices << std::endl;
    //std::cout << "Faces: " << scene.facets << std::endl;
}

void generate_mesh_box(Eigen::Vector3d min, Eigen::Vector3d max, Eigen::MatrixXd &V_box, Eigen::MatrixXi &F_box) {
    // Corners of the bounding box
    V_box.resize(4,3);
    V_box <<
          min(0), min(1), 0.0,
          max(0), min(1), 0.0,
          max(0), max(1), 0.0,
          min(0), max(1), 0.0;

    // Faces of the bounding box
    F_box.resize(2,3);
    F_box <<
          0, 1, 2,
          0, 2, 3;
}

void generate_box(Eigen::Vector3d min, Eigen::Vector3d max, Eigen::MatrixXd &V_box, Eigen::MatrixXi &E_box) {
    // Corners of the bounding box
    V_box.resize(4,3);
    V_box <<
          min(0), min(1), 0.0,
          max(0), min(1), 0.0,
          max(0), max(1), 0.0,
          min(0), max(1), 0.0;

    // Edges of the bounding box
    E_box.resize(4,2);
    E_box <<
          0, 1,
          1, 2,
          2, 3,
          3, 0;
}

std::string convert_mesh_to_off(std::string input_file) {
    auto tmp_dir = boost::filesystem::temp_directory_path();
    auto f_off = tmp_dir / boost::filesystem::unique_path("so_%%%%-%%%%-%%%%-%%%%.off");

    std::string app(PYTHON_DIR "convert_mesh.py");
    std::string cmd = app + " " + input_file + " " + f_off.string();
    if (::system(cmd.c_str()) == 0) {
        return f_off.string();
    }

    return "";
}

// Split a string into tokens
std::vector<std::string> split(const std::string &str, const std::string &delimiters = " ") {
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    std::vector<std::string> tokens;
    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }

    return tokens;
}

inline bool startswith(const std::string &str, const std::string &prefix) {
    return (str.compare(0, prefix.size(), prefix) == 0);
}

template<int NUM_SIDES = 3>
bool load_obj(const std::string &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    std::string line;
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("failed to open file " + filename);
    }

    std::vector<Eigen::RowVector3d> VV;
    std::vector<Eigen::Matrix<int, 1, NUM_SIDES>> FF;
    while (std::getline(in, line)) {
        if (startswith(line, "# Vertices:")) {
            int n;
            std::sscanf(line.c_str(), "# Vertices: %d", &n);
            VV.reserve(n);

            std::cout << "# Vertices: " << n << std::endl;
            continue;
        }
        if (startswith(line, "# Faces:")) {
            int n;
            std::sscanf(line.c_str(), "# Faces: %d", &n);
            FF.reserve(n);

            std::cout << "# Faces: " << n << std::endl;
            continue;
        }
        std::istringstream iss(line);
        std::string key;
        if (iss >> key) {
            if (startswith(key, "#")) {
                continue;
            } else if (key == "v") {
                double x, y, z;
                iss >> x >> y >> z;
                VV.emplace_back();
                VV.back() << x, y, z;
            } else if (key == "f" || key == "l") {
                auto tokens = split(line.substr(1));
                if (tokens.size() != NUM_SIDES) {
                    std::cerr << "Facet has incorrect size: " << line << std::endl;
                    return false;
                } else {
                    FF.emplace_back();
                    for (int lv = 0; lv < NUM_SIDES; ++lv) {
                        std::string str = tokens[lv];
                        int v;
                        if (str.find('/') != std::string::npos) {
                            v = std::stoi(split(str, "/").front());
                        } else {
                            v = std::stoi(str);
                        }
                        FF.back()[lv] = v - 1; // Shift indices by 1 (start from 0)
                    }
                }
            }
        }
    }

    V.resize(VV.size(), 3);
    F.resize(FF.size(), NUM_SIDES);
    for (size_t v = 0; v < VV.size(); ++v) {
        V.row(v) = VV[v];
    }
    for (size_t f = 0; f < FF.size(); ++f) {
        F.row(f) = FF[f];
    }

    // std::cout << "Read a mesh with " << VV.size() << " vertices, " << FF.size() << " elements." << std::endl;

    return true;
}

#endif //EXAMPLE_MESH_TOOLS_H
