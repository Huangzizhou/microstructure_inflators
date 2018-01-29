//
// Created by Davi Colli Tozoni on 1/28/18.
//

#ifndef EXAMPLE_MESH_TOOLS_H
#define EXAMPLE_MESH_TOOLS_H

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

#endif //EXAMPLE_MESH_TOOLS_H
