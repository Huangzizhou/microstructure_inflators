#include <igl/viewer/Viewer.h>
#include <igl/cat.h>
#include <nanogui/slider.h>
#include <nanogui/combobox.h>
#include "simulation_tools.h"
#include "mesh_tools.h"
#include <vector>
#include <string>

using namespace std;

struct BoundaryCondition {
    string config;
    string type;
    Mesh mesh;

    Eigen::Vector3d min_corner;
    Eigen::Vector3d max_corner;

    Eigen::Vector3d value;
};

Mesh mesh;
vector<BoundaryCondition> boundary_conditions;
igl::viewer::Viewer viewer;
vector<string> boundary_condition_configs;

bool endswith(std::string const &str, std::string const &ending) {
    if (str.length() >= ending.length()) {
        return (0 == str.compare(str.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

string relative_path(std::string str) {
    std::size_t pos = str.find_last_of("/");
    if (pos != string::npos)
        return str.substr(pos+1);
    else
        return str;
}

void read_boundary_condition_config(string json_path) {
    vector<BoundaryCondition> result;

    // reading JSON file
    string prefix = BOUNDARY_CONDITIONS_DIR;
    std::ifstream input(prefix + string("/") + json_path);
    json j;
    input >> j;

    json regions = j["regions"];
    // iterate the array
    for (json::iterator it = regions.begin(); it != regions.end(); ++it) {
        BoundaryCondition bc;
        bc.config = json_path;
        json current_region = *it;

        json box = current_region["box"];
        vector<float> min_corner = box["minCorner"];
        vector<float> max_corner = box["maxCorner"];

        bc.min_corner << min_corner[0], min_corner[1], min_corner[2];
        bc.max_corner << max_corner[0], max_corner[1], max_corner[2];

        bc.type = current_region["type"];
        Eigen::RowVector3d color;
        color << 0.0, 1.0, 0.0;
        if (bc.type.compare("dirichlet") == 0) {
            color << 0.0, 0.0, 0.0;
        }

        vector<float> value = current_region["value"];
        bc.value << value[0], value[1], value[2];

        generate_mesh_box(bc.min_corner, bc.max_corner, bc.mesh.vertices, bc.mesh.facets);

        bc.mesh.colors.resize(bc.mesh.facets.rows(), 3);
        bc.mesh.colors << color.replicate(bc.mesh.facets.rows(), 1);

        result.push_back(bc);
    }

    boundary_conditions = result;
}

void populate_boundary_condition_configs() {
    boost::filesystem::path dir = BOUNDARY_CONDITIONS_DIR;

    for (const auto & p : boost::filesystem::directory_iterator(dir)) {
        std::string s = p.path().filename().string();
        if (endswith(s, ".json")) {
            string filename = relative_path(p.path().string());
            boundary_condition_configs.push_back(filename);
            cout << filename << endl;
        }
    }

    read_boundary_condition_config(boundary_condition_configs[0]);
}

void display_mesh() {

    // Set displacements
    Mesh mesh_copy = mesh;
    if (mesh.vertexAttrs.count("u")) {
        mesh_copy.vertices = mesh.vertices + mesh._displacements * mesh.vertexAttrs.at("u");
    }

    // Show facet attribute
    mesh._attr_range[0] = 0.0f;
    mesh._attr_range[1] = 1.0f;
    if (mesh.facetAttrs.count(mesh._selected_fa)) {
        const Eigen::MatrixXd & A = mesh.facetAttrs.at(mesh._selected_fa);
        if (A.cols() == 9) {
            Eigen::VectorXd S;
            Eigen::MatrixXd C;
            von_mises_stress(A, S);
            double vmin = mesh._color_range[0];
            double vmax = mesh._color_range[1];
            if (mesh._color_normalize) {
                vmin = S.minCoeff();
                vmax = S.maxCoeff();
            }
            mesh._attr_range[0] = (float) S.minCoeff();
            mesh._attr_range[1] = (float) S.maxCoeff();

            //std::cout << "Mesh vertices: " << mesh_copy.vertices << std::endl;
            //std::cout << "Selected fa: " << mesh._selected_fa << std::endl;
            //std::cout << "Mesh attr range: " << mesh._attr_range[0] << ", " << mesh._attr_range[1] << std::endl;
            //std::cout << "Mesh color range: " << mesh._color_range[0]<< ", " << mesh._color_range[1] << std::endl;

            igl::jet(S, vmin, vmax, C);
            mesh.colors = C;
        } else if (A.cols() == 3) {
            mesh.colors = A;
        }
    }

    mesh_copy.colors = mesh.colors;
    Mesh scene_mesh_temp, scene_mesh;
    scene_mesh = boundary_conditions[0].mesh;
    for (unsigned i = 1; i < boundary_conditions.size(); i++) {
        append_meshes(scene_mesh, boundary_conditions[i].mesh, scene_mesh_temp);
        scene_mesh = scene_mesh_temp;
    }
    append_meshes(scene_mesh, mesh_copy, scene_mesh_temp);
    scene_mesh = scene_mesh_temp;
    viewer.data.set_mesh(scene_mesh.vertices, scene_mesh.facets);
    viewer.data.set_colors(scene_mesh.colors);
}

int main(int argc, char *argv[])
{
    populate_boundary_condition_configs();

    // Inline mesh of a cube
    Eigen::MatrixXd V(4,3);
    Eigen::MatrixXi F(2,3);

    V <<
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0;

    F <<
        0, 1, 2,
        0, 2, 3;

    //Mesh mesh(V, F);
    mesh.vertices = V;
    mesh.facets = F;
    mesh.colors.resize(F.rows(), 3);
    mesh.colors << Eigen::RowVector3d(1.0,0.7,0.2).replicate(F.rows(), 1);

    // Find the bounding box
    Eigen::Vector3d m = V.colwise().minCoeff();
    Eigen::Vector3d M = V.colwise().maxCoeff();
    Eigen::MatrixXd V_box;
    Eigen::MatrixXi E_box;
    generate_box(m, M, V_box, E_box);

    /*BoundaryCondition bc;
    Eigen::Vector3d min_dirichlet, max_dirichlet;
    min_dirichlet << -0.01, -0.01, 0;
    max_dirichlet << 1.01, 0.01, 0;
    Eigen::MatrixXd V_dirichlet;
    Eigen::MatrixXi E_dirichlet, F_dirichlet;
    generate_box(min_dirichlet, max_dirichlet, V_dirichlet, E_dirichlet);
    generate_mesh_box(min_dirichlet, max_dirichlet, V_dirichlet, F_dirichlet);
    bc.mesh.vertices = V_dirichlet;
    bc.mesh.facets = F_dirichlet;
    bc.mesh.colors.resize(F_dirichlet.rows(), 3);
    bc.mesh.colors << Eigen::RowVector3d(0.0, 1.0, 0.0).replicate(F_dirichlet.rows(), 1);
    boundary_conditions.push_back(bc);

    Eigen::Vector3d min_neumann, max_neumann;
    min_neumann << -0.01, 0.99, 0;
    max_neumann << 1.01, 1.01, 0;
    Eigen::MatrixXd V_neumann;
    Eigen::MatrixXi E_neumann, F_neumann;
    generate_box(min_neumann, max_neumann, V_neumann, E_neumann);
    generate_mesh_box(min_neumann, max_neumann, V_neumann, F_neumann);
    bc.mesh.vertices = V_neumann;
    bc.mesh.facets = F_neumann;
    bc.mesh.colors.resize(F_neumann.rows(), 3);
    bc.mesh.colors << Eigen::RowVector3d(0.0, 0.0, 1.0).replicate(F_neumann.rows(), 1);
    boundary_conditions.push_back(bc);*/

    // Plot the mesh
    // Extend viewer menu
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        // Add new group
        viewer.ngui->addGroup("Simulation tools:");

        nanogui::ComboBox *combo_box = new nanogui::ComboBox(viewer.ngui->window(), boundary_condition_configs);
        combo_box->setFontSize(16);
        combo_box->setFixedSize(Eigen::Vector2i(120,20));
        combo_box->setCallback([](int value) {
            cout << "Read: " << boundary_condition_configs[value] << endl;
            read_boundary_condition_config(boundary_condition_configs[value]);
            display_mesh();
        });
        viewer.ngui->addWidget("BCs", combo_box);

        viewer.ngui->addButton("Run simulation",[](){
            std::cout << "Simulating...\n";
            simulate(mesh, 200.0, 0.35, boundary_conditions[0].config);

            display_mesh();
        });

        // Expose variable directly ...
        viewer.ngui->addVariable("displacement", mesh._displacements);

        nanogui::Slider *slider = new nanogui::Slider(viewer.ngui->window());
        slider->setValue(0.0f);
        slider->setRange(std::pair<float, float>(-200.0f, 200.f));
        slider->setFixedWidth(200);
        /* Propagate slider changes to the text box */
        slider->setCallback([](float value) {
            mesh._displacements = value;
            display_mesh();
        });
        viewer.ngui->addWidget("", slider);

        // Generate menu
        viewer.screen->performLayout();

        return false;
    };

    // initializing camera and rotation mode
    using RT = igl::viewer::ViewerCore::RotationType;
    viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
    viewer.core.orthographic = true;
    viewer.core.set_rotation_type(RT::ROTATION_TYPE_NO_ROTATION);
    viewer.core.camera_zoom = 1.0;

    // plot min and max coordinates, to make it easier to understand the drawing
    std::stringstream l1;
    l1 << m(0) << ", " << m(1);
    viewer.data.add_label(m,l1.str());
    std::stringstream l2;
    l2 << M(0) << ", " << M(1);
    viewer.data.add_label(M,l2.str());

    viewer.core.background_color << 1.0, 1.0, 1.0, 1.0;
    viewer.core.shininess = 1e8;
    display_mesh();

    viewer.launch();
}
