#include <igl/viewer/Viewer.h>
#include <igl/copyleft/cgal/convex_hull.h>
#include <igl/cat.h>
#include <nanogui/slider.h>
#include <nanogui/combobox.h>
#include "simulation_tools.h"
#include "mesh_tools.h"
#include "file_dialog.h"
#include <vector>
#include <string>

using namespace std;

enum EditorState{
    ADDING_BOUNDARY_CONDITION, NONE
};
typedef EditorState EditorState;

EditorState current_state = NONE; // current state
Eigen::MatrixX3d current_points;

nanogui::Window * pop_up;

struct BoundaryCondition {
    string config;
    string type;
    Mesh mesh;

    Eigen::Vector3d min_corner;
    Eigen::Vector3d max_corner;

    Eigen::MatrixX3d area_points;

    Eigen::Vector3d value;
};

Mesh mesh;
vector<BoundaryCondition> boundary_conditions;
igl::viewer::Viewer viewer;
vector<string> boundary_condition_configs;

void display_mesh();
void clear_mesh();

// This function is called every time a keyboard button is pressed
bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (button == 0) {
        Eigen::Vector3f bc;
        double x = viewer.current_mouse_x;
        double y = viewer.core.viewport(3) - viewer.current_mouse_y;

        double left_x = viewer.core.viewport(0);
        double top_y = viewer.core.viewport(1);
        double width = viewer.core.viewport(2);
        double height = viewer.core.viewport(3);
        double canonical_x = 2/width * (x - left_x) - 1.;
        double canonical_y = 2/height * (y - top_y) - 1.;

        Eigen::Vector4f viewport_coord;
        viewport_coord <<  canonical_x, canonical_y, 0.0, 1.0;
        Eigen::Vector4f world_coord = (viewer.core.proj * viewer.core.view * viewer.core.model).inverse() * viewport_coord;

        cout << "(x,y): " << world_coord(0) << ", " << world_coord(1) << endl;

        Eigen::RowVector3d new_point(world_coord(0), world_coord(1), 0.0);

        if (current_state == ADDING_BOUNDARY_CONDITION) {
            current_points.conservativeResize(current_points.rows() + 1, 3);
            current_points.row(current_points.rows() - 1) = new_point;
            viewer.data.add_points(new_point, Eigen::RowVector3d(1, 0, 0));
        }
    }
    else {
        viewer.data.points.resize(0,0);
        current_points.resize(0,0);
        viewer.data.dirty |= igl::viewer::ViewerData::DIRTY_OVERLAY_POINTS;
    }

    return false;
}

// This function is called every time a keyboard button is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
    if (key == 1)
    {
        cout << "Enter!" << endl;
        if (current_points.rows() == 2) {
            // TODO: box!
        }
        else if (current_points.rows() > 2) {
            BoundaryCondition bc;

            cout << "Drawing polygon..." << endl;
            // draw polygon
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            cout << "Current points " << current_points << endl;

            igl::copyleft::cgal::convex_hull(current_points, V, F);

            int origin = 0;
            //Eigen::RowVector3d origin = V.row(0);
            F.resize(V.rows()-2, 3);
            for (unsigned v = 2; v < V.rows(); v++) {
                F.row(v-2) << origin, v-1, v;
            }

            bc.area_points = current_points;
            bc.config = "extra";
            bc.mesh.vertices = V;
            bc.mesh.facets = F;

            cout << "V " << V << endl;
            cout << "F " << F << endl;

            Eigen::RowVector3d color;
            color << 1.0, 0.0, 0.0;
            bc.mesh.colors.resize(bc.mesh.facets.rows(), 3);
            bc.mesh.colors << color.replicate(bc.mesh.facets.rows(), 1);
            boundary_conditions.push_back(bc);

            current_points.resize(0,0);
            viewer.data.points.resize(0,3);
            clear_mesh();
            display_mesh();

            // Add an additional menu window
            nanogui::Window * new_window = viewer.ngui->addWindow(Eigen::Vector2i(220,10),"New Boundary Condition");
            pop_up = new_window;
            // Expose the same variable directly ...
            viewer.ngui->addButton("Ok", [&](){
                cout << "Clicked Ok" << endl;

                viewer.screen->removeChild(pop_up);

                //pop_up->setVisible(false);
            });

            // Generate menu
            viewer.screen->performLayout();
        }
    }

    return false;
}

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

    cout << "BC: " << json_path << endl;

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
        json box_percentage = current_region["box%"];

        if (box.size() > 0) {
            cout << "boxes: " << box.size() << endl;

            vector<float> min_corner = box["minCorner"];
            vector<float> max_corner = box["maxCorner"];

            bc.min_corner << min_corner[0], min_corner[1], min_corner[2];
            bc.max_corner << max_corner[0], max_corner[1], max_corner[2];
        }
        if (box_percentage.size() > 0){
            cout << "box%es: " << box_percentage.size() << endl;

            // Find the bounding box
            Eigen::Vector3d m;
            Eigen::Vector3d M;
            if (mesh.vertices.rows() > 0) {
                m = mesh.vertices.colwise().minCoeff();
                M = mesh.vertices.colwise().maxCoeff();
            }
            else {
                m << -1.0 , -1.0, 0.0;
                M << 1.0 , 1.0, 0.0;
            }

            vector<float> min_corner = box_percentage["minCorner"];
            vector<float> max_corner = box_percentage["maxCorner"];

            Eigen::Array3d min_corner_array;
            Eigen::Array3d max_corner_array;
            min_corner_array << min_corner[0], min_corner[1], min_corner[2];
            max_corner_array << max_corner[0], max_corner[1], max_corner[2];

            // translate to 0,0,0
            //min_corner_array += 1.0;
            //max_corner_array += 1.0;

            // scale (w/2, h/2)
            double w = M(0) - m(0);
            double h = M(1) - m(1);
            double d = M(2) - m(2);
            min_corner_array(0) *= w;///2;
            min_corner_array(1) *= h;///2;
            min_corner_array(2) *= d;///2;
            max_corner_array(0) *= w;///2;
            max_corner_array(1) *= h;///2;
            max_corner_array(2) *= d;///2;

            // translate to min corner
            min_corner_array += m.array();
            max_corner_array += m.array();

            bc.min_corner = min_corner_array.matrix();
            bc.max_corner = max_corner_array.matrix();
        }

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

void clear_mesh() {
    viewer.data.V = Eigen::MatrixXd (0,3);
    viewer.data.F = Eigen::MatrixXi (0,3);
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
    Mesh scene_mesh;
    scene_mesh = boundary_conditions[0].mesh;
    for (unsigned i = 1; i < boundary_conditions.size(); i++) {
        Mesh mesh_temp;
        append_meshes(scene_mesh, boundary_conditions[i].mesh, mesh_temp);
        scene_mesh = mesh_temp;
    }

    Mesh mesh_temp;
    append_meshes(scene_mesh, mesh_copy, mesh_temp);
    scene_mesh = mesh_temp;

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

    // Plot the mesh
    // Extend viewer menu
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        // Add new group
        viewer.ngui->addGroup("Simulation tools:");

        viewer.ngui->addButton("Load mesh", []() {
            std::cout << "Loading mesh...\n";

            string filename = FileDialog::openFileName();

            cout << "Converting msh file " << filename << endl;

            string off_path = convert_mesh_to_off(filename);

            cout << "Opening obj file " << off_path << endl;

            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            //load_obj(obj_path, V, F);
            igl::readOFF(off_path, V, F);

            if (V.cols() == 3 && F.cols() == 3) {
                clear_mesh();

                mesh.vertices = V;
                mesh.facets = F;
                Eigen::RowVector3d color;
                color << 1.0, 1.0, 0.0;
                mesh.colors.resize(F.rows(), 3);
                mesh.colors << color.replicate(F.rows(), 1);

                mesh.vertexAttrs.clear();
                mesh.facetAttrs.clear();

                cout << "Vertices: " << mesh.vertices.rows() << endl;
                cout << "Faces: " << mesh.facets.rows() << endl;
                cout << "Colors: " << mesh.colors.rows() << endl;
            }

            //viewer.load_mesh_from_file(obj_path);
            if (boundary_conditions.size() > 0) {
                read_boundary_condition_config(boundary_conditions[0].config);
            }

            display_mesh();
        });

        nanogui::ComboBox *combo_box = new nanogui::ComboBox(viewer.ngui->window(), boundary_condition_configs);
        combo_box->setFontSize(16);
        combo_box->setFixedSize(Eigen::Vector2i(120,20));
        combo_box->setCallback([](int value) {
            cout << "Read: " << boundary_condition_configs[value] << endl;
            read_boundary_condition_config(boundary_condition_configs[value]);
            display_mesh();
        });
        viewer.ngui->addWidget("BCs", combo_box);

        viewer.ngui->addButton("Add BC",[](){
            std::cout << "Adding boundary condition...\n";
            current_state = ADDING_BOUNDARY_CONDITION;
        });

        viewer.ngui->addButton("Run simulation",[](){
            std::cout << "Simulating...\n";
            simulate(mesh, 200.0, 0.35, boundary_conditions[0].config);

            display_mesh();
        });

        // Expose variable directly ...
        viewer.ngui->addVariable("displacement", mesh._displacements);

        nanogui::Slider *slider = new nanogui::Slider(viewer.ngui->window());
        slider->setValue(0.0f);
        slider->setRange(std::pair<float, float>(-20.0f, 20.f));
        slider->setFixedWidth(200);
        /* Propagate slider changes to the text box */
        slider->setCallback([](float value) {
            cout << "Sliding..." << endl;
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

    if (boundary_conditions.size() > 0) {
        cout << "Reading conditions " << boundary_conditions[0].config << endl;
        read_boundary_condition_config(boundary_conditions[0].config);
    }
    display_mesh();

    viewer.callback_mouse_down = &mouse_down;
    viewer.callback_key_down = &key_down;

    viewer.launch();
}
