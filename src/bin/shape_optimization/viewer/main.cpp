#include <igl/viewer/Viewer.h>
#include <igl/copyleft/cgal/convex_hull.h>
#include <igl/cat.h>
#include <igl/unproject_onto_mesh.h>
#include <nanogui/slider.h>
#include <nanogui/combobox.h>
#include <vector>
#include <string>

#include "simulation_tools.h"
#include "mesh_tools.h"
#include "file_dialog.h"
#include "boundary_condition_tools.h"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
using namespace std;

enum EditorState{
    ADDING_BOUNDARY_CONDITION, EDITING_NEW_BOUNDARY_CONDITION, NONE
};
typedef EditorState EditorState;

EditorState current_state = NONE; // current state
Eigen::MatrixX3d current_points;
nanogui::Window * pop_up;
nanogui::Window * selected_pop_up;
Mesh mesh;
bool is_boundary_condition_selected = false;
vector<BoundaryCondition> boundary_conditions;
BoundaryCondition * selected_boundary_condition;
int selected_boundary_condition_index = -1;
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

            if (is_boundary_condition_selected) {
                is_boundary_condition_selected = false;
                viewer.screen->removeChild(selected_pop_up);
                display_mesh();
            }

        }
        else if (current_state == EDITING_NEW_BOUNDARY_CONDITION) {

        }
        else {

            if (is_boundary_condition_selected) {
                cout << "removing previously selected bc" << endl;

                Eigen::RowVector3d color;
                color << 0.0, 1.0, 0.0;
                if (selected_boundary_condition->type.compare("dirichlet") == 0) {
                    color << 0.0, 0.0, 0.0;
                }

                selected_boundary_condition->mesh.colors.resize(selected_boundary_condition->mesh.facets.rows(), 3);
                selected_boundary_condition->mesh.colors << color.replicate(selected_boundary_condition->mesh.facets.rows(), 1);

                is_boundary_condition_selected = false;
                viewer.screen->removeChild(selected_pop_up);
                display_mesh();
            }

            for (unsigned b=0; b<boundary_conditions.size(); b++) {
                BoundaryCondition &bc = boundary_conditions[b];
                int fid;
                Eigen::RowVector3d intersection;

                if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
                                             viewer.core.proj, viewer.core.viewport, bc.mesh.vertices, bc.mesh.facets, fid, intersection)) {
                    cout << "hit" << endl;

                    // paint hit red
                    Eigen::RowVector3d color;
                    color << 1.0, 0.0, 0.0;
                    bc.mesh.colors << color.replicate(bc.mesh.facets.rows(), 1);

                    cout << "before displaying mesh" << endl;
                    display_mesh();
                    cout << "after displaying mesh" << endl;

                    is_boundary_condition_selected = true;
                    selected_boundary_condition = &(boundary_conditions[b]);
                    selected_boundary_condition_index = b;

                    cout << "finish hit" << endl;

                    break;
                }
                else {
                    cout << "NOT a hit" << endl;
                }
            }

            if (is_boundary_condition_selected) {
                cout << "drawing new menu" << endl;

                // Add an additional menu window
                nanogui::Window * new_window = viewer.ngui->addWindow(Eigen::Vector2i(470,10),"Boundary Condition");
                selected_pop_up = new_window;

                cout << "drawing new menu 2" << endl;

                nanogui::ComboBox *combo_box = new nanogui::ComboBox(new_window, {"dirichlet", "force"});
                combo_box->setCallback((const function<void(int)> &) [&](int value) {
                    vector<string> options = {"dirichlet", "force"};
                    cout << "type: " << options[value] << endl;
                    selected_boundary_condition->type = options[value];
                });
                if (selected_boundary_condition->type.compare("dirichlet") == 0)
                    combo_box->setSelectedIndex(0);
                else
                    combo_box->setSelectedIndex(1);

                viewer.ngui->addWidget("type", combo_box);

                cout << "drawing new menu 3" << endl;

                nanogui::FloatBox<float> *float_box_x = new nanogui::FloatBox<float>(new_window, selected_boundary_condition->value(0));
                float_box_x->setEditable(true);
                float_box_x->setCallback((const function<void(float)> &) [&](float value) {
                    selected_boundary_condition->value(0) = value;
                });
                nanogui::FloatBox<float> *float_box_y = new nanogui::FloatBox<float>(new_window, selected_boundary_condition->value(1));
                float_box_y->setEditable(true);
                float_box_y->setCallback((const function<void(float)> &) [&](float value) {
                    selected_boundary_condition->value(1) = value;
                });
                nanogui::FloatBox<float> *float_box_z = new nanogui::FloatBox<float>(new_window, selected_boundary_condition->value(2));
                float_box_z->setEditable(true);
                float_box_z->setCallback((const function<void(float)> &) [&](float value) {
                    selected_boundary_condition->value(2) = value;
                });
                viewer.ngui->addWidget("x", float_box_x);
                viewer.ngui->addWidget("y", float_box_y);
                viewer.ngui->addWidget("z", float_box_z);

                // Expose delete option
                viewer.ngui->addButton("Delete", [&](){
                    is_boundary_condition_selected = false;
                    viewer.screen->removeChild(selected_pop_up);

                    cout << "Deleting " << selected_boundary_condition_index << endl;
                    boundary_conditions.erase(boundary_conditions.begin() + selected_boundary_condition_index);

                    clear_mesh();
                    display_mesh();
                });

                viewer.ngui->addButton("Ok", [&](){
                    Eigen::RowVector3d color;
                    color << 0.0, 1.0, 0.0;
                    if (selected_boundary_condition->type.compare("dirichlet") == 0) {
                        color << 0.0, 0.0, 0.0;
                    }

                    selected_boundary_condition->mesh.colors.resize(selected_boundary_condition->mesh.facets.rows(), 3);
                    selected_boundary_condition->mesh.colors << color.replicate(selected_boundary_condition->mesh.facets.rows(), 1);

                    is_boundary_condition_selected = false;
                    viewer.screen->removeChild(selected_pop_up);
                    display_mesh();

                });

                cout << "drawing new menu 3" << endl;

                // Generate menu
                viewer.screen->performLayout();
                display_mesh();
            }
        }
    }
    else {
        viewer.data.points.resize(0,0);
        current_points.resize(0,3);
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
        BoundaryCondition bc;
        if (current_points.rows() == 2) {
            cout << "Drawing box..." << endl;
            // draw polygon
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            Eigen::RowVector3d m, M;
            if (current_points(0,0) < current_points(1,0)) {
                m(0) = current_points(0,0);
                M(0) = current_points(1,0);
            }
            else {
                m(0) = current_points(1,0);
                M(0) = current_points(0,0);
            }

            if (current_points(0,1) < current_points(1,1)) {
                m(1) = current_points(0,1);
                M(1) = current_points(1,1);
            }
            else {
                m(1) = current_points(1,1);
                M(1) = current_points(0,1);
            }

            if (current_points(0,2) < current_points(1,2)) {
                m(2) = current_points(0,2);
                M(2) = current_points(1,2);
            }
            else {
                m(2) = current_points(1,2);
                M(2) = current_points(0,2);
            }

            bc.min_corner = m;
            bc.max_corner = M;

            generate_mesh_box(bc.min_corner, bc.max_corner, bc.mesh.vertices, bc.mesh.facets);

            cout << "Generated mesh box: " << endl;
            cout << bc.mesh.vertices << endl;
            cout << bc.mesh.facets << endl;

            bc.area_points.resize(2, 3);
            bc.area_points.row(0) = m;
            bc.area_points.row(1) = M;
            bc.config = "extra";
        }
        else if (current_points.rows() > 2) {
            cout << "Drawing polygon..." << endl;
            // draw polygon
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            cout << "Current points " << current_points << endl;

            igl::copyleft::cgal::convex_hull(current_points, V, F);

            int origin = 0;
            //Eigen::RowVector3d origin = V.row(0);
            F.resize(V.rows() - 2, 3);
            for (unsigned v = 2; v < V.rows(); v++) {
                F.row(v - 2) << origin, v - 1, v;
            }

            bc.area_points = current_points;
            bc.config = "extra";
            bc.mesh.vertices = V;
            bc.mesh.facets = F;
        }

        Eigen::RowVector3d color;
        color << 1.0, 0.0, 0.0;
        bc.mesh.colors.resize(bc.mesh.facets.rows(), 3);
        bc.mesh.colors << color.replicate(bc.mesh.facets.rows(), 1);

        cout << "Generated mesh box: " << endl;
        cout << bc.mesh.vertices << endl;
        cout << bc.mesh.facets << endl;

        boundary_conditions.push_back(bc);

        current_points.resize(0,3);
        viewer.data.points.resize(0,3);
        clear_mesh();
        display_mesh();

        // Add an additional menu window
        nanogui::Window * new_window = viewer.ngui->addWindow(Eigen::Vector2i(470,10),"New Boundary Condition");
        pop_up = new_window;

        boundary_conditions.back().type = "dirichlet";
        nanogui::ComboBox *combo_box = new nanogui::ComboBox(new_window, {"dirichlet", "force"});
        combo_box->setCallback((const function<void(int)> &) [](int value) {
            vector<string> options = {"dirichlet", "force"};
            cout << "option: " << options[value] << endl;
            boundary_conditions.back().type = options[value];
        });
        viewer.ngui->addWidget("type", combo_box);

        boundary_conditions.back().value << 0.0, 0.0, 0.0;
        nanogui::FloatBox<float> *float_box_x = new nanogui::FloatBox<float>(new_window, 0.0);
        float_box_x->setEditable(true);
        float_box_x->setCallback((const function<void(float)> &) [&](float value) {
            boundary_conditions.back().value(0) = value;
        });
        nanogui::FloatBox<float> *float_box_y = new nanogui::FloatBox<float>(new_window, 0.0);
        float_box_y->setEditable(true);
        float_box_y->setCallback((const function<void(float)> &) [&](float value) {
            boundary_conditions.back().value(1) = value;
        });
        nanogui::FloatBox<float> *float_box_z = new nanogui::FloatBox<float>(new_window, 0.0);
        float_box_z->setEditable(true);
        float_box_z->setCallback((const function<void(float)> &) [&](float value) {
            boundary_conditions.back().value(2) = value;
        });
        viewer.ngui->addWidget("x", float_box_x);
        viewer.ngui->addWidget("y", float_box_y);
        viewer.ngui->addWidget("z", float_box_z);


        // Expose the same variable directly ...
        viewer.ngui->addButton("Ok", [&](){
            cout << "Clicked Ok" << endl;

            cout << "Value: " << boundary_conditions.back().value << endl;

            Eigen::RowVector3d color;
            color << 0.0, 1.0, 0.0;
            if (boundary_conditions.back().type.compare("dirichlet") == 0) {
                color << 0.0, 0.0, 0.0;
            }

            boundary_conditions.back().mesh.colors.resize(boundary_conditions.back().mesh.facets.rows(), 3);
            boundary_conditions.back().mesh.colors << color.replicate(boundary_conditions.back().mesh.facets.rows(), 1);

            viewer.screen->removeChild(pop_up);
            display_mesh();
        });

        // Generate menu
        viewer.screen->performLayout();
        display_mesh();

        current_state = NONE;

    }

    return false;
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

    boundary_conditions = read_boundary_condition_config(boundary_condition_configs[0], mesh);
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
    if (boundary_conditions.size() > 0)
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

void load_mesh() {
    cout << "Loading mesh...\n";

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

        // Add new window
        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Simulation tools");

        viewer.ngui->addButton("Load mesh", (const function<void()> &) [](){
            load_mesh();

            //viewer.load_mesh_from_file(obj_path);
            if (boundary_conditions.size() > 0) {
                boundary_conditions = read_boundary_condition_config(boundary_conditions[0].config, mesh);
            }

            display_mesh();
        });

        nanogui::ComboBox *combo_box = new nanogui::ComboBox(viewer.ngui->window(), boundary_condition_configs);
        //combo_box->setFontSize(16);
        combo_box->setFixedSize(Eigen::Vector2i(100,20));
        combo_box->setCallback((const function<void(int)> &) [](int value) {
                    cout << "Read: " << boundary_condition_configs[value] << endl;
                    boundary_conditions = read_boundary_condition_config(boundary_condition_configs[value], mesh);
                    clear_mesh();
                    display_mesh();
                });
        viewer.ngui->addWidget("BCs", combo_box);

        viewer.ngui->addButton("Save BC condiguration", (const function<void()> &) [](){
            std::cout << "Saving configuration...\n";
            json bcs_json = boundary_conditions_to_json(boundary_conditions);
            save_json_to_config_file(bcs_json);
        });

        viewer.ngui->addButton("Add BC", (const function<void()> &) [](){
                    std::cout << "Adding boundary condition...\n";
                    current_state = ADDING_BOUNDARY_CONDITION;
                });

        viewer.ngui->addButton("Run simulation", (const function<void()> &) [](){
                    std::cout << "Simulating...\n";
                    json bcs_json = boundary_conditions_to_json(boundary_conditions);
                    simulate(mesh, 200.0, 0.35, bcs_json);

                    display_mesh();
                });

        // Expose variable directly ...
        viewer.ngui->addVariable("displacement", mesh._displacements);

        nanogui::Slider *slider = new nanogui::Slider(viewer.ngui->window());
        slider->setValue(0.0f);
        slider->setRange(std::pair<float, float>(-20.0f, 20.f));
        slider->setFixedWidth(200);
        /* Propagate slider changes to the text box */
        slider->setCallback((const function<void(float)> &) [](float value) {
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
    viewer.core.point_size = 15;

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
        boundary_conditions = read_boundary_condition_config(boundary_conditions[0].config, mesh);
    }
    display_mesh();

    viewer.callback_mouse_down = &mouse_down;
    viewer.callback_key_down = &key_down;


    viewer.launch();
}

#pragma clang diagnostic pop