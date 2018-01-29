#include <igl/viewer/Viewer.h>
#include <igl/cat.h>
#include <nanogui/slider.h>
#include "simulation_tools.h"
#include "mesh_tools.h"

Mesh mesh, mesh_bc1, mesh_bc2;
igl::viewer::Viewer viewer;

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
    append_meshes(mesh_bc1, mesh_bc2, scene_mesh_temp);
    append_meshes(scene_mesh_temp, mesh_copy, scene_mesh);
    viewer.data.set_mesh(scene_mesh.vertices, scene_mesh.facets);
    viewer.data.set_colors(scene_mesh.colors);
}

int main(int argc, char *argv[])
{
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

    Eigen::Vector3d min_dirichlet, max_dirichlet;
    min_dirichlet << -0.01, -0.01, 0;
    max_dirichlet << 1.01, 0.01, 0;
    Eigen::MatrixXd V_dirichlet;
    Eigen::MatrixXi E_dirichlet, F_dirichlet;
    generate_box(min_dirichlet, max_dirichlet, V_dirichlet, E_dirichlet);
    generate_mesh_box(min_dirichlet, max_dirichlet, V_dirichlet, F_dirichlet);
    mesh_bc1.vertices = V_dirichlet;
    mesh_bc1.facets = F_dirichlet;
    mesh_bc1.colors.resize(F_dirichlet.rows(), 3);
    mesh_bc1.colors << Eigen::RowVector3d(0.0, 1.0, 0.0).replicate(F_dirichlet.rows(), 1);

    Eigen::Vector3d min_neumann, max_neumann;
    min_neumann << -0.01, 0.99, 0;
    max_neumann << 1.01, 1.01, 0;
    Eigen::MatrixXd V_neumann;
    Eigen::MatrixXi E_neumann, F_neumann;
    generate_box(min_neumann, max_neumann, V_neumann, E_neumann);
    generate_mesh_box(min_neumann, max_neumann, V_neumann, F_neumann);
    mesh_bc2.vertices = V_neumann;
    mesh_bc2.facets = F_neumann;
    mesh_bc2.colors.resize(F_neumann.rows(), 3);
    mesh_bc2.colors << Eigen::RowVector3d(0.0, 0.0, 1.0).replicate(F_neumann.rows(), 1);

    // Plot the mesh
    // Extend viewer menu
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        // Add new group
        viewer.ngui->addGroup("Simulation tools:");
        viewer.ngui->addButton("Run simulation",[](){
            std::cout << "Simulating...\n";
            simulate(mesh, 200.0, 0.35);

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
