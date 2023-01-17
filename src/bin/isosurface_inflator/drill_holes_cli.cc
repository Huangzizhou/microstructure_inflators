////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>
#include <isosurface_inflator/VDBTools.hh>

#include <openvdb/tools/Diagnostics.h>
#include <openvdb/tools/LevelSetSphere.h>

#include <igl/point_mesh_squared_distance.h>

#include <fstream>
////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

FloatGrid::Ptr erode(const std::vector<Vec3s> &V, const std::vector<Vec3I> &F, const double volume_ratio, double tol, std::vector<Vec3s> &Ve, std::vector<Vec3I> &Tri) 
{
    // std::vector<Vec3s> V;
    // std::vector<Vec3I> F;
    // read_mesh(input, V, F);

    double half_diag;
    {
        Eigen::MatrixXd V_mat;
        openvdb_to_eigen(V, V_mat);
        half_diag = igl::bounding_box_diagonal(V_mat) / 2;
    }

    const double initial_vol = abs(compute_surface_mesh_volume(V, F));
    const double hole_vol = 1 - volume_ratio;

    // std::cout << "target volume ratio: " << hole_vol << "\n";
    // std::cout << "half diag: " << half_diag << "\n";

    // eroded mesh
    // std::vector<Vec3s> Ve;
    // std::vector<Vec3I> Tri;
    double voxel = 1e-1 * half_diag;
    int halfWidth = 5;// std::max(5, (int)(min_width / voxel));
    double current;
    double current_vol;
    FloatGrid::Ptr grid;
    while (voxel > 1e-3 * half_diag)
    {
        math::Transform::Ptr xform(nullptr);
        xform = math::Transform::createLinearTransform(voxel);
        grid = tools::meshToLevelSet<FloatGrid>(*xform, V, F, halfWidth);

        double upper = 0.2 * half_diag;
        double upper_vol = volume_after_offset(grid, upper, Ve, Tri) / initial_vol;
        double lower = 0;
        double lower_vol = volume_after_offset(grid, lower, Ve, Tri) / initial_vol;

        while (upper_vol > hole_vol && upper < half_diag * 2)
        {
            upper *= 1.5;
            upper_vol = volume_after_offset(grid, upper, Ve, Tri) / initial_vol;

            // cout << "current radius interval: " << upper << "\n";
            // cout << "current volume interval: " << upper_vol << "\n";
        }

        assert(upper_vol <= hole_vol && lower_vol >= hole_vol);

        current = (lower + upper) / 2;
        current_vol = volume_after_offset(grid, current, Ve, Tri) / initial_vol;
        while (abs(current_vol - hole_vol) > tol && abs(upper - lower) > 1e-6 * half_diag)
        {
            if (current_vol > hole_vol)
            {
                lower = current;
                lower_vol = current_vol;
            }
            else
            {
                upper = current;
                upper_vol = current_vol;
            }

            current = (lower + upper) / 2;
            current_vol = volume_after_offset(grid, current, Ve, Tri) / initial_vol;

            // cout << "current radius interval: " << lower << " " << current << " " << upper << "\n";
            // cout << "current volume interval: " << lower_vol << " " << current_vol << " " << upper_vol << "\n";

            if (current_vol > lower_vol || current_vol < upper_vol)
                break;
        }

        if (abs(current_vol - hole_vol) < tol)
            break;
        else
        {
            voxel /= 2;
            halfWidth *= 1.5;
            std::cout << current_vol << " doesn't satisfy volume requirement, refine voxel to " << voxel / half_diag << "\n";
        }
    }

    if (abs(current_vol - hole_vol) > tol)
        throw std::runtime_error("Doesn't satisfy volume requirement");

    std::unique_ptr<tools::LevelSetFilter<FloatGrid>> filter = std::make_unique<tools::LevelSetFilter<FloatGrid>>(*gridPtrCast<FloatGrid>(grid));
    filter->offset(current);
    
    return grid;
}

/*
patch json format
[
    {
        "index": [2,2,3],
        "ratio": 0.5
    },
    ...
]
*/

int main(int argc, char * argv[]) {
    // Default arguments
    struct {
        std::string volume, surface;
        std::string patch_config;
        std::string output = "out.obj";
        double gridSize = 0.1;
        int resolution = 50;
        double tunnel = 0;
        bool preserve_original_surface = false;
    } args;

    // Parse arguments
    CLI::App app{"drill_holes_cli"};
    app.add_option("volume,--vol", args.volume, "Volume mesh.")->required();
    app.add_option("patch,-p,--patch", args.patch_config, "Patch description (json file).")->required();
    app.add_option("--gridSize", args.gridSize, "Grid size.")->required();
    app.add_option("-t,--tunnel", args.tunnel, "Tunnel size.")->required();
    app.add_option("--surface", args.surface, "Surface mesh.");
    app.add_option("-o,--output", args.output, "Output triangle mesh.");
    app.add_option("-r,--resolution", args.resolution, "Density field resolution.");
    app.add_flag("--preserve-surface", args.preserve_original_surface, "Preserve original object surface");
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    if (args.tunnel > 0)
        args.preserve_original_surface = false;

    // Load patch config
    json patch;
    std::ifstream patchFile(args.patch_config);
    try {
        patchFile >> patch;
    } catch (...) {
        std::cerr << "Error parsing the json file" << std::endl;
        return 0;
    }

    std::map<Eigen::Vector3i, json, myless> material_patterns;
    int i = 0;
    for (auto entry : patch)
    {
        Eigen::Vector3i x;// = entry["index"];
        std::copy_n(entry["index"].begin(), 3, x.data());
        entry["i"] = i++;
        material_patterns[x] = entry;
    }

    Eigen::MatrixXd V, Vsurf;
    Eigen::MatrixXi T, Fsurf;
    igl::readMSH(args.volume, V, T);

    if (is_file_exist(args.surface.c_str())) {
        Eigen::MatrixXi tmpF;
        igl::readOBJ(args.surface, Vsurf, Fsurf);
    }
    else {
        Eigen::MatrixXi tmpF;
        igl::boundary_facets(T, tmpF);
        Eigen::VectorXi tmpI;
        igl::remove_unreferenced(V, tmpF, Vsurf, Fsurf, tmpI);
    }

    std::vector<std::vector<int>> cell_tets(material_patterns.size());
    {
        Eigen::MatrixXd centers;
        igl::barycenter(V, T, centers);
        centers /= args.gridSize;
        
        const int n_tets = T.rows();
        for (i = 0; i < n_tets; i++)
        {
            Eigen::Vector3i center;
            center << floor(centers(i, 0)), floor(centers(i, 1)), floor(centers(i, 2));
            if (auto search = material_patterns.find(center); search != material_patterns.end())
                cell_tets[search->second["i"]].push_back(i);
        }
    }

    const float voxel = args.gridSize / args.resolution;
    math::Transform::Ptr xform = math::Transform::createLinearTransform(voxel);
    FloatGrid::Ptr whole_grid;
    {
        std::vector<Vec3s> Vsurf_tmp;
        std::vector<Vec3I> Fsurf_tmp;
        eigen_to_openvdb(Vsurf, Vsurf_tmp);
        eigen_to_openvdb(Fsurf, Fsurf_tmp);
        whole_grid = tools::meshToLevelSet<FloatGrid>(*xform, Vsurf_tmp, Fsurf_tmp, 3);
    }

    // center of internal voids
    Eigen::MatrixXd centers;
    centers.setZero(cell_tets.size(), 3);
    Eigen::Matrix<bool, -1, 1> void_flags;
    void_flags.setZero(cell_tets.size());

    std::vector<FloatGrid::Ptr> void_grids(material_patterns.size());
    for (auto const& it : material_patterns)
    {
        const int i = it.second["i"];
        const auto& sub_ids = cell_tets[i];
        if (sub_ids.size() == 0 || it.second["internal"].get<bool>())
            continue;

        std::vector<Vec3s> subV_openvdb, outV;
        std::vector<Vec3I> subF_openvdb, outF;
        // intersection of cube and volume
        {
            Vec3f center(it.first(0) + 0.5, it.first(1) + 0.5, it.first(2) + 0.5);
            FloatGrid::Ptr tmp_grid = openvdb::tools::createLevelSetCube<FloatGrid>(args.gridSize, args.gridSize * center, voxel);
            FloatGrid::Ptr cell_grid = openvdb::tools::csgIntersectionCopy(*whole_grid, *tmp_grid);
            std::vector<Vec4I> Quad;
            tools::volumeToMesh(*gridPtrCast<FloatGrid>(cell_grid), subV_openvdb, subF_openvdb, Quad, 0, 0, true);
            clean_quads(subF_openvdb, Quad);
        }

        if (subF_openvdb.size() == 0)
            continue;
        
        if (it.second["ratio"] >= 1 - 1e-5)
            continue;

        erode(subV_openvdb, subF_openvdb, it.second["ratio"], 1e-2, outV, outF);
        
        for (const auto &p : outV)
            for (int d = 0; d < 3; d++)
                centers(i, d) += p(d);
        centers.row(i) /= outV.size();

        if (args.tunnel > 0)
        {
            // check if the cell is large enough to hold a sphere, if not, ignore the void inside
            Eigen::MatrixXd outV_;
            Eigen::MatrixXi outF_;
            openvdb_to_eigen(outV, outV_);
            openvdb_to_eigen(outF, outF_);
            Eigen::VectorXd sqrD;
            Eigen::VectorXi I;
            Eigen::MatrixXd C;
            Eigen::MatrixXd center = centers.row(i);
            igl::point_mesh_squared_distance(center, outV_, outF_, sqrD, I, C);

            if (sqrt(sqrD.minCoeff()) <= args.tunnel + voxel)
            {
                void_flags[i] = false;
                continue;
            }

            void_flags[i] = true;
            void_grids[i] = tools::meshToLevelSet<FloatGrid>(*xform, outV, outF, 3);

            Vec3f center_(center(0),center(1),center(2));
            auto grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>((float)args.tunnel, center_, (float)voxel);
            openvdb::tools::csgUnion(*(void_grids[i]), *grid);
        }
        else
        {
            void_flags[i] = true;
            void_grids[i] = tools::meshToLevelSet<FloatGrid>(*xform, outV, outF, 3);
        }
        openvdb::tools::csgDifference(*whole_grid, *(void_grids[i]));

        if (args.tunnel > 0 && void_flags[i])
        {
            Eigen::Vector3i index;
            for (int k = -1; k <= 1; k += 2)
            {
                for (int d = 0; d < 3; d++)
                {
                    index = it.first;
                    index(d) += k;
                    if (auto search = material_patterns.find(index); search != material_patterns.end())
                    {
                        const int j = search->second["i"];
                        if (void_flags[j])
                        {
                            // check whether the tunnel break any object surface
                            {
                                bool intersect_flag = false;
                                for (int v = 0; v < Vsurf.rows(); v++)
                                {
                                    if (cylinderSDF(centers.row(i), centers.row(j), args.tunnel, Vsurf.row(v)) < voxel)
                                    {
                                        intersect_flag = true;
                                        break;
                                    }
                                }
                                
                                if (intersect_flag)
                                    continue;
                            }
                            std::cout << it.first.transpose() << ", " << search->first.transpose() << "\n";
                            
                            auto grid = createLevelSetCylinder(centers.row(i), centers.row(j), args.tunnel, voxel); 
                            openvdb::tools::csgDifference(*whole_grid, *grid);

                            // tools::Diagnose<FloatGrid> d(*whole_grid);
                            // tools::CheckNan<FloatGrid> c;
                            // std::string str = d.check(c);
                            // if (!str.empty())
                            // {
                            //     grid = createLevelSetCylinder(centers.row(i), centers.row(search->second["i"]), args.tunnel, voxel);
                            //     tools::Diagnose<FloatGrid> d(*grid);
                            //     tools::CheckNan<FloatGrid> c;
                            //     std::string str = d.check(c);
                            //     throw std::runtime_error("cylinder grid error between " + std::to_string(i.get<int>()) + " " + std::to_string(search->second["i"].get<int>()) + " : " + str);
                            // }
                        }
                        else if (search->second["internal"].get<bool>())
                        {
                            std::cout << it.first.transpose() << ", " << search->first.transpose() << "\n";

                            Eigen::Vector3d top = centers.row(i);
                            top(d) = (it.first(d) + 0.5 + 0.5 * k) * args.gridSize;
                            
                            auto grid = createLevelSetCylinder(centers.row(i), top, args.tunnel, voxel); 
                            openvdb::tools::csgDifference(*whole_grid, *grid);
                        }
                    }
                }
            }
        
            // tunnel to the object surface
            Eigen::MatrixXd center = centers.row(i);
            Eigen::VectorXd sqrD;
            Eigen::VectorXi I;
            Eigen::MatrixXd C;
            igl::point_mesh_squared_distance(center, Vsurf, Fsurf, sqrD, I, C);

            Eigen::Index minRow;
            double min = sqrD.minCoeff(&minRow);

            auto grid = createLevelSetCylinder(center, C.row(minRow), args.tunnel, voxel); 
            openvdb::tools::csgDifference(*whole_grid, *grid);

            Vec3f sphere_center(C(minRow,0),C(minRow,1),C(minRow,2));
            grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>((float)args.tunnel, sphere_center, (float)voxel);
            openvdb::tools::csgDifference(*whole_grid, *grid);
        }
    }

    std::vector<Vec3s> Ve;
    std::vector<Vec3I> Tri;
    sdf2mesh(whole_grid, Ve, Tri);

    if (args.preserve_original_surface)
    {
        Eigen::MatrixXi F;
        Eigen::MatrixXd V;
        openvdb_to_eigen(Tri, F);
        openvdb_to_eigen(Ve, V);

        Eigen::SparseMatrix<int> adj;
        igl::adjacency_matrix(F, adj);

        Eigen::VectorXi C, K;
        int nc = connected_components(adj, C, K);

        Eigen::Index maxv;
        V.rowwise().squaredNorm().maxCoeff(&maxv);

        std::vector<int> hole_faces;
        for (int i = 0; i < F.rows(); i++)
            if (C[maxv] != C[F(i, 0)])
                hole_faces.push_back(i);
        
        Eigen::VectorXi tmpI;
        Eigen::MatrixXd Vholes;
        Eigen::MatrixXi Fholes;
        igl::remove_unreferenced(V, F(hole_faces, Eigen::all), Vholes, Fholes, tmpI);

        Fholes = (Fholes.array() + Vsurf.rows()).matrix().eval();

        F.setZero(Fholes.rows() + Fsurf.rows(), Fholes.cols());
        F << Fholes, Fsurf;

        V.setZero(Vsurf.rows() + Vholes.rows(), Vsurf.cols());
        V << Vsurf, Vholes;

        eigen_to_openvdb(V, Ve);
        eigen_to_openvdb(F, Tri);
    }

    write_mesh(args.output, Ve, Tri);

    return 0;
}