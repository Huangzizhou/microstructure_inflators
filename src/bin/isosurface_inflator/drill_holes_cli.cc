////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>
#include <isosurface_inflator/VDBTools.hh>

////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

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

    std::cout << "target volume ratio: " << hole_vol << "\n";
    std::cout << "half diag: " << half_diag << "\n";

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

            cout << "current radius interval: " << upper << "\n";
            cout << "current volume interval: " << upper_vol << "\n";
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

            cout << "current radius interval: " << lower << " " << current << " " << upper << "\n";
            cout << "current volume interval: " << lower_vol << " " << current_vol << " " << upper_vol << "\n";

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
        std::string volume;
        std::string patch_config;
        std::string output = "out.obj";
        double gridSize = 0.1;
        int resolution = 50;
    } args;

    // Parse arguments
    CLI::App app{"drill_holes_cli"};
    app.add_option("volume,--vol", args.volume, "Volume mesh.")->required();
    app.add_option("patch,-p,--patch", args.patch_config, "Patch description (json file).")->required();
    app.add_option("--gridSize", args.gridSize, "Grid size.")->required();
    app.add_option("-o,--output", args.output, "Output triangle mesh.");
    app.add_option("-r,--resolution", args.resolution, "Density field resolution.");
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

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

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    std::vector<std::vector<int>> cell_tets(material_patterns.size());
    {
        igl::readMSH(args.volume, V, T);
        
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

    FloatGrid::Ptr whole_grid = volmesh2sdf(V, T, 1. / args.resolution);

    for (auto const& it : material_patterns)
    {
        const auto& sub_ids = cell_tets[it.second["i"]];
        if (sub_ids.size() == 0)
            continue;
        
        Eigen::MatrixXi tmpT(sub_ids.size(), 4);
        for (i = 0; i < sub_ids.size(); i++)
            tmpT.row(i) = T.row(sub_ids[i]);

        std::vector<Vec3s> subV_openvdb, outV;
        std::vector<Vec3I> subF_openvdb, outF;
        volmesh2surface(V, tmpT, subV_openvdb, subF_openvdb);

        if (subF_openvdb.size() == 0)
            continue;

        auto grid = erode(subV_openvdb, subF_openvdb, it.second["ratio"], 1e-2, outV, outF);
        math::Transform::Ptr xform = math::Transform::createLinearTransform(1. / args.resolution);
        grid = tools::meshToLevelSet<FloatGrid>(*xform, outV, outF, 3);
        openvdb::tools::csgDifference(*whole_grid, *grid);
    }

    std::vector<Vec3s> Ve;
    std::vector<Vec3I> Tri;
    sdf2mesh(whole_grid, Ve, Tri);
    write_mesh(args.output, Ve, Tri);

    return 0;
}