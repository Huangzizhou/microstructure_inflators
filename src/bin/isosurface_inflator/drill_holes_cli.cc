////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>

#include <igl/readOBJ.h>
#include <igl/readMSH.h>
#include <igl/barycenter.h>
#include <igl/adjacency_matrix.h>
#include <igl/writeOBJ.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/remove_unreferenced.h>
#include <igl/boundary_facets.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/SignedFloodFill.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetMeasure.h>

////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;
using WireMeshBasePtr = std::shared_ptr<WireMeshBase>;
using namespace openvdb;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

struct myless {
    bool operator()(const Eigen::Vector3i &lhs, const Eigen::Vector3i &rhs) const 
    {
        if (lhs[0] != rhs[0])
            return lhs[0] < rhs[0];
        else if (lhs[1] != rhs[1])
            return lhs[1] < rhs[1];
        else if (lhs[2] != rhs[2])
            return lhs[2] < rhs[2];
        else
            return false;
    }
};

void read_mesh(const std::string &input, std::vector<Vec3s> &V, std::vector<Vec3I> &F)
{
    Eigen::MatrixXd v;
    Eigen::MatrixXi f;
    igl::readOBJ(input, v, f);

    V.reserve(v.rows());
    F.reserve(f.rows());

    for (int i = 0; i < v.rows(); i++)
        V.emplace_back(v(i, 0), v(i, 1), v(i, 2));
    
    for (int j = 0; j < f.rows(); j++)
        F.emplace_back(f(j, 0), f(j, 1), f(j, 2));
}

template < typename DerivedC, typename DerivedK>
int connected_components(
  const Eigen::SparseMatrix<int> & A,
  Eigen::PlainObjectBase<DerivedC> & C,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  typedef typename Eigen::SparseMatrix<int>::Index Index;
  const auto m = A.rows();
  assert(A.cols() == A.rows() && "A should be square");
  // 1.1 sec
  // m  means not yet visited
  C.setConstant(m,1,m);
  // Could use amortized dynamic array but didn't see real win.
  K.setZero(m,1);
  typename DerivedC::Scalar c = 0;
  for(Eigen::Index f = 0;f<m;f++)
  {
    // already seen
    if(C(f)<m) continue;
    // start bfs
    std::queue<Index> Q;
    Q.push(f);
    while(!Q.empty())
    {
      const Index g = Q.front();
      Q.pop();
      // already seen
      if(C(g)<m) continue;
      // see it
      C(g) = c;
      K(c)++;
      for(typename Eigen::SparseMatrix<int>::InnerIterator it (A,g); it; ++it)
      {
        const Index n = it.index();
        // already seen
        if(C(n)<m) continue;
        Q.push(n);
      }
    }
    c++;
  }
  K.conservativeResize(c,1);
  return c;
}

template < typename DerivedC, typename DerivedD >
void eigen_to_openvdb(const Eigen::MatrixBase<DerivedC> &mat, std::vector<openvdb::v10_0::math::Vec3<DerivedD>> &out)
{
    assert(mat.cols() == 3);
    out.reserve(mat.rows());
    for (int i = 0; i < mat.rows(); i++)
    {
        out.emplace_back(mat(i, 0), mat(i, 1), mat(i, 2));
    }
}

FloatGrid::Ptr mesh2sdf(const std::string &path, const double voxel)
{
    Eigen::MatrixXd v;
    Eigen::MatrixXi f;
    igl::readOBJ(path, v, f);

    Eigen::SparseMatrix<int> adj;
    igl::adjacency_matrix(f, adj);

    Eigen::VectorXi C, K;
    int nc = connected_components(adj, C, K);

    Eigen::Index maxv;
    v.rowwise().squaredNorm().maxCoeff(&maxv);

    std::vector<Vec3s> SV;
    std::vector<Vec3I> SF;

    for (int i = 0; i < v.rows(); i++)
        SV.emplace_back(v(i,0),v(i,1),v(i,2));

    for (int i = 0; i < f.rows(); i++)
        if (C[f(i,0)] == C[maxv])
            SF.emplace_back(f(i,0),f(i,1),f(i,2));

    math::Transform::Ptr xform = math::Transform::createLinearTransform(voxel);
    FloatGrid::Ptr grid = tools::meshToLevelSet<FloatGrid>(*xform, SV, SF, 3);

    for (int k = 0; k < K.size(); k++)
    {
        if (k == C[maxv])
            continue;
        
        SF.clear();
        for (int i = 0; i < f.rows(); i++)
            if (C[f(i,0)] == k)
                SF.emplace_back(f(i,0),f(i,1),f(i,2));

        math::Transform::Ptr xform = math::Transform::createLinearTransform(voxel);
        FloatGrid::Ptr tmp_grid = tools::meshToLevelSet<FloatGrid>(*xform, SV, SF, 3);

        openvdb::tools::csgDifference(*grid, *tmp_grid);
    }

    return grid;
}

void clean_quads(std::vector<Vec3I> &Tri, std::vector<Vec4I> &Quad)
{
    for (const auto &f : Quad)
    {
        Tri.emplace_back(f(0), f(1), f(2));
        Tri.emplace_back(f(0), f(2), f(3));
    }

    Quad.clear();
}

void write_mesh(const std::string &output, const std::vector<Vec3s> &V, const std::vector<Vec3I> &Tri, const std::vector<Vec4I> &Quad)
{
    FILE * obj_file = fopen(output.c_str(),"w");
    
    for(int i = 0;i<(int)V.size();i++)
    {
        fprintf(obj_file,"v");
        for(int j = 0;j<3;++j)
        {
        fprintf(obj_file," %0.17g", V[i](j));
        }
        fprintf(obj_file,"\n");
    }

    for(int i = 0;i<(int)Tri.size();++i)
    {
        fprintf(obj_file,"f");
        for(int j = 0; j<3;++j)
        {
        // OBJ is 1-indexed
        fprintf(obj_file," %u",Tri[i](2 - j)+1);
        }
        fprintf(obj_file,"\n");
    }

    for(int i = 0;i<(int)Quad.size();++i)
    {
        fprintf(obj_file,"f");
        for(int j = 0; j<4;++j)
        {
        // OBJ is 1-indexed
        fprintf(obj_file," %u",Quad[i](j)+1);
        }
        fprintf(obj_file,"\n");
    }

    fclose(obj_file);
}

double SignedVolumeOfTriangle(const Vec3s &p1, const Vec3s &p2, const Vec3s &p3) {
    double v321 = p3(0)*p2(1)*p1(2);
    double v231 = p2(0)*p3(1)*p1(2);
    double v312 = p3(0)*p1(1)*p2(2);
    double v132 = p1(0)*p3(1)*p2(2);
    double v213 = p2(0)*p1(1)*p3(2);
    double v123 = p1(0)*p2(1)*p3(2);
    return (1 / 6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

double compute_surface_mesh_volume(const std::vector<Vec3s> &V, const std::vector<Vec3I> &Tri)
{
    Vec3s center(0, 0, 0);
    for (const auto &v : V)
        center += v;
    center /= V.size();
    
    double vol = 0;
    for (const auto &f : Tri)
        vol += SignedVolumeOfTriangle(V[f(0)] - center, V[f(1)] - center, V[f(2)] - center);

    return vol;
}

double volume_after_offset(const FloatGrid::Ptr &grid, const double radius, std::vector<Vec3s> &Ve, std::vector<Vec3I> &Tri)
{
    // erode
    auto tmp_grid = grid->deepCopyGrid();
    std::unique_ptr<tools::LevelSetFilter<FloatGrid>> filter = std::make_unique<tools::LevelSetFilter<FloatGrid>>(*gridPtrCast<FloatGrid>(tmp_grid));
    filter->offset(radius);

    // ls2mesh
    Ve.clear();
    Tri.clear();
    std::vector<Vec4I> Quad;
    tools::volumeToMesh(*gridPtrCast<FloatGrid>(tmp_grid), Ve, Tri, Quad, 0, 0, true);
    clean_quads(Tri, Quad);

    double vol = abs(compute_surface_mesh_volume(Ve, Tri));
        
    return vol;
}

FloatGrid::Ptr erode(const std::vector<Vec3s> &V, const std::vector<Vec3I> &F, const double volume_ratio, double tol, std::vector<Vec3s> &Ve, std::vector<Vec3I> &Tri) 
{
    // std::vector<Vec3s> V;
    // std::vector<Vec3I> F;
    // read_mesh(input, V, F);

    double half_diag;
    {
        Eigen::MatrixXd V_mat(V.size(), 3);
        for (int i = 0; i < V.size(); i++)
            for (int d = 0; d < 3; d++)
                V_mat(i, d) = V[i](d);
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

    Eigen::MatrixXd V, centers;
    Eigen::MatrixXi T;
    igl::readMSH(args.volume, V, T);
    igl::barycenter(V, T, centers);
    const int n_tets = T.rows();

    centers /= args.gridSize;
    std::vector<std::vector<int>> cell_tets(material_patterns.size());
    for (i = 0; i < n_tets; i++)
    {
        Eigen::Vector3i center;
        center << floor(centers(i, 0)), floor(centers(i, 1)), floor(centers(i, 2));
        if (auto search = material_patterns.find(center); search != material_patterns.end())
        {
            // std::cout << search->second << "\n";
            cell_tets[search->second["i"]].push_back(i);
        }
    }

    FloatGrid::Ptr whole_grid;
    {
        Eigen::MatrixXi F;
        igl::boundary_facets(T, F);

        Eigen::MatrixXi subF;
        Eigen::MatrixXd subV;
        Eigen::VectorXi tmpI;
        igl::remove_unreferenced(V, F, subV, subF, tmpI);

        std::vector<Vec3s> subV_openvdb;
        std::vector<Vec3I> subF_openvdb;
        eigen_to_openvdb(subV, subV_openvdb);
        eigen_to_openvdb(subF, subF_openvdb);

        math::Transform::Ptr xform = math::Transform::createLinearTransform(1. / args.resolution);
        whole_grid = tools::meshToLevelSet<FloatGrid>(*xform, subV_openvdb, subF_openvdb, 3);
    }

    for (auto const& it : material_patterns)
    {
        const auto& sub_ids = cell_tets[it.second["i"]];
        if (sub_ids.size() == 0)
            continue;
        Eigen::MatrixXi tmpT(sub_ids.size(), 4);
        Eigen::MatrixXi subF;
        Eigen::MatrixXd subV;

        for (i = 0; i < sub_ids.size(); i++)
            tmpT.row(i) = T.row(sub_ids[i]);

        Eigen::MatrixXi F;
        igl::boundary_facets(tmpT, F);

        Eigen::VectorXi tmpI;
        igl::remove_unreferenced(V, F, subV, subF, tmpI);

        std::vector<Vec3s> subV_openvdb, outV;
        std::vector<Vec3I> subF_openvdb, outF;
        eigen_to_openvdb(subV, subV_openvdb);
        eigen_to_openvdb(subF, subF_openvdb);

        if (subF_openvdb.size() == 0)
            continue;

        auto grid = erode(subV_openvdb, subF_openvdb, it.second["ratio"], 1e-2, outV, outF);
        math::Transform::Ptr xform = math::Transform::createLinearTransform(1. / args.resolution);
        grid = tools::meshToLevelSet<FloatGrid>(*xform, outV, outF, 3);
        openvdb::tools::csgDifference(*whole_grid, *grid);
    }

    std::vector<Vec3s> Ve;
    std::vector<Vec3I> Tri;
    std::vector<Vec4I> Quad;
    tools::volumeToMesh(*whole_grid, Ve, Tri, Quad, 0, 0, true);
    clean_quads(Tri, Quad);
    write_mesh(args.output, Ve, Tri, std::vector<Vec4I>());

    return 0;
}