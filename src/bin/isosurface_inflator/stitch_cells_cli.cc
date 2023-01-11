////////////////////////////////////////////////////////////////////////////////
// Example run:
// ./isosurface_inflator/stitch_cells_cli $MICRO_DIR/isosurface_inflator/tests/patch.json -o out.obj
////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>

#include <igl/readOBJ.h>
#include <igl/adjacency_matrix.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/SignedFloodFill.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Composite.h>

////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;
using WireMeshBasePtr = std::shared_ptr<WireMeshBase>;
using namespace openvdb;

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

std::string lowercase(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
    return data;
}

#define TRY_SYMMETRY(s, x, p)                                  \
    if (lowercase(x) == lowercase(#s))                         \
    {                                                          \
        return std::make_shared<WireMesh<Symmetry::s<>>>((p)); \
    }

#define TRY_KEY_VAL(s, a, x, p)                                \
    if (lowercase(x) == lowercase(#a))                         \
    {                                                          \
        return std::make_shared<WireMesh<Symmetry::s<>>>((p)); \
    }

WireMeshBasePtr load_wire_mesh(const std::string &sym, const std::string &path) {
    TRY_SYMMETRY(Square, sym, path);
    TRY_SYMMETRY(Cubic, sym, path);
    TRY_SYMMETRY(Orthotropic, sym, path);
    TRY_SYMMETRY(Diagonal, sym, path);
    TRY_KEY_VAL(DoublyPeriodic, Doubly_Periodic, sym, path);
    TRY_KEY_VAL(TriplyPeriodic, Triply_Periodic, sym, path);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

/*
patch json format
[
    {
        "params": [
            0.5,
            0.333333,
            0.666667,
            0.333333,
            ...
        ],
        "symmetry": "Orthotropic",
        "pattern": "./data/patterns/3D/reference_wires/pattern0646.wire",
        "index": [2,2,3]
    },
    ...
]
*/

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
    // Default arguments
    struct {
        std::string patch_config;
        std::string object_surface = "";
        std::string output = "out.obj";
        double gridSize = 0.1;
        int resolution = 50;
    } args;

    // Parse arguments
    CLI::App app{"stitch_cells_cli"};
    app.add_option("patch,-p,--patch", args.patch_config, "Patch description (json file).")->required();
    app.add_option("--gridSize", args.gridSize, "Grid size.")->required();
    app.add_option("--surface", args.object_surface, "Object surface.");
    app.add_option("-o,--output", args.output, "Output triangle mesh.");
    app.add_option("-r,--resolution", args.resolution, "Density field resolution.");
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    // execute<3>(args.patch_config, args.output, args.resolution);
    const int resolution = args.resolution;

    // Load patch config
    json patch;
    std::ifstream patchFile(args.patch_config);
    try {
        patchFile >> patch;
    } catch (...) {
        std::cerr << "Error parsing the json file" << std::endl;
        return 0;
    }

    openvdb::initialize();


    std::map<Eigen::Vector3i, json, myless> material_patterns;
    for (auto entry : patch)
    {
        Eigen::Vector3i x;// = entry["index"];
        std::copy_n(entry["index"].begin(), 3, x.data());
        material_patterns[x] = entry;
    }

    /* create sdf with internal microstructure cells empty */

    FloatGrid::Ptr surf_grid;
    {
        // math::Transform::Ptr xform = math::Transform::createLinearTransform(args.gridSize / (resolution - 1));
        // surf_grid = tools::meshToLevelSet<FloatGrid>(*xform, SV, SF, 3);
        surf_grid = mesh2sdf(args.object_surface, args.gridSize / (resolution - 1));

        for (auto const& it : material_patterns)
        {
            Eigen::Vector3i id = it.first;
            Vec3f center(id(0) + 0.5,id(1) + 0.5,id(2) + 0.5);
            auto tmp_grid = openvdb::tools::createLevelSetCube<FloatGrid>((resolution - 1), center * (resolution - 1), 1);
            openvdb::tools::csgDifference(*surf_grid, *tmp_grid);
        }
    }

    /* create sdf for internal microstructure cells */

    const double bg_val = 3.0 / resolution;
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(bg_val);
    grid->setTransform(math::Transform::createLinearTransform(args.gridSize / (resolution - 1)));
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

    int cur = 0;
    for (auto const& it : material_patterns)
    {
        // Assign topologies and parameters to each cell
        NDCubeArray<WireMeshBasePtr, 3, 3> topologyGrid;
        NDCubeArray<std::vector<double>, 3, 3> parameterGrid;

        for (int i = -1; i <= 1; i+=1)
        for (int d = 0; d < 3; d++)
        {
            if (d != 0 && i == 0)
                continue;
            Eigen::Vector3i id = it.first;
            id[d] += i;
            json entry;
            if (auto search = material_patterns.find(id); search != material_patterns.end())
                entry = search->second;
            else
                entry = it.second;
            
            NDArrayIndex<3> index;
            Eigen::Vector3i local_id = (id - it.first).array() + 1;
            std::copy_n(local_id.data(), 3, index.idxs.begin());

            topologyGrid(index) = load_wire_mesh(entry["symmetry"], entry["pattern"]);
            parameterGrid(index) = entry["params"].get<std::vector<double>>();
        }

        auto swm = make_stitched_wire_mesh<3, true>(topologyGrid);

        auto params = swm.paramsFromParamGrid(parameterGrid);

        PatternSignedDistance<double, StitchedWireMesh<3, true>> sdf(swm);
        sdf.setUseAabbTree(true);
        std::unique_ptr<MesherBase> mesher = std::make_unique<IGLSurfaceMesherMC>();
        sdf.setParameters(params, mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
        
        const auto &bbox = sdf.boundingBox();
        Point3d minCorner = bbox.minCorner;
        Point3d maxCorner = bbox.maxCorner;

        const size_t nsamples = resolution * resolution * resolution;

        // Evaluation point locations;
        // flattened to be accessed as:
        // xi + resolution * (yi + resolution * zi)
        Eigen::MatrixXd sampleLocations(nsamples, 3);
        Eigen::MatrixXd ijks(nsamples, 3);
        {
            size_t i = 0;
            for (size_t zi = 0; zi < resolution; ++zi) {
                for (size_t yi = 0; yi < resolution; ++yi) {
                    for (size_t xi = 0; xi < resolution; ++xi) {
                        sampleLocations.row(i) = bbox.interpolatePoint(
                                Point3D(xi / Real(resolution - 1.0),
                                        yi / Real(resolution - 1.0),
                                        zi / Real(resolution - 1.0)));
                        ijks.row(i) << xi, yi, zi;
                        ++i;
                    }
                }
            }
        }

        // Evaluate signed distances at each grid point
        Eigen::VectorXd signedDistances(nsamples);
    #if MICRO_WITH_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
                [&](const tbb::blocked_range<size_t> &r) {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                        signedDistances(i) = sdf.signedDistance(sampleLocations.row(i));
                }
            );
    #else
        for (size_t i = 0; i < nsamples; ++i)
            signedDistances(i) = sdf.signedDistance(sampleLocations.row(i));
    #endif

        int compressed = 0;
        Eigen::Vector3i offset = it.first * (resolution - 1);
        for (int i = 0; i < nsamples; i++)
        {
            if (signedDistances(i) < bg_val)
                accessor.setValue(openvdb::Coord(ijks(i, 0) + offset(0), ijks(i, 1) + offset(1), ijks(i, 2) + offset(2)), signedDistances(i));
            else
                compressed++;
        }

        std::cout << "Saved space " << compressed / (double) nsamples * 100 << " %\n";
        std::cout << "Completed " << ((++cur) * 100.0) / material_patterns.size() << " %\n";
    }

    openvdb::tools::csgUnion(*grid, *surf_grid);

    /* create tunnels to remove internal materials after 3d printing */
    
    double tunnel_size = 0.2;
    for (auto const& it : material_patterns)
    {
        for (int i = -1; i <= 1; i+=1)
        for (int d = 0; d < 3; d++)
        {
            if (d != 0 && i == 0)
                continue;
            Eigen::Vector3i id = it.first;
            id[d] += i;
            if (material_patterns.find(id) == material_patterns.end()) // need to create tunnel
            {
                Vec3f center(id(0) + 0.5,id(1) + 0.5,id(2) + 0.5);
                Vec3f corner1 = center - tunnel_size / 2;
                Vec3f corner2 = center + tunnel_size / 2;
                corner1(d) = id(d);
                corner2(d) = id(d) + 1;
                openvdb::math::BBox<Vec3f> bbox(corner1 * (resolution - 1), corner2 * (resolution - 1));
                math::Transform::Ptr xform = math::Transform::createLinearTransform(1);
                auto tmp_grid = openvdb::tools::createLevelSetBox<FloatGrid>(bbox, *xform);
                openvdb::tools::csgDifference(*grid, *tmp_grid);
            }
        }
    }

    /* remove top */
{
    Vec3f corner1(-5,-5,1);
    Vec3f corner2(5,5,2);
    openvdb::math::BBox<Vec3f> bbox(corner1 * (resolution - 1), corner2 * (resolution - 1));
    math::Transform::Ptr xform = math::Transform::createLinearTransform(1);
    auto tmp_grid = openvdb::tools::createLevelSetBox<FloatGrid>(bbox, *xform);
    openvdb::tools::csgDifference(*grid, *tmp_grid);
}{
    Vec3f corner1(-5,-5,-5);
    Vec3f corner2(0,5,5);
    openvdb::math::BBox<Vec3f> bbox(corner1 * (resolution - 1), corner2 * (resolution - 1));
    math::Transform::Ptr xform = math::Transform::createLinearTransform(1);
    auto tmp_grid = openvdb::tools::createLevelSetBox<FloatGrid>(bbox, *xform);
    openvdb::tools::csgDifference(*grid, *tmp_grid);
}{
    Vec3f corner1(-5,-5,-5);
    Vec3f corner2(5,0,5);
    openvdb::math::BBox<Vec3f> bbox(corner1 * (resolution - 1), corner2 * (resolution - 1));
    math::Transform::Ptr xform = math::Transform::createLinearTransform(1);
    auto tmp_grid = openvdb::tools::createLevelSetBox<FloatGrid>(bbox, *xform);
    openvdb::tools::csgDifference(*grid, *tmp_grid);
}
    /* sdf to mesh */

    // openvdb::tools::signedFloodFill(grid->tree());
    grid->setName("density");
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    std::vector<Vec3s> Ve;
    std::vector<Vec3I> Tri;
    std::vector<Vec4I> Quad;
    tools::volumeToMesh(*grid, Ve, Tri, Quad, 0, 0, true);
    clean_quads(Tri, Quad);
    write_mesh(args.output, Ve, Tri, std::vector<Vec4I>());

    // openvdb::io::File file(args.output);
    // openvdb::GridPtrVec(grids);
    // grids.push_back(grid);
    // file.write(grids);
    // file.close();

    return 0;
}
