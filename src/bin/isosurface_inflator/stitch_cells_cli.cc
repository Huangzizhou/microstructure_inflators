////////////////////////////////////////////////////////////////////////////////
// Example run:
// ./isosurface_inflator/stitch_cells_cli $MICRO_DIR/isosurface_inflator/tests/patch.json -o out.obj
////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>

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
        fprintf(obj_file," %u",Tri[i](j)+1);
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

openvdb::FloatGrid::Ptr sdf2vdb(const SignedDistanceRegion<3> &sdf, const int resolution, const double background)
{
    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(background);
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    
    const auto &bbox = sdf.boundingBox();
    Point3d minCorner = bbox.minCorner;
    Point3d maxCorner = bbox.maxCorner;
    std::cout << "minCorner " << minCorner.transpose() << "\n";
    std::cout << "maxCorner " << maxCorner.transpose() << "\n";

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
    for (int i = 0; i < nsamples; i++)
    {
        if (abs(signedDistances(i)) < background)
            accessor.setValue(openvdb::Coord(ijks(i, 0), ijks(i, 1), ijks(i, 2)), signedDistances(i));
        else
            compressed++;
    }
    openvdb::tools::signedFloodFill(grid->tree());

    std::cout << "Saved space " << compressed / (double) nsamples * 100 << " %\n";

    grid->setName("density");
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    return grid;
}

////////////////////////////////////////////////////////////////////////////////

template<size_t N>
void execute(const std::string &patchFilename, const std::string &meshingOptions, const std::string &outname) {
    // Load patch config
    json patch;
    std::ifstream patchFile(patchFilename);
    try {
        patchFile >> patch;
    } catch (...) {
        std::cerr << "Error parsing the json file" << std::endl;
        return;
    }

    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }
    if (!meshingOptions.empty()) { mesher->meshingOptions.load(meshingOptions); }

    // Assign topologies and parameters to each cell
    NDCubeArray<WireMeshBasePtr, N, 3> topologyGrid;
    NDCubeArray<std::vector<double>, N, 3> parameterGrid;

    for (auto entry : patch) {
        NDArrayIndex<N> index;
        std::copy_n(entry["index"].begin(), N, index.idxs.begin());
        topologyGrid(index) = load_wire_mesh(entry["symmetry"], entry["pattern"]);
        parameterGrid(index) = entry["params"].get<std::vector<double>>();
    }

    auto swm = make_stitched_wire_mesh<N, true>(topologyGrid);

    auto params = swm.paramsFromParamGrid(parameterGrid);

    PatternSignedDistance<double, StitchedWireMesh<N, true>> sdf(swm);
    sdf.setUseAabbTree(true);
    sdf.setParameters(params, mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);

    if (N == 3) {
        const int resolution = 100;
        const double bg_val = 5.0 / resolution;
        auto grid = sdf2vdb(sdf, resolution, bg_val);

        openvdb::io::File file(outname);

        openvdb::GridPtrVec(grids);
        grids.push_back(grid);
        file.write(grids);
        file.close();
    }
}

template<size_t N>
void execute(const std::string &patchFilename, const std::string &outname, const int resolution) {
    // Load patch config
    json patch;
    std::ifstream patchFile(patchFilename);
    try {
        patchFile >> patch;
    } catch (...) {
        std::cerr << "Error parsing the json file" << std::endl;
        return;
    }

    assert(N == 3);

    std::map<Eigen::Vector3i, json, myless> material_patterns;
    for (auto entry : patch)
    {
        Eigen::Vector3i x;// = entry["index"];
        std::copy_n(entry["index"].begin(), 3, x.data());
        material_patterns[x] = entry;
    }

    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }

    // const int resolution = 50;
    const double bg_val = 3.0 / resolution;
    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(bg_val);
    // grid->setTransform(math::Transform::createLinearTransform(1. / (resolution - 1)));
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

    int cur = 0;
    for (auto const& it : material_patterns)
    {
        // Assign topologies and parameters to each cell
        NDCubeArray<WireMeshBasePtr, N, 3> topologyGrid;
        NDCubeArray<std::vector<double>, N, 3> parameterGrid;

        bool store_all = false;
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
            {
                entry = it.second;
                // store_all = true;
            }
            
            NDArrayIndex<N> index;
            Eigen::Vector3i local_id = (id - it.first).array() + 1;
            std::copy_n(local_id.data(), N, index.idxs.begin());

            topologyGrid(index) = load_wire_mesh(entry["symmetry"], entry["pattern"]);
            parameterGrid(index) = entry["params"].get<std::vector<double>>();
        }

        auto swm = make_stitched_wire_mesh<N, true>(topologyGrid);

        auto params = swm.paramsFromParamGrid(parameterGrid);

        PatternSignedDistance<double, StitchedWireMesh<N, true>> sdf(swm);
        sdf.setUseAabbTree(true);
        sdf.setParameters(params, mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
        
        const auto &bbox = sdf.boundingBox();
        Point3d minCorner = bbox.minCorner;
        Point3d maxCorner = bbox.maxCorner;
        std::cout << "minCorner " << minCorner.transpose() << "\n";
        std::cout << "maxCorner " << maxCorner.transpose() << "\n";

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
            if (store_all || signedDistances(i) < bg_val)
                accessor.setValue(openvdb::Coord(ijks(i, 0) + offset(0), ijks(i, 1) + offset(1), ijks(i, 2) + offset(2)), signedDistances(i));
            else
                compressed++;
        }

        std::cout << "Saved space " << compressed / (double) nsamples * 100 << " %\n";
        std::cout << "Completed " << ((++cur) * 100.0) / material_patterns.size() << " %\n";

        for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++)
        for (int k = -1; k <= 1; k++)
        {
            Eigen::Vector3i id = it.first;
            id[0] += i;
            id[1] += j;
            id[2] += k;
            if (material_patterns.count(id) == 0)
            {
                Vec3f center(id(0) + 0.5,id(1) + 0.5,id(2) + 0.5);
                auto tmp_grid = openvdb::tools::createLevelSetCube<FloatGrid>((resolution - 1), center * (resolution - 1), 1);
                // tmp_grid->setTransform(math::Transform::createLinearTransform(1. / (resolution - 1)));
                openvdb::tools::csgUnion(*grid, *tmp_grid);
            }
        }
    }

    // openvdb::tools::signedFloodFill(grid->tree());
    grid->setName("density");
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    std::vector<Vec3s> Ve;
    std::vector<Vec3I> Tri;
    std::vector<Vec4I> Quad;
    tools::volumeToMesh(*grid, Ve, Tri, Quad, 0, 0, true);
    for (auto &v : Ve)
        v /= resolution;
    clean_quads(Tri, Quad);
    write_mesh(outname, Ve, Tri, std::vector<Vec4I>());

    // openvdb::io::File file(outname);
    // openvdb::GridPtrVec(grids);
    // grids.push_back(grid);
    // file.write(grids);
    // file.close();
}
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
    // Default arguments
    struct {
        std::string patch_config;
        std::string output = "out.obj";
        int resolution = 50;
    } args;

    // Parse arguments
    CLI::App app{"stitch_cells_cli"};
    app.add_option("patch,-p,--patch", args.patch_config, "3x3 patch description (json file).")->required();
    app.add_option("output,-o,--output", args.output, "Output triangle mesh.");
    app.add_option("-r,--resolution", args.resolution, "Density field resolution.");
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    execute<3>(args.patch_config, args.output, args.resolution);

    return 0;
}
