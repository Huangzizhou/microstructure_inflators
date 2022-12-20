#include <isosurface_inflator/WireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/MultiResGrid.h>
#include <openvdb/tools/Mask.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include <iostream>
#include <cstdio>
#include <functional>

namespace po = boost::program_options;
using namespace std;
using namespace openvdb;

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
    double vol = 0;
    for (const auto &f : Tri)
    {
        vol += SignedVolumeOfTriangle(V[f(0)], V[f(1)], V[f(2)]);
    }

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

    return abs(compute_surface_mesh_volume(Ve, Tri));
}

void execute(const string &input, const string &output, const double volume_ratio) 
{
    std::vector<Vec3s> V;
    std::vector<Vec3I> F;
    read_mesh(input, V, F);
    const double initial_vol = abs(compute_surface_mesh_volume(V, F));
    const double hole_vol = 1 - volume_ratio;

    // mesh2ls
    double voxel = 1e-2;
    math::Transform::Ptr xform(nullptr);
    xform = math::Transform::createLinearTransform(voxel);
    FloatGrid::Ptr grid = tools::meshToLevelSet<FloatGrid>(*xform, V, F);

    // eroded mesh
    std::vector<Vec3s> Ve;
    std::vector<Vec3I> Tri;

    double upper = 1;
    double upper_vol = volume_after_offset(grid, upper, Ve, Tri) / initial_vol;
    double lower = 0;
    double lower_vol = volume_after_offset(grid, lower, Ve, Tri) / initial_vol;

    while (upper_vol > hole_vol)
    {
        upper *= 2;
        upper_vol = volume_after_offset(grid, upper, Ve, Tri) / initial_vol;
    }

    assert(upper_vol <= hole_vol && lower_vol >= hole_vol);

    double current = (lower + upper) / 2;
    double current_vol = volume_after_offset(grid, current, Ve, Tri) / initial_vol;
    while (abs(current_vol - hole_vol) > 1e-3 && abs(upper - lower) > 1e-8)
    {
        cout << "current radius interval: " << lower << " " << current << " " << upper << "\n";
        cout << "current volume interval: " << lower_vol << " " << current_vol << " " << upper_vol << "\n";
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
    }

    // merge with original mesh
    if (compute_surface_mesh_volume(V, F) * compute_surface_mesh_volume(Ve, Tri) >= 0)
    {
        std::cout << "Wrong orientation of eroded mesh!\n";
    }
    
    for (auto &f : F)
        for (int i = 0; i < 3; i++)
            f(i) += Ve.size();
    Ve.insert(Ve.end(), V.begin(), V.end());
    Tri.insert(Tri.end(), F.begin(), F.end());

    // save polygon mesh
    write_mesh(output, Ve, Tri, std::vector<Vec4I>());
    
    // igl::writeOBJ(output, V, F);
}

void usage(int status, const po::options_description &visible_opts) {
    cerr << "Usage: ./erode_cli input_mesh output_mesh volume_ratio" << endl;
    cerr << "eg: ./erode_cli mesh.obj out.obj 0.4" << endl;
    cout << visible_opts << endl;
    exit(status);
}

po::variables_map parseCmdLine(int argc, char *argv[]) {
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("input" ,  po::value<string>(), "input mesh file")
        ("output",  po::value<string>(), "output mesh file")
        ("ratio" ,  po::value<double>(), "solid volume")
        ;
    po::positional_options_description p;
    p.add("input", 1);
    p.add("output", 1);
    p.add("ratio", 1);

    // Options visible in the help message.
    po::options_description visible_opts;
    visible_opts.add_options()("help,h", "Produce this help message");

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cerr << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;
    {
        auto ratio = vm["ratio"].as<double>();
        if (ratio <= 0 || ratio > 1) {
            cerr << "Error: invalid ratio " << ratio << endl;
            fail = true;
        }
    }

    if (fail)
        usage(1, visible_opts);

    return vm;
}

int main(int argc, char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    auto &config = IsosurfaceInflatorConfig::get();

    execute(args["input"].as<string>(), args["output"].as<string>(), args["ratio"].as<double>());
}
