#include "VDBTools.hh"
#include <queue>

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

void volmesh2surface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, std::vector<Vec3s> &v, std::vector<Vec3I> &f)
{
    Eigen::MatrixXi F;
    igl::boundary_facets(T, F);

    Eigen::MatrixXi subF;
    Eigen::MatrixXd subV;
    Eigen::VectorXi tmpI;
    igl::remove_unreferenced(V, F, subV, subF, tmpI);

    eigen_to_openvdb(subV, v);
    eigen_to_openvdb(subF, f);
}

FloatGrid::Ptr volmesh2sdf(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const double voxel)
{
    std::vector<Vec3s> subV_openvdb;
    std::vector<Vec3I> subF_openvdb;
    volmesh2surface(V, T, subV_openvdb, subF_openvdb);

    math::Transform::Ptr xform = math::Transform::createLinearTransform(voxel);
    FloatGrid::Ptr grid = tools::meshToLevelSet<FloatGrid>(*xform, subV_openvdb, subF_openvdb, 3);

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

void sdf2mesh(const FloatGrid::Ptr &grid, std::vector<Vec3s> &V, std::vector<Vec3I> &F)
{
    std::vector<Vec4I> Quad;
    tools::volumeToMesh(*grid, V, F, Quad, 0, 0, true);
    clean_quads(F, Quad);
}

void write_mesh(const std::string &output, const std::vector<Vec3s> &V, const std::vector<Vec3I> &Tri)
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

    // for(int i = 0;i<(int)Quad.size();++i)
    // {
    //     fprintf(obj_file,"f");
    //     for(int j = 0; j<4;++j)
    //     {
    //     // OBJ is 1-indexed
    //     fprintf(obj_file," %u",Quad[i](j)+1);
    //     }
    //     fprintf(obj_file,"\n");
    // }

    fclose(obj_file);
}

double signed_volume_of_triangle(const Vec3s &p1, const Vec3s &p2, const Vec3s &p3) {
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
        vol += signed_volume_of_triangle(V[f(0)] - center, V[f(1)] - center, V[f(2)] - center);

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
    sdf2mesh(gridPtrCast<FloatGrid>(tmp_grid), Ve, Tri);

    double vol = abs(compute_surface_mesh_volume(Ve, Tri));
        
    return vol;
}
