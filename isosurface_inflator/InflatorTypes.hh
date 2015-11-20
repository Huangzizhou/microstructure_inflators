#ifndef INFLATORTYPES_HH
#define INFLATORTYPES_HH

#include <Eigen/Dense>
#include <vector>
#include <Fields.hh>
#include <Geometry.hh>

// Inflator needs templated vector types for automatic differentiation
template<typename Real> using Point3  = Eigen::Matrix<Real, 3, 1>;
template<typename Real> using Point2  = Eigen::Matrix<Real, 2, 1>;
template<typename Real> using Vector3 = Eigen::Matrix<Real, 3, 1>;
template<typename Real> using Vector2 = Eigen::Matrix<Real, 2, 1>;
typedef Point3<double>  Point3d;
typedef Point3< float>  Point3f;
typedef Vector3<double> Vector3d;
typedef Vector3< float> Vector3f;
typedef Vector2<double> Vector2d;
typedef Vector2< float> Vector2f;

#endif /* end of include guard: INFLATORTYPES_HH */
