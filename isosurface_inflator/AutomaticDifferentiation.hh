////////////////////////////////////////////////////////////////////////////////
// AutomaticDifferentiation.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Includes and functions needed for automatic differentiation
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 15:26:43
////////////////////////////////////////////////////////////////////////////////
#ifndef AUTOMATICDIFFERENTIATION_HH
#define AUTOMATICDIFFERENTIATION_HH

// #include <adept.h>
#include <unsupported/Eigen/AutoDiff>
#include <algorithm>
#include <type_traits>

// using ADReal = adept::adouble;

template<typename T> struct IsAutoDiffType : public std::false_type { };
template<typename T> struct IsAutoDiffType<Eigen::AutoDiffScalar<T>> : public std::true_type { };
// template<>           struct IsAutoDiffType<          adept::adouble> : public std::true_type { };

// A note on Eigen's norm() vs squaredNorm():
// Adept's sqrt overload is not visible to Eigen, so we must use
// sqrt(x.squaredNorm()) (or include Adept before Eigen).

// Wrapper to get the underlying value of an autodiff type (or do nothing for
// primitive types)
template<typename T>
struct StripAutoDiffImpl {
    using result_type = T;
    static result_type run(const T &v) { return v; }
};

template<typename _DerType>
struct StripAutoDiffImpl<Eigen::AutoDiffScalar<_DerType>> {
    using result_type = typename Eigen::internal::traits<typename Eigen::internal::remove_all<_DerType>::type>::Scalar;
    static result_type run(const Eigen::AutoDiffScalar<_DerType> &v) { return v.value(); }
};

// Cast autodiff vectors/matrices to plain vectors/matrices.
template<typename _DerType, int... I>
struct StripAutoDiffImpl<Eigen::Matrix<Eigen::AutoDiffScalar<_DerType>, I...>> {
    using autodiff_type = Eigen::Matrix<Eigen::AutoDiffScalar<_DerType>, I...>;
    using scalar_type = typename Eigen::internal::traits<typename Eigen::internal::remove_all<_DerType>::type>::Scalar;
    using result_type = Eigen::Matrix<scalar_type, I...>;

    static result_type run(const autodiff_type &v) {
        result_type r(v.rows(), v.cols());
        for (size_t i = 0; i < v.rows(); ++i) {
            for (size_t j = 0; j < v.cols(); ++j)
                r(i, j) = v(i, j).value();
        }
        return r;
    }
};

// template<>           struct StripAutoDiffImpl<adept::adouble> { static double run(const adept::adouble &v) { return v.value(); } };

// Strip automatic differentiation wrapper from a scalar value type (does
// nothing when applied to a non-autodiff type).
template<typename T>
typename StripAutoDiffImpl<T>::result_type
stripAutoDiff(const T &val) {
    return StripAutoDiffImpl<T>::run(val);
}

template<typename T>
constexpr bool isAutodiffType() {
    return !std::is_same<typename StripAutoDiffImpl<T>::result_type, T>::value;
}

template<typename T>
bool isAutodiffType(const T &val) { return isAutodiffType<T>(); }

// For casting to non autodiff types, we must strip
template<bool IsAutodiffTarget>
struct AutodiffCastImpl {
    template<typename TNew, typename TOrig>
    static TNew run(const TOrig &val) { return TNew(stripAutoDiff(val)); }
};

// Casting to autodiff type just works
template<>
struct AutodiffCastImpl<true> {
    template<typename TNew, typename TOrig>
    static TNew run(const TOrig &val) { return TNew(val); }
};

template<typename TNew, typename TOrig>
TNew autodiffCast(const TOrig &orig) {
    return AutodiffCastImpl<isAutodiffType<TNew>()>::template run<TNew>(orig);
}

// Implement tanh for Eigen autodiff (argument dependent lookup)
namespace Eigen {
#define EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(FUNC,CODE) \
  template<typename DerType> \
  inline const Eigen::AutoDiffScalar<Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<typename Eigen::internal::traits<typename Eigen::internal::remove_all<DerType>::type>::Scalar>, const typename Eigen::internal::remove_all<DerType>::type> > \
  FUNC(const Eigen::AutoDiffScalar<DerType>& x) { \
    using namespace Eigen; \
    typedef typename Eigen::internal::traits<typename Eigen::internal::remove_all<DerType>::type>::Scalar Scalar; \
    typedef AutoDiffScalar<CwiseUnaryOp<Eigen::internal::scalar_multiple_op<Scalar>, const typename Eigen::internal::remove_all<DerType>::type> > ReturnType; \
    CODE; \
  }

    // Implement tanh for eigen
    EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(tanh,
      using std::tanh;
      using std::cosh;
      using numext::abs2;
      return ReturnType(tanh(x.value()),
                        x.derivatives() * (Scalar(1)/abs2(cosh(x.value()))));
    )

#undef EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY

}
#endif /* end of include guard: AUTOMATICDIFFERENTIATION_HH */
