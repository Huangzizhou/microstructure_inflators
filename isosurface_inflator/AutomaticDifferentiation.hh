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
#include <cmath>
#include <algorithm>
#include <limits>
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

    // Implement tanh for eigen autodiff
    EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(tanh,
      using std::tanh;
      using std::cosh;
      using numext::abs2;
      return ReturnType(tanh(x.value()),
                        x.derivatives() * (Scalar(1)/abs2(cosh(x.value()))));
    )

    // Implement log(cosh(x)) with derivative; useful for stable exp_smin computation
    EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(log_cosh,
      return ReturnType(std::log(std::cosh(x.value())),
                        x.derivatives() * std::tanh(x.value()));
    )

#undef EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY
}

// Implement log_cosh for non-autodiff types.
template<typename T>
typename std::enable_if<!isAutodiffType<T>(), T>::type
log_cosh(const T val) {
    return log(cosh(val));
}

// std::numeric_limits is dangerous! If you use it on Eigen's autodiff types you
// will get undefined behavior.
template<typename T>
struct safe_numeric_limits : public std::numeric_limits<typename StripAutoDiffImpl<T>::result_type>
{
    using NonADType = typename StripAutoDiffImpl<T>::result_type;
    static_assert(std::is_arithmetic<NonADType>::value,
                  "std::numeric_limits is broken for non-arithmetic types!");
};

////////////////////////////////////////////////////////////////////////////////
// Derivative debugging
////////////////////////////////////////////////////////////////////////////////
// Check for Inf/NaN in derivative fields
template<typename T>
typename std::enable_if<isAutodiffType<T>(), bool>::type
hasInvalidDerivatives(const T &val) {
    const auto &der = val.derivatives();
    for (size_t i = 0; i < der.rows(); ++i)
        if (std::isnan(der[i]) || std::isinf(der[i])) return true;
    return false;
}

// Return false for non-autodiff types.
template<typename T>
typename std::enable_if<!isAutodiffType<T>(), bool>::type
hasInvalidDerivatives(const T &val) { return false; }

template<typename T>
typename std::enable_if<isAutodiffType<T>(), void>::type
reportDerivatives(std::ostream &os, const T &val) {
    auto prec = os.precision(5);
    const auto &der = val.derivatives();
    for (size_t i = 0; i < der.rows(); ++i)
        os << "\t" << der[i];
    os.precision(prec);
}

// do nothing for non-autodiff types.
template<typename T>
typename std::enable_if<!isAutodiffType<T>(), void>::type
reportDerivatives(std::ostream &os, const T &val) { }

#endif /* end of include guard: AUTOMATICDIFFERENTIATION_HH */
