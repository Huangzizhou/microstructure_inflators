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

// Wrapper to get the underlying value of an adept double (or do nothing for
// primitive types)
template<typename T> struct StripAutoDiffImpl                 { static double run(const              T &v) { return v;         } };
// template<>           struct StripAutoDiffImpl<adept::adouble> { static double run(const adept::adouble &v) { return v.value(); } };

template<typename T> struct StripAutoDiffImpl<Eigen::AutoDiffScalar<T>> {
    static double run(const Eigen::AutoDiffScalar<T> &v) { return v.value(); }
};

template<typename T> double stripAutoDiff(const T &val) { return StripAutoDiffImpl<T>::run(val); }

#endif /* end of include guard: AUTOMATICDIFFERENTIATION_HH */
