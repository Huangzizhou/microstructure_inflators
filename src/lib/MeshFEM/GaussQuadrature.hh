////////////////////////////////////////////////////////////////////////////////
// GaussQuadrature.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Gaussian quadrature rules for edges, triangles, and tetrahedra for
//      degrees up to 4.
//
//      These routines work both on functions with K + 1 Real parameters (where
//      K + 1 is the number of nodes of the K simplex) and functions with a
//      single EvalPt parameter.
//
//      SFINAE is used to "overload" the integration routines to work in both of
//      these cases.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/10/2014 17:13:25
////////////////////////////////////////////////////////////////////////////////
#ifndef GAUSSQUADRATURE_HH
#define GAUSSQUADRATURE_HH
#include <MeshFEM/Types.hh>
#include <MeshFEM/Functions.hh>
#include <MeshFEM/function_traits.hh>

// Edge function (1D)
// 1 point quadrature for const and linear, 2 point for quadratic and cubic, 3 for quartic and quintic
template<size_t _Deg, typename F, typename std::enable_if<(function_traits<F>::arity == 2) && (_Deg <= 13), int>::type = 0>
typename function_traits<F>::result_type integrate_edge(const F &f, Real vol = 1.0) {
    if (_Deg <= 1) { return vol * f(0.5, 0.5); }
    if ((_Deg == 2) || (_Deg == 3)) {
        constexpr double c0 = 0.78867513459481288225; // (3 + sqrt(3)) / 6
        constexpr double c1 = 0.21132486540518711775; // (3 - sqrt(3)) / 6
        typename function_traits<F>::result_type result(f(c0, c1));
        result += f(c1, c0);
        result *= vol / 2.0;
        return result;
    }
    if ((_Deg == 4) || (_Deg == 5)) {
        constexpr double c0 = 0.11270166537925831148; // (1 - sqrt(3/5)) / 2
        constexpr double c1 = 0.88729833462074168852; // (1 + sqrt(3/5)) / 2
        typename function_traits<F>::result_type result(f(c0, c1));
        result += f(c1, c0);
        result *= 5.0 / 18.0;
        result += (4.0 / 9.0) * f(0.5, 0.5);
        result *= vol;
        return result;
    }
    if ((_Deg == 6) || (_Deg == 7)) {
        constexpr double c0 = 0.33000947820757186760; // ( 1 - sqrt(3/7 - 2/7 sqrt(6/5)) ) / 2
        constexpr double c1 = 0.66999052179242813240; // ( 1 + sqrt(3/7 - 2/7 sqrt(6/5)) ) / 2

        constexpr double c2 = 0.06943184420297371239; // ( 1 - sqrt(3/7 + 2/7 sqrt(6/5)) ) / 2
        constexpr double c3 = 0.93056815579702628761; // ( 1 + sqrt(3/7 + 2/7 sqrt(6/5)) ) / 2

        constexpr double w0 = 0.32607257743127307131; // ( 18 + sqrt(30) ) / 72
        constexpr double w1 = 0.17392742256872692869; // ( 18 - sqrt(30) ) / 72

        typename function_traits<F>::result_type result(f(c0, c1));
        result += f(c1, c0);
        result *= w0;
        result += w1 * (f(c2, c3) + f(c3,c2)) ;
        result *= vol;
        return result;
    }
    if ((_Deg == 8) || (_Deg == 9)) {
        constexpr double c0 = 0.5;                    // 0.5

        constexpr double c1 = 0.23076534494715845448; // ( 1 - 1/3 sqrt(5 - 2 sqrt(10/7)) ) / 2
        constexpr double c2 = 0.76923465505284154552; // ( 1 + 1/3 sqrt(5 - 2 sqrt(10/7)) ) / 2

        constexpr double c3 = 0.04691007703066800360; // ( 1 - 1/3 sqrt(5 + 2 sqrt(10/7)) ) / 2
        constexpr double c4 = 0.95308992296933199640; // ( 1 + 1/3 sqrt(5 + 2 sqrt(10/7)) ) / 2

        constexpr double w0 = 0.28444444444444444444; // 128 / 450
        constexpr double w1 = 0.23931433524968323402; // (322 + 13 sqrt(70)) / 1800
        constexpr double w2 = 0.11846344252809454376; // (322 - 13 sqrt(70)) / 1800

        typename function_traits<F>::result_type result(f(c0, c0));
        result *= w0;

        result += w1 * (f(c1, c2) + f(c2,c1)) ;

        result += w2 * (f(c3, c4) + f(c4,c3)) ;

        result *= vol;
        return result;
    }
    if ((_Deg == 10) || (_Deg == 11)) {
        constexpr double c0 = 0.03376524289842398609385;  // (1 - 0.9324695142031520278123) / 2
        constexpr double c1 = 0.96623475710157601390615;  // (1 + 0.9324695142031520278123) / 2

        constexpr double c2 = 0.1693953067668677431695; // (1 - 0.661209386466264513661) / 2
        constexpr double c3 = 0.8306046932331322568305; // (1 + 0.661209386466264513661) / 2

        constexpr double c4 = 0.3806904069584015456845; // (1 - 0.238619186083196908631) / 2
        constexpr double c5 = 0.6193095930415984543155; // (1 + 0.238619186083196908631) / 2

        constexpr double w0 = 0.08566224618958517252015; // 0.1713244923791703450403 / 2
        constexpr double w1 = 0.1803807865240693037849; // 0.3607615730481386075698 / 2
        constexpr double w2 = 0.23395696728634552369495; // 0.4679139345726910473899 / 2

        typename function_traits<F>::result_type result(f(c0, c1));
        result += f(c1, c0);
        result *= w0;

        result += w1 * (f(c2, c3) + f(c3,c2)) ;

        result += w2 * (f(c4, c5) + f(c5,c4)) ;

        result *= vol;
        return result;
    }
    if ((_Deg == 12) || (_Deg == 13)) {
        constexpr double c0 = 0.5;                    // 0.5

        constexpr double c1 = 0.0254460438286207377369; // ( 1 - 0.9491079123427585245262 ) / 2
        constexpr double c2 = 0.9745539561713792622631; // ( 1 + 0.9491079123427585245262 ) / 2

        constexpr double c3 = 0.12923440720030278006805; // ( 1 - 0.7415311855993944398639 ) / 2
        constexpr double c4 = 0.87076559279969721993195; // ( 1 + 0.7415311855993944398639 ) / 2

        constexpr double c5 = 0.2970774243113014165467; // ( 1 - 0.4058451513773971669066 ) / 2
        constexpr double c6 = 0.7029225756886985834533; // ( 1 + 0.4058451513773971669066 ) / 2

        constexpr double w0 = 0.20897959183673469388; // 0.417959183673469387755 / 2
        constexpr double w1 = 0.06474248308443484663; // 0.1294849661688696932706 / 2
        constexpr double w2 = 0.13985269574463833395; // 0.2797053914892766679015 / 2
        constexpr double w3 = 0.19091502525255947248; // 0.38183005050511894495 / 2

        typename function_traits<F>::result_type result(f(c0, c0));
        result *= w0;

        result += w1 * (f(c1, c2) + f(c2,c1)) ;

        result += w2 * (f(c3, c4) + f(c4,c3)) ;

        result += w3 * (f(c5, c6) + f(c6,c5)) ;

        result *= vol;

        return result;
    }

    assert(false);
}

template<size_t _Deg, typename F, typename std::enable_if<function_traits<F>::arity == 1, int>::type = 0>
typename function_traits<F>::result_type integrate_edge(const F &f, Real vol = 1.0) {
    return integrate_edge<_Deg>([&](Real p0, Real p1) { return f(std::make_tuple(p0, p1)); }, vol); }

// Triangle function (2D)
// 1 point quadrature for const and linear, 3 for quadratic, 4 for cubic, and 6
// for quartic
// For efficiency, a negative weight rule is used for cubic
// integrals (the nonnegative weight rule would use 6 points)
// This means that the rule should not be used for stiffness matrix construction
// to avoid ruining positive semidefiniteness (This is only a problem for FEM
// degree 3+, which is not currently implemented)
template<size_t _Deg, typename F, typename std::enable_if<(function_traits<F>::arity == 3) && (_Deg <= 13), int>::type = 0>
typename function_traits<F>::result_type integrate_tri(const F &f, Real vol = 1.0) {
    if (_Deg <= 1) { return vol * f(1 / 3.0, 1 / 3.0, 1 / 3.0); }
    if (_Deg == 2) {
        // More accurate than the simpler midpoint rule
        constexpr double c0 = 2 / 3.0;
        constexpr double c1 = 1 / 6.0;
        typename function_traits<F>::result_type result(f(c0, c1, c1));
        result += f(c1, c0, c1);
        result += f(c1, c1, c0);
        result *= vol / 3.0;
        return result;
    }
    if (_Deg == 3) {
        constexpr double c0 = 3 / 5.0;
        constexpr double c1 = 1 / 5.0;
        typename function_traits<F>::result_type result(f(c0, c1, c1));
        result += f(c1, c0, c1);
        result += f(c1, c1, c0);
        result *= (25.0 / 48);
        result += (-9.0 / 16) * f(1 / 3.0, 1 / 3.0, 1 / 3.0); // NEGATIVE WEIGHT
        result *= vol;
        return result;
    }
    if (_Deg == 4) {
        // The analytic expressions of these weights are complicated...
        // See Derivations/TriangleGaussFelippa.nb
        // (From the Mathematica code in:
        // http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf )
        constexpr double w_0 =  0.22338158967801146570;
        constexpr double c0_0 = 0.10810301816807022736;
        constexpr double c1_0 = 0.44594849091596488632;
        typename function_traits<F>::result_type tmp(f(c0_0, c1_0, c1_0));
        tmp += f(c1_0, c0_0, c1_0);
        tmp += f(c1_0, c1_0, c0_0);
        tmp *= w_0;

        constexpr double w_1 =  0.10995174365532186764;
        constexpr double c0_1 = 0.81684757298045851308;
        constexpr double c1_1 = 0.09157621350977074346;
        typename function_traits<F>::result_type result(f(c0_1, c1_1, c1_1));
        result += f(c1_1, c0_1, c1_1);
        result += f(c1_1, c1_1, c0_1);
        result *= w_1;

        result += tmp;
        result *= vol;
        return result;
    }
    if (_Deg == 5) {
        // The analytic expressions of these weights are complicated...
        // See Derivations/TriangleGaussFelippa.nb
        // (From the Mathematica code in:
        // http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf )
        constexpr double w_0 =  0.12593918054482715260;
        constexpr double c0_0 = 0.79742698535308732240;
        constexpr double c1_0 = 0.10128650732345633880;
        typename function_traits<F>::result_type tmp(f(c0_0, c1_0, c1_0));
        tmp += f(c1_0, c0_0, c1_0);
        tmp += f(c1_0, c1_0, c0_0);
        tmp *= w_0;

        constexpr double w_1 =  0.13239415278850618074;
        constexpr double c0_1 = 0.059715871789769820459;
        constexpr double c1_1 = 0.47014206410511508977;
        typename function_traits<F>::result_type result(f(c0_1, c1_1, c1_1));
        result += f(c1_1, c0_1, c1_1);
        result += f(c1_1, c1_1, c0_1);
        result *= w_1;

        result += tmp;
        result += (9.0 / 40) * f(1.0 / 3, 1.0 / 3, 1.0 / 3);

        result *= vol;
        return result;
    }
    if (_Deg >= 6 && _Deg <= 13 ) {
        throw std::runtime_error("Triangle quadrature for degree is not implemented!");
    }
    assert(false);
}
template<size_t _Deg, typename F, typename std::enable_if<function_traits<F>::arity == 1, int>::type = 0>
typename function_traits<F>::result_type integrate_tri(const F &f, Real vol = 1.0) {
    return integrate_tri<_Deg>([&](Real p0, Real p1, Real p2) { return f(std::make_tuple(p0, p1, p2)); }, vol);
}

// Tet function (3D)
// 1 point quadrature for const and linear, 4 point for quadratic, 5 for cubic,
// and 11 for quartic.
// For efficiency, negative weight rules are used for cubic and quartic
// integrals (the nonnegative weight rules use 8 and 16 points respectively).
// This means that those rules should not be used for stiffness matrix
// construction to avoid ruining positive semidefiniteness (This is only a
// problem for FEM degree 3+, which is not currently implemented)
template<size_t _Deg, typename F, typename std::enable_if<(function_traits<F>::arity == 4) && (_Deg <= 4), int>::type = 0>
typename function_traits<F>::result_type integrate_tet(const F &f, Real vol = 1.0) {
    if (_Deg <= 1) { return vol * f(1 / 4.0, 1 / 4.0, 1 / 4.0, 1 / 4.0); }
    if (_Deg == 2) {
        constexpr double c0 = 0.58541019662496845446; // (5 + 3 sqrt(5)) / 20
        constexpr double c1 = 0.13819660112501051518; // (5 -   sqrt(5)) / 20
        typename function_traits<F>::result_type result(f(c0, c1, c1, c1));
        result += f(c1, c0, c1, c1);
        result += f(c1, c1, c0, c1);
        result += f(c1, c1, c1, c0);
        result *= vol / 4;
        return result;
    }
    if (_Deg == 3) {
        // http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
        constexpr double c0 = 0.5;
        constexpr double c1 = 1 / 6.0;
        typename function_traits<F>::result_type result(f(c0, c1, c1, c1));
        result += f(c1, c0, c1, c1);
        result += f(c1, c1, c0, c1);
        result += f(c1, c1, c1, c0);
        result *= 0.45;
        result += (-0.8) * f(1 / 4.0, 1 / 4.0, 1 / 4.0, 1 / 4.0); // NEGATIVE WEIGHT
        result *= vol;
        return result;
    }
    if (_Deg == 4) {
        // This rule is from
        // http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
        // but the weights there are off by a factor of 6!
        typename function_traits<F>::result_type result(f(0.25, 0.25, 0.25, 0.25));
        result *= -148.0 / 1875.0; // NEGATIVE WEIGHT

        constexpr double c0_0 = 11.0 / 14.0;
        constexpr double c1_0 =  1.0 / 14.0;
        typename function_traits<F>::result_type tmp(f(c0_0, c1_0, c1_0, c1_0));
        tmp += f(c1_0, c0_0, c1_0, c1_0);
        tmp += f(c1_0, c1_0, c0_0, c1_0);
        tmp += f(c1_0, c1_0, c1_0, c0_0);
        tmp *= 343.0 / 7500.0;
        result += tmp;

        constexpr double c0_1 = 0.39940357616679920500; // (14 + sqrt(70)) / 56
        constexpr double c1_1 = 0.10059642383320079500; // (14 - sqrt(70)) / 56
        tmp  = f(c0_1, c0_1, c1_1, c1_1);
        tmp += f(c0_1, c1_1, c0_1, c1_1);
        tmp += f(c0_1, c1_1, c1_1, c0_1);
        tmp += f(c1_1, c0_1, c0_1, c1_1);
        tmp += f(c1_1, c0_1, c1_1, c0_1);
        tmp += f(c1_1, c1_1, c0_1, c0_1);
        tmp *= 56.0 / 375.0;
        result += tmp;

        result *= vol;
        return result;
    }
    assert(false);
}
template<size_t _Deg, typename F, typename std::enable_if<function_traits<F>::arity == 1, int>::type = 0>
typename function_traits<F>::result_type integrate_tet(const F &f, Real vol = 1.0) {
    return integrate_tet<_Deg>([&](Real p0, Real p1, Real p2, Real p3) { return f(std::make_tuple(p0, p1, p2, p3)); }, vol);
}

// Integration on a _K simplex (runs the implementations above).
// Usage:
// Quadrature<Simplex::{Edge,Triangle,Tetrahedron}, Degree>::integrate(f);
template<size_t _K, size_t _Deg>
class Quadrature { };

template<size_t _Deg> class Quadrature<Simplex::Edge,        _Deg> { public: template<typename F> static auto integrate(const F &f, Real vol = 1.0) -> decltype(integrate_edge<_Deg>(f)) { return integrate_edge<_Deg>(f, vol); } };
template<size_t _Deg> class Quadrature<Simplex::Triangle,    _Deg> { public: template<typename F> static auto integrate(const F &f, Real vol = 1.0) -> decltype(integrate_tri <_Deg>(f)) { return integrate_tri< _Deg>(f, vol); } };
template<size_t _Deg> class Quadrature<Simplex::Tetrahedron, _Deg> { public: template<typename F> static auto integrate(const F &f, Real vol = 1.0) -> decltype(integrate_tet <_Deg>(f)) { return integrate_tet< _Deg>(f, vol); } };

#endif /* end of include guard: GAUSSQUADRATURE_HH */
