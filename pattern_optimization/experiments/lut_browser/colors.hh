////////////////////////////////////////////////////////////////////////////////
// colors.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Some inline color utilities.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/06/2010 14:38:33
////////////////////////////////////////////////////////////////////////////////
#ifndef COLORS_HH
#define COLORS_HH

#include <cstdlib>
#include <cassert>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
// Forward declarations
////////////////////////////////////////////////////////////////////////////////
template<typename Real>
struct RGBColor;
template<typename Real>
struct HSVColor;

template<typename Real>
Real clamp(Real val, Real min_val = 0.0, Real max_val = 1.0)
{
    return std::max<Real>((Real) min_val, std::min<Real>((Real) max_val, val));
}

////////////////////////////////////////////////////////////////////////////////
/* Implements a color in the HSV space
// @tparam  Real    floating point type
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
struct HSVColor {
    /** HSVA components (all should stay in range [0, 1]) */
    Real hsva[4];
    Real &h, &s, &v, &a;

    HSVColor()
        : h(hsva[0]), s(hsva[1]), v(hsva[2]), a(hsva[3])
    {
        set(0, 0, 0);
    }

    HSVColor(Real h, Real s, Real v, Real a = 1.0)
        : h(hsva[0]), s(hsva[1]), v(hsva[2]), a(hsva[3])
    {
        set(h, s, v, a);
    }

    HSVColor(const HSVColor &c)
        : h(hsva[0]), s(hsva[1]), v(hsva[2]), a(hsva[3])
    {
        set(c.h, c.s, c.v, c.a);
    }

    void set(Real h_new, Real s_new, Real v_new, Real a_new = 1.0)
    {
        hsva[0] = h_new; hsva[1] = s_new; hsva[2] = v_new; hsva[3] = a_new;
    }

    void clamp()
    {
        for (int i = 0; i < 4; ++i)
            hsva[i] = ::clamp(hsva[i]);
    }

    HSVColor &operator=(const HSVColor &c)
    {
        set(c.h, c.s, c.v, c.a);
        return *this;
    }

    operator const Real *() const {
        return hsva;
    }

    operator Real *()  {
        return hsva;
    }

    operator RGBColor<Real>() {
        return hsvToRGB(*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
/* Implements a color in the RGB space
// @tparam  Real    floating point type
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
struct RGBColor {
    /** RGBA components (all should stay in range [0, 1]) */
    Real rgba[4];

    Real &r, &g, &b, &a;

    RGBColor()
        : r(rgba[0]), g(rgba[1]), b(rgba[2]), a(rgba[3])
    {
        set(0, 0, 0);
    }

    RGBColor(const RGBColor &c)
        : r(rgba[0]), g(rgba[1]), b(rgba[2]), a(rgba[3])
    {
        set(c.r, c.g, c.b);
    }

    RGBColor(Real r, Real g, Real b, Real a = 1.0)
        : r(rgba[0]), g(rgba[1]), b(rgba[2]), a(rgba[3])
    {
        set(r, g, b, a);
    }

    void set(Real r_new, Real g_new, Real b_new, Real a_new = 1.0)
    {
        rgba[0] = r_new; rgba[1] = g_new; rgba[2] = b_new; rgba[3] = a_new;
    }

    void clamp()
    {
        for (int i = 0; i < 4; ++i)
            rgba[i] = ::clamp(rgba[i]);
    }

    RGBColor &operator=(const RGBColor &c)
    {
        set(c.r, c.g, c.b, c.a);
        return *this;
    }

    operator const Real *() const {
        return rgba;
    }

    operator Real *() {
        return rgba;
    }
};

////////////////////////////////////////////////////////////////////////////////
/* Implements a linear color gradient (from a to b)
// @tparam  Color    color type
*///////////////////////////////////////////////////////////////////////////////
template<typename Color>
class ColorGradient
{
private:
    Color   m_start, m_end;
    int     m_varComponent;
public:
    ColorGradient(const Color &startColor, const Color &endColor)
        : m_start(startColor), m_end(endColor)
    { }

    Color operator()(float s) const
    {
        s = clamp(s);
        Color result;
        for (int i = 0; i < 4; ++i)
            result[i] = (1 - s) * m_start[i]  + s * m_end[i];
        return result;
    }

    void setAlpha(float alpha)
    {
        m_start[3] = m_end[3] = alpha;
    }

};

////////////////////////////////////////////////////////////////////////////////
/*! Convert a color from the hsv colorspace (all components in [0, 1]), to the
//  rgb colorspace (all components in [0, 1])
//  @param[in]  h       hue
//  @param[in]  s       saturation
//  @param[in]  v       value
//  @param[in]  a       alpha
//  @return     color in RGB space
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
RGBColor<Real> hsvToRGB(Real h, Real s, Real v, Real a = 1.0f)
{
    h *= 6;

    int i = int(floor(h));
    Real f = (i & 1) ? h - i : 1.0 - (h - i);
    Real m = v * (1.0 - s);
    Real n = v * (1.0 - s * f);

    Real r = v;
    Real g = n;
    Real b = m;

    r = (i == 2 || i == 3) ? m : ((i == 1 || i == 4) ? n : r);
    g = (i == 1 || i == 2) ? v : ((i == 4 || i == 5) ? m : g);
    b = (i == 2 || i == 5) ? n : ((i == 3 || i == 4) ? v : b);

    return RGBColor<Real>(r, g, b, a);
}

////////////////////////////////////////////////////////////////////////////////
/*! Convert a color from the hsv colorspace (all components in [0, 1]), to the
//  rgb colorspace (all components in [0, 1])
//  @param[in]  hsv     the color in hsv space
//  @return     color in RGB space
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
RGBColor<Real> hsvToRGB(const HSVColor<Real> &hsv)
{
    return hsvToRGB(hsv.h, hsv.s, hsv.v, hsv.a);
}

////////////////////////////////////////////////////////////////////////////////
/*! Generates a random RGB color.
//  @return     random RGB color object
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
RGBColor<Real> RandomRGB()
{
    return RGBColor<Real>(rand() / ((Real) RAND_MAX)
                , rand() / ((Real) RAND_MAX), rand() / ((Real) RAND_MAX));
}

////////////////////////////////////////////////////////////////////////////////
// Common template instantiations
////////////////////////////////////////////////////////////////////////////////
typedef HSVColor<float> HSVColorf;
typedef RGBColor<float> RGBColorf;

#endif  // COLORS_HH
