////////////////////////////////////////////////////////////////////////////////
// data_viewer.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Fast OpenGL-based lookup table viewer for visualizing mapping from
//      pattern parameters to (E, nu)
//      Useful for exploring pattern "foldover"
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/29/2015 17:56:46
////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <set>
#include <functional>
#include <iomanip>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#include <colors.hh>
#include <Future.hh>
#include "../../LookupTable.hh"
#include "../../../isosurface_inflator/PatternSignedDistance.hh"
typedef IsotropicLookupTable<double> LUT;

// TODO: make this configurable
using WMesh = WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>;
using PSD = PatternSignedDistance<float, WMesh>;

#include <boost/filesystem/operations.hpp>

using namespace std;
// #include <Eigen/Dense>
// using Eigen::Vector2f;

// views' render methods assume that the view transforms have been set up so
// that they should render in [0, 1]^2
// Also, expects event x, y to be in [0, 1]^2
class View;
typedef shared_ptr<View> ViewPtr;

// Stores a 2D view coordinate system transform (to map event locations to
// subview coordinates)
template<typename Real>
class ViewTransform {
public:
    // Identity transform by default
    ViewTransform() : m_cf1([](Real x, Real y) -> Real { return x; }),
                   m_cf2([](Real x, Real y) -> Real { return y; }) { }

    typedef std::function<Real(Real, Real)> CoordinateFunction2;
    template<typename F1, typename F2>
    ViewTransform(const F1 &f1, const F2 &f2) : m_cf1(f1), m_cf2(f2) { }

    float transformedX(float x, float y) const { return m_cf1(x, y); }
    float transformedY(float x, float y) const { return m_cf2(x, y); }
    void  apply (float x, float y, float &xout, float &yout) const { xout = transformedX(x, y); yout = transformedY(x, y); }
    void  apply (float &xinout, float &yinout) const { apply(xinout, yinout, xinout, yinout); }

    // Compose with b on the left (apply this transform after b).
    ViewTransform composeLeft(const ViewTransform &b) const {
        // WARNING: capturing member variables by value is impossible! A
        // by-value default capture really just captures the "this" pointer and
        // accesses the members through it!
        ViewTransform copy(*this);
        return ViewTransform(
            [=](Real x, Real y) -> Real { return copy.transformedX(b.transformedX(x, y), b.transformedY(x, y)); },
            [=](Real x, Real y) -> Real { return copy.transformedY(b.transformedX(x, y), b.transformedY(x, y)); });
    }

    // Compose with b on the right (apply this transform before b).
    ViewTransform composeRight(const ViewTransform &b) const {
        // WARNING: capturing member variables by value is impossible! A
        // by-value default capture really just captures the "this" pointer and
        // accesses the members through it!
        ViewTransform copy(*this);
        return ViewTransform(
            [=](Real x, Real y) -> Real { return b.transformedX(copy.transformedX(x, y), copy.transformedY(x, y)); },
            [=](Real x, Real y) -> Real { return b.transformedY(copy.transformedX(x, y), copy.transformedY(x, y)); });
    }
private:
    CoordinateFunction2 m_cf1, m_cf2;
};

typedef ViewTransform<float> ViewTransformf;

void drawString(const std::string &str) {
    for (char c : str)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, c);
}

// Empty, transparent view by default.
class View : public std::enable_shared_from_this<View> {
public:
    virtual void render() const { };

    virtual ViewPtr viewForPoint(float x, float y) { return shared_from_this(); }
    // Gets the transformation mapping coordinates of a point in the
    // top-level-view's system to coordinates in the system of the subview at (x, y)
    virtual ViewTransformf viewTransformForPoint(float x, float y) { return ViewTransformf(); }

    virtual bool mousePress(int button, float x, float y)   { return false; }
    virtual bool mouseRelease(int button, float x, float y) { return false; }
    virtual bool mouseDrag(float x, float y)                { return false; }
    virtual bool keyPressed(unsigned char c)                { return false; }

    typedef std::function<void(std::shared_ptr<View>)> CallBack;
};

// A range that supports auto-min and auto-max
// The interval is [min, max)
struct DataRange {
    // Empty range supporting union
    DataRange() : m_min(std::numeric_limits<float>::infinity()), m_max(-std::numeric_limits<float>::infinity()) { }
    DataRange(const std::vector<float> &data) { setAutoMin(data); setAutoMax(data); }
    DataRange(float min, float max) : m_min(min), m_max(max) { }

    // With a padding so the extreme values are included in the range.
    void setAutoMin(const std::vector<float> &data, float padding=1e-6) { m_min = (data.size() == 0) ? 0.0 : *std::min_element(data.begin(), data.end()) - padding; }
    void setAutoMax(const std::vector<float> &data, float padding=1e-6) { m_max = (data.size() == 0) ? 0.0 : *std::max_element(data.begin(), data.end()) + padding; }
    void setMin(float m) { m_min = m; }
    void setMax(float m) { m_max = m; }

    float min() const { return m_min; }
    float max() const { return m_max; }
    float size() const { return max() - min(); }

    // n: number of subintervals
    // i: which subinterval (in {0..n-1})
    DataRange subrange(size_t n, size_t i) const {
        float fullMin = min(), fullSize = size();
        float subsize = fullSize / n;
        return DataRange(fullMin + subsize * i, fullMin + subsize * (i + 1));
    }

    // Query if in [min, max)
    bool inRange(float f) const { return (f >= min()) && (f < max()); }

    DataRange unionRange(const DataRange &r) const {
        return DataRange (std::min(min(), r.min()),
                          std::max(max(), r.max()));
    }

    DataRange intersectRange(const DataRange &r) const {
        return DataRange (std::max(min(), r.min()),
                          std::min(max(), r.max()));
    }

    friend std::ostream &operator<<(std::ostream &os, const DataRange &dr) {
        os << "[" << dr.min() << ", " << dr.max() << ")";
        return os;
    }

private:
    float m_min, m_max;
};

struct DataBins {
    DataBins() { } // No bins
    DataBins(const std::vector<float> &data, const DataRange &dr, size_t numBins = 50) { m_computeBins(data, dr, numBins); }
    DataBins(const std::vector<float> &data, size_t numBins = 50) { m_computeBins(data, DataRange(data), numBins); }
    size_t numBins() const { return m_bins.size(); }

    size_t count(size_t i) const { return m_bins.at(i); }
    size_t maxCount() const {
        if (numBins() == 0) return 0;
        return *std::max_element(m_bins.begin(), m_bins.end());
    }

private:
    void m_computeBins(const std::vector<float> &data, const DataRange &dr, size_t numBins) {
        m_bins.assign(numBins, 0);
        float range = dr.size();
        for (auto f : data) {
            int b = floor(numBins * ((f - dr.min()) / range));
            if ((b >= 0) && (size_t(b) < numBins))
                ++m_bins[b];
        }
    }
    std::vector<size_t> m_bins;
};

class Plot : public View {
public:
    typedef ColorMap<RGBColorf, float> CMap;
    typedef std::shared_ptr<CMap> CMapPointer;

    void setForegroundColor(const RGBColorf &c) { m_fgcolor = c; }
    void setBackgroundColor(const RGBColorf &c) { m_bgcolor = c; }
    void setSelectionColor (const RGBColorf &c) { m_slcolor = c; }
    void setForegroundColor(float r, float g, float b, float a = 1.0) { setForegroundColor(RGBColorf(r, g, b, a)); }
    void setBackgroundColor(float r, float g, float b, float a = 1.0) { setBackgroundColor(RGBColorf(r, g, b, a)); }
    void setSelectionColor (float r, float g, float b, float a = 1.0) { setSelectionColor (RGBColorf(r, g, b, a)); }

    const RGBColorf &foregroundColor() const { return m_fgcolor; }
    const RGBColorf &backgroundColor() const { return m_bgcolor; }
    const RGBColorf &selectionColor () const { return m_slcolor; }

    // Set the colormap used to color each data series independently
    // (otherwise foregroundColor() is used for all series)
    void setSeriesColormap(CMapName name) { m_seriesCMap = make_shared<CMap>(name); }
    CMapPointer seriesColormap() { return m_seriesCMap; }

    // By default, plots support only one series.
    virtual size_t numSeries() const { return 1; }

    RGBColorf seriesColor(size_t s) const {
        if (m_seriesCMap) return (*m_seriesCMap)(float(s) / float(numSeries()));
        else              return foregroundColor();
    }

    void setPointSize(float ps) { m_pointSize = ps; }

    // Renders background (Axes in future?)
    virtual void render() const {
        glColor4fv(m_bgcolor);
        glBegin(GL_QUADS); {
            glVertex3f(0, 0, 0);
            glVertex3f(1, 0, 0);
            glVertex3f(1, 1, 0);
            glVertex3f(0, 1, 0);
        }
        glEnd();
    }

protected:
    RGBColorf m_fgcolor, m_bgcolor, m_slcolor;
    CMapPointer m_seriesCMap;
    float m_pointSize = 1.0;
};


class HistogramPlot : public Plot {
public:
    HistogramPlot() { }
    HistogramPlot(const std::vector<float> &data) { setData(data); }
    HistogramPlot(const std::vector<float> &data, size_t numBins) { setData(data, numBins); }
    HistogramPlot(const std::vector<float> &data, size_t numBins, const DataRange &dr) { setData(data, dr, numBins); }
    HistogramPlot(const std::vector<float> &data, const DataRange &dr) { setData(data, dr); }

    virtual void render() const {
        Plot::render();
        float maxCount = m_bins.maxCount();
        glColor4fv(foregroundColor());
        glBegin(GL_QUADS); {
            float width = 1.0 / m_bins.numBins();
            float xoffset = 0;
            float yoffset = 0;
            for (size_t i = 0; i < m_bins.numBins(); ++i) {
                float height = m_bins.count(i) / maxCount;
                glVertex3f(xoffset        , yoffset, 0);
                glVertex3f(xoffset + width, yoffset, 0);
                glVertex3f(xoffset + width, yoffset + height, 0);
                glVertex3f(xoffset        , yoffset + height, 0);
                xoffset += width;
            }
        }
        glEnd();
    }

    void setData(const std::vector<float> &data, size_t numBins = 50) {
        setData(data, DataRange(data), numBins);
    }
    void setData(const std::vector<float> &data, const DataRange &dr, size_t numBins = 50) {
        m_bins = DataBins(data, dr, numBins);
        m_dataRange = dr;
    }

    const DataRange &dataRange() const { return m_dataRange; }

    virtual bool mousePress  (int button, float x, float y) { return false; }
    virtual bool mouseRelease(int button, float x, float y) { return false; }

private:
    DataBins m_bins;
    DataRange m_dataRange;
};

class LogscaleYPlot : public Plot {
public:
    LogscaleYPlot() { }
    void setData(const std::vector<std::vector<float>> &series_xdata,
                 const std::vector<std::vector<float>> &series_ydata) {
        if (series_xdata.size() != series_ydata.size())
            throw std::runtime_error("Unmatched data sizes");
        for (size_t s = 0; s < series_xdata.size(); ++s) {
            if (series_xdata[s].size() != series_ydata[s].size())
                throw std::runtime_error("Unmatched data sizes");
        }
        m_series_xdata = series_xdata;
        m_series_ydata = series_ydata;
        m_updateViewCoordCache();
        updateSelection();
    }

    void setData(const std::vector<float> &xdata, const std::vector<float> &ydata) {
        setData(std::vector<std::vector<float>>(1, xdata), std::vector<std::vector<float>>(1, ydata));
    }

    void setPlotRange(DataRange xRange, DataRange yRange) {
        m_xRange = xRange;
        m_yRange = yRange;
        m_updateViewCoordCache();
    }

    void autoPlotRange() {
        DataRange xrange, yrange;
        for (size_t s = 0; s < numSeries(); ++s) {
            xrange = xrange.unionRange(DataRange(m_series_xdata[s]));
            yrange = yrange.unionRange(DataRange(m_series_ydata[s]));
        }
        setPlotRange(xrange, yrange);
    }

    DataRange xRange() const { return m_xRange; }
    DataRange yRange() const { return m_yRange; }

    virtual void render() const {
        Plot::render();
        glPointSize(m_pointSize);
        size_t index = 0;
        // TODO: this is now a bottleneck--should be sped up by using a VBO.
        glBegin(GL_POINTS); {
            for (size_t s = 0; s < numSeries(); ++s) {
                // Get data points (already transformed into view coordinates!)
                const auto &xdata = m_series_viewx_cache[s];
                const auto &ydata = m_series_viewy_cache[s];
                auto fgcolor = seriesColor(s);
                glColor4fv(fgcolor);
                for (size_t i = 0; i < xdata.size(); ++i) {
                    float x = xdata[i], y = ydata[i];
                    if (m_pointSelected.at(index)) glColor3fv(selectionColor());
                    glVertex3f(x, y, 0.0);
                    if (m_pointSelected.at(index)) glColor3fv(fgcolor);
                    ++index;
                }
            }
        }
        glEnd();

        size_t selSize = selectionSize();
        if (selSize > 0) {
            glColor3fv(selectionColor());
            glRasterPos2i(0, 0);
            float sx, sy;
            m_viewToPlot(m_lastSelectCenterX, m_lastSelectCenterY, sx, sy);
            drawString(std::to_string(selSize) + " points selected around ("
                    + std::to_string(sx) + ", " + std::to_string(sy) + ")");
        }
    }

    virtual size_t numSeries() const { return m_series_xdata.size(); }
    size_t numPoints() const {
        size_t count = 0;
        for (size_t s = 0; s < numSeries(); ++s) count += numPointsInSeries(s);
        return count;
    }
    size_t numPointsInSeries(size_t s) const { return m_series_xdata.at(s).size(); }

    void setSelectionRadius(float radius) {
        m_selectionRadius = radius;
        m_selectPoint(m_lastSelectCenterX, m_lastSelectCenterY);
    }
    float selectionRadius() const { return m_selectionRadius; }

    // The aspect ratio of the selection ellipse (width / height)
    void setSelectionAspect(float aspect) { m_selectionAspect = aspect; }

    void resetSelection() { m_pointSelected.assign(numPoints(), false); }
    void updateSelection() { m_selectPoint(m_lastSelectCenterX, m_lastSelectCenterY); }

    size_t selectionSize() const {
        size_t count = 0;
        for (bool b : m_pointSelected) if (b) ++count;
        return count;
    }

    // Index in concatenated series data!
    std::vector<size_t> selectedPointIndices() const {
        std::vector<size_t> result;
        result.reserve(selectionSize());
        for (size_t i = 0; i < numPoints(); ++i) if (m_pointSelected[i]) result.push_back(i);
        return result;
    }

    std::pair<size_t, size_t> getSeriesAndOffset(size_t index) const {
        if (index >= numPoints()) throw std::runtime_error("Index out of bounds");
        size_t series = 0;
        while (index >= numPointsInSeries(series)) index -= numPointsInSeries(series++);
        return std::make_pair(series, index);
    }

    void setSelectionChangedCallback(const CallBack &cb) { m_selChangedCallback = cb; }

    virtual bool mousePress(  int button, float x, float y) { m_selectPoint(x, y); return true; }
    virtual bool mouseRelease(int button, float x, float y) { m_selectPoint(x, y); return true; }
    virtual bool mouseDrag(               float x, float y) { m_selectPoint(x, y); return true; }
    virtual bool keyPressed(unsigned char c) {
        bool handled = false;
        if (c == '+' || c == '=') {
            setSelectionRadius(selectionRadius() * 1.25);
            handled = true;
        }
        if (c == '-' || c == '_') {
            setSelectionRadius(selectionRadius() * 0.8);
            handled = true;
        }
        return handled;
    }

private:
    // Select the data points around a click location (in view coordinates)
    // The distance is measured in view space, taking aspect ratio into account.
    void m_selectPoint(float x, float y) {
        m_lastSelectCenterX = x, m_lastSelectCenterY = y;

        std::vector<bool> oldSelection(std::move(m_pointSelected));
        bool selectionChanged = false;
        if (oldSelection.size() != numPoints()) {
            selectionChanged = true;
            oldSelection.assign(numPoints(), false);
        }
        resetSelection();
        size_t index = 0;
        for (size_t s = 0; s < numSeries(); ++s) {
            for (size_t si = 0; si < numPointsInSeries(s); ++si) {
                float xd = m_series_viewx_cache[s][si], yd = m_series_viewy_cache[s][si];
                float distx = xd - x, disty = yd - y;
                distx /= m_selectionAspect;
                float dist = sqrt(distx * distx + disty * disty);
                if (dist < m_selectionRadius)
                    m_pointSelected[index] = true;
                selectionChanged  = selectionChanged || (m_pointSelected[index] != oldSelection[index]);
                ++index;
            }
        }
        if (selectionChanged && m_selChangedCallback)
            m_selChangedCallback(shared_from_this());
    }

    // in-place coordinate transforms
    void m_plotToView(float &inoutx, float &inouty) const {
        inoutx = (    inoutx  -     m_xRange.min() ) / (    m_xRange.max()  -     m_xRange.min() );
        inouty = (log(inouty) - log(m_yRange.min())) / (log(m_yRange.max()) - log(m_yRange.min()));
    }
    void m_viewToPlot(float &inoutx, float &inouty) const {
        inoutx =     inoutx * (    m_xRange.max()  -     m_xRange.min() ) +     m_xRange.min()  ;
        inouty = exp(inouty * (log(m_yRange.max()) - log(m_yRange.min())) + log(m_yRange.min()));
    }

    // out-of-place coordinate transforms
    void m_plotToView(float inx, float iny, float &outx, float &outy) const { outx = inx, outy = iny; m_plotToView(outx, outy); }
    void m_viewToPlot(float inx, float iny, float &outx, float &outy) const { outx = inx, outy = iny; m_viewToPlot(outx, outy); }

    void m_updateViewCoordCache() {
        m_series_viewx_cache.resize(numSeries());
        m_series_viewy_cache.resize(numSeries());
        for (size_t s = 0; s < numSeries(); ++s) {
            size_t np = numPointsInSeries(s);
            m_series_viewx_cache[s].assign(np, 0.0);
            m_series_viewy_cache[s].assign(np, 0.0);
            for (size_t p = 0; p < np; ++p) {
                m_plotToView(m_series_xdata[s][p], m_series_ydata[s][p],
                             m_series_viewx_cache[s][p], m_series_viewy_cache[s][p]);
            }
        }
    }

    // Data members
    float m_selectionRadius = 5.0 / 512, m_selectionAspect = 1.0;
    DataRange m_xRange, m_yRange;
    std::vector<bool> m_pointSelected;
    CallBack m_selChangedCallback;
    float m_lastSelectCenterX = std::numeric_limits<float>::infinity(),
          m_lastSelectCenterY = std::numeric_limits<float>::infinity();
    vector<vector<float>> m_series_xdata, m_series_ydata;
    // cached view coordinates of the data--must be updated with
    // m_updateViewCoordCache when the data or the axes change
    vector<vector<float>> m_series_viewx_cache, m_series_viewy_cache;
};

// TODO: make this configurable
unique_ptr<PSD> g_patternSignedDistance;
class PatternView : public View {
public:
    PatternView(size_t resolution) : m_resolution(resolution) { }

    virtual void render() const {
        glColor3f(0.8, 0.8, 0.8);
        glBegin(GL_QUADS); {
            glVertex3f(0, 0, 0);
            glVertex3f(1, 0, 0);
            glVertex3f(1, 1, 0);
            glVertex3f(0, 1, 0);
        }

        glColor3f(0.5, 0.06, 0.5);
        float width = 1.0 / m_resolution;
        float yoffset = 0;
        for (size_t r = 0; r < m_resolution; ++r) {
            float xoffset = 0;
            for (size_t c = 0; c < m_resolution; ++c) {
                // Map [0, 1] to [-1, 1]
                Point3<float> p(2 * (xoffset + 0.5 * width) - 1,
                                2 * (yoffset + 0.5 * width) - 1,
                                0);
                if (g_patternSignedDistance->isInside(p)) {
                    glVertex3f(xoffset        , yoffset, 0);
                    glVertex3f(xoffset + width, yoffset, 0);
                    glVertex3f(xoffset + width, yoffset + width, 0);
                    glVertex3f(xoffset        , yoffset + width, 0);
                }
                xoffset += width;
            }
            yoffset += width;
        }

        glEnd();

#if 0
        glBegin(GL_QUADS); {
            float width = 1.0 / m_bins.numBins();
            float xoffset = 0;
            float yoffset = 0;
            for (size_t i = 0; i < m_bins.numBins(); ++i) {
                float height = m_bins.count(i) / maxCount;
                glVertex3f(xoffset        , yoffset, 0);
                glVertex3f(xoffset + width, yoffset, 0);
                glVertex3f(xoffset + width, yoffset + height, 0);
                glVertex3f(xoffset        , yoffset + height, 0);
                xoffset += width;
            }
        }
        glEnd();
#endif
    }

    void setResolution(size_t resolution) { m_resolution = resolution; }

private:
    size_t m_resolution;
};

// Wrap a view, providing a proportional margin around it.
class Margin : public View {
public:
    Margin(ViewPtr view, float margin = 0.02)
        : m_subview(view), m_margin(margin) { }

    void setColor(const RGBColorf &c) { m_color = c; }

    virtual void render() const {
        glPushMatrix();

        glTranslatef(m_margin, m_margin, 0);
        glScalef(scale(), scale(), 1);

        glColor4fv(m_color);
        glBegin(GL_QUADS); {
            glVertex3f(0, 0, 0);
            glVertex3f(1, 0, 0);
            glVertex3f(1, 1, 0);
            glVertex3f(0, 1, 0);
        }
        glEnd();

        m_subview->render();

        glPopMatrix();
    }

    float scale() const { return 1 - 2 * m_margin; }

    virtual ViewPtr viewForPoint(float x, float y) {
        ViewTransformf xf = m_xf();
        return m_subview->viewForPoint(xf.transformedX(x, y), xf.transformedY(x, y));
    }

    virtual ViewTransformf viewTransformForPoint(float x, float y) {
        // cout << "Getting xform for margined view" << endl;
        ViewTransformf xf = m_xf();
        auto subviewTransform = m_subview->viewTransformForPoint(xf.transformedX(x, y),
                                                                 xf.transformedY(x, y));
        return xf.composeRight(subviewTransform);
    }

private:
    ViewTransformf m_xf() const {
        // Capture margin and scale by value, not this ptr
        float m = m_margin;
        float s = scale();
        return ViewTransformf([=](float xx, float yy) { return (xx - m) / s; },
                              [=](float xx, float yy) { return (yy - m) / s; });
    }

    float m_margin;
    ViewPtr m_subview;
    RGBColorf m_color;
};

class Layout : public View {
public:
    // Stacked views are on top of each other and render in the order added
    // (last added on top).
    enum class Direction { Horizontal, Vertical, Stacked };
    Layout(Direction dir) : m_layoutDir(dir) { }

    // Warning: resizes all subviews to equal percentages
    ViewPtr addSubview(ViewPtr view) {
        m_subviews.push_back(view);
        m_subviewPercentages.assign(m_subviews.size(), 1.0 / m_subviews.size());
        return view;
    }

    void clearSubviews() {
        m_subviews.clear();
        m_subviewPercentages.clear();
    }

    ViewPtr subview(size_t i) { return m_subviews.at(i); }

    void setSubviewPercentages(const vector<float> &pct) {
        if (pct.size() != m_subviews.size()) throw runtime_error("Must specify every subview's layout pct");
        float tot = 0;
        for (float p : pct) tot += p;
        if (std::abs(tot - 1.0) > 1e-8) throw runtime_error("Layout pct must add to 100");
        m_subviewPercentages = pct;
    }

    virtual void render() const {
        // Offset and size determine the 1D "bounding box" of the subview.
        float offset = 0;
        for (size_t i = 0; i < m_subviews.size(); ++i) {
            glPushMatrix();
            float size = m_subviewPercentages.at(i);
            if (m_layoutDir == Direction::Horizontal) {
                glTranslatef(offset, 0, 0);
                glScalef(size, 1, 1);
            }
            else if (m_layoutDir == Direction::Vertical) {
                glTranslatef(0, offset, 0);
                glScalef(1, size, 1);
            }
            else if (m_layoutDir == Direction::Stacked) {
                // Each stacked view gets the full current canvas
            }
            else assert(false);

            m_subviews[i]->render();

            glPopMatrix();
            offset += size;
        }
    }

    virtual ViewPtr viewForPoint(float x, float y) {
        float size, offset;
        auto xf = m_xf(size, offset);
        return m_findSubview(x, y, size, offset)->viewForPoint(xf.transformedX(x, y), xf.transformedY(x, y));
    }

    virtual ViewTransformf viewTransformForPoint(float x, float y) {
        // cout << "Getting xform for layout" << endl;
        float size, offset;
        ViewPtr subview = m_findSubview(x, y, size, offset);
        if (!subview) return ViewTransformf();

        // left-compose the transform from a point in this view's coordinates to
        // the subview's coordinates
        auto xf = m_xf(size, offset);
        ViewTransformf subviewTransform = subview->viewTransformForPoint(xf.transformedX(x, y), xf.transformedY(x, y));

        return xf.composeRight(subviewTransform);
    }

protected:
    // Find the subview at a particular point in this view's coordinate
    // and get its size/offset in this view.
    ViewPtr m_findSubview(float x, float y, float &size, float &offset) {
        ViewPtr responder;
        offset = 0;
        for (size_t i = 0; i < m_subviews.size(); ++i) {
            size = m_subviewPercentages[i];
            if (m_layoutDir == Direction::Horizontal) {
                if (x < size) { responder = m_subviews[i]; break; }
                x -= size;
                offset += size;
            }
            else if (m_layoutDir == Direction::Vertical) {
                if (y < size) { responder = m_subviews[i]; break; }
                y -= size;
                offset += size;
            }
            else if (m_layoutDir == Direction::Stacked) {
                // Return the topmost (*last*) subview that responds
                size = 1.0; // Stacked views keep the same coordinate system.
                auto rsv = m_subviews[i]->viewForPoint(x, y);
                if (rsv) responder = rsv;
            }
            else assert(false);
        }
        return responder;
    }

    // Get the view transform into the subview with given size and offset
    ViewTransformf m_xf(float size, float offset) const {
        ViewTransformf xf;
        if (m_layoutDir == Direction::Horizontal) {
            xf = ViewTransformf([=](float xx, float yy) { return (xx - offset) / size; },
                                [=](float xx, float yy) { return yy; });
        }
        else if (m_layoutDir == Direction::Vertical) {
            xf = ViewTransformf([=](float xx, float yy) { return xx; },
                                [=](float xx, float yy) { return (yy - offset) / size; });
        }
        else if (m_layoutDir == Direction::Stacked) {
            // xf = leave as identity
        }
        else assert(false);
        return xf;
    }

    Direction m_layoutDir;
    vector<ViewPtr> m_subviews;
    vector<float> m_subviewPercentages;
};


// A clicked review continues to respond to events until all mouse buttons are
// released.
// TODO: notify if views are cleared to disable active view?
class UIHandler {
public:
    UIHandler(ViewPtr view) : m_mainView(view) { }

    int decodeGLUTButton(int button) {
        switch(button) {
            case GLUT_LEFT_BUTTON:   return 0;
            case GLUT_MIDDLE_BUTTON: return 1;
            case GLUT_RIGHT_BUTTON:  return 2;
        }
        return 3;
    }

    // Return true if handled
    // TODO: check if the view we clicked on accepted the event before making it
    // active!
    bool mouseFunc(int button, int state, float x, float y) {
        bool handled = false;
        button = decodeGLUTButton(button);
        if (button > 2) return false; // Only handle left/middle/right button
        if (state == GLUT_DOWN) {
            assert(m_buttonsHeld[button] == false);
            if (m_anyHeld()) assert(m_activeView);
            else {
                m_activeView = m_mainView->viewForPoint(x, y);
                m_activeViewTransform = m_mainView->viewTransformForPoint(x, y);
            }
            m_buttonsHeld[button] = true;
            handled = m_activeView->mousePress(button, m_activeViewTransform.transformedX(x, y),
                                                       m_activeViewTransform.transformedY(x, y));
        }
        else if (state == GLUT_UP) {
            assert(m_activeView);
            m_buttonsHeld[button] = false;
            handled = m_activeView->mouseRelease(button, m_activeViewTransform.transformedX(x, y),
                                                         m_activeViewTransform.transformedY(x, y));
            if (!m_anyHeld()) { m_activeView.reset(); }
        }

        if (m_activeView) { m_firstResponder = m_activeView; }
        return handled;
    }

    bool motionFunc(float x, float y) {
        // motion should only be called between a mouse down and a mouse up
        // (where m_activeView should be set), but GLUT appears to sometimes
        // reorder things :(
        // if (!m_activeView) cout << "WARNING: drag on inactive view" << endl;
        if (!m_activeView) return false;
        return m_activeView->mouseDrag(m_activeViewTransform.transformedX(x, y),
                                       m_activeViewTransform.transformedY(x, y));
    }

    // Forward keypress to most recently clicked view (first responder)
    bool keyboardFunc(unsigned char c) {
        if (!m_firstResponder) return false;
        return m_firstResponder->keyPressed(c);
    }

    void setFirstResponder(ViewPtr view) { m_firstResponder = view; }

private:
    ViewPtr m_mainView;
    bool m_anyHeld() { return m_buttonsHeld[0] || m_buttonsHeld[1] || m_buttonsHeld[2]; }
    std::vector<bool> m_buttonsHeld = std::vector<bool>(4, false); // GLUT_LEFT_BUTTON, MIDDLE, RIGHT, other
    ViewPtr m_activeView, m_firstResponder;
    ViewTransformf m_activeViewTransform;
};

////////////////////////////////////////////////////////////////////////////////
// Globals
////////////////////////////////////////////////////////////////////////////////
ViewPtr g_mainView;
shared_ptr<UIHandler> g_uiHandler;
shared_ptr<Layout> g_histLayoutLeft, g_histLayoutRight;

vector<shared_ptr<HistogramPlot>> g_selectedParamHistograms;
// Ranges of each parameter over all tables.
vector<DataRange> g_paramRanges;
shared_ptr<LogscaleYPlot> g_allTablesLogscalePlot, g_currentTableLogscalePlot;

vector<LUT> g_dataTables, g_dataTablesCurrPattern;
vector<string> g_dataTablePaths;

set<size_t> g_patterns;
size_t g_currPattern, g_currTableIdx;

////////////////////////////////////////////////////////////////////////////////
// View data for a particular pattern
////////////////////////////////////////////////////////////////////////////////
void selectTableAndPattern(size_t table, size_t pat) {
    g_currTableIdx = table;
    g_currPattern = pat;
    glutSetWindowTitle(("LUT Viewer: " + g_dataTablePaths.at(g_currTableIdx) +
                        ", Pattern " + to_string(g_currPattern)).c_str());

    // Get the data for the current pattern
    g_dataTablesCurrPattern.clear();
    g_dataTablesCurrPattern.reserve(g_dataTables.size());
    for (size_t i = 0; i < g_dataTables.size(); ++i)
        g_dataTablesCurrPattern.push_back(g_dataTables[i].selectPattern(g_currPattern));

    vector<vector<vector<float>>> series_paramValues; // series->param->pt_idx
    size_t nParams = g_dataTablesCurrPattern.front().numParams();

    ////////////////////////////////////////////////////////////////////////////
    // Extract the series data and update the all-series logscale plot of all
    // series
    ////////////////////////////////////////////////////////////////////////////
    vector<vector<float>> series_xdata, series_ydata;
    series_xdata.reserve(g_dataTablesCurrPattern.size());
    series_ydata.reserve(g_dataTablesCurrPattern.size());
    for (const auto &tab : g_dataTablesCurrPattern) {
        series_xdata.emplace_back(tab.getNus<float>());
        series_ydata.emplace_back(tab.getEs<float>());
    }

    g_allTablesLogscalePlot->setData(series_xdata, series_ydata);
    g_allTablesLogscalePlot->autoPlotRange();

    ////////////////////////////////////////////////////////////////////////////
    // Create new parameter histograms
    // WARNING: This must happen before changing the *current* logscale plot data,
    // otherwise the selection changed callback will access histograms that
    // don't exist yet.
    ////////////////////////////////////////////////////////////////////////////
    // Get parameter values and compute their ranges
    g_paramRanges.assign(nParams, DataRange());
    for (const auto &tab : g_dataTablesCurrPattern) {
        series_paramValues.emplace_back(nParams, vector<float>());
        for (size_t p = 0; p < nParams; ++p) {
            series_paramValues.back()[p] = tab.getParamValues<float>(p);
            g_paramRanges[p] = g_paramRanges[p].unionRange(series_paramValues.back()[p]);
        }
    }

    g_histLayoutLeft->clearSubviews();
    g_histLayoutRight->clearSubviews();
    g_selectedParamHistograms.clear();

    for (size_t p = 0; p < nParams; ++p) {
        auto stacker = make_shared<Layout>(Layout::Direction::Stacked);
        bool bottom = true;
        for (size_t s = 0; s < series_paramValues.size(); ++s) {
            if (s == g_currTableIdx) continue; // histogram for all tables but the current.
            auto paramHistSeries = make_shared<HistogramPlot>(series_paramValues.at(s).at(p), g_paramRanges[p]);
            if (bottom) paramHistSeries->setBackgroundColor(1, 1, 1, 1);
            else        paramHistSeries->setBackgroundColor(1, 1, 1, 0);
            bottom = false;
            paramHistSeries->setForegroundColor(g_allTablesLogscalePlot->seriesColor(s));
            stacker->addSubview(paramHistSeries);
        }

        // Put current histogram on top.
        {
            size_t s = g_currTableIdx; // histogram for current table
            auto currTableParamHist = make_shared<HistogramPlot>(series_paramValues.at(s).at(p), g_paramRanges[p]);
            if (bottom) currTableParamHist->setBackgroundColor(1, 1, 1, 1);
            else        currTableParamHist->setBackgroundColor(1, 1, 1, 0.25);
            currTableParamHist->setForegroundColor(0, 0, 0, 0.75);
            stacker->addSubview(currTableParamHist);
        }

        auto paramHistSelected = make_shared<HistogramPlot>();
        paramHistSelected->setBackgroundColor(1, 1, 1, 0);
        paramHistSelected->setForegroundColor(1, 0, 0, 0.95);
        stacker->addSubview(paramHistSelected);
        auto marginedStacker = make_shared<Margin>(stacker);
        if (p & 1) g_histLayoutRight->addSubview(marginedStacker);
        else       g_histLayoutLeft ->addSubview(marginedStacker);

        g_selectedParamHistograms.push_back(paramHistSelected);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Update the *current* logscale plot
    ////////////////////////////////////////////////////////////////////////////
    g_currentTableLogscalePlot->setData(series_xdata.at(g_currTableIdx), series_ydata.at(g_currTableIdx));
    g_currentTableLogscalePlot->setPlotRange(g_allTablesLogscalePlot->xRange(), g_allTablesLogscalePlot->yRange());

    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////
// IO Callback Routines
////////////////////////////////////////////////////////////////////////////////
void mouseToViewCoordinates(int x, int y, float &xv, float &yv) {
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    float width  = viewport[2];
    float height = viewport[3];
    xv = x / width;
    yv = (height - y) / height;
}

void PassiveMotionFunc(int x, int y) { ; }

void MotionFunc(int x, int y) {
    float xv, yv;
    mouseToViewCoordinates(x, y, xv, yv);
    g_uiHandler->motionFunc(xv, yv);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
    float xv, yv;
    mouseToViewCoordinates(x, y, xv, yv);
    g_uiHandler->mouseFunc(button, state, xv, yv);
    glutPostRedisplay();
}

void SpecialKeyboardFunc(int k, int x, int y) { ; }
void KeyboardFunc(unsigned char c, int x, int y) {
    // x, y don't seem to be set on os x--ignore them
    bool handled = false;
    if (c == 'w') {
        size_t i = 0;
        auto filePath = [](size_t i) { return "selection" + std::to_string(i) + ".txt"; };
        while (boost::filesystem::exists(filePath(i))) ++i;

        // Write selected subset.
        LUT selectedLUT = g_dataTablesCurrPattern.at(g_currTableIdx).subset(g_currentTableLogscalePlot->selectedPointIndices());
        if (selectedLUT.size() > 0) {
            selectedLUT.write(filePath(i));
            cout << "Wrote " << selectedLUT.size() << " points to " << filePath(i) << endl;
        }
    }
    if ((tolower(c) == 'p') && (g_patterns.size() > 1)) {
        // Cycle through patterns
        auto p = find(g_patterns.begin(), g_patterns.end(), g_currPattern);
        if (p == g_patterns.end()) p = g_patterns.begin();
        else {
            // p: increment, P: decrement
            if (c == 'p') { if (++p == g_patterns.end()) p = g_patterns.begin(); }
            else {
                if (p == g_patterns.begin()) p = g_patterns.end();
                --p;
            }
        }
        selectTableAndPattern(g_currTableIdx, *p);
        handled = true;
    }
    if ((tolower(c) == 't') && (g_dataTables.size() > 1)) {
        // cycle through tables
        size_t newIdx = g_currTableIdx + (c == 't' ? 1 : -1);
        if (newIdx == g_dataTables.size()) newIdx = 0;
        if (newIdx > g_dataTables.size()) newIdx = g_dataTables.size() - 1;
        selectTableAndPattern(newIdx, g_currPattern);
        handled = true;
    }
    if (!handled) handled = g_uiHandler->keyboardFunc(c);
    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////
// View setup
////////////////////////////////////////////////////////////////////////////////
void Reshape(int width, int height) {
    // Set OpenGL viewport and camera
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////
/*! Called by GLUT when redisplay needed
*///////////////////////////////////////////////////////////////////////////////
void Display()
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    int width  = viewport[2];
    int height = viewport[3];

    glShadeModel(GL_FLAT);
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    // Use antialiasing
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-1, -1, 0.0);
    glScalef(2.0, 2.0, 1.0);

    g_mainView->render();

    glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on sucess)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    if (argc < 2) {
        cerr << "Usage: data_viewer dataTable [dataTable2...]" << endl;
        exit(-1);
    }

    for (int i = 1; i < argc; ++i)
        g_dataTablePaths.emplace_back(argv[i]);

    // Initialize GLUT
    glutInit(&argc, argv);

    int width = 1024, height = 768;
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("Lookup Table Viewer");

    cout << "Help" << endl;
    cout << "+/-:\tchange selection radius" << endl;
    cout << "p/P:\tcycle through patterns" << endl;
    cout << "t/T:\tcycle through data tables" << endl;

    // Set GLUT event callbacks
    glutMouseFunc(MouseFunc);
    glutMotionFunc(MotionFunc);
    glutPassiveMotionFunc(PassiveMotionFunc);
    glutKeyboardFunc(KeyboardFunc);
    glutSpecialFunc(SpecialKeyboardFunc);

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);

    // TODO: make pattern configurable.
    WMesh wmesh("../../../Luigi/wireinflator2D/meshes/octa_cell.obj");
    g_patternSignedDistance = Future::make_unique<PSD>(wmesh);

    // try {
        auto mainLayout = make_shared<Layout>(Layout::Direction::Horizontal);
        g_mainView = mainLayout;
        g_allTablesLogscalePlot = make_shared<LogscaleYPlot>();
        g_allTablesLogscalePlot->setBackgroundColor(1, 1, 1, 1.0); // white
        g_allTablesLogscalePlot->setSeriesColormap(COLORMAP_JET);
        g_allTablesLogscalePlot->seriesColormap()->setAlpha(0.5);
        g_allTablesLogscalePlot->setPointSize(2.0);

        g_currentTableLogscalePlot = make_shared<LogscaleYPlot>();
        g_currentTableLogscalePlot->setBackgroundColor(1, 1, 1, 0); // transparent
        g_currentTableLogscalePlot->setForegroundColor(0.0,0.0,0.0);
        g_currentTableLogscalePlot->setSelectionColor(1.0,0.0,0.0);
        g_currentTableLogscalePlot->setPointSize(4.0);

        // stack g_currentTableLogscalePlot atop g_allTablesLogscalePlot
        auto plotStacker = make_shared<Layout>(Layout::Direction::Stacked);
        plotStacker->addSubview(g_allTablesLogscalePlot);
        plotStacker->addSubview(g_currentTableLogscalePlot);


        auto sideBar = make_shared<Layout>(Layout::Direction::Vertical);
        g_histLayoutLeft  = make_shared<Layout>(Layout::Direction::Vertical);
        g_histLayoutRight = make_shared<Layout>(Layout::Direction::Vertical);

        auto plotWithMargin = make_shared<Margin>(plotStacker);

        auto histLayouts = make_shared<Layout>(Layout::Direction::Horizontal);
        histLayouts->addSubview(g_histLayoutLeft);
        histLayouts->addSubview(g_histLayoutRight);
        sideBar->addSubview(histLayouts);
        sideBar->addSubview(make_shared<PatternView>(256));
        sideBar->setSubviewPercentages({0.67, 0.33});

        mainLayout->addSubview(plotWithMargin);
        mainLayout->addSubview(sideBar);

        mainLayout->setSubviewPercentages({0.67, 0.33});

        g_uiHandler = make_shared<UIHandler>(g_mainView);
        g_uiHandler->setFirstResponder(g_currentTableLogscalePlot);

        for (const string &s : g_dataTablePaths) {
            g_dataTables.emplace_back(s);
            g_dataTables.back().patterns(g_patterns);
        }

        cout << "Patterns:";
        for (size_t pat : g_patterns)
            cout << " " << pat;
        cout << endl;

        if (g_patterns.size() == 0) throw runtime_error("No patterns read.");

        selectTableAndPattern(0, *g_patterns.begin());

        // Must set the callback after creating the histograms (in
        // selectTableAndPattern)--otherwise it may be called before the
        // histogram views exist.
        g_currentTableLogscalePlot->setSelectionChangedCallback(
            [&](ViewPtr v) {
                auto lv = dynamic_pointer_cast<LogscaleYPlot>(v);
                assert(lv);
                auto selSubset = g_dataTablesCurrPattern.at(g_currTableIdx).subset(lv->selectedPointIndices());
                size_t nParams = g_dataTablesCurrPattern.at(g_currTableIdx).numParams();
                std::vector<float> firstPatternParams(nParams), paramValues;
                for (size_t p = 0; p < nParams; ++p) {
                    if (selSubset.size()) {
                        paramValues = selSubset.getParamValues<float>(p);
                        firstPatternParams[p] = paramValues.front();
                    }
                    else paramValues.clear();


                    g_selectedParamHistograms[p]->setData(paramValues, g_paramRanges[p]);
                }

                // Render the first selected pattern
                g_patternSignedDistance->setParameters(firstPatternParams);
            }
        );

        // Call the GLUT main loop
        glutMainLoop();
    // }
    // catch (const exception &e) {
    //     cerr << "Caught exception: " << e.what() << endl;
    // }

    return 0;
}
