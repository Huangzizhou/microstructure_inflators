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

#include "colors.hh"

#include <boost/filesystem/operations.hpp>

using namespace std;
// #include <Eigen/Dense>
// using Eigen::Vector2f;

// views' render methods assume that the view transforms have been set up so
// that they should render in [0, 1]^2
// Also, expects event x, y to be in [0, 1]^2
class View;
typedef shared_ptr<View> ViewPtr;

#include <sstream>
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
        elems.push_back(item);
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void drawString(const std::string &str) {
    for (char c : str)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, c);
}

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

    // Compose b on the right (apply this transform after b).
    ViewTransform composeRight(const ViewTransform &b) const {
        // WARNING: capturing member variables by value is impossible! A
        // by-value default capture really just captures the "this" pointer and
        // accesses the members through it!
        ViewTransform copy(*this);
        return ViewTransform(
            [=](Real x, Real y) -> Real { return copy.transformedX(b.transformedX(x, y), b.transformedY(x, y)); },
            [=](Real x, Real y) -> Real { return copy.transformedY(b.transformedX(x, y), b.transformedY(x, y)); });
    }

    // Compose b on the left (apply this transform before b).
    ViewTransform composeLeft(const ViewTransform &b) const {
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

class View : public std::enable_shared_from_this<View> {
public:
    virtual void render() const = 0;

    virtual ViewPtr viewForPoint(float x, float y) = 0;
    virtual ViewTransformf viewTransformForPoint(float x, float y) = 0;

    virtual bool mousePress(int button, float x, float y)   { return false; }
    virtual bool mouseRelease(int button, float x, float y) { return false; }
    virtual bool mouseDrag(float x, float y)                { return false; }
    virtual bool keyPressed(unsigned char c)                { return false; }

    typedef std::function<void(std::shared_ptr<View>)> CallBack;
};

// A range that supports auto-min and auto-max
// The interval is [min, max)
struct DataRange {
    DataRange() : m_min(0), m_max(0) { }
    DataRange(const std::vector<float> &data) { setAutoMin(data); setAutoMax(data); }
    DataRange(float min, float max) : m_min(min), m_max(max) { }

    void setAutoMin(const std::vector<float> &data) { m_min = (data.size() == 0) ? 0.0 : *std::min_element(data.begin(), data.end()); }
    void setAutoMax(const std::vector<float> &data) { m_max = (data.size() == 0) ? 0.0 : *std::max_element(data.begin(), data.end()); }
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

private:
    float m_min, m_max;
};

struct DataBins {
    DataBins() { } // No bins
    DataBins(const std::vector<float> &data, const DataRange &dr, size_t numBins = 20) { m_computeBins(data, dr, numBins); }
    DataBins(const std::vector<float> &data, size_t numBins = 20) { m_computeBins(data, DataRange(data), numBins); }
    size_t numBins() const { return m_bins.size(); }

    size_t count(size_t i) const { return m_bins.at(i); }
    size_t maxCount() const {
        if (numBins() == 0) return 0;
        return *std::max_element(m_bins.begin(), m_bins.end());
    }

private:
    void m_computeBins(const std::vector<float> &data, const DataRange &dr, size_t numBins) {
        m_bins.assign(numBins, 0);
        for (auto f : data) {
            // could optimize
            for (size_t i = 0; i < numBins; ++i) {
                if (dr.subrange(numBins, i).inRange(f)) {
                    ++m_bins[i];
                    break;
                }
            }
        }
    }
    std::vector<size_t> m_bins;
};

class Plot : public View {
public:
    void setForegroundColor(const RGBColorf &c) { m_fgcolor = c; }
    void setBackgroundColor(const RGBColorf &c) { m_bgcolor = c; }
    void setSelectionColor (const RGBColorf &c) { m_slcolor = c; }
    void setForegroundColor(float r, float g, float b, float a = 1.0) { setForegroundColor(RGBColorf(r, g, b, a)); }
    void setBackgroundColor(float r, float g, float b, float a = 1.0) { setBackgroundColor(RGBColorf(r, g, b, a)); }
    void setSelectionColor (float r, float g, float b, float a = 1.0) { setSelectionColor (RGBColorf(r, g, b, a)); }

    const RGBColorf &foregroundColor() const { return m_fgcolor; }
    const RGBColorf &backgroundColor() const { return m_bgcolor; }
    const RGBColorf &selectionColor () const { return m_slcolor; }

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

    virtual ViewPtr viewForPoint(float x, float y) { return shared_from_this(); }
    // Gets the transformation mapping coordinates of a point in the
    // top-level-view's system to coordinates in the system of the subview at (x, y)
    virtual ViewTransformf viewTransformForPoint(float x, float y) { return ViewTransformf(); }
    
protected:
    RGBColorf m_fgcolor, m_bgcolor, m_slcolor;
    float m_pointSize = 1.0;
};


class HistogramPlot : public Plot {
public:
    HistogramPlot() { }
    HistogramPlot(const std::vector<float> &data) { setData(data); }
    HistogramPlot(const std::vector<float> &data, size_t numBins) { setData(data, numBins); }
    HistogramPlot(const std::vector<float> &data, size_t numBins, const DataRange &dr) { setData(data, dr, numBins); }

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
                xoffset += width;
                glVertex3f(xoffset        , yoffset, 0);
                glVertex3f(xoffset + width, yoffset, 0);
                glVertex3f(xoffset + width, yoffset + height, 0);
                glVertex3f(xoffset        , yoffset + height, 0);
            }
        }
        glEnd();
    }

    void setData(const std::vector<float> &data, size_t numBins = 20) {
        setData(data, DataRange(data), numBins);
    }
    void setData(const std::vector<float> &data, const DataRange &dr, size_t numBins = 20) {
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

class LogscalePlot : public Plot {
public:
    LogscalePlot() { }
    void setData(const std::vector<float> &xdata, const std::vector<float> &ydata) {
        if (xdata.size() != ydata.size())
            throw std::runtime_error("Unmatched data sizes");
        m_xdata = xdata;
        m_ydata = ydata;
        resetSelection();
    }

    void setPlotRange(DataRange xRange, DataRange yRange) {
        m_xRange = xRange;
        m_yRange = yRange;
    }
    void autoPlotRange() {
        m_xRange = DataRange(m_xdata);
        m_yRange = DataRange(m_ydata);
    }

    DataRange xRange() const { return m_xRange; }
    DataRange yRange() const { return m_yRange; }

    virtual void render() const {
        Plot::render();
        glColor3fv(foregroundColor());
        glPointSize(m_pointSize);
        glBegin(GL_POINTS); {
            for (size_t i = 0; i < numPoints(); ++i) {
                float x = m_xdata[i], y = m_ydata[i];
                m_plotToView(x, y);
                if (m_pointSelected[i]) glColor3fv(selectionColor());
                glVertex3f(x, y, 0.0);
                if (m_pointSelected[i]) glColor3fv(foregroundColor());
            }
        }
        glEnd();

        size_t selSize = selectionSize();
        if (selSize > 0) {
            glColor3fv(selectionColor());
            glRasterPos2i(0, 0);
            drawString(std::to_string(selSize) + " points selected around ("
                    + std::to_string(m_lastSelectCenterX) + ", "
                    + std::to_string(m_lastSelectCenterY) + ")");
        }
    }

    size_t numPoints() const { return m_xdata.size(); }

    void setSelectionRadius(float radius) {
        m_selectionRadius = radius;
        m_selectPoint(m_lastSelectCenterX, m_lastSelectCenterY);
    }
    float selectionRadius() const { return m_selectionRadius; }

    // The aspect ratio of the selection ellipse (width / height)
    void setSelectionAspect(float aspect) { m_selectionAspect = aspect; }

    void resetSelection() { m_pointSelected.assign(numPoints(), false); }
    size_t selectionSize() const {
        size_t count = 0;
        for (bool b : m_pointSelected) if (b) ++count;
        return count;
    }
    std::vector<size_t> selectedPointIndices() const {
        std::vector<size_t> result;
        result.reserve(selectionSize());
        for (size_t i = 0; i < numPoints(); ++i) if (m_pointSelected[i]) result.push_back(i);
        return result;
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
    float m_selectionRadius = 5.0 / 512, m_selectionAspect = 1.0;
    DataRange m_xRange, m_yRange;
    std::vector<float> m_xdata, m_ydata;
    std::vector<bool> m_pointSelected;
    CallBack m_selChangedCallback;
    float m_lastSelectCenterX = std::numeric_limits<float>::min(),
          m_lastSelectCenterY = std::numeric_limits<float>::min();

    // Select the data points around a click location (in view coordinates)
    // The distance is measured in view space, taking aspect ratio into account.
    void m_selectPoint(float x, float y) {
        m_lastSelectCenterX = x, m_lastSelectCenterY = y;

        std::vector<bool> oldSelection(m_pointSelected);
        bool selectionChanged = false;
        resetSelection();
        for (size_t i = 0; i < numPoints(); ++i) {
            float xd = m_xdata[i], yd = m_ydata[i];
            m_plotToView(xd, yd);
            float distx = xd - x, disty = yd - y;
            distx /= m_selectionAspect;
            float dist = sqrt(distx * distx + disty * disty);
            if (dist < m_selectionRadius) 
                m_pointSelected[i] = true;
            selectionChanged  = selectionChanged || (m_pointSelected[i] != oldSelection[i]);
        }
        if (selectionChanged && m_selChangedCallback)
            m_selChangedCallback(shared_from_this());
    }

    // in-place coordinate transforms
    void m_plotToView(float &inoutx, float &inouty) const {
        inoutx = (log(inoutx) - log(m_xRange.min())) / (log(m_xRange.max()) - log(m_xRange.min()));
        inouty = (log(inouty) - log(m_yRange.min())) / (log(m_yRange.max()) - log(m_yRange.min()));
    }
    void m_viewToPlot(float &inoutx, float &inouty) const {
        inoutx = exp(inoutx * (log(m_xRange.max()) - log(m_xRange.min())) + log(m_xRange.min()));
        inouty = exp(inouty * (log(m_xRange.max()) - log(m_xRange.min())) + log(m_xRange.min()));
    }

    // out-of-place coordinate transforms
    void m_plotToView(float inx, float iny, float &outx, float &outy) const { outx = inx, outy = iny; m_plotToView(outx, outy); }
    void m_viewToPlot(float inx, float iny, float &outx, float &outy) const { outx = inx, outy = iny; m_viewToPlot(outx, outy); }
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
        return m_findSubview(x, y, size, offset);
    }

    virtual ViewTransformf viewTransformForPoint(float x, float y) { 
        float size, offset;
        ViewPtr subview = m_findSubview(x, y, size, offset);
        if (!subview) return ViewTransformf();

        // left-compose the transform from a point in this view's coordinates to
        // the subview's coordinates
        ViewTransformf subviewTransform = subview->viewTransformForPoint(x, y);
        if (m_layoutDir == Direction::Horizontal) {
            return subviewTransform.composeLeft(
                    ViewTransformf([=](float xx, float yy) { return xx / size - offset; },
                                   [=](float xx, float yy) { return yy; }));
        }
        else if (m_layoutDir == Direction::Vertical) {
            return subviewTransform.composeLeft(
                    ViewTransformf([=](float xx, float yy) { return xx; },
                                   [=](float xx, float yy) { return yy / size - offset; }));
        }
        else if (m_layoutDir == Direction::Stacked) {
            return subviewTransform;
        }
        else assert(false);
        return ViewTransformf();
    }

protected:
    // Find the subview at a particular point in this view's coordinate
    // and get its size/offset in this view.
    ViewPtr m_findSubview(float x, float y, float &size, float &offset) {
        for (size_t i = 0; i < m_subviews.size(); ++i) {
            size = m_subviewPercentages[i];
            offset = 0;
            if (m_layoutDir == Direction::Horizontal) {
                if (x < size)
                    return m_subviews[i]->viewForPoint(x / size, y);
                x -= size;
                offset += size;
            }
            else if (m_layoutDir == Direction::Vertical) {
                if (y < size)
                    return m_subviews[i]->viewForPoint(x, y / size);
                y -= size;
                offset += size;
            }
            else if (m_layoutDir == Direction::Stacked) {
                // Return the first subview that responds
                size = 1.0; // Stacked views keep the same coordinate system.
                auto respondedSubview = m_subviews[i]->viewForPoint(x, y); 
                if (respondedSubview) return respondedSubview;
            }
            else assert(false);
        }
        return ViewPtr();
    }

    Direction m_layoutDir;
    vector<ViewPtr> m_subviews;
    vector<float> m_subviewPercentages;
};


// A clicked review continues to respond to events until all mouse buttons are
// released.
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

        if (m_activeView) m_firstResponder = m_activeView;
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

struct LookupTable {
    LookupTable() { }
    LookupTable(const std::string &path) { load(path); }

    void load(const std::string &path) {
        std::ifstream lutStream(path);
        if (!lutStream.is_open())
            throw std::runtime_error("Couldn't open " + path);

        std::string line;
        while (std::getline(lutStream, line)) {
            records.emplace_back(line);
        }
    }

    void write(const std::string &path) const {
        std::ofstream outFile(path);
        if (!outFile.is_open())
            throw std::runtime_error("Couldn't open " + path + " for writing");
        write(outFile);
    }

    void write(std::ostream &os) const {
        for (auto r : records)
            os << r << endl;
    }

    struct Record {
        Record(const std::string &lutLine) {
            std::runtime_error invalidColSize("Invalid number of columns in LUT.");
            auto tokens = split(lutLine, '\t');
            if (tokens.size() < 5) throw invalidColSize;
            size_t nParams = std::stoi(tokens[4]);
            if (tokens.size() != 5 + nParams) throw invalidColSize;
            pattern = std::stoi(tokens[0]);
            E = std::stod(tokens[1]);
            nu = std::stod(tokens[2]);
            anisotropy = std::stod(tokens[3]);

            for (size_t i = 0; i < nParams; ++i)
                patternParams.push_back(std::stod(tokens[5 + i]));
        }

        size_t numParams() const { return patternParams.size(); }
        size_t pattern;
        double E, nu, anisotropy;
        std::vector<double> patternParams;

        // Write record without trailing newline.
        friend std::ostream &operator<<(std::ostream &os, const Record &r) {
            std::ios state(NULL);
            state.copyfmt(os);
            os << std::setfill('0') << std::setw(4) << r.pattern << '\t';
            os.copyfmt(state);
            os << r.E << '\t' << r.nu << '\t' << r.anisotropy << '\t' << r.numParams(); 
            for (auto p : r.patternParams)
                os << '\t' << p;
            return os;
        }
    };

    std::set<size_t> patterns() const {
        std::set<size_t> result;
        for (const auto &r : records) result.insert(r.pattern);
        return result;
    }

    size_t  pattern(size_t r) const { return records.at(r).pattern; }
    double        E(size_t r) const { return records.at(r).E; }
    double       nu(size_t r) const { return records.at(r).nu; }
    double        A(size_t r) const { return records.at(r).anisotropy; }
    size_t numParams() const { 
        std::set<size_t> result;
        for (const auto &r : records) result.insert(r.numParams());
        if (result.size() > 1) throw std::runtime_error("All LUT entries must have the same number of paramters for this operation.");
        return *result.begin();
    }

    // All of this is extremely inefficient for now!
    LookupTable selectPattern(size_t pattern) {
        LookupTable result;
        for (const auto &r : records) {
            if (r.pattern == pattern)
                result.records.push_back(r);
        }
        return result;
    }
    LookupTable subset(const std::vector<size_t> idxs) {
        LookupTable result;
        for (size_t i : idxs) result.records.push_back(records.at(i));
        return result;
    }

    size_t size() const { return records.size(); }

    // Allow type conversion of fields
    template<typename T> void  getEs(                   std::vector<T> &out) const { out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.E); }
    template<typename T> void getNus(                   std::vector<T> &out) const { out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.nu); }
    template<typename T> void  getAs(                   std::vector<T> &out) const { out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.anisotropy); }
    template<typename T> void  getParamValues(size_t p, std::vector<T> &out) const {
        if (p >= numParams()) { throw std::runtime_error("Invaild parameter"); }     out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.patternParams.at(p));
    }

    template<typename T> std::vector<T>  getEs()                  const { std::vector<T> result;              getEs(result); return result; }
    template<typename T> std::vector<T> getNus()                  const { std::vector<T> result;             getNus(result); return result; }
    template<typename T> std::vector<T>  getAs()                  const { std::vector<T> result;              getAs(result); return result; }
    template<typename T> std::vector<T>  getParamValues(size_t p) const { std::vector<T> result;  getParamValues(p, result); return result; }

    std::vector<Record> records;
};

////////////////////////////////////////////////////////////////////////////////
// Globals
////////////////////////////////////////////////////////////////////////////////
ViewPtr g_mainView;
shared_ptr<UIHandler> g_uiHandler;
LookupTable g_lut;
vector<shared_ptr<HistogramPlot>> g_selectedParamHistograms, g_fullParamHistograms;
shared_ptr<LogscalePlot> g_logscalePlot;

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
        LookupTable selectedLUT = g_lut.subset(g_logscalePlot->selectedPointIndices());
        if (selectedLUT.size() > 0) {
            selectedLUT.write(filePath(i));
            cout << "Wrote " << selectedLUT.size() << " points to " << filePath(i) << endl;
        }
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
    glClearColor(0, 0, 0, 1);
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
        cerr << "Usage: data_viewer dataTable [pattern]" << endl;
        exit(-1);
    }

    // Initialize GLUT
    glutInit(&argc, argv);

    int width = 1024, height = 768;
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("Lookup Table Viewer");

    // Set GLUT event callbacks
    glutMouseFunc(MouseFunc);
    glutMotionFunc(MotionFunc);
    glutPassiveMotionFunc(PassiveMotionFunc);
    glutKeyboardFunc(KeyboardFunc);
    glutSpecialFunc(SpecialKeyboardFunc);

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);

    try {
    g_lut.load(argv[1]);
    cout << "Patterns:";
    for (size_t pat : g_lut.patterns())
        cout << " " <<pat;
    cout << endl;
    if (argc < 3) {
        cout << "Selecting " << *g_lut.patterns().begin() << endl;
        g_lut = g_lut.selectPattern(*g_lut.patterns().begin());
    }
    else {
        g_lut = g_lut.selectPattern(stoi(argv[2]));
    }

    auto mainLayout = make_shared<Layout>(Layout::Direction::Horizontal);
    g_mainView = mainLayout;
    g_logscalePlot = make_shared<LogscalePlot>();
    g_logscalePlot->setData(g_lut.getNus<float>(), g_lut.getEs<float>());
    g_logscalePlot->setBackgroundColor(1, 1, 1);
    g_logscalePlot->setForegroundColor(0.0,0.0,0.0);
    g_logscalePlot->setSelectionColor(1.0,0.0,0.0);
    g_logscalePlot->setPointSize(2.0);
    g_logscalePlot->autoPlotRange();
    g_logscalePlot->setSelectionChangedCallback(
        [&](ViewPtr v) {
            auto lv = dynamic_pointer_cast<LogscalePlot>(v);
            assert(lv);
            auto selSubset = g_lut.subset(lv->selectedPointIndices());
            size_t nParams = g_lut.numParams();
            for (size_t p = 0; p < nParams; ++p) {
                if (selSubset.size()) g_selectedParamHistograms[p]->setData(selSubset.getParamValues<float>(p), g_fullParamHistograms[p]->dataRange());
                else g_selectedParamHistograms[p]->setData(std::vector<float>());
            }
        }
    );

    mainLayout->addSubview(g_logscalePlot);
    auto histLayout = make_shared<Layout>(Layout::Direction::Vertical);
    mainLayout->addSubview(histLayout);

    size_t nParams = g_lut.numParams();
    for (size_t p = 0; p < nParams; ++p) {
        auto paramHistFull = make_shared<HistogramPlot>(g_lut.getParamValues<float>(p));
        auto paramHistSelected = make_shared<HistogramPlot>();
        paramHistFull->setBackgroundColor(1, 1, 1);
        paramHistFull->setForegroundColor(0.75,0.75,0.75);
        paramHistSelected->setBackgroundColor(0, 0, 0, 0);
        paramHistSelected->setForegroundColor(1.0,0.0,0.0);

        auto stacker = make_shared<Layout>(Layout::Direction::Stacked);
        stacker->addSubview(paramHistFull);
        stacker->addSubview(paramHistSelected);
        histLayout->addSubview(stacker);

        g_selectedParamHistograms.push_back(paramHistSelected);
        g_fullParamHistograms.push_back(paramHistFull);
    }

    mainLayout->setSubviewPercentages({0.75, 0.25});

    g_uiHandler = make_shared<UIHandler>(g_mainView);
    g_uiHandler->setFirstResponder(g_logscalePlot);

    // Call the GLUT main loop
    glutMainLoop();
    }
    catch (const exception &e) {
        cerr << "Caught exception: " << e.what() << endl;
    }

    return 0;
}
