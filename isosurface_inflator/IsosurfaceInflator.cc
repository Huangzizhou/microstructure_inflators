#include "IsosurfaceInflator.hh"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <MeshIO.hh>
#include <SimplicialMesh.hh>
#include "AutomaticDifferentiation.hh"
#include "WireMesh.hh"
#include "PatternSignedDistance.hh"
#include "SnapAndReflect.hh"
#include "Isometries.hh"

#if 0
#include "CGALClippedVolumeMesher.hh"
#include "VCGSurfaceMesher.hh"
#endif

#include "BoxIntersectionMesher.hh"
#include "MidplaneMesher.hh"

#include "IsosurfaceInflatorConfig.hh"

using namespace std;

template<class WMesh, template<class> class Mesher>
class IsosurfaceInflatorImpl;

// Postprocess:
//    Snap to base cell and then reflect if necessary
//    Compute vertex normals and normal shape velocities
template<size_t N, class Point>
void postProcess(vector<MeshIO::IOVertex>  &vertices,
                 vector<MeshIO::IOElement> &elements,
                 vector<vector<Real>>      &normalShapeVelocities,
                 vector<Point>             &vertexNormals,
                 const IsosurfaceInflator::Impl &inflator,
                 bool needsReflecting);

class IsosurfaceInflator::Impl {
public:
    virtual vector<Real>  defaultParameters() const = 0;
    virtual size_t                numParams() const = 0;
    virtual bool   isThicknessParam(size_t p) const = 0;
    virtual bool    isPositionParam(size_t p) const = 0;
    virtual bool    isBlendingParam(size_t p) const = 0;

    virtual MeshingOptions &meshingOptions() = 0;

    // Delegates to derived IsosurfaceInflatorImpl for WMesh-dependent stuff
    // (via the virtual functions below).
    void inflate(const vector<Real> &params) {
        inflatedParams = params;
        meshPattern(params);
        // Terminate after initial meshing, for debugging.
        if (m_noPostprocess) return;
        
        bool needsReflecting = generateFullPeriodCell && _mesherGeneratesOrthoCell();

        // Determine if meshed domain is 2D or 3D and postprocess accordingly
        BBox<Point> bbox(vertices);
        if (std::abs(bbox.dimensions()[2]) < 1e-8) postProcess<2>(vertices, elements, normalShapeVelocities, vertexNormals, *this, needsReflecting);
        else                                       postProcess<3>(vertices, elements, normalShapeVelocities, vertexNormals, *this, needsReflecting);
    }

    // Mesh the param (fills vertices, elements member arrays)
    virtual void meshPattern(const vector<Real> &params) = 0;

    // Derivative of signed distance function with respect to evaluation point
    virtual vector<Point> signedDistanceGradient(const vector<Point> &evalPoints) const = 0;

    // Derivative of signed distance function wrt each pattern parameter
    virtual vector<vector<Real>> signedDistanceParamPartials(const vector<Point> &evalPoints) const = 0;

    // Whether the inflator generates only the orthotropic symmetry base cell
    virtual bool _mesherGeneratesOrthoCell() const = 0;

    virtual ~Impl() { }

    ////////////////////////////////////////////////////////////////////////////
    // "Private" data memebers
    ////////////////////////////////////////////////////////////////////////////
    // Chooses whether the full period cell is created in cases where the
    // representative cell is smaller (e.g. orthotropic patterns, whose base
    // cell is 1/8th the period cell).
    bool generateFullPeriodCell = true;
    bool m_noPostprocess = false; // for debugging
    vector<MeshIO::IOVertex>  vertices;
    vector<MeshIO::IOElement> elements;
    vector<vector<Real>>      normalShapeVelocities;
    vector<Point>             vertexNormals;

    // The params that were most recently inflated (to which
    // (vertices, elements) correspond).
    vector<Real> inflatedParams;
};

// The WMesh-dependent implementation details.
// E.g.: WMesh = WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>,
//       Mesher = CGALClippedVolumeMesher 
template<class WMesh, template<class> class Mesher>
class IsosurfaceInflatorImpl : public IsosurfaceInflator::Impl {
public:
    using Point = IsosurfaceInflator::Point;
    typedef PatternSignedDistance<Real, WMesh> PSD;
    typedef typename WMesh::PatternSymmetry PatternSymmetry;
    IsosurfaceInflatorImpl(const string &wireMeshPath)
        : wmesh(wireMeshPath), pattern(wmesh) { }

    virtual void meshPattern(const vector<Real> &params) {
        // std::cout << "Meshing parameters:";
        // for (auto p : params)
        //     std::cout << "\t" << p;
        // std::cout << std::endl;

        // Optional debugging graph output.
        const auto &config = IsosurfaceInflatorConfig::get();
        if (config.dumpInflationGraph()) { wmesh.saveInflationGraph(config.inflationGraphPath, params); }
        if (config.dumpReplicatedGraph()) { wmesh.saveReplicatedBaseUnit(config.replicatedGraphPath); }

        pattern.setParameters(params);
        mesher.mesh(pattern, this->vertices, this->elements);
        cout << vertices.size() << " vertices, " << elements.size() << " elements" << endl;
    }

    virtual bool _mesherGeneratesOrthoCell() const {
        return is_base_of<Symmetry::Orthotropic<typename PatternSymmetry::Tolerance>,
                          PatternSymmetry>::value;
    }

    // Derivative of signed distance function with respect to evaluation point
    // (autodiff-based).
    virtual vector<Point> signedDistanceGradient(const vector<Point> &evalPoints) const {
        using adept::adouble;
        typedef Point3<adouble> Pt;
        adept::Stack stack;

        vector<Point> distGradX(evalPoints.size());
        for (size_t p = 0; p < evalPoints.size(); ++p) {
            Pt x(evalPoints[p].template cast<adouble>());
            stack.new_recording();
            adouble dist = pattern.signedDistance(x);
            dist.set_gradient(1);
            stack.reverse();
            for (size_t i = 0; i < 3; ++i)
                 distGradX[p][i] =  x[i].get_gradient();
        }
        return distGradX;
    }

    // Derivative of signed distance function with respect to each pattern
    // parameter (autodiff-based).
    virtual vector<vector<Real>> signedDistanceParamPartials(const vector<Point> &evalPoints) const {
        size_t nEvals = evalPoints.size(), nParams = pattern.numParams();
        vector<vector<Real>> partials(nParams, vector<Real>(nEvals));

        // Brute force version: rebuild full autodiff pattern signed distance point
        // for each evaluation point. The obvious optimization is to
        // only rerun the signed distance evaluation for each evaluation, but
        // this may be difficult with Adept and probably gives negligible
        // speedup.
        for (size_t e = 0; e < nEvals; ++e) {
            using adept::adouble;
            adept::Stack stack;
            PatternSignedDistance<adouble, WMesh> patternAutodiff(wmesh);

            vector<adouble> params(inflatedParams.size());
            for (size_t i = 0; i < params.size(); ++i)
                params[i] = inflatedParams[i];

            stack.new_recording();
            patternAutodiff.setParameters(params);
            adouble sd = patternAutodiff.signedDistance(evalPoints[e].template cast<adouble>().eval());
            sd.set_gradient(1.0);
            stack.reverse();
            for (size_t p = 0; p < nParams; ++p)
                partials[p][e] = params[p].get_gradient();
        }

        return partials;
    }

    virtual ~IsosurfaceInflatorImpl() { }

    virtual vector<Real> defaultParameters() const { return wmesh.defaultParameters(); }
    virtual size_t               numParams() const { return wmesh.numParams(); }
    virtual bool  isThicknessParam(size_t p) const { return wmesh.isThicknessParam(p); }
    virtual bool   isPositionParam(size_t p) const { return wmesh.isPositionParam(p); }
    virtual bool   isBlendingParam(size_t p) const { return wmesh.isBlendingParam(p); }

    virtual MeshingOptions &meshingOptions() { return mesher.meshingOptions; }

    WMesh wmesh;
    PSD pattern;
    Mesher<PSD> mesher;
};

void IsosurfaceInflator::inflate(const vector<Real> &params) { m_imp->inflate(params); }

vector<Real> IsosurfaceInflator::defaultParameters()        const { return m_imp->defaultParameters(); }
size_t       IsosurfaceInflator::numParams()                const { return m_imp->numParams(); }
bool         IsosurfaceInflator::isThicknessParam(size_t p) const { return m_imp->isThicknessParam(p); }
bool         IsosurfaceInflator:: isPositionParam(size_t p) const { return m_imp->isPositionParam(p); }
bool         IsosurfaceInflator:: isBlendingParam(size_t p) const { return m_imp->isBlendingParam(p); }

void IsosurfaceInflator::setGenerateFullPeriodCell(bool onoff) { m_imp->generateFullPeriodCell = onoff; }
bool IsosurfaceInflator::shouldGenerateFullPeriodCell() const { return m_imp->generateFullPeriodCell; }

const vector<MeshIO::IOVertex > &IsosurfaceInflator::vertices()              const { return m_imp->vertices; }
const vector<MeshIO::IOElement> &IsosurfaceInflator::elements()              const { return m_imp->elements; }
const vector<vector<Real>>      &IsosurfaceInflator::normalShapeVelocities() const { return m_imp->normalShapeVelocities; }

auto IsosurfaceInflator::vertexNormals() const -> const vector<Point> & { return m_imp->vertexNormals; }

void IsosurfaceInflator::disablePostprocess() { m_imp->m_noPostprocess = true; }
void IsosurfaceInflator:: enablePostprocess() { m_imp->m_noPostprocess = false; }

MeshingOptions &IsosurfaceInflator::meshingOptions() { return m_imp->meshingOptions(); }

IsosurfaceInflator::IsosurfaceInflator(const string &type, bool vertexThickenss, const string &wireMeshPath) {
    if (!vertexThickenss) throw runtime_error("Only per-vertex thickness is currently supported.");
    if (type == "cubic")                throw std::runtime_error("Disabled."); /* m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>         , CGALClippedVolumeMesher>(wireMeshPath); */
    else if (type == "orthotropic")     throw std::runtime_error("Disabled."); /* m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>   , CGALClippedVolumeMesher>(wireMeshPath); */
    else if (type == "triply_periodic") throw std::runtime_error("Disabled."); /* m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>, CGALClippedVolumeMesher>(wireMeshPath); */
    else if (type == "cubic_preview")   {
        throw std::runtime_error("Disabled");
        // m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>, VCGSurfaceMesher>(wireMeshPath);
        // Triangle mesh doesn't support our post-processing
        disablePostprocess();
    }
    else if (type == "cubic_features")   {
        // Output the sharp 1D features created by bounding box intersection.
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>, BoxIntersectionMesher>(wireMeshPath);
        // Line mesh doesn't support our post-processing
        disablePostprocess();
    }
    else if (type == "2D_square") {
        throw std::runtime_error("2D square symmetry unimplemented."); 
        // m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>, MidplaneMesher>(wireMeshPath);
        // Line mesh doesn't support our post-processing
        // disablePostprocess();
    }
    else if (type == "2D_orthotropic") {
        // throw std::runtime_error("Disabled."); 
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>, MidplaneMesher>(wireMeshPath);
    }
    else throw runtime_error("Unknown inflator type: " + type);
}

IsosurfaceInflator::~IsosurfaceInflator() {
    delete m_imp;
}

// Postprocess:
//    Snap to base cell and then reflect if necessary
//    Compute vertex normals and normal shape velocities
template<size_t N, class Point>
void postProcess(vector<MeshIO::IOVertex>  &vertices,
                 vector<MeshIO::IOElement> &elements,
                 vector<vector<Real>>      &normalShapeVelocities,
                 vector<Point>             &vertexNormals,
                 const IsosurfaceInflator::Impl &inflator,
                 bool needsReflecting)
{
    vector<vector<bool>> onMinFace, onMaxFace;
    // WARNING: for non-reflecting inflators, this should snap to bbmin and bbmax!!!
    // TODO: change to pass the meshing cell
    snapVerticesToUnitCell<>(vertices, onMinFace, onMaxFace);

    SimplicialMesh<N> symBaseCellMesh(elements, vertices.size());
    // Mark internal cell-face vertices: vertices on the meshing cell
    // boundary that actually lie inside the object (i.e. they are only mesh
    // boundary vertices because of the intersection of the periodic pattern with
    // the meshing box).
    // This this is not the case if any non-cell-face triangle is incident
    vector<bool>internalCellFaceVertex(symBaseCellMesh.numBoundaryVertices(), true);
    for (auto bs : symBaseCellMesh.boundarySimplices()) {
        bool isCellFace = false;
        for (size_t d = 0; d < N; ++d) {
            bool onMin = true, onMax = true;
            for (size_t bvi = 0; bvi < bs.numVertices(); ++bvi) {
                onMin &= onMinFace[d].at(bs.vertex(bvi).volumeVertex().index());
                onMax &= onMinFace[d].at(bs.vertex(bvi).volumeVertex().index());
            }
            isCellFace |= (onMin | onMax);
        }
        if (isCellFace) continue;

        internalCellFaceVertex.at(bs.vertex(0).index()) = false;
        internalCellFaceVertex.at(bs.vertex(1).index()) = false;
        internalCellFaceVertex.at(bs.vertex(2).index()) = false;
    }

    // Compute parameter shape velocities on all (true) boundary vertices
    vector<Point> evaluationPoints;
    // Tie the boundary vertices to the associated data evaluation point
    vector<size_t> evalPointIndex(symBaseCellMesh.numBoundaryVertices(),
                                  numeric_limits<size_t>::max());
    for (auto bv : symBaseCellMesh.boundaryVertices()) {
        if (internalCellFaceVertex.at(bv.index())) continue;
        evalPointIndex[bv.index()] = evaluationPoints.size();
        evaluationPoints.push_back(vertices.at(bv.volumeVertex().index()));
    }

    // sd(x, p) = 0
    // grad_x(sd) . dx/dp + d(sd)/dp = 0,   grad_x(sd) = n |grad_x(sd)|
    // ==>  v . n = -[ d(sd)/dp ] / |grad_x(sd)|
    vector<vector<Real>> vnp = inflator.signedDistanceParamPartials(evaluationPoints);
    vector<Point>    sdGradX = inflator.signedDistanceGradient(evaluationPoints);
    vector<Real> sdGradNorms(evaluationPoints.size());
    for (size_t i = 0; i < evaluationPoints.size(); ++i)
        sdGradNorms[i] = sdGradX[i].norm();

    for (auto &vn : vnp) {
        for (size_t i = 0; i < vn.size(); ++i)
            vn[i] *= -1.0 / sdGradNorms[i];
    }
    
    // If the mesher only generates the orthotropic base cell, the mesh must
    // be reflected to get the full period cell (if requested).
    if (needsReflecting) {
        vector<MeshIO::IOVertex>  reflectedVertices;
        vector<MeshIO::IOElement> reflectedElements;
        vector<size_t>   vertexOrigin;
        vector<Isometry> vertexIsometry;
        reflectXYZ(N, vertices, elements, onMinFace,
                   reflectedVertices, reflectedElements,
                   vertexOrigin, vertexIsometry);

        // Copy over normal velocity and vertex normal data.
        // So that we don't rely on consistent boundary vertex numbering
        // (which is nevertheless guaranteed by our mesh datastructure),
        // we store this data per-vertex (setting the data on internal
        // vertices to 0).
        normalShapeVelocities.assign(vnp.size(), vector<Real>(reflectedVertices.size()));
        vertexNormals.assign(reflectedVertices.size(), Point::Zero());

        for (size_t i = 0; i < reflectedVertices.size(); ++i) {
            auto v = symBaseCellMesh.vertex(vertexOrigin[i]);
            assert(v);
            auto bv = v.boundaryVertex();
            if (!bv || internalCellFaceVertex[bv.index()]) continue;
            size_t evalIdx = evalPointIndex.at(bv.index());
            assert(evalIdx < evaluationPoints.size());
            // Transform normals by the reflection isometry and normalize
            vertexNormals[i] = vertexIsometry[i].apply(sdGradX[evalIdx]) / sdGradNorms[evalIdx];
            for (size_t p = 0; p < vnp.size(); ++p)
                normalShapeVelocities[p][i] = vnp[p][evalIdx];
        }

        reflectedVertices.swap(vertices);
        reflectedElements.swap(elements);
    }
    else {
        // Otherwise, we still need to copy over the normal shape velocities
        normalShapeVelocities.assign(vnp.size(), vector<Real>(vertices.size()));
        vertexNormals.assign(vertices.size(), Point::Zero());
        for (auto bv : symBaseCellMesh.boundaryVertices()) {
            size_t evalIdx = evalPointIndex[bv.index()];
            if (evalIdx >= evalPointIndex.size()) continue;
            size_t vi = bv.volumeVertex().index();
            vertexNormals[vi] = sdGradX[evalIdx] / sdGradNorms[evalIdx];
            for (size_t p = 0; p < vnp.size(); ++p)
                normalShapeVelocities[p][vi] = vnp[p][evalIdx];
        }
    }

    // static int _run_num = 0;
    // MeshIO::save("debug_" + std::to_string(_run_num++) + ".msh", vertices, elements);
}
