#include "IsosurfaceInflator.hh"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <MeshIO.hh>
#include <MSHFieldWriter.hh>
#include <SimplicialMesh.hh>
#include <Future.hh>
#include <filters/remove_dangling_vertices.hh>
#include <GlobalBenchmark.hh>

#include "AutomaticDifferentiation.hh"
#include "WireMesh.hh"
#include "PatternSignedDistance.hh"
#include "SnapAndReflect.hh"
#include "Isometries.hh"

#define DEBUG_EVALPTS 0

#if 1
#include "CGALClippedVolumeMesher.hh"
#include "VCGSurfaceMesher.hh"
#endif

#include "BoxIntersectionMesher.hh"
#include "MidplaneMesher.hh"

#include "IsosurfaceInflatorConfig.hh"

#include "parallel_for.h"

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
                 bool                      needsReflecting,
                 const BBox<Point>         &meshCell,
                 const MeshingOptions      &opts);

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

        // std::cout << "Meshed params:";
        // for (Real p : params) std::cout << "\t" << p;
        // std::cout << std::endl;

        // Determine if meshed domain is 2D or 3D and postprocess accordingly
        BBox<Point> bbox(vertices);
        if (std::abs(bbox.dimensions()[2]) < 1e-8) postProcess<2>(vertices, elements, normalShapeVelocities, vertexNormals, *this, needsReflecting, meshingCell(), meshingOptions());
        else                                       postProcess<3>(vertices, elements, normalShapeVelocities, vertexNormals, *this, needsReflecting, meshingCell(), meshingOptions());
    }

    // Mesh the param (fills vertices, elements member arrays)
    virtual void meshPattern(const vector<Real> &params) = 0;

    // Derivative of signed distance function with respect to evaluation point
    virtual vector<Point> signedDistanceGradient(const vector<Point> &evalPoints) const = 0;

    // Derivative of signed distance function wrt each pattern parameter
    virtual vector<vector<Real>> signedDistanceParamPartials(const vector<Point> &evalPoints) const = 0;

    // Whether the inflator generates only the orthotropic symmetry base cell
    // (by default--reflectiveInflator will override this)
    virtual bool _mesherGeneratesOrthoCell() const = 0;

    // The meshing cell box.
    virtual BBox<Point> meshingCell() const = 0;

    void clear() {
        vertices.clear(), elements.clear();
        normalShapeVelocities.clear();
        vertexNormals.clear();
        inflatedParams.clear();
    }

    virtual ~Impl() { }

    ////////////////////////////////////////////////////////////////////////////
    // "Private" data memebers
    ////////////////////////////////////////////////////////////////////////////
    // Chooses whether the full period cell is created in cases where the
    // representative cell is smaller (e.g. orthotropic patterns, whose base
    // cell is 1/8th the period cell).
    bool generateFullPeriodCell = true;
    // Chooses whether we mesh only the symmetry base cell and reflect to the
    // full period cell if necessary (default), or we mesh the full period cell
    // when generateFullPeriodCell == true
    bool reflectiveInflator = true;
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

        BENCHMARK_START_TIMER_SECTION("meshPattern");
        // Optional debugging graph output.
        const auto &config = IsosurfaceInflatorConfig::get();
        if (config.dumpInflationGraph())  { wmesh.saveInflationGraph(    config.inflationGraphPath, params); }
        if (config.dumpReplicatedGraph()) { wmesh.saveReplicatedBaseUnit(config.replicatedGraphPath); }
        if (config.dumpBaseUnitGraph())   { wmesh.saveBaseUnit(          config.baseUnitGraphPath); }

        pattern.setParameters(params);

        // Change the pattern's meshing domain if we're forcing meshing of the
        // full TriplyPeriodic base cell.
        if (generateFullPeriodCell && !reflectiveInflator)
            pattern.setBoundingBox(Symmetry::TriplyPeriodic<>::representativeMeshCell<Real>());

        mesher.mesh(pattern, this->vertices, this->elements);
        BENCHMARK_STOP_TIMER_SECTION("meshPattern");
        // cout << vertices.size() << " vertices, " << elements.size() << " elements" << endl;
    }

    // Note: when reflectiveInflator = false, the mesher generates the full
    // period cell.
    virtual bool _mesherGeneratesOrthoCell() const {
        return is_base_of<Symmetry::Orthotropic<typename PatternSymmetry::Tolerance>,
                          PatternSymmetry>::value
               && (reflectiveInflator || !generateFullPeriodCell);
    }

    // The meshing cell box.
    virtual BBox<Point> meshingCell() const { return PatternSymmetry::template representativeMeshCell<Real>(); }

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
        const size_t nEvals = evalPoints.size(), nParams = pattern.numParams();
        vector<vector<Real>> partials(nParams, vector<Real>(nEvals));

        // Brute force version: rebuild full autodiff pattern signed distance point
        // for each evaluation point. The obvious optimization is to
        // only rerun the signed distance evaluation for each evaluation, but
        // this may be difficult with Adept and probably gives negligible
        // speedup.
        // auto paramPartialLoop = [&partials, &evalPoints, nParams, this](size_t e) {
        using adept::adouble;
        adept::Stack stack;
        for (size_t e = 0; e <= nEvals; ++ e) {
            // assert(!stack.is_thread_unsafe());
            PatternSignedDistance<adouble, WMesh> patternAutodiff(wmesh);

            vector<adouble> params(inflatedParams.size());
            for (size_t i = 0; i < params.size(); ++i)
                params[i] = inflatedParams[i];

            stack.new_recording();
            patternAutodiff.setParameters(params);
            adouble sd = patternAutodiff.signedDistance(evalPoints[e].template cast<adouble>().eval());
            sd.set_gradient(1.0);
            stack.reverse();
            for (size_t p = 0; p < nParams; ++p) {
                partials[p][e] = params[p].get_gradient();
                if (std::isnan(partials[p][e])) {
                    std::cerr << "nan sd partial " << p << " at evaluation point "
                              << evalPoints[e] << std::endl;
                    std::cerr << "sd at pt:\t" << pattern.signedDistance(evalPoints[e]) << std::endl;
                    std::cerr << "autodiff sd at pt:\t" << sd << std::endl;
                    // throw std::runtime_error("nan sd");
                }
            }
        };
        // igl::parallel_for(nEvals, paramPartialLoop);

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
void IsosurfaceInflator::setReflectiveInflator(bool onoff) { m_imp->reflectiveInflator = onoff; }
BaseCellType IsosurfaceInflator::baseCellType() const {
    if      (m_imp->generateFullPeriodCell     ) return BaseCellType::TriplyPeriodic;
    else if (m_imp->_mesherGeneratesOrthoCell()) return BaseCellType::Orthotropic;
    else throw std::runtime_error("Unknown/incompatible base cell type.");
}

const vector<MeshIO::IOVertex > &IsosurfaceInflator::vertices()              const { return m_imp->vertices; }
const vector<MeshIO::IOElement> &IsosurfaceInflator::elements()              const { return m_imp->elements; }
const vector<vector<Real>>      &IsosurfaceInflator::normalShapeVelocities() const { return m_imp->normalShapeVelocities; }

void                             IsosurfaceInflator::clear()                       { m_imp->clear(); }

auto IsosurfaceInflator::vertexNormals() const -> const vector<Point> & { return m_imp->vertexNormals; }

void IsosurfaceInflator::disablePostprocess() { m_imp->m_noPostprocess = true; }
void IsosurfaceInflator:: enablePostprocess() { m_imp->m_noPostprocess = false; }

MeshingOptions &IsosurfaceInflator::meshingOptions() { return m_imp->meshingOptions(); }

IsosurfaceInflator::IsosurfaceInflator(const string &type, bool vertexThickness, const string &wireMeshPath) {
    if (!vertexThickness) throw runtime_error("Only per-vertex thickness is currently supported.");
    if (type == "cubic")                                                          m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>         , CGALClippedVolumeMesher>(wireMeshPath);
    else if (type == "orthotropic")                                               m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>   , CGALClippedVolumeMesher>(wireMeshPath);
    else if (type == "triply_periodic")                                           m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>, CGALClippedVolumeMesher>(wireMeshPath);
    else if (type == "cubic_preview")   {
        // throw std::runtime_error("Disabled");
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>, VCGSurfaceMesher>(wireMeshPath);
        // Triangle mesh doesn't support our post-processing
        disablePostprocess();
    }
    else if (type == "cubic_features")   {
        // Output the sharp 1D features created by bounding box intersection.
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>, BoxIntersectionMesher>(wireMeshPath);
        // Line mesh doesn't support our post-processing
        disablePostprocess();
    }
    else if (type == "orthotropic_features")   {
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>, BoxIntersectionMesher>(wireMeshPath);
        disablePostprocess();
    }
    else if (type == "triply_periodic_features")   {
        // Output the sharp 1D features created by bounding box intersection.
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>, BoxIntersectionMesher>(wireMeshPath);
        // Line mesh doesn't support our post-processing
        disablePostprocess();
    }
    else if (type == "triply_periodic_preview")   {
        // Output the sharp 1D features created by bounding box intersection.
        m_imp = new IsosurfaceInflatorImpl<WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>, VCGSurfaceMesher>(wireMeshPath);
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
                 bool                      needsReflecting,
                 const BBox<Point>         &meshCell,
                 const MeshingOptions      &opts)
{
    BENCHMARK_START_TIMER_SECTION("postProcess");
    {
        const size_t origSize = vertices.size();
        remove_dangling_vertices(vertices, elements);
        if (vertices.size() != origSize)
            std::cerr << "WARNING: removed " << origSize - vertices.size() << " dangling vertices." << std::endl;
    }

    // MeshIO::save("pre_snap.msh", vertices, elements);

    BENCHMARK_START_TIMER("SnapAndDetermineEvaluationPts");
    vector<vector<bool>> onMinFace, onMaxFace;

    using SMesh = SimplicialMesh_t<N>;
    std::unique_ptr<SMesh> bcm;
    try {
        bcm = Future::make_unique<SMesh>(elements, vertices.size());
    }
    catch (...) {
        std::cerr << "Exception while building mesh" << std::endl;
        std::cerr << "Dumping debug.msh" << std::endl;
        MeshIO::save("debug.msh", vertices, elements);
        throw;
    }
    SMesh &symBaseCellMesh = *bcm;

    try {
        if (N == 3) smartSnap3D(vertices, symBaseCellMesh, meshCell);
    }
    catch (std::exception &e) {
        // For bad meshes caused, e.g, by self intersections, the feature 
        std::cerr << "WARNING: smartSnap3D failed--probably nonsmooth geometry at cell interface. Resorting to dumbSnap3D." << std::endl;
        std::cerr << "(" << e.what() << ")" << std::endl;
        dumbSnap3D(vertices, meshCell, opts.facetDistance);
    }
    BBox<Point3D> snappedBB(vertices);
    if (snappedBB != meshCell) {
        std::cerr << "Failed to snap mesh. Dumping debug.msh" << std::endl;
        MeshIO::save("debug.msh", vertices, elements);
        throw std::runtime_error("Snapping failed.");
    }

    using FM = PeriodicBoundaryMatcher::FaceMembership<3>;
    onMinFace.assign(3, std::vector<bool>(vertices.size(), false));
    onMaxFace.assign(3, std::vector<bool>(vertices.size(), false));
    std::vector<FM> faceMembership;
    for (size_t i = 0; i < vertices.size(); ++i) {
        FM fm(vertices[i], meshCell, 0); // zero tolerance!
        onMinFace[0][i] = fm.onMinFace(0), onMaxFace[0][i] = fm.onMaxFace(0);
        onMinFace[1][i] = fm.onMinFace(1), onMaxFace[1][i] = fm.onMaxFace(1);
        onMinFace[2][i] = fm.onMinFace(2), onMaxFace[2][i] = fm.onMaxFace(2);
    }

    // MeshIO::save("post_snap.msh", vertices, elements);

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
            for (auto bv : bs.vertices()) {
                onMin &= onMinFace[d].at(bv.volumeVertex().index());
                onMax &= onMaxFace[d].at(bv.volumeVertex().index());
            }
            isCellFace |= (onMin | onMax);
        }
        if (isCellFace) continue;

        for (auto bv : bs.vertices())
            internalCellFaceVertex.at(bv.index()) = false;
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

#if DEBUG_EVALPTS
    {
        MSHFieldWriter debug("debug_evalpts.msh", vertices, elements);
        for (size_t d = 0; d < N; ++d) {
            ScalarField<Real> minFaceIndicator(vertices.size()), maxFaceIndicator(vertices.size());
            for (size_t i = 0; i < vertices.size(); ++i) {
                minFaceIndicator[i] = onMinFace[d].at(i) ? 1.0 : 0.0;
                maxFaceIndicator[i] = onMaxFace[d].at(i) ? 1.0 : 0.0;
            }
            debug.addField("onMinFace[" + std::to_string(d) + "]", minFaceIndicator);
            debug.addField("onMaxFace[" + std::to_string(d) + "]", maxFaceIndicator);
        }
        ScalarField<Real> evalPtIndicator(vertices.size());
        ScalarField<Real> icfv(vertices.size());
        evalPtIndicator.clear();
        icfv.clear();
        for (auto bv : symBaseCellMesh.boundaryVertices()) {
            if (evalPointIndex[bv.index()] < evaluationPoints.size())
                evalPtIndicator[bv.volumeVertex().index()] = 1.0;
            icfv[bv.volumeVertex().index()] = internalCellFaceVertex.at(bv.index()) ? 1.0 : 0.0;
        }
        debug.addField("isEvalPoint", evalPtIndicator);
        debug.addField("internalCellFaceVertex", icfv);
    }
#endif

    BENCHMARK_STOP_TIMER("SnapAndDetermineEvaluationPts");

    BENCHMARK_START_TIMER("SignedDistanceGradientsAndPartials");
    // sd(x, p) = 0
    // grad_x(sd) . dx/dp + d(sd)/dp = 0,   grad_x(sd) = n |grad_x(sd)|
    // ==>  v . n = -[ d(sd)/dp ] / |grad_x(sd)|
    vector<vector<Real>> vnp;
    vector<Point> sdGradX;
    // try {
        vnp = inflator.signedDistanceParamPartials(evaluationPoints);
        sdGradX = inflator.signedDistanceGradient(evaluationPoints);
    // }
    // catch(...) {
    //     MSHFieldWriter debug("debug.msh", vertices, elements);
    //     BENCHMARK_STOP_TIMER_SECTION("SignedDistanceGradientsAndPartials");
    //     BENCHMARK_STOP_TIMER_SECTION("postProcess");
    //     throw;
    // }
    vector<Real> sdGradNorms(evaluationPoints.size());
    for (size_t i = 0; i < evaluationPoints.size(); ++i) {
        sdGradNorms[i] = sdGradX[i].norm();
        // We evaluate on the boundary--there should be a well-defined normal
        if (std::abs(sdGradNorms[i]) < 1e-8) {
            std::cerr << "Undefined normal at evaluation point "
                      << evaluationPoints[i] << std::endl;
            BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");
            BENCHMARK_STOP_TIMER_SECTION("postProcess");
            throw std::runtime_error("Normal undefined.");
        }
    }

    for (auto &vn : vnp) {
        for (size_t i = 0; i < vn.size(); ++i) {
            vn[i] *= -1.0 / sdGradNorms[i];
            if (std::isnan(vn[i])) {
                ScalarField<Real> fail(vertices.size());
                fail.clear();
                for (const auto bv : symBaseCellMesh.boundaryVertices()) {
                    size_t e = evalPointIndex.at(bv.index());
                    if (e < evaluationPoints.size()) {
                        if (std::isnan(vn.at(e))) fail[bv.volumeVertex().index()] = 1.0;
                    }
                }

                MSHFieldWriter debug("debug.msh", vertices, elements);
                debug.addField("fail", fail);
                BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");
                BENCHMARK_STOP_TIMER_SECTION("postProcess");
                throw std::runtime_error("nan vn");
                // assert(false);
            }
        }
    }
    BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");

    BENCHMARK_START_TIMER("Reflecting");
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
    BENCHMARK_STOP_TIMER("Reflecting");

    // static int _run_num = 0;
    // MeshIO::save("debug_" + std::to_string(_run_num++) + ".msh", vertices, elements);
    BENCHMARK_STOP_TIMER_SECTION("postProcess");
}
