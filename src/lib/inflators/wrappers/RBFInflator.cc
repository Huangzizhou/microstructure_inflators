#include "RBFInflator.hh"

using namespace std;

RBFInflator::RBFInflator(Real epsilon, size_t dim) {
    m_dim = dim;
    m_epsilon = epsilon;

    // dim1 relates to x, so number of columns
    // dim2 relates to y, so number of rows
    m_coeffMatrix.resize(dim);
    for (size_t i = 0; i < dim; i++) {
        m_coeffMatrix[i].resize(dim);
    }

    Point2D min = Point2D(-1.0, -1.0);
    Point2D max = Point2D( 1.0,  1.0);

    m_bbox = BBox<Point2D>(min, max);
}

RBFInflator::RBFInflator(std::string png_path, Real epsilon, size_t dim) {
    m_dim = dim;
    m_epsilon = epsilon;

    RBF<Real> levelSet = RBF<Real>(png_path, epsilon, dim);

    m_coeffMatrix = levelSet.coefficients();

    Point2D min = Point2D(-1.0, -1.0);
    Point2D max = Point2D( 1.0,  1.0);

    m_bbox = BBox<Point2D>(min, max);
}

// VoxelsInflator methods
void
RBFInflator::m_inflate(const std::vector<Real> &params) {
    // Makes sure we use the last meshing options set
    m_mesher.meshingOptions = m_meshingOptions;

    // Makes sure we use periodic mesher
    m_mesher.meshInterfaceConsistently = true;

    std::vector<std::vector<Real>> coeffMatrix = vecToMat(params, m_dim, m_dim);

    RBF<Real> levelSet = RBF<Real>(coeffMatrix, m_epsilon);

    m_mesher.meshSlice(levelSet, m_vertices, m_elements);
}

std::vector<VectorField<Real, 2>>
RBFInflator::volumeShapeVelocities() const {
    std::vector<VectorField<Real, N>> result(numParameters());
    SimplicialMesh<N> mesh = SimplicialMesh<N>(m_elements, m_vertices.size());

    RBF<Real> levelSet = RBF<Real>(m_coeffMatrix, m_epsilon);

    std::vector<Vector2D> gradients(mesh.numBoundaryVertices());
    std::vector<Real> gradientNorms(mesh.numBoundaryVertices());
    std::vector<Vector2D> normals(mesh.numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices()) {
        auto vv = bv.volumeVertex();
        Vector2D p = m_vertices[vv.index()];
        Vector2D gradient = levelSet.gradient(p);
        Real gradientNorm = gradient.norm();
        Vector2D normal = gradient / gradientNorm;

        gradients[bv.index()] = gradient;
        gradientNorms[bv.index()] = gradientNorm;
        normals[bv.index()] = normal;
    }

    for (unsigned param = 0; param < numParameters(); param++) {
        VectorField<Real, N> velocityForP(mesh.numVertices());
        velocityForP.clear();

        // Param to dimension indices
        size_t i = param / m_dim;
        size_t j = param -  m_dim * i;

        for (auto bv : mesh.boundaryVertices()) {
            auto vv = bv.volumeVertex();
            Vector2D p = m_vertices[vv.index()];
            Vector2D gradient = gradients[bv.index()];
            Real gradientNorm = gradientNorms[bv.index()];
            Vector2D normal = normals[bv.index()];
            //std::cout << "Normal: " << normal << std::endl;

            Real partial = levelSet.partialDerivative(i, j, p);

            velocityForP(vv.index()) = - 1.0 / gradientNorm * partial * normal;
        }

        result[param] = velocityForP;
    }

    // Post processing. Based on what is done in Isoinflator.
    using FM = PeriodicBoundaryMatcher::FaceMembership<2>;
    vector<vector<bool>> onMinFace, onMaxFace;
    onMinFace.assign(2, std::vector<bool>(m_vertices.size(), false));
    onMaxFace.assign(2, std::vector<bool>(m_vertices.size(), false));
    std::vector<FM> faceMembership;
    for (size_t i = 0; i < m_vertices.size(); ++i) {
        FM fm(m_vertices[i], m_bbox, 0); // zero tolerance!
        onMinFace[0][i] = fm.onMinFace(0), onMaxFace[0][i] = fm.onMaxFace(0);
        onMinFace[1][i] = fm.onMinFace(1), onMaxFace[1][i] = fm.onMaxFace(1);
    }

    // Mark internal cell-face vertices: vertices on the meshing cell
    // boundary that actually lie inside the object (i.e. they are only mesh
    // boundary vertices because of the intersection of the periodic pattern with
    // the meshing box).
    // This this is not the case if any non-cell-face triangle is incident
    vector<bool> internalCellFaceVertex(mesh.numBoundaryVertices(), true);
    for (auto bs : mesh.boundarySimplices()) {
        bool isCellFace = false;
        for (size_t d = 0; d < 2; ++d) {
            bool onMin = true, onMax = true;
            for (auto bv : bs.vertices()) {
                onMin &= onMinFace[d].at(bv.volumeVertex().index());
                onMax &= onMaxFace[d].at(bv.volumeVertex().index());
            }

            if (onMin || onMax) {
                for (auto bv : bs.vertices()) {
                    auto vv = bv.volumeVertex();
                    for (unsigned param = 0; param < numParameters(); param++) {
                        result[param](vv.index())[d] = 0.0;
                    }
                }
                isCellFace = true;
            }
        }

        if (isCellFace) continue;

        for (auto bv : bs.vertices())
            internalCellFaceVertex.at(bv.index()) = false;
    }

    for (auto bv : mesh.boundaryVertices()) {
        if (internalCellFaceVertex.at(bv.index())) {
            auto vv = bv.volumeVertex();
            for (unsigned param = 0; param < numParameters(); param++) {
                for (size_t d = 0; d < 2; d++) {
                    result[param](vv.index())[d] = 0.0;
                }
            }
        }
    }

    return result;
}

size_t
RBFInflator::numParameters() const {
    return m_coeffMatrix.size() * m_coeffMatrix[0].size();
}

void
RBFInflator::loadMeshingOptions(const std::string &moptsPath) {
    m_meshingOptions.load(moptsPath);
}

void
RBFInflator::setMaxElementVolume(Real maxElementVol) {
    m_meshingOptions.maxArea = maxElementVol;
}

Real
RBFInflator::getMaxElementVolume() const {
    return m_meshingOptions.maxArea;
}


// Implementation of auxiliary functions
std::vector<std::vector<Real>>
RBFInflator::vecToMat(const std::vector<Real> &vec, size_t nRows, size_t nCols) const {
    std::vector<std::vector<Real>> result(nRows);

    for (size_t i = 0 ; i < nRows; i++) {
        vector<Real>::const_iterator first = vec.begin() + i * nCols;
        vector<Real>::const_iterator last = vec.begin() + (i+1) * nCols;
        vector<Real> newVec(first, last);
        result[i] = newVec;
    }

    return result;
}

std::vector<Real>
RBFInflator::matToVec(const std::vector<std::vector<Real>> &mat) const {
    std::vector<Real> result;

    for (size_t i = 0; i < mat.size(); i++) {
        result.insert(result.end(), mat[i].begin(), mat[i].end() );
    }

    return result;
}