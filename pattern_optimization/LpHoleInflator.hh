////////////////////////////////////////////////////////////////////////////////
// LpHoleInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Simple 2D inflator for debugging: an p-norm hole inflator.
//      There are two parameters: the radius, r in (0, 1), and the Lp-norm
//      parameter, p in (0, inf). The base cell is the square [-1, 1]^2 with a
//      hole of radius r and centered at (0, 0) subtracted.
//
//      To ensure a periodic boundary with minimal effort, the base cell is then
//      reflected into a 4x4 grid of cells.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/10/2015 22:09:49
////////////////////////////////////////////////////////////////////////////////
#ifndef LPHOLEINFLATOR_HH
#define LPHOLEINFLATOR_HH

#include "InflatorBase.hh"
#include <Triangulate.h>
#include <filters/reflect.hh>
#include <Fields.hh>
#include <cmath>

struct LpHoleInflator : InflatorBase<LpHoleInflator> {
    using Base = InflatorBase<LpHoleInflator>;
    using NSV = NormalShapeVelocity<2>;
    
    LpHoleInflator() { }

    // Parameters: (radius, p)
    size_t numParameters() const { return 2; }
    ParameterType parameterType(size_t p) const {
        if (p == 0) return ParameterType::Thickness;
        if (p == 1) return ParameterType::Blending; // Sorta...
        assert(false);
    }

    void setNumSubdiv(size_t ns) { m_nsubdiv = ns; }

    void inflate(const std::vector<Real> &params) {
        assert(params.size() == numParameters());
        m_radius = params[0];
        m_p = params[1];

        // Create the square
        std::vector<MeshIO::IOVertex> inVertices = { {-1, -1, 0},
                                                     { 1, -1, 0},
                                                     { 1,  1, 0},
                                                     {-1,  1, 0} };
        std::vector<std::pair<size_t, size_t>> inEdges = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };

        size_t firstHoleVertex = inVertices.size();
        inVertices.emplace_back(m_radius, 0);

        double degreesPerSubdiv = 2 * M_PI / m_nsubdiv;
        // Create all hole boundary segments except the last
        for (size_t i = 1; i < m_nsubdiv; ++i) {
            double theta = degreesPerSubdiv * i;
            // https://www.mathworks.com/matlabcentral/newsreader/view_thread/279050
            inVertices.emplace_back(
                m_radius * copysign(pow(fabs(cos(theta)), 2 / m_p), cos(theta)),
                m_radius * copysign(pow(fabs(sin(theta)), 2 / m_p), sin(theta)));
            inEdges.push_back({inVertices.size() - 2,
                               inVertices.size() - 1});
        }

        // Close the hole path.
        inEdges.push_back({inVertices.size() - 1, firstHoleVertex});

        // Pick a point in the hole: the first segment forms a triangle with the
        // origin that lies entirely within the hole. Choose its barycenter.
        std::vector<MeshIO::IOVertex> holes;
        holes.emplace_back(((1 / 3.0) * (inVertices.at(firstHoleVertex    ).point +
                                         inVertices.at(firstHoleVertex + 1).point)).eval());

        triangulatePSLC(inVertices, inEdges, holes, m_vertices, m_elements,
                        m_triangleArea, "Q");

        reflect(2, m_vertices, m_elements, m_vertices, m_elements);
    }

    void setMaxElementVolume(Real maxElementVol) { m_triangleArea = maxElementVol; }

    // Isosurface normal velocity:
    // phi(x, p) = 0 = r^p - |x|^p - |y|^p
    // grad phi . dx/dp + dphi/dp = 0
    // |grad phi| n . dx/dp = -dphi/dp
    // n . dx/dp = -1/|grad phi| dphi/dp
    template<class _FEMMesh>
    std::vector<NSV> computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        size_t numBE = mesh.numBoundaryElements();
        size_t nParams = numParameters();
        std::vector<NSV> vn_p(nParams, NSV(numBE));

        for (auto be : mesh.boundaryElements()) {
            assert(be.numVertices() == 2);
            if (be->isPeriodic) {
                vn_p[0][be.index()] = 0;
                vn_p[1][be.index()] = 0;
                continue;
            }
            for (size_t v = 0; v < 2; ++v) {
                Real x = be.vertex(v).volumeVertex().node()->p[0];
                Real y = be.vertex(v).volumeVertex().node()->p[1];

                // map back to [-1, 1]^2 base cell
                if (x < -1) { x += 2; x *= -1; }
                if (y < -1) { y += 2; y *= -1; }
                assert(x >= -1 && x <= 1.0);
                assert(y >= -1 && y <= 1.0);

                Real absX = fabs(x), absY = fabs(y);

                Real phi_x = -copysign(m_p * pow(absX, m_p - 1.0), x);
                Real phi_y = -copysign(m_p * pow(absY, m_p - 1.0), y);
                Real inv_grad_norm = 1.0 / sqrt(phi_x * phi_x + phi_y * phi_y);

                Real phi_r = m_p * pow(m_radius, m_p - 1.0);
                Real phi_p = pow(m_radius, m_p) * log(m_radius)
                           - ((absX > 1e-9) ? pow(absX, m_p) * log(absX) : 0.0)
                           - ((absY > 1e-9) ? pow(absY, m_p) * log(absY) : 0.0);
                
                Real nsv_r = -inv_grad_norm * phi_r;
                Real nsv_p = -inv_grad_norm * phi_p;
                if (std::isnan(nsv_r) || std::isnan(nsv_p))
                    throw std::runtime_error("BLARGH");
                vn_p[0][be.index()][v] = nsv_r;
                vn_p[1][be.index()][v] = nsv_p;
            }
        }

        return vn_p;
    }

    // Boundary vector field corresponding to boundary normal shape velocity
    // scalar field nsv.
    template<class _FEMMesh>
    VectorField<Real, 2> shapeVelocity(const _FEMMesh &mesh, const NSV &nsv) {
        VectorField<Real, 2> result = analyticNormals(mesh);
        assert(nsv.size() == mesh.numBoundaryElements());
        std::vector<bool> isSet(mesh.numBoundaryVertices(), false);
        for (auto be : mesh.boundaryElements()) {
            const auto &nsvbe = nsv.at(be.index());
            assert(nsvbe.size() == be.numVertices());
            for (size_t i = 0; i < be.numVertices(); ++i) {
                size_t bvi = be.vertex(i).index();
                // Zero velocity overrides
                if (!isSet[bvi] || (nsvbe[i] == 0.0)) {
                    result(bvi) *= nsvbe[i];
                    isSet[bvi] = true;
                }
            }
        }
        return result;
    }

    // Get a per-boundary-vertex perturbation vector field induced by changing
    // each parameter.
    // This can be used with the discrete shape derivative.
    template<class _FEMMesh>
    std::vector<VectorField<Real, 2>>
    shapeVelocities(const _FEMMesh &mesh) const {
        size_t numBV = mesh.numBoundaryVertices();

        std::vector<NSV> nsv   = computeShapeNormalVelocities(mesh);
        VectorField<Real, 2> n = analyticNormals(mesh);

        // TODO: it'd be faster to just recompute shape velocity at each vertex...
        std::vector<VectorField<Real, 2>> svel(numParameters());
        for (size_t p = 0; p < svel.size(); ++p) {
            VectorField<Real, 2> &svel_p = svel.at(p);
            const auto &nsv_p = nsv.at(p);
            assert(nsv_p.size() == mesh.numBoundaryElements());
            svel_p.resizeDomain(numBV);
            for (auto be : mesh.boundaryElements()) {
                const auto &nsv_be = nsv_p.at(be.index());
                for (size_t v = 0; v < 2; ++v) {
                    size_t bvi = be.vertex(v).index();
                    svel_p(bvi)  = n(bvi);
                    svel_p(bvi) *= nsv_be[v];
                }
            }
        }

        return svel;
    }

    // Isosurface normal:
    // phi(x, p) = 0 = r^p - x^p - y^p
    // grad phi / |grad phi|
    template<class _FEMMesh>
    VectorField<Real, 2> analyticNormals(const _FEMMesh &mesh) const {
        VectorField<Real, 2> normals;
        normals.resizeDomain(mesh.numBoundaryVertices());

        for (auto bv : mesh.boundaryVertices()) {
            auto v = bv.volumeVertex().node()->p;
            Vector2D n;
            for (size_t c = 0; c < 2; ++c) {
                bool flip = false;
                // map back to base cell
                if (v[c] < -1) { v[c] += 2; v[c] *= -1; flip = true; }
                n[c] = -copysign(m_p * pow(fabs(v[c]), m_p - 1.0), v[c]); // grad_c phi
                if (flip) n[c] *= -1;
            }
            normals(bv.index()) = n.normalized();
        }

        // Clear normals on the periodic boundary
        for (auto be : mesh.boundaryElements()) {
            if (be->isPeriodic) {
                normals(be.vertex(0).index()).setZero();
                normals(be.vertex(1).index()).setZero();
            }
        }

        return normals;
    }

    // 2D is always printable.
    bool isPrintable(const std::vector<Real> &/* params */) const { return true; }

private:
    size_t m_nsubdiv = 64;
    Real   m_triangleArea = 0.001;

    Real   m_radius, m_p; // Set by inflation
};

#endif /* end of include guard: LPHOLEINFLATOR_HH */
