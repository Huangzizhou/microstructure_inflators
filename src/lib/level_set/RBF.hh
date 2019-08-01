////////////////////////////////////////////////////////////////////////////////
// RBF.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes signed distance represented by voxel densities
*/
////////////////////////////////////////////////////////////////////////////////
#ifndef RBF_HH
#define RBF_HH

#include <isosurface_inflator/SignedDistanceRegion.hh>
#include <isosurface_inflator/SignedDistance.hh>

#include <MeshFEM/Future.hh>
#include <unordered_map>
#include <unordered_set>
#include <math.h>


namespace SD = SignedDistance;

template<typename _Real>
class RBF : public SignedDistanceRegion<2> {

public:

    using Real = _Real;
    RBF(const std::vector<std::vector<Real>> &densityMatrix, Real epsilon);

    // Always support double type for compatibility with SignedDistanceRegion
    virtual double signedDistance(const Point2D &p) const override {
        return stripAutoDiff(m_signedDistanceImpl(autodiffCast<Point2<Real>>(p)));
    }

    // Also support automatic differentiation types
    template<typename Real2, bool DebugDerivatives = false>
    Real2 signedDistance(const Point2<Real2> &p) const { return m_signedDistanceImpl(p); }

    size_t numParams() const { return m_coeffMatrix.size() * m_coeffMatrix[0].size(); }

    // Representative cell bounding box (region to be meshed)
    virtual const BBox<Point2D> &boundingBox() const override { return m_bbox; }

    virtual ~RBF() override = default;


private:

    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 m_signedDistanceImpl(const Point2<Real2> p) const {
        /*
        results = np.zeros(len(x1))

        coeff_idx = 0
        for i in range(0, self.d1):
            for j in range(0, self.d2):
                # Find rbf point correspond to i and j
                xi = self.min1 + i * self.dt1
                xj = self.min2 + j * self.dt1

                r = np.sqrt(np.power(x1 - xi, 2) + np.power(x2 - xj, 2))

                results += self.coeffs[coeff_idx] * np.exp(- np.power((self.epsilon * r), 2))

                coeff_idx += 1
        return results
        */

        Real2 result = 0.0;

        for (size_t i=0; i<m_d1; i++) {
            for (size_t j=0; j<m_d2; j++) {
                // Find rbf point correspond to i and j
                Real xi = m_min1 + i * m_dt1;
                Real xj = m_min2 + j * m_dt1;

                Real2 r = sqrt((p[0] - xi) * (p[0] - xi) + (p[1] - xj) * (p[1] - xj));

                result += m_coeffMatrix[i][j] * exp(- (m_epsilon * r) * (m_epsilon * r));
            }
        }

        return result;
    }

    // Bounding box for the meshing cell.
    BBox<Point2D> m_bbox;

    // density matrix
    std::vector<std::vector<Real>> m_coeffMatrix;
    size_t m_d1, m_d2;
    Real m_epsilon;
    Real m_dt1, m_dt2;
    Real m_min1, m_min2;
    Real m_max1, m_max2;
};

#endif /* end of include guard: VOXELSIGNEDDISTANCE_HH */
