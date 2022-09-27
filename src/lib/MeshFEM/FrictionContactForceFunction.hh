//
// Created by Davi Colli Tozoni on 10/7/18.
//

#ifndef FRICTIONCONTACTFORCEFUNCTION_H
#define FRICTIONCONTACTFORCEFUNCTION_H

#define DEBUG_QUADRATURE 0
#define QUADRATURE_ORDER 13

#include "NonLinearElasticityFunction.hh"

template<size_t N>
using MatrixND = Eigen::Matrix<Real, N, N>;

template<size_t N>
class SmoothNorm {
public:
    SmoothNorm(double eta = 1e-2) : m_eta(eta) {
    }

    double value(PointND<N> x) const {
        double result;

        double xNorm = x.norm();
        if (xNorm >= m_eta) {
            result = xNorm;
        }
        else {
            result = xNorm * innerValueZ(x/m_eta)   +   3.0/8.0*m_eta;

            double directComputation = -1.0/(8.0*m_eta*m_eta*m_eta) * xNorm*xNorm*xNorm*xNorm   +   3.0/(4.0*m_eta) * xNorm*xNorm   +   3.0/8.0*m_eta;
            if (abs(result - directComputation) > 1e-12) {
                std::cerr << "Difference in computing norm: " << abs(result - directComputation) << std::endl;
            }
        }

        return result;
    }

    PointND<N> derivative(PointND<N> x) const {
        PointND<N> result;

        double xNorm = x.norm();
        if (xNorm >= m_eta) {
            result = 1.0/xNorm * x;
        }
        else {
            result = innerDerivativeZ(x/m_eta);

            PointND<N> directComputation = -(xNorm*xNorm)/(2.0*m_eta*m_eta*m_eta) * x   +   3.0/(2.0*m_eta) * x;
            if ((result - directComputation).norm() > 1e-12) {
                std::cerr << "Difference in computing norm derivative: " << (result - directComputation).norm() << std::endl;
            }
        }

        return result;
    }

    MatrixND<N> secondDerivative(PointND<N> x) const {
        MatrixND<N> result;

        MatrixND<N> identity = Eigen::Matrix<double, N, N>::Identity();
        MatrixND<N> xxt = x * x.transpose();
        double xNorm = x.norm();

        if (xNorm >= m_eta) {
            result = 1.0/xNorm * identity   -   1./(xNorm*xNorm*xNorm) * xxt;
        }
        else {
            result = 1.0/m_eta * innerSecondDerivativeZ(x/m_eta);

            MatrixND<N> directComputation = -1.0/(m_eta*m_eta*m_eta) * xxt   -   (xNorm*xNorm)/(2.0*m_eta*m_eta*m_eta) * identity   +   3.0/(2.0*m_eta) * identity;
            if ((result - directComputation).norm() > 1e-10) {
                std::cerr << "Difference in computing norm seco nd derivative: " << (result - directComputation).norm() << std::endl;
            }
        }

        return result;
    }

private:
    double m_eta;

    // z = x / m_eta
    double innerValueZ(PointND<N> z) const {
        double zNorm = z.norm();

        double result = zNorm * (- zNorm*zNorm/8.0  +   3.0/4.0);

        return result;
    }

    // z = x / m_eta
    PointND<N> innerDerivativeZ(PointND<N> z) const {
        double zNorm = z.norm();

        PointND<N> result = - zNorm*zNorm/2.0 * z   +   3.0/2.0 * z;

        return result;
    }

    // z = x / m_eta
    MatrixND<N> innerSecondDerivativeZ(PointND<N> z) const {
        double zNorm = z.norm();
        MatrixND<N> zzt = z * z.transpose();
        MatrixND<N> identity = Eigen::Matrix<double, N, N>::Identity();

        MatrixND<N> result = - zzt   +  (- zNorm*zNorm/2.0  +  3.0/2.0) * identity;

        return result;
    }
};

template<typename Real, class _Mesh>
class FrictionContactForceFunction : public NonLinearElasticityFunction<Real> {
public:
    typedef _Mesh    Mesh;
    static constexpr size_t N = Mesh::FEMData::N;
    static constexpr size_t Degree = Mesh::FEMData::Degree;
    typedef VectorField<Real, N> VField;

    FrictionContactForceFunction(_Mesh &mesh, Real alpha, Real etaRelu, Real etaNorm, Real mu) : m_mesh(mesh), m_alpha(alpha), m_h(etaRelu), m_norm(etaNorm), m_mu(mu) { }

    virtual ~FrictionContactForceFunction() = default;

    virtual std::vector<Real> evaluate(std::vector<Real> u) const override {
        std::vector<Real> result(N * m_mesh.numNodes(), 0.0);
        VField uVectorField = dofToNodeField(u);

        MSHFieldWriter writer("frictionContactForce.msh", m_mesh);

        for (auto ce : m_mesh.boundaryElements()) {
            if (!ce->isInContactRegion || ce->contactElement > 0)
                continue;

            VectorND<N> normal = ce->normal();
            Real area = ce->volume();

            for (size_t A = 0; A < ce.numNodes(); ++A) {
                size_t globalA = ce.node(A).volumeNode().index();

                for (size_t i = 0; i < N; i++) {

                    // Create interpolants and use it for integrating on u and phi
                    Interpolant<VectorND<N>, N - 1, Degree> uInterpolant;
                    size_t uIndex = 0;
                    for (auto n : ce.nodes()) {
                        uInterpolant[uIndex] = uVectorField(n.volumeNode().index());
                        uIndex++;
                    }

                    Interpolant<Real, N - 1, Degree> phiInterpolant;
                    phiInterpolant = 0;
                    phiInterpolant[A] = 1.0;

                    // Compute according to shape optimization documentation.
                    Real value = m_mu/m_alpha;
                    value *= Quadrature<N-1, QUADRATURE_ORDER>::integrate([&](const EvalPt<N-1> &pt) {
                        VectorND<N> uT = uInterpolant(pt) - (uInterpolant(pt).dot(normal)) * normal;
                        VectorND<N> NetaDerivative = m_norm.derivative(uT);
                        return phiInterpolant(pt) * m_h.value(uInterpolant(pt).dot(normal)) * NetaDerivative(i);
                    }, area);

                    // save result in corresponding
                    result[N * globalA + i] += value;
                }
            }
        }

        writer.addField("u" , dofToNodeField(u), DomainType::PER_NODE);
        writer.addField("friction force" , -1.0 * dofToNodeField(result), DomainType::PER_NODE);

        return result;
    }

    virtual TripletMatrix<Triplet<Real>> jacobian(std::vector<Real> u) const override {
        TripletMatrix<Triplet<Real>> result;
        VField uVectorField = dofToNodeField(u);

        MatrixND<N> identity = Eigen::Matrix<double, N, N>::Identity();

        // Initialize matrix, trying to guess final size
        constexpr size_t perElementSize = Mesh::ElementData::PerElementStiffness::RowsAtCompileTime;
        result.init(N * m_mesh.numNodes(), N * m_mesh.numNodes());
        result.reserve(perElementSize * perElementSize * m_mesh.numBoundaryElements());

        for (auto ce : m_mesh.boundaryElements()) {
            if (!ce->isInContactRegion || ce->contactElement > 0)
                continue;

            VectorND<N> normal = ce->normal();
            Real area = ce->volume();

            for (size_t A = 0; A < ce.numNodes(); ++A) {
                size_t globalA = ce.node(A).volumeNode().index();
                for (size_t B = 0; B < ce.numNodes(); ++B) {
                    size_t globalB = ce.node(B).volumeNode().index();

                    for (size_t k = 0; k < N; k++) {
                        for (size_t l = 0; l < N; l++) {

                            // Create interpolants and use it for integrating on u and phi
                            Interpolant < VectorND < N > , N - 1, Degree > uInterpolant;
                            size_t uIndex = 0;
                            for (auto n : ce.nodes()) {
                                uInterpolant[uIndex] = uVectorField(n.volumeNode().index());
                                uIndex++;
                            }

                            Interpolant < Real, N - 1, Degree > phiA;
                            phiA = 0;
                            phiA[A] = 1.0;

                            Interpolant < Real, N - 1, Degree > phiB;
                            phiB = 0;
                            phiB[B] = 1.0;

                            Real value = m_mu / m_alpha;
                            value *= Quadrature<N - 1, QUADRATURE_ORDER>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> uT = uInterpolant(pt) - (uInterpolant(pt).dot(normal)) * normal;
                                VectorND<N> NetaDerivative = m_norm.derivative(uT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uT);

                                Real firstInnerPart = m_h.derivative(uInterpolant(pt).dot(normal)) * NetaDerivative(k) * normal(l);
                                Real secondInnerPart = m_h.value(uInterpolant(pt).dot(normal)) * NetaSecondDerivative.row(k) * (identity.col(l) - normal(l) * normal);

                                return phiA(pt) * phiB(pt) * (firstInnerPart +  secondInnerPart);

                            }, area);

#if DEBUG_QUADRATURE
                            Real firstTerm5 = m_mu / m_alpha;
                            firstTerm5 *= Quadrature<N - 1, 5>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> uT = uInterpolant(pt) - (uInterpolant(pt).dot(normal)) * normal;
                                VectorND<N> NetaDerivative = m_norm.derivative(uT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uT);

                                Real firstInnerPart = m_h.derivative(uInterpolant(pt).dot(normal)) * NetaDerivative(k) * normal(l);
                                Real secondInnerPart = m_h.value(uInterpolant(pt).dot(normal)) * NetaSecondDerivative.row(k) * (identity.col(l) - normal(l) * normal);

                                return phiA(pt) * phiB(pt) * (firstInnerPart +  secondInnerPart);

                            }, area);

                            Real secondTerm5 = m_mu / m_alpha;
                            secondTerm5 *= Quadrature<N - 1, 5>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> uT = uInterpolant(pt) - (uInterpolant(pt).dot(normal)) * normal;
                                VectorND<N> NetaDerivative = m_norm.derivative(uT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uT);

                                Real firstInnerPart = m_h.derivative(uInterpolant(pt).dot(normal)) * NetaDerivative.dot(normal) * normal(l);
                                Real secondInnerPart = m_h.value(uInterpolant(pt).dot(normal)) * (NetaSecondDerivative * (identity.col(l) - normal(l) * normal)).dot(normal);
                                return phiA(pt) * phiB(pt) * normal(k) * (firstInnerPart + secondInnerPart);

                            }, area);

                            Real value5 = firstTerm5 - secondTerm5;

                            Real firstTerm3 = m_mu / m_alpha;
                            firstTerm3 *= Quadrature<N - 1, 3>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> uT = uInterpolant(pt) - (uInterpolant(pt).dot(normal)) * normal;
                                VectorND<N> NetaDerivative = m_norm.derivative(uT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uT);

                                Real firstInnerPart = m_h.derivative(uInterpolant(pt).dot(normal)) * NetaDerivative(k) * normal(l);
                                Real secondInnerPart = m_h.value(uInterpolant(pt).dot(normal)) * NetaSecondDerivative.row(k) * (identity.col(l) - normal(l) * normal);

                                return phiA(pt) * phiB(pt) * (firstInnerPart +  secondInnerPart);

                            }, area);

                            Real secondTerm3 = m_mu / m_alpha;
                            secondTerm3 *= Quadrature<N - 1, 3>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> uT = uInterpolant(pt) - (uInterpolant(pt).dot(normal)) * normal;
                                VectorND<N> NetaDerivative = m_norm.derivative(uT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uT);

                                Real firstInnerPart = m_h.derivative(uInterpolant(pt).dot(normal)) * NetaDerivative.dot(normal) * normal(l);
                                Real secondInnerPart = m_h.value(uInterpolant(pt).dot(normal)) * (NetaSecondDerivative * (identity.col(l) - normal(l) * normal)).dot(normal);
                                return phiA(pt) * phiB(pt) * normal(k) * (firstInnerPart + secondInnerPart);

                            }, area);


                            Real value3 = firstTerm3 - secondTerm3;

                            if (abs(value - value5) > 1e-10) {
                                std::cerr << "Difference between using 5th order quadrature is " << abs(value - value5) <<  std::endl;
                            }

                            if (abs(value - value3) > 1e-10) {
                                std::cerr << "Difference between using 3th order quadrature is " << abs(value - value3) <<  std::endl;
                            }
#endif
                            // save result in corresponding
                            result.addNZ(N * globalA + k, N * globalB + l, value);
                        }
                    }
                }
            }

            // Ended adding all element information
        }

        // Sort and sum corresponding matrix positions (to have only one value per matrix position)
        result.sumRepeated();
        return result;
    }

    virtual size_t size() const override {
        return N * m_mesh.numNodes();
    }

private:

    template<class _Vec>
    VField dofToNodeField(const _Vec &x) const {
        assert(x.size() == size());

        VField f(m_mesh.numNodes());
        for (size_t d = 0; d < m_mesh.numNodes(); d++) {
            for (size_t c = 0; c < N; ++c)
                f(d)[c] = x[N * d + c];
        }

        return f;
    }

    _Mesh &m_mesh;
    Real m_alpha;
    SmoothPositiveMax m_h;
    SmoothNorm<N> m_norm;
    Real m_mu; // friction coefficient
};


#endif //FRICTIONCONTACTFORCEFUNCTION_H
