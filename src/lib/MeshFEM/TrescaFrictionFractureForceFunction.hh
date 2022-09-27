//
// Created by Davi Colli Tozoni on 11/08/18.
//

#ifndef TRESCAFRICTIONFRACTUREFORCEFUNCTION_H
#define TRESCAFRICTIONFRACTUREFORCEFUNCTION_H

#include "MSHFieldWriter.hh"
#include "FrictionContactForceFunction.hh"

template<typename Real, class _Mesh>
class TrescaFrictionFractureForceFunction : public NonLinearElasticityFunction<Real> {
public:
    typedef _Mesh    Mesh;
    static constexpr size_t N = Mesh::FEMData::N;
    static constexpr size_t Degree = Mesh::FEMData::Degree;
    typedef VectorField<Real, N> VField;

    TrescaFrictionFractureForceFunction(_Mesh &mesh, Real alpha, Real etaRelu, Real etaNorm, Real mu, std::vector<Real> z) : m_mesh(mesh), m_alpha(alpha), m_h(etaRelu), m_norm(etaNorm), m_mu(mu), m_z(z) { }

    virtual ~TrescaFrictionFractureForceFunction() = default;

    virtual std::vector<Real> evaluate(std::vector<Real> u) const override {
        std::vector<Real> result(N * m_mesh.numNodes(), 0.0);
        VField uVectorField = dofToNodeField(u);
        VField zVectorField = dofToNodeField(m_z);

        MSHFieldWriter writer("trescaFrictionFractureForce.msh", m_mesh);

        for (auto ce : m_mesh.boundaryElements()) {
            if (!ce->isInContactRegion || ce->contactElement < 0)
                continue;

            VectorND<N> normalNegative = ce->normal();
            Real area = ce->volume();

            // Create interpolants and use it for integrating on z
            Interpolant<VectorND<N>, N - 1, Degree> zNegativeInterpolant;
            Interpolant<VectorND<N>, N - 1, Degree> zPositiveInterpolant;
            size_t zIndex = 0;
            for (auto n : ce.nodes()) {
                VectorND<N> negativePoint = n.volumeNode()->p;
                size_t negativeIndex = n.volumeNode().index();

                // Search for opposite corresponding node. To make sure it works always, check geometrically!
                size_t positiveIndex = 0;
                bool found = false;
                for (auto on : m_mesh.boundaryElement(ce->contactElement).nodes()) {
                    VectorND<N> positivePoint = on.volumeNode()->p;

                    if ((negativePoint - positivePoint).norm() < 1e-12) {
                        positiveIndex = on.volumeNode().index();
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw std::runtime_error("Did not find corresponding node on other fracture side");
                }

                zNegativeInterpolant[zIndex] = zVectorField(negativeIndex);
                zPositiveInterpolant[zIndex] = zVectorField(positiveIndex);
                zIndex++;
            }

            // Create interpolants and use it for integrating on u
            Interpolant<VectorND<N>, N - 1, Degree> uNegativeInterpolant;
            Interpolant<VectorND<N>, N - 1, Degree> uPositiveInterpolant;
            size_t uIndex = 0;
            for (auto n : ce.nodes()) {
                VectorND<N> negativePoint = n.volumeNode()->p;
                size_t negativeIndex = n.volumeNode().index();

                // Search for opposite corresponding node. To make sure it works always, check geometrically!
                size_t positiveIndex = 0;
                bool found = false;
                for (auto on : m_mesh.boundaryElement(ce->contactElement).nodes()) {
                    VectorND<N> positivePoint = on.volumeNode()->p;

                    if ((negativePoint - positivePoint).norm() < 1e-12) {
                        positiveIndex = on.volumeNode().index();
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw std::runtime_error("Did not find corresponding node on other fracture side");
                }

                uNegativeInterpolant[uIndex] = uVectorField(negativeIndex);
                uPositiveInterpolant[uIndex] = uVectorField(positiveIndex);
                uIndex++;
            }

            for (size_t A = 0; A < ce.numNodes(); ++A) {
                size_t globalA = ce.node(A).volumeNode().index();

                Interpolant<Real, N - 1, Degree> phiInterpolant;
                phiInterpolant = 0;
                phiInterpolant[A] = 1.0;

                for (size_t i = 0; i < N; i++) {

                    // Compute according to shape optimization documentation.
                    Real value = m_mu/m_alpha;
                    value *= Quadrature<N-1, 9>::integrate([&](const EvalPt<N-1> &pt) {
                        VectorND<N> zDiff = zNegativeInterpolant(pt) - zPositiveInterpolant(pt);
                        VectorND<N> uDiff = uNegativeInterpolant(pt) - uPositiveInterpolant(pt);
                        VectorND<N> uDiffT = uDiff - (uDiff.dot(normalNegative)) * normalNegative;
                        VectorND<N> NetaDerivative = m_norm.derivative(uDiffT);

                        return phiInterpolant(pt) * m_h.value(zDiff.dot(normalNegative)) * NetaDerivative(i);
                    }, area);

                    // save result in corresponding
                    result[N * globalA + i] += value;
                }
            }
        }

        writer.addField("u" , dofToNodeField(u), DomainType::PER_NODE);
        writer.addField("z" , dofToNodeField(m_z), DomainType::PER_NODE);
        writer.addField("tresca friction force" , -1.0 * dofToNodeField(result), DomainType::PER_NODE);

        return result;
    }

    virtual TripletMatrix<Triplet<Real>> jacobian(std::vector<Real> u) const override {
        TripletMatrix<Triplet<Real>> result;
        VField uVectorField = dofToNodeField(u);
        VField zVectorField = dofToNodeField(m_z);

        MatrixND<N> identity = Eigen::Matrix<double, N, N>::Identity();

        // Initialize matrix, trying to guess final size
        constexpr size_t perElementSize = Mesh::ElementData::PerElementStiffness::RowsAtCompileTime;
        result.init(N * m_mesh.numNodes(), N * m_mesh.numNodes());
        result.reserve(perElementSize * perElementSize * m_mesh.numBoundaryElements());

        for (auto ce : m_mesh.boundaryElements()) {
            if (!ce->isInContactRegion || ce->contactElement < 0)
                continue;

            VectorND<N> normalNegative = ce->normal();
            Real area = ce->volume();

            // Create interpolants and use it for integrating on z
            Interpolant<VectorND<N>, N - 1, Degree> zNegativeInterpolant;
            Interpolant<VectorND<N>, N - 1, Degree> zPositiveInterpolant;
            size_t zIndex = 0;
            for (auto n : ce.nodes()) {
                VectorND<N> negativePoint = n.volumeNode()->p;
                size_t negativeIndex = n.volumeNode().index();

                // Search for opposite corresponding node. To make sure it works always, check geometrically!
                size_t positiveIndex = 0;
                bool found = false;
                for (auto on : m_mesh.boundaryElement(ce->contactElement).nodes()) {
                    VectorND<N> positivePoint = on.volumeNode()->p;

                    if ((negativePoint - positivePoint).norm() < 1e-12) {
                        positiveIndex = on.volumeNode().index();
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw std::runtime_error("Did not find corresponding node on other fracture side");
                }

                zNegativeInterpolant[zIndex] = zVectorField(negativeIndex);
                zPositiveInterpolant[zIndex] = zVectorField(positiveIndex);
                zIndex++;
            }

            // Create interpolants and use it for integrating on u
            Interpolant<VectorND<N>, N - 1, Degree> uNegativeInterpolant;
            Interpolant<VectorND<N>, N - 1, Degree> uPositiveInterpolant;
            std::vector<size_t> oppositeElementOrderOfNodes;
            size_t uIndex = 0;
            for (auto n : ce.nodes()) {
                VectorND<N> negativePoint = n.volumeNode()->p;
                size_t negativeIndex = n.volumeNode().index();

                // Search for opposite corresponding node. To make sure it works always, check geometrically!
                size_t positiveIndex = 0;
                size_t oIndex = 0;
                bool found = false;
                for (auto on : m_mesh.boundaryElement(ce->contactElement).nodes()) {
                    VectorND<N> positivePoint = on.volumeNode()->p;

                    if ((negativePoint - positivePoint).norm() < 1e-12) {
                        positiveIndex = on.volumeNode().index();
                        oppositeElementOrderOfNodes.push_back(oIndex);
                        found = true;
                        break;
                    }

                    oIndex++;
                }
                if (!found) {
                    throw std::runtime_error("Did not find corresponding node on other fracture side");
                }

                uNegativeInterpolant[uIndex] = uVectorField(negativeIndex);
                uPositiveInterpolant[uIndex] = uVectorField(positiveIndex);
                uIndex++;
            }

            for (size_t A = 0; A < ce.numNodes(); ++A) {
                size_t globalA = ce.node(A).volumeNode().index();

                Interpolant < Real, N - 1, Degree > phiA;
                phiA = 0;
                phiA[A] = 1.0;

                for (size_t B = 0; B < ce.numNodes(); ++B) {
                    size_t globalB = ce.node(B).volumeNode().index();

                    Interpolant<Real, N - 1, Degree> phiB;
                    phiB = 0;
                    phiB[B] = 1.0;

                    for (size_t k = 0; k < N; k++) {
                        for (size_t l = 0; l < N; l++) {

                            Real value = m_mu / m_alpha;
                            value *= Quadrature<N - 1, 9>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> zDiff = zNegativeInterpolant(pt) - zPositiveInterpolant(pt);
                                VectorND<N> uDiff = uNegativeInterpolant(pt) - uPositiveInterpolant(pt);
                                VectorND<N> uDiffT = uDiff - (uDiff.dot(normalNegative)) * normalNegative;
                                VectorND<N> NetaDerivative = m_norm.derivative(uDiffT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uDiffT);

                                Real firstInnerPart = 0.0;
                                Real secondInnerPart = m_h.value(zDiff.dot(normalNegative)) * NetaSecondDerivative.row(k) * (identity.col(l) - normalNegative(l) * normalNegative);

                                return phiA(pt) * phiB(pt) * (firstInnerPart +  secondInnerPart);

                            }, area);

                            // save result in corresponding
                            result.addNZ(N * globalA + k, N * globalB + l, value);
                        }
                    }
                }

                // Now, also need to add derivative with respect to opposite edge
                for (size_t B = 0; B < ce.numNodes(); ++B) {
                    // Notice that global node B should be selected according to order of corresponding nodes in the opposite elements
                    size_t globalB = m_mesh.boundaryElement(ce->contactElement).node(oppositeElementOrderOfNodes[B]).volumeNode().index();

                    Interpolant<Real, N - 1, Degree> phiB;
                    phiB = 0;
                    phiB[B] = 1.0;

                    for (size_t k = 0; k < N; k++) {
                        for (size_t l = 0; l < N; l++) {
                            Real value = m_mu / m_alpha;
                            value *= Quadrature<N - 1, 9>::integrate([&](const EvalPt<N - 1> &pt) {
                                VectorND<N> zDiff = zNegativeInterpolant(pt) - zPositiveInterpolant(pt);
                                VectorND<N> uDiff = uNegativeInterpolant(pt) - uPositiveInterpolant(pt);
                                VectorND<N> uDiffT = uDiff - (uDiff.dot(normalNegative)) * normalNegative;
                                VectorND<N> NetaDerivative = m_norm.derivative(uDiffT);
                                MatrixND<N> NetaSecondDerivative = m_norm.secondDerivative(uDiffT);

                                Real firstInnerPart = 0.0;
                                Real secondInnerPart = m_h.value(zDiff.dot(normalNegative)) * NetaSecondDerivative.row(k) * (identity.col(l) - normalNegative(l) * normalNegative);

                                return - phiA(pt) * phiB(pt) * (firstInnerPart +  secondInnerPart);

                            }, area);

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

    void setTrescaDisplacements(std::vector<Real> newZ) {
        m_z = newZ;
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
    std::vector<Real> m_z; // fake displacement for Tresca model
};


#endif //TRESCAFRICTIONFRACTUREFORCEFUNCTION_H
