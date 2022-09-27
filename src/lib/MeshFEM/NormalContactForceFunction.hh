//
// Created by Davi Colli Tozoni on 7/18/18.
//

#ifndef NORMALCONTACTFORCEFUNCTION_H
#define NORMALCONTACTFORCEFUNCTION_H

#include "NonLinearElasticityFunction.hh"

class SmoothPositiveMax {
public:
    SmoothPositiveMax(double eta = 1e-4) : m_eta(eta) {
    }

    double value(double x) const {
        double result;

        if (x <= -m_eta) {
            result = 0.0;
        }
        else if (x < m_eta) {
            result = x * innerValueY(x/m_eta)  +  m_eta/4.0;

            //double directComputation = 1.0/(4.0*m_eta) * x*x + 1.0/2.0 * x + m_eta/4.0;
            //if (abs(result - directComputation) > 1e-10) {
            //    std::cerr << "Difference in computing relu: " << abs(result - directComputation) << std::endl;
            //}
        }
        else {
            result = x;
        }

        return result;
    }

    double derivative(double x) const {
        double result;

        if (x <= -m_eta) {
            result = 0.0;
        }
        else if (x < m_eta) {
            result = innerDerivativeY(x/m_eta);

            //result = x/(2.0*m_eta) + 1.0/2.0;
            //double directComputation = x/(2.0*m_eta) + 1.0/2.0;
            //if (abs(result - directComputation) > 1e-10) {
            //    std::cerr << "Difference in computing relu derivative: " << abs(result - directComputation) << std::endl;
            //}
        }
        else {
            result = 1.0;
        }

        return result;
    }

private:

    // y = x / m_eta
    double innerValueY(double y) const {
        return (y + 2.0)/4.0;
    }

    // y = x / m_eta
    double innerDerivativeY(double y) const {
        return (y + 1)/2.0;
    }

    double m_eta;
};

template<typename Real, class _Mesh>
class NormalContactForceFunction : public NonLinearElasticityFunction<Real> {
public:
    typedef _Mesh    Mesh;
    static constexpr size_t N = Mesh::FEMData::N;
    static constexpr size_t Degree = Mesh::FEMData::Degree;
    typedef VectorField<Real, N> VField;

    NormalContactForceFunction(_Mesh &mesh, Real alpha, Real etaRelu) : m_mesh(mesh), m_alpha(alpha), m_h(etaRelu) { }

    virtual ~NormalContactForceFunction() = default;

    virtual std::vector<Real> evaluate(std::vector<Real> u) const override {
        std::vector<Real> result(N * m_mesh.numNodes(), 0.0);
        VField uVectorField = dofToNodeField(u);

        MSHFieldWriter writer("normalContactForce.msh", m_mesh);

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
                    // TODO: Take a better look at the degree we need for this quadrature, cause phi is has degree
                    // "Degree" and h_eta has maximum degree of 2*Degree on (u.n)
                    Real value = 1.0/m_alpha;
                    value *= Quadrature<N-1, 7>::integrate([&](const EvalPt<N-1> &pt) {
                        return phiInterpolant(pt) * normal(i) * m_h.value(uInterpolant(pt).dot(normal));
                    }, area);

                    // save result in corresponding
                    result[N * globalA + i] += value;
                }
            }
        }

        writer.addField("u" , dofToNodeField(u), DomainType::PER_NODE);
        writer.addField("normal force" , -1.0 * dofToNodeField(result), DomainType::PER_NODE);

        return result;
    }

    virtual TripletMatrix<Triplet<Real>> jacobian(std::vector<Real> u) const override {
        TripletMatrix<Triplet<Real>> result;
        VField uVectorField = dofToNodeField(u);

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

                            // Compute according to equation (10) in shape optimization documentation.
                            Real value = 1.0 / m_alpha;
                            value *= Quadrature<N - 1, 7>::integrate([&](const EvalPt<N - 1> &pt) {
                                return phiA(pt) * normal(k) * m_h.derivative(uInterpolant(pt).dot(normal)) * phiB(pt) * normal(l);
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
};




#endif //NORMALCONTACTFORCEFUNCTION_H
