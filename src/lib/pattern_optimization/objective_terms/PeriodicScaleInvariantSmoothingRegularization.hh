//
// Attempt to make regularization not dependend on "area" of vertices
//

#ifndef SCALEINVARIANTSMOOTHINGREGULARIZATION_H
#define SCALEINVARIANTSMOOTHINGREGULARIZATION_H

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include <MeshFEM/OneForm.hh>
#include <pattern_optimization/SDConversions.hh>
#include <Eigen/Sparse>

namespace PatternOptimization {
    namespace ObjectiveTerms {

        template<class _Sim>
        struct PeriodicScaleInvariantSmoothingRegularizationTerm : public ObjectiveTerm<_Sim::N> {
            using SField = ScalarField<Real>;
            using VField = VectorField<Real, _Sim::N>;
            using  OForm = ScalarOneForm<_Sim::N>;

            static constexpr size_t N = _Sim::N;

            static SField computeSmoothingField(_Sim &simulator) {
                std::unique_ptr<PeriodicCondition<N>> pc = Future::make_unique<PeriodicCondition<N>>(simulator.mesh());
                std::vector<std::vector<size_t>> correspondingVertices(simulator.mesh().numVertices());

                for (auto v : simulator.mesh().boundaryVertices()) {
                    size_t vNodeIndex = v.node().volumeNode().index();
                    std::vector<size_t> vNodeIndices = pc->identifiedNodes(vNodeIndex);
                    std::vector<size_t> vIndices;

                    for (size_t i = 0; i < vNodeIndices.size(); i++) {
                        size_t vIndex = simulator.mesh().node(vNodeIndices[i]).vertex().index();
                        vIndices.push_back(vIndex);
                    }

                    // Sort corresponding vertices. This way, always first in the list has lowest index
                    std::sort(vIndices.begin(), vIndices.end());
                    correspondingVertices[v.node().volumeNode().vertex().index()] = vIndices;
                }

                // Computing adjacency list for each vertex
                std::vector<std::set<size_t>> adj(simulator.mesh().numVertices());
                for (auto be : simulator.mesh().boundaryElements()) {
                    if (be->isInternal) {
                        continue;
                    }

                    for (auto v1 : be.vertices()) {
                        size_t v1Index = v1.node().volumeNode().vertex().index();

                        for (auto v2 : be.vertices()) {
                            size_t v2Index = v2.node().volumeNode().vertex().index();

                            if (v2 == v1)
                                continue;
                            else if (!adj[v1Index].count(v2Index)) {
                                adj[v1Index].insert(v2Index);
                            }
                        }
                    }
                }

                int numVertices = simulator.mesh().numVertices();
                std::vector<Real> result(numVertices, 0.0);
                for (auto ve : simulator.mesh().boundaryVertices()) {
                    size_t veIndex = ve.node().volumeNode().vertex().index();

                    size_t number_neighbors = 0;
                    for (auto idx : correspondingVertices[veIndex]) {
                        for (auto neighbor : adj[idx]) {
                            number_neighbors++;
                        }
                    }

                    if (number_neighbors == 0) {
                        continue;
                    }

                    // For each neighbor, add the corresponding terms
                    Eigen::Matrix<Real, N, 1> numerator = Eigen::MatrixXd::Zero(N, 1);
                    Real denominator = 0.0;
                    for (auto idx : correspondingVertices[veIndex]) {
                        Eigen::Matrix<Real, N, 1> v = simulator.mesh().vertex(idx).node()->p;
                        for (auto neighbor : adj[idx]) {
                            Eigen::Matrix<Real, N, 1> u = simulator.mesh().vertex(neighbor).node()->p;
                            Eigen::Matrix<Real, N, 1> diff = u - v;

                            numerator += diff;
                            denominator += diff.norm();
                        }
                    }

                    /*if (correspondingVertices[veIndex].size() > 1) {
                        std::cout << "UHUULLL" << std::endl;
                        std::cout << "ve.node().volumeNode()->p: " << ve.node().volumeNode()->p << std::endl;
                        std::cout << "numerator: " << numerator << std::endl;
                        std::cout << "denominator: " << denominator << std::endl;
                    }*/

                    Eigen::Matrix<Real, N, 1> smoothness = numerator / denominator;
                    result[veIndex] = smoothness.squaredNorm();
                }

                return result;
            }


            static Real computeCost(_Sim &simulator) {
                SField field = computeSmoothingField(simulator);
                Real result = 0;

                for (size_t i=0; i<field.size(); i++) {
                    //std::cout << i << ": " << field[i] << std::endl;
                    result += field[i];
                }

                return result;
            }

            static OForm computeVolumeDifferential(_Sim &simulator) {
                size_t numVertices = simulator.mesh().numVertices();
                Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> differential = Eigen::MatrixXd::Zero(N, numVertices);

                std::unique_ptr<PeriodicCondition<N>> pc = Future::make_unique<PeriodicCondition<N>>(simulator.mesh());
                std::vector<std::vector<size_t>> correspondingVertices(simulator.mesh().numVertices());

                for (auto v : simulator.mesh().boundaryVertices()) {
                    size_t vNodeIndex = v.node().volumeNode().index();
                    std::vector<size_t> vNodeIndices = pc->identifiedNodes(vNodeIndex);
                    std::vector<size_t> vIndices;

                    for (size_t i = 0; i < vNodeIndices.size(); i++) {
                        size_t vIndex = simulator.mesh().node(vNodeIndices[i]).vertex().index();
                        vIndices.push_back(vIndex);
                    }

                    // Sort corresponding vertices. This way, always first in the list has lowest index
                    std::sort(vIndices.begin(), vIndices.end());
                    correspondingVertices[v.node().volumeNode().vertex().index()] = vIndices;
                }

                std::vector<std::set<size_t>> adj(simulator.mesh().numVertices());
                for (auto be : simulator.mesh().boundaryElements()) {
                    if (be->isInternal) {
                        continue;
                    }

                    for (auto v1 : be.vertices()) {
                        size_t v1Index = v1.node().volumeNode().vertex().index();

                        for (auto v2 : be.vertices()) {
                            size_t v2Index = v2.node().volumeNode().vertex().index();

                            if (v2 == v1)
                                continue;
                            else if (!adj[v1Index].count(v2Index)) {
                                adj[v1Index].insert(v2Index);
                            }
                        }
                    }
                }

                for (auto ve : simulator.mesh().boundaryVertices()) {
                    size_t veIndex = ve.node().volumeNode().vertex().index();
                    //Eigen::Matrix<Real, N, 1> v = ve.node().volumeNode()->p;

                    size_t number_neighbors = 0;
                    for (auto idx : correspondingVertices[veIndex]) {
                        for (auto neighbor : adj[idx]) {
                            number_neighbors++;
                        }
                    }

                    if (number_neighbors == 0)
                        continue;

                    Eigen::Matrix<Real, N, 1> laplacian = Eigen::MatrixXd::Zero(N, 1);
                    Real denominator = 0.0;

                    std::vector<Eigen::Matrix<Real, N, 1>> diffs;
                    std::vector<Real> diffs_norms;
                    Real sum_norms = 0.0;
                    for (auto idx : correspondingVertices[veIndex]) {
                        Eigen::Matrix<Real, N, 1> v = simulator.mesh().vertex(idx).node()->p;
                        for (auto neighbor : adj[idx]) {
                            Eigen::Matrix<Real, N, 1> u = simulator.mesh().vertex(neighbor).node()->p;
                            Eigen::Matrix<Real, N, 1> diff = u - v;
                            Real diff_norm = diff.norm();

                            laplacian += diff;
                            sum_norms += diff_norm;

                            diffs.push_back(diff);
                            diffs_norms.push_back(diff_norm);
                        }
                    }

                    // - 4 (u+w - 2v) / (|u-v| + |w-v|)^2
                    Eigen::Matrix<Real, N, 1> first_part = - 4 * laplacian / (sum_norms * sum_norms);

                    // 2 |u+w-2v|^2 / (|u-v| + |w-v|)^3
                    Eigen::Matrix<Real, N, 1> second_part = Eigen::MatrixXd::Zero(N, 1);;
                    Real second_part_multiplier = 2 * laplacian.squaredNorm() / (sum_norms * sum_norms * sum_norms);
                    size_t i = 0;
                    for (auto idx : correspondingVertices[veIndex]) {
                        for (auto neighbor : adj[idx]) {
                            second_part += diffs[i] / diffs_norms[i];
                            i++;
                        }
                    }
                    second_part = second_part_multiplier * second_part;

                    Eigen::Matrix<Real, N, 1> dRdv = first_part + second_part;
                    differential.col(veIndex) += dRdv;

                    // We also need to add the term of the regularization with respect to the neighbors.
                    for (auto idx : correspondingVertices[veIndex]) {
                        Eigen::Matrix<Real, N, 1> v = simulator.mesh().vertex(idx).node()->p;
                        for (auto neighbor : adj[idx]) {
                            Eigen::Matrix<Real, N, 1> u = simulator.mesh().vertex(neighbor).node()->p;
                            Eigen::Matrix<Real, N, 1> diff = u - v;

                            Eigen::Matrix<Real, N, 1> dRdu;

                            // 2 / (\u-v\ + \w-v\)^2 * (u+w-2v)
                            Eigen::Matrix<Real, N, 1> first_part_neighbor = 2.0 / (sum_norms * sum_norms) * laplacian;

                            // -2 |u+w-2v|^2 / (\u-v\ + \w-v\)^3 * (u-v) / |u-v|
                            Eigen::Matrix<Real, N, 1> second_part_neighbor =
                                    -2.0 * laplacian.squaredNorm() / (sum_norms * sum_norms * sum_norms) * diff /
                                    diff.norm();

                            dRdu = first_part_neighbor + second_part_neighbor;

                            differential.col(neighbor) += dRdu;
                        }
                    }
                }

                differential.resize(differential.cols() * differential.rows(), 1);
                VField vectorField(differential);
                OForm oneform(vectorField);

                return oneform;
            }

            static OForm computeDifferential(_Sim &simulator) {
                OForm volumeDifferential = computeVolumeDifferential(simulator);

                return SDConversions::diff_bdry_from_diff_vol(volumeDifferential, simulator);
            }

            template<class _Iterate>
            PeriodicScaleInvariantSmoothingRegularizationTerm(_Iterate &it) : m_sim(it.simulator()) {
                m_smoothingField = computeSmoothingField(it.simulator());
                this->m_differential = computeDifferential(it.simulator());
            }

            // Evaluate objective function (without weight)
            virtual Real evaluate() const override {
                Real result = computeCost(m_sim);
                return result;
            }

            virtual void writeFields(MSHFieldWriter &writer) const override {
                writer.addField("Scale Invariant Smoothing", m_smoothingField);

                auto bdryVel = SDConversions::descent_from_diff_bdry(this->m_differential, m_sim);
                VField xferBdryVel(m_sim.mesh().numVertices());
                xferBdryVel.clear();
                for (auto v : m_sim.mesh().vertices()) {
                    auto bv = v.boundaryVertex();
                    if (!bv) continue;
                    xferBdryVel(v.index()) = bdryVel(bv.index());
                }

                writer.addField("SI Smoothing Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
            }

            virtual ~PeriodicScaleInvariantSmoothingRegularizationTerm() { }
        private:
            SField m_smoothingField;
            _Sim &m_sim;
        };

        // Configuration to be applied by iterate factory
        template<class _Sim>
        struct IFConfigPeriodicSISmoothingRegularization: public IFConfig {
            template<class _Iterate>
            void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
                auto sr = Future::make_unique<PeriodicScaleInvariantSmoothingRegularizationTerm<_Sim>>(*it);
                sr->setWeight(weight);

                if (!normalizations.isSet("PeriodicScaleInvariantSmoothingRegularization"))
                    normalizations.set("PeriodicScaleInvariantSmoothingRegularization", 1.0);

                sr->setNormalization(normalizations["PeriodicScaleInvariantSmoothingRegularization"]);
                it->addObjectiveTerm("PeriodicScaleInvariantSmoothingRegularization", std::move(sr));
            }
            Real weight = 0.0;
        };

    }
}

#endif //SCALEINVARIANTSMOOTHINGREGULARIZATION_H