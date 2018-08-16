//
// Created by Davi Colli Tozoni on 8/16/18.
//

#include "SolutionManager.hh"

#define EQ_TOL 1e-2
#define INEQ_TOL 1e-14

bool SolutionManager::areNewParameters(const std::vector<Real> &params) {
    bool differ = true;
    if (m_prevParams.size() == params.size()) {
        differ = false;
        for (size_t i = 0; i < params.size(); ++i) {
            if (std::abs(m_prevParams[i] - params[i]) > 1e-9) {
                differ = true;
                break;
            }
        }
    }
    m_prevParams = params;
    if (!differ) return false;

    ++m_niters;
    return true;
}

bool SolutionManager::hasViableSolution(const std::vector<Real> &params) {
    auto &it = m_im.get(params.size(), params.data());

    return it.hasViableSolution(EQ_TOL, INEQ_TOL);
}

bool SolutionManager::isImprovement(Real val, bool respectConstraints) {
    if (val < m_best && respectConstraints) {
        m_best = val;
        m_bestIter = m_niters;
        return true;
    }
    return false;
}

void SolutionManager::updateAndReport(const std::vector<Real> &x) {
    auto &it = m_im.get();
    double val = it.evaluate();

    if (areNewParameters(x) && it.shouldReport()) {
        std::cout << "Iteration: " << m_niters << std::endl;
        it.writeDescription(std::cout);

        bool successful = isImprovement(val, it.hasViableSolution(EQ_TOL, INEQ_TOL));
        if (successful) {
            std::cout << "\\o/ GREAT! Found a good solution that is also the best so far..." << std::endl;

            // Uses the opportunity to update the inflators, if needed
            m_im.update();
        }
        else {
            std::cout << "=( Not the best solution so far..." << std::endl;
        }

        std::cout << "Write mesh" << std::endl;
        if (m_outPath != "") {
            it.writeMeshAndFields(m_outPath + "_" + std::to_string(m_niters));

            // If succeeded, then also print two extra files with information about current best solution
            if (successful) {
                std::string infoPath = m_outPath + "_sol.txt";
                std::ofstream out(infoPath);
                out << "Iteration: " << m_niters << std::endl;
                it.writeDescription(out);
                it.writeMeshAndFields(m_outPath + "_best");
                out.close();
            }
        }

        std::cout << std::endl;
    }
}
