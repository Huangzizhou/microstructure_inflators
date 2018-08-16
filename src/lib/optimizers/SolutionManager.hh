//
// Created by Davi Colli Tozoni on 6/6/18.
//

#ifndef BESTSOLUTIONMANAGER_H
#define BESTSOLUTIONMANAGER_H

#include "IterateManagerBase.hh"
#include "BoundConstraints.hh"
#include "OptimizerConfig.hh"

// Class that manages best solution for all solvers.
class SolutionManager {
public:
    SolutionManager(PatternOptimization::IterateManagerBase &im) : m_im(im), m_best(std::numeric_limits<Real>::max()) { }

    bool areNewParameters(const std::vector<Real> &params);

    bool hasViableSolution(const std::vector<Real> &params);

    bool isImprovement(Real val, bool respectConstraints);

    void updateAndReport(const std::vector<Real> &x);

    std::vector<Real> m_prevParams;
    size_t m_niters = 0;
    size_t m_bestIter = 0;
    PatternOptimization::IterateManagerBase &m_im;
    Real m_best;
    std::string m_outPath;
};

#endif //BESTSOLUTIONMANAGER_H
