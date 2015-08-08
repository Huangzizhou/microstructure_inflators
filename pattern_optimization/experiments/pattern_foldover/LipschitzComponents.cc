////////////////////////////////////////////////////////////////////////////////
// LipschitzComponents.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Extract the Lipschitz-continuous components of an lookup table.
//      Entries are considered in the same component if they satisfy the
//      Lipschitz continuity criterion:
//          d_param(p_1, p_2) <= lipschitz_constant * d_mat(Ch(p1), Ch(p2))
//      *OR* if their parameter distance is below the specified threshold:
//          d_param(p_1, p_2) < threshold
//      where d_param and d_mat are some norm-induced metrics.
//      We nondimensionalize by dividing each side of the inequalities by the
//      norm of p1/Ch(p1).
//
//      We start by putting the "midpoint" of the material property space in its
//      own component. We sort the rest of the points by material-space
//      distance and attempt to add them to a component one at a time. If the
//      point belongs to none of the existing components, a new component is
//      created.
//
//      When deciding if a point belongs to a component, we apply the continuity
//      criteria above against every point in the component.
//
//      The material property midpoint and d_param are computed in (log E, nu)
//      space to reflect the sampling strategy.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/30/2015 19:01:04
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "../../LookupTable.hh"
#include "utils.hh"

using namespace std;

typedef IsotropicLookupTable<double> LUT;

double nuWeight = 10.0;
double materialNorm(double logE, double nu) { return sqrt(logE * logE + (nuWeight * nu) * (nuWeight * nu)); }
double materialDist(double logE1, double nu1, double logE2, double nu2) { return materialNorm(logE1 - logE2, nu1 - nu2); }

// Compute the logspace dist to indexed points of logE, nu
vector<double> distToPoints(double ptLogE, double ptNu,
        const vector<double> &logEs, const vector<double> &nus,
        const vector<size_t> &indices = vector<size_t>())
{
    std::vector<double> dists;
    if (indices.size() == 0) {
        // All points
        dists.reserve(nus.size());
        for (size_t i = 0; i < nus.size(); ++i)
            dists.push_back(materialDist(logEs[i], nus[i], ptLogE, ptNu));
    }
    else {
        dists.reserve(indices.size());
        for (size_t i : indices)
            dists.push_back(materialDist(logEs[i], nus[i], ptLogE, ptNu));
    }

    return dists;
}

int main(int argc, char *argv[])
{
    if (argc != 5) {
        cerr << "Usage: ./LipschitzComponents dataTable pattern lipschitz_constant threshold" << endl;
        exit(-1);
    }

    LUT lut(argv[1]);
    LUT patternLut = lut.selectPattern(stoi(argv[2]));
    double lipschitz_constant = stod(argv[3]);
    double threshold = stod(argv[4]);

    auto Es = patternLut.getEs<double>();
    auto nus = patternLut.getNus<double>();
    auto params = patternLut.getParamValues<double>();
    auto logEs = Es;
    for (double &e : logEs) e = log(e);

    double logEAvg = 0, nuAvg = 0;
    for (double &le : logEs) logEAvg += le;
    for (double &nu : nus) nuAvg += nu;

    logEAvg /= logEs.size();
    nuAvg /= nus.size();

    std::vector<double> dists = distToPoints(logEAvg, nuAvg, logEs, nus);

    vector<vector<size_t>> components;

    // Visit points in increasing distance from midpoint
    vector<size_t> perm;
    sortPermutation(dists, perm);
    size_t count = 0;
    for (size_t i : perm) {
        vector<size_t> containingComponents;
        // vector<double> componentMatDists;
        // vector<double> componentParamDists;
        for (size_t c = 0; c < components.size(); ++c) {
            // double ptLogE = logEs[i], ptNu = nus[i];
            // dists = distToPoints(ptLogE, ptNu, logEs, nus, components[c]);
            // auto it = min_element(dists.begin(), dists.end());
            // size_t closestCompIdx = std::distance(dists.begin(), it);
            // size_t closestIdx = components[c].at(closestCompIdx);


            // // Continuity criteria
            // double matDist = materialDist(logEs[closestIdx], nus[closestIdx], ptLogE, ptNu);
            // double paramDist = (params.col(closestIdx) - params.col(i)).norm();
            // if ((paramDist < threshold) || (paramDist < lipschitz_constant * matDist)) {
            //     containingComponents.push_back(c);
            // }

            bool containedInC = true;
            for (size_t ci : components[c]) {
                double matDist = materialDist(logEs[ci], nus[ci], logEs[i], nus[i]);
                double paramDist = (params.col(ci) - params.col(i)).norm();
                if ((ci == 36839 || ci == 36844) && (i == 36839 || i == 36844)) {
                    cout << "matDist: " << matDist << ", paramDist: " << paramDist << endl;
                    cout << "line 1: " <<  Es[ci] << ", " << nus[ci] << ", " << params.col(ci)[0] << ", " << params.col(ci)[1] << ", " << params.col(ci)[2] << ", " << params.col(ci)[3] << ", " << params.col(ci)[4] << ", " << params.col(ci)[5] << endl;
                    cout << "line 2: " <<  Es[i ] << ", " << nus[i ] << ", " << params.col(i )[0] << ", " << params.col(i )[1] << ", " << params.col(i )[2] << ", " << params.col(i )[3] << ", " << params.col(i )[4] << ", " << params.col(i )[5] << endl;
                }
                if (!((paramDist < threshold) || (paramDist < lipschitz_constant * matDist))) {
                    containedInC = false;
                    break;
                }
            }
            if (containedInC)
                containingComponents.push_back(c);
        }
        if (containingComponents.size() == 0) {
            // New component
            components.emplace_back(1, i);
            cout << "New component: " << components.size() << endl;
            // cout << "mat dists:";
            // for (double md : componentMatDists) cout << " " << md;
            // cout << endl << "param dists:";
            // for (double pd : componentParamDists) cout << " " << pd;
            // cout << endl;
        }
        else {
            if (containingComponents.size() > 1) {
                // std::cerr << "Warning: " << Es[i] << ", " << nus[i]
                //           << " belongs to multiple components." << endl;
            }
            components.at(containingComponents.at(0)).push_back(i);
        }
        if (count++ % (perm.size() / 10) == 1) {
            cout << "Processed " << count << "/" << perm.size() << endl;
        }
    }

    for (size_t c = 0; c < components.size(); ++c) {
        patternLut.subset(components[c]).write("component" + to_string(c) + ".txt");
    }

    return 0;
}
