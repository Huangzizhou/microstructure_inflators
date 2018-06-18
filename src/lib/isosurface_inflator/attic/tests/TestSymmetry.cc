#include "Symmetry.hh"

#include <iostream>
using namespace std;

template<class SGroup>
vector<size_t> getOperationCounts(const SGroup &g) {
    vector<size_t> counts;
    for (auto &isometry : g) {
        size_t length = isometry.operations.size();
        if (length >= counts.size())
            counts.resize(length + 1);
        counts[length] += 1;
    }
    return counts;
}

template<class SGroup>
void report(const SGroup &g) {
    cout << g.size() << ":";
    for (size_t c : getOperationCounts(g)) {
        cout << " " << c;
    }
    cout << endl;
}

void merge_duplicates(vector<vector<Vector3d>> &nodeSets) {
    for (auto &nodeSet : nodeSets) {
        vector<bool> isDuplicate(nodeSet.size(), false);
        // Mark the the second (and later) copies of duplicated vertices
        for (size_t i = 0; i < nodeSet.size(); ++i) {
            for (size_t j = i + 1; j < nodeSet.size(); ++j) {
                if ((nodeSet[i] - nodeSet[j]).norm() < 1e-6)
                    isDuplicate[j] = true;
            }
        }

        // Replace node set with only the non-duplicates
        std::vector<Vector3d> fullSet;
        fullSet.swap(nodeSet);
        for (size_t i = 0; i < isDuplicate.size(); ++i)
            if (!isDuplicate[i]) nodeSet.push_back(fullSet[i]);
    }
}

void report(const string &name, const vector<vector<Vector3d>> &nodeSets) {
    cout << name << " node copies:";
    for (const auto &ns : nodeSets)
        cout << " " << ns.size();
    cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    report(Symmetry::TriplyPeriodic<>::symmetryGroup());
    report(Symmetry::Orthotropic<>::symmetryGroup());
    report(Symmetry::Cubic<>::symmetryGroup());

    size_t count = 0;
    auto g = Symmetry::Cubic<>::symmetryGroup();
    vector<Vector3d> vertexNodes = { Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0), Vector3d(1.0, 1.0, 0.0), Vector3d(1.0, 1.0, 1.0) };
    vector<Vector3d> edgeNodes = {
        Vector3d(0.5, 0.0, 0.0), Vector3d(0.5, 0.5, 0.0), Vector3d(0.5, 0.5, 0.5),  // 0-1 0-2 0-3
        Vector3d(1.0, 0.5, 0.0), Vector3d(1.0, 0.5, 0.5), Vector3d(1.0, 1.0, 0.5),  // 1-2 1-3 2-3
    };
    vector<Vector3d> faceNodes = { // 012 013 023 123
        Vector3d(2.0 / 3.0, 1.0 / 3.0, 0.0 / 3.0), // 012
        Vector3d(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), // 013
        Vector3d(2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0), // 023
        Vector3d(3.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0), // 123
    };

    vector<Vector3d> volumeNodes = { Vector3D(0.75, 0.5, 0.25) };

    vector<vector<Vector3d>> vertexCopies(vertexNodes.size());
    vector<vector<Vector3d>> edgeCopies(edgeNodes.size());
    vector<vector<Vector3d>> faceCopies(faceNodes.size());
    vector<vector<Vector3d>> volumeCopies(volumeNodes.size());

    for (const auto &isometry : g) {
        bool hasTranslation = false;
        for (auto op : isometry.operations) {
            if (dynamic_pointer_cast<const Isometry::Translation>(op)) hasTranslation = true;
        }
        if (hasTranslation) continue;
        ++count;

        for (size_t i = 0; i < vertexNodes.size(); ++i)
            vertexCopies[i].push_back(isometry.apply(vertexNodes[i]));
        for (size_t i = 0; i < edgeNodes.size(); ++i)
            edgeCopies[i].push_back(isometry.apply(edgeNodes[i]));
        for (size_t i = 0; i < faceNodes.size(); ++i)
            faceCopies[i].push_back(isometry.apply(faceNodes[i]));
        for (size_t i = 0; i < volumeNodes.size(); ++i)
            volumeCopies[i].push_back(isometry.apply(volumeNodes[i]));
    }

    merge_duplicates(vertexCopies);
    merge_duplicates(edgeCopies);
    merge_duplicates(faceCopies);
    merge_duplicates(volumeCopies);

    report("vertex", vertexCopies);
    report("edge", edgeCopies);
    report("face", faceCopies);
    report("volume", volumeCopies);

    cout << "non-translational cubic isometries: " << count << endl;

    return 0;
}
