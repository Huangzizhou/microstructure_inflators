#include "WireMesh.hh"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>> wmesh("pattern0746.wire");
    cout << wmesh.numVertices() << endl;
    cout << wmesh.numEdges() << endl;
    cout << wmesh.numBaseVertices() << endl;
    cout << wmesh.numBaseEdges() << endl;

    wmesh.saveBaseUnit("unit.wire");
    wmesh.save("scaled.wire");
    wmesh.saveReplicatedBaseUnit("replicated.wire");
    wmesh.saveInflationGraph("igraph.wire");
    return 0;
}
