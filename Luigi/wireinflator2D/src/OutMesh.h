#ifndef OUTMESH_H
#define OUTMESH_H

#include <array>
#include <vector>
#include <unordered_map>

template <size_t Dim = 3, size_t NodesPerElement = 3>
struct OutMesh
{
	typedef int                                    IndexType;
	typedef double                                 ScalarType;
	typedef std::array<ScalarType, Dim>            CoordType;
	typedef std::array<ScalarType, Dim>            NormalType;
	typedef std::array<IndexType, NodesPerElement> FaceType;
	typedef std::pair<IndexType, IndexType>        EdgeType;

	typedef std::vector<CoordType>  NodeVector;
	typedef std::vector<FaceType>   ElementVector;

	typedef std::vector<ScalarType>              Fields;
	typedef std::unordered_map<EdgeType, Fields> EdgeFields;

	NodeVector    nodes;
	ElementVector elements;
	EdgeFields    edge_fields;
};

#endif // OUTMESH_H
