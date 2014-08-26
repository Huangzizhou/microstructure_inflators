#ifndef EDGEMESHTYPE_H
#define EDGEMESHTYPE_H

#include <vcg/complex/complex.h>

class EdgeType;
class VertexType;

struct EUsedTypes : public vcg::UsedTypes< vcg::Use<VertexType>::AsVertexType, vcg::Use<EdgeType>::AsEdgeType > {};

class VertexType : public vcg::Vertex< EUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::Qualityd, vcg::vertex::BitFlags, vcg::vertex::VEAdj> {};
class EdgeType   : public vcg::Edge< EUsedTypes, vcg::edge::VertexRef, vcg::edge::Qualityd, vcg::edge::BitFlags, vcg::edge::EEAdj, vcg::edge::VEAdj> {};
class EMesh      : public vcg::tri::TriMesh< std::vector<VertexType>, std::vector<EdgeType> > {};

#endif // EDGEMESHTYPE_H
