#ifndef EDGEMESHTYPE_H
#define EDGEMESHTYPE_H

#include <vcg/complex/complex.h>

class EEdgeType;
class EVertexType;

struct EUsedTypes : public vcg::UsedTypes< vcg::Use<EVertexType>::AsVertexType,
                                           vcg::Use<EEdgeType>::AsEdgeType> {};

class EVertexType : public vcg::Vertex< EUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::Qualityd, vcg::vertex::BitFlags, vcg::vertex::VEAdj> {};
class EEdgeType   : public vcg::Edge< EUsedTypes, vcg::edge::VertexRef, vcg::edge::Qualityd, vcg::edge::BitFlags, vcg::edge::EEAdj, vcg::edge::VEAdj> {};
class EMesh       : public vcg::tri::TriMesh< std::vector<EVertexType>, std::vector<EEdgeType> > {};

#endif // EDGEMESHTYPE_H
