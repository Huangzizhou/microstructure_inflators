#ifndef OCTACELLPATTERN_H
#define OCTACELLPATTERN_H

#include "EdgeMeshPattern.h"
#include <vcg/space/index/index2D/grid_static_ptr_2D.h>

#include <utility>
#include <array>
#include <vector>
#include <unordered_map>
#include <map>

namespace std
{
template <>
struct hash<pair<int, int> >
{
	typedef pair<int, int> argument_type;
	typedef size_t result_type;

	result_type operator()(const argument_type& a) const
	{
		hash<int> hasher;
		argument_type e = a;
		if (e.first > e.second)
			std::swap(e.first, e.second);

		return hasher(e.first) * 31 + hasher(e.second);
	}
};
}

// PATTERN CONFIGURATION
//
//    (Y)
//     ^
// 1.0 |             F
//     |             |
//     |             |
// 0.75|      A______|______B
//     |      |      J      |
//     |      |             |
//     |      |             |
//  0.5|E-----I             K-----G
//     |      |             |
//     |      |             |
//     |      D_____ L _____C
// 0.25|             |
//     |             |
//     |             |
//     +-------------H-------------> (X)
//   0,0    0.25    0.5   0.75   1.0

struct OctaCellParameters : public CellParameters<8>
{

	// Thickness parameters
	// 0 - thicknes JL
	// 1 - thicknes FH
	// 2 - thicknes IK
	// 3 - thicknes EG
	// 4 - thicknes ABCD

	// Offset parameters
	// 5 - offset JL
	// 6 - offset IK
	// 7 - offset ABCD

	double params[NumberOfParameters];

	OctaCellParameters(void)
	{
		for (int i=0; i<NumberOfParameters; ++i)
		{
			params[i] = (i<5) ? 0.1 : 0;
		}
	}

	double & parameter(int i)
	{
		assert(i>=0 && i<NumberOfParameters);
		return params[i];
	}

	const double & cParameter(int i) const
	{
		assert(i>=0 && i<NumberOfParameters);
		return params[i];
	}

	std::pair<double,double> parameterRange(int i) const
	{
		assert(i>=0 && i<NumberOfParameters);

		if (i<5)
		{
			return std::make_pair(0.01, 0.1);
		}
		else
		{
			return std::make_pair(-0.15, 0.15);
		}
	}
};

template <class TriMesh>
class OctaCellPattern : public EdgeMeshPattern<TriMesh, OctaCellParameters>
{
public:
	typedef OctaCellPattern<TriMesh>                     ThisType;
	typedef EdgeMeshPattern<TriMesh, OctaCellParameters> BaseType;

	typedef typename BaseType::PatternParameters PatternParameters;
	typedef typename BaseType::ScalarType        ScalarType;
	typedef typename BaseType::Coord2Type        Coord2Type;

	void computeVelocityField(TriMesh & mesh,
	                          std::unordered_map<std::pair<int, int>, std::array<ScalarType, PatternParameters::NumberOfParameters> > & fields)
	{
		typedef typename EMesh::CoordType CoordType;

		// get the generating edge mesh
		EMesh em;
		this->getEdgeMesh(em);

		// update border topology
		vcg::tri::UpdateTopology<TriMesh>::FaceFace(mesh);
		vcg::tri::UpdateFlags<TriMesh>::FaceBorderFromFF(mesh);
		vcg::tri::UpdateFlags<TriMesh>::VertexBorderFromFace(mesh);

		// get all border edges and their normals
		std::unordered_map<std::pair<int,int>, vcg::Point2d> edges;
		for (size_t i=0; i<mesh.face.size(); ++i)
		{
			typename TriMesh::FaceType & f = mesh.face[i];
			for (char j=0; j<3; ++j)
			{
				if (vcg::face::IsBorder(f, j))
				{
					const typename TriMesh::VertexType * v0 = f.cV(j);
					const typename TriMesh::VertexType * v1 = f.cV((j+1)%3);
					int idx0 = vcg::tri::Index(mesh, v0);
					int idx1 = vcg::tri::Index(mesh, v1);
					std::pair<int, int> e = { idx0, idx1 };

					if (edges.count(e) == 0)
					{
						// compute normal of edge
						CoordType n = CoordType::Construct((v1->cP() - v0->cP()).Normalize());
						edges[e] = vcg::Point2d(n[1], -n[0]); // rotated by -90 degrees
					}
				}
			}
		}

		// collect all border vertex pointers
		std::vector<typename TriMesh::VertexPointer> vtx;
		for (size_t i=0; i<mesh.vert.size(); ++i)
		{
			typename TriMesh::VertexType & v = mesh.vert[i];
			if (v.IsB())
				vtx.push_back(&v);
		}

		// compute closest segment acceleration structure
		std::vector<SegmentType> segments;
		for (size_t i=0; i<em.edge.size(); ++i)
		{
			// generate the segments to the edges
			for (size_t i=0; i<em.edge.size(); ++i)
			{
				// get the edge
				const EMesh::EdgeType & e = em.edge[i];
				if (e.IsD())
					continue;

				// get their vertices
				const EMesh::VertexType & v0 = *e.cV(0);
				const EMesh::VertexType & v1 = *e.cV(1);

				// get the vertices coordinate
				Coord2Type p0(ScalarType(v0.cP()[0]), ScalarType(v0.cP()[1]));
				Coord2Type p1(ScalarType(v1.cP()[0]), ScalarType(v1.cP()[1]));

				// get radius per vertex
				double r0 = v0.cQ();
				double r1 = v1.cQ();
				assert(r0 > 0 && r1 > 0);

				// rotated vector orthogonal to edge direction
				Coord2Type edge = (p0-p1);
				double d = edge.Norm();
				Coord2Type ortho_dir = edge.Normalize();
				ortho_dir = Coord2Type(-ortho_dir.Y(), ortho_dir.X());

				double rotation = asin((r0 - r1)/d);

				Coord2Type leftRot, rightRot;
				leftRot = rightRot = ortho_dir;
				leftRot.Rotate(rotation);
				rightRot.Rotate(-rotation);

				segments.push_back(SegmentType(p0 + (leftRot * r0),  p1 + (leftRot * r1),  vcg::tri::Index(em, v0), vcg::tri::Index(em, v1)));
				segments.push_back(SegmentType(p0 - (rightRot * r0), p1 - (rightRot * r1), vcg::tri::Index(em, v0), vcg::tri::Index(em, v1)));
			}

			// generate vertex circles
			for (size_t i=0; i<em.vert.size(); ++i)
			{
				const EMesh::VertexType & v = em.vert[i];
				if (v.IsD())
					continue;

				std::vector<Coord2Type> c;
				this->generateCircle(c, Coord2Type(v.cP()[0], v.cP()[1]), v.cQ());

				for (size_t k=0; k<c.size(); k++)
				{
					segments.push_back(SegmentType(c[k], c[(k+1)%c.size()], vcg::tri::Index(em, &v)));
				}
			}
		}

		// compute map between border vertices and vertex-/edge-segments
		vcg::GridStaticPtr2D<SegmentType, ScalarType> grid2D;
		grid2D.Set(segments.begin(), segments.end());
		std::vector<std::pair<const SegmentType *, Coord2Type> > vtxToSegs;
		vtxToSegs.resize(vtx.size());
		SegmentMarker marker;
		for (size_t i=0; i<vtx.size(); ++i)
		{
			const typename TriMesh::VertexPointer v = vtx[i];
			Coord2Type & closest = vtxToSegs[i].second;
			vtxToSegs[i].first = GetClosestSegment(grid2D, Coord2Type(v->cP()[0], v->cP()[1]), marker, closest);
//			vtxToSegs[i].first = GetClosestSegmentBruteF(segments, Coord2Type(v->cP()[0], v->cP()[1]), closest);
			assert(vtxToSegs[i].first != NULL);
		}

		// compute displacement per vertex changing each parameter
		std::vector<ParameterChange> params_op;
		this->getParameterOperations(params_op);
		assert(params_op.size() == PatternParameters::NumberOfParameters);

		fields.clear();
		typedef typename TriMesh::VertexType::NormalType NormalType;
		for (size_t p=0; p<PatternParameters::NumberOfParameters; ++p)
		{
			this->m_params.patternParams.cParameter(p);
			ParameterChange & par = params_op[p];

			// for each parameter compute the displacement per vertex
			if (par.type == ParameterChange::Radius)
			{
				for (size_t i=0; i<vtxToSegs.size(); ++i)
				{
					typename TriMesh::VertexPointer v = vtx[i];
					const SegmentType * s       = vtxToSegs[i].first;
					const Coord2Type  & closest = vtxToSegs[i].second;

					switch (s->belonging) {
					case SegmentType::EdgeSegment :
					{
						NormalType n0 = (par.node_ops.count(NodeID(s->index0)) == 0) ?
						                    NormalType(0,0,0) :
						                    NormalType::Construct(CoordType(s->P0()[0], s->P0()[1], 0) - em.vert[s->index0].cP()).Normalize();

						NormalType n1 = (par.node_ops.count(NodeID(s->index1)) == 0) ?
						                    NormalType(0,0,0) :
						                    NormalType::Construct(CoordType(s->P1()[0], s->P1()[1], 0) - em.vert[s->index1].cP()).Normalize();

						ScalarType t =  s->interpolationParameter(closest);
						v->N() = n0 * (1-t) + n1 * t;

						break;
					}
					case SegmentType::VertexSegment :
					{
						v->N() = (par.node_ops.count(NodeID(s->index0)) == 0) ?
						             NormalType(0,0,0) :
						             (NormalType::Construct(v->cP()) - NormalType::Construct(em.vert[s->index0].cP())).Normalize();
						break;
					}
					default:
						assert(0);
						break;
					}
				}
			}
			else if (par.type == ParameterChange::Translation)
			{
				// for each parameter compute the displacement per vertex
				for (size_t i=0; i<vtxToSegs.size(); ++i)
				{
					typename TriMesh::VertexPointer v = vtx[i];
					const SegmentType * s       = vtxToSegs[i].first;
					const Coord2Type  & closest = vtxToSegs[i].second;

					switch (s->belonging) {
					case SegmentType::EdgeSegment :
					{
						NormalType n0 = (par.node_ops.count(NodeID(s->index0)) == 0) ?
						                    NormalType(0,0,0) :
						                    NormalType(par.node_ops[NodeID(s->index0)][0], par.node_ops[NodeID(s->index0)][1], 0);

						NormalType n1 = (par.node_ops.count(NodeID(s->index1)) == 0) ?
						                    NormalType(0,0,0) :
						                    NormalType(par.node_ops[NodeID(s->index1)][0], par.node_ops[NodeID(s->index1)][1], 0);

						ScalarType t =  s->interpolationParameter(closest);
						v->N() = n0 * (1-t) + n1 * t;
						break;
					}
					case SegmentType::VertexSegment :
					{
						v->N() = (par.node_ops.count(NodeID(s->index0)) == 0) ?
						             NormalType(0,0,0) :
						             NormalType(par.node_ops[NodeID(s->index1)][0], par.node_ops[NodeID(s->index1)][1], 0);
						break;
					}
					default:
						assert(0);
						break;
					}
				}
			}
			else
			{
				assert(0);
			}

			// compute edge velocity field
			for (auto & e : edges)
			{
				NormalType n0 = NormalType::Construct(mesh.vert[e.first.first].cN());
				NormalType n1 = NormalType::Construct(mesh.vert[e.first.second].cN());
				ScalarType val = (e.second.dot(vcg::Point2d(n0[0], n0[1])) + e.second.dot(vcg::Point2d(n1[0], n1[1]))) / 2;
				std::array<ScalarType, PatternParameters::NumberOfParameters> & f = fields[e.first];
				f[p] = val;

//				std::cout << e.first.first << " " << e.first.second << ": " << val << std::endl << std::flush;
			}
		}
	}

protected:
	typedef typename BaseType::ParameterChange ParameterChange;

	typedef enum
	{
		A = 0,
		B,
		C,
		D,
		E,
		F,
		G,
		H,
		I,
		J,
		K,
		L,
	} NodeID;

	void getParameterOperations(std::vector<ParameterChange> & params_op) const
	{
		params_op.clear();
		params_op.resize(PatternParameters::NumberOfParameters);
		{
			ParameterChange & p = params_op[0];
			p.type = ParameterChange::Radius;
			p.node_ops[J] = p.node_ops[L] = Coord2Type();
		}
		{
			ParameterChange & p = params_op[1];
			p.type = ParameterChange::Radius;
			p.node_ops[F] = p.node_ops[H] = Coord2Type();
		}
		{
			ParameterChange & p = params_op[2];
			p.type = ParameterChange::Radius;
			p.node_ops[I] = p.node_ops[K] = Coord2Type();
		}
		{
			ParameterChange & p = params_op[3];
			p.type = ParameterChange::Radius;
			p.node_ops[E] = p.node_ops[G] = Coord2Type();
		}
		{
			ParameterChange & p = params_op[4];
			p.type = ParameterChange::Radius;
			p.node_ops[A] = p.node_ops[B] = p.node_ops[C] = p.node_ops[D] = Coord2Type();
		}
		{
			ParameterChange & p = params_op[5];
			p.type = ParameterChange::Translation;
			Coord2Type JLDispl(0,1);
			p.node_ops[J] =  JLDispl;
			p.node_ops[L] = -JLDispl;
		}
		{
			ParameterChange & p = params_op[6];
			p.type = ParameterChange::Translation;
			Coord2Type IKDispl(1,0);
			p.node_ops[I] = -IKDispl;
			p.node_ops[K] =  IKDispl;
		}
		{
			ParameterChange & p = params_op[7];
			p.type = ParameterChange::Translation;
			p.node_ops[A] = Coord2Type(-M_SQRT1_2,  M_SQRT1_2);
			p.node_ops[B] = Coord2Type( M_SQRT1_2,  M_SQRT1_2);
			p.node_ops[C] = Coord2Type( M_SQRT1_2, -M_SQRT1_2);
			p.node_ops[D] = Coord2Type(-M_SQRT1_2, -M_SQRT1_2);
		}
	}

	void generateUnmodifiedEdgeMesh(EMesh & em)
	{
		typedef typename EMesh::CoordType CoordType;

		em.Clear();

		// Generate nodes
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(0.25, 0.75, 0)); // A
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(0.75, 0.75, 0)); // B
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(0.75, 0.25, 0)); // C
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(0.25, 0.25, 0)); // D

		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(   0,  0.5, 0)); // E
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType( 0.5,  1.0, 0)); // F
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType( 1.0,  0.5, 0)); // G
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType( 0.5,    0, 0)); // H

		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(0.25,  0.5, 0)); // I
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType( 0.5, 0.75, 0)); // J
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType(0.75,  0.5, 0)); // K
		vcg::tri::Allocator<EMesh>::AddVertex(em, CoordType( 0.5, 0.25, 0)); // L

		// Generate edges
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[A], &em.vert[J]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[J], &em.vert[B]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[B], &em.vert[K]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[K], &em.vert[C]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[C], &em.vert[L]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[L], &em.vert[D]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[D], &em.vert[I]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[I], &em.vert[A]);

		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[E], &em.vert[I]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[F], &em.vert[J]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[G], &em.vert[K]);
		vcg::tri::Allocator<EMesh>::AddEdge(em, &em.vert[H], &em.vert[L]);
	}

	// Segment type for accelerated query
	class SegmentType : public vcg::Segment2<ScalarType>
	{
	public:
		typedef typename vcg::Segment2<ScalarType> BaseType;

		typedef enum
		{
			VertexSegment,
			EdgeSegment
		} Type;

		int mark;
		size_t index0;
		size_t index1;
		Type belonging;

		bool IsD() { return false; }

		SegmentType(void) { ; }

		SegmentType(const vcg::Point2<ScalarType> & _P0,
		            const vcg::Point2<ScalarType> & _P1,
		            size_t index)
		    : BaseType(_P0, _P1)
		    , mark(0)
		    , index0(index)
		    , index1(index)
		    , belonging(VertexSegment)
		{
			;
		}

		SegmentType(const vcg::Point2<ScalarType> & _P0,
		            const vcg::Point2<ScalarType> & _P1,
		            size_t idx0,
		            size_t idx1)
		    : BaseType(_P0, _P1)
		    , mark(0)
		    , index0(idx0)
		    , index1(idx1)
		    , belonging(EdgeSegment)
		{
			;
		}

		SegmentType(const SegmentType & s1)
		    : BaseType(s1)
		{
			this->index0 = s1.index0;
			this->index1 = s1.index1;
			this->belonging = s1.belonging;
			this->mark = s1.mark;
		}

		int & TMark() { return mark; }

		ScalarType interpolationParameter(const vcg::Point2<ScalarType> & p) const
		{
			vcg::Point2<ScalarType> v(p - this->P0());
			vcg::Point2<ScalarType> v_base(this->P1() - this->P0());

			return ScalarType(v.dot(v_base))/v_base.SquaredNorm();
		}
	};

	// marker type
	class SegmentMarker
	{
	public:
		int mark;

		SegmentMarker() { mark=0; }

		void UnMarkAll() { mark++; }

		bool IsMarked(SegmentType * obj)
		{
			int markObj = obj->TMark();
			return (markObj == mark);
		}

		void Mark(SegmentType * obj)
		{
			obj->TMark() = mark;
		}
	};

	// query
	static SegmentType * GetClosestSegment(vcg::GridStaticPtr2D<SegmentType, ScalarType> & grid,
	                                       const Coord2Type & _p,
	                                       SegmentMarker & marker,
	                                       Coord2Type &_closestPt)
	{
		vcg::PointSegment2DEPFunctor<ScalarType> PDistFunct;

		ScalarType _minDist = 0;
		ScalarType _maxDist = std::numeric_limits<ScalarType>::max();
		return (grid.GetClosest(PDistFunct, marker, _p, _maxDist, _minDist, _closestPt));
	}

	static const SegmentType * GetClosestSegmentBruteF(const std::vector<SegmentType> & v,
	                                                   const Coord2Type & _p,
	                                                   Coord2Type &_closestPt)
	{
		ScalarType _minDist = std::numeric_limits<ScalarType>::max();
		const SegmentType * ret = NULL;
		for (size_t i=0; i<v.size(); i++)
		{
			vcg::Point2<ScalarType> test;
			test = vcg::ClosestPoint(v[i], _p);
			ScalarType currD = (test-_p).Norm();
			if (currD < _minDist)
			{
				_closestPt = test;
				_minDist = currD;
				ret = &v[i];
			}
		}
		return ret;
	}
};

#endif // OCTACELLPATTERN_H
