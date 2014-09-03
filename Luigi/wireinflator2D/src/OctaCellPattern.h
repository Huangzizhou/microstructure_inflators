#ifndef OCTACELLPATTERN_H
#define OCTACELLPATTERN_H

#include "EdgeMeshPattern.h"

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

protected:


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
			p.node_ops[J] = p.node_ops[L] = vcg::Point2d();
		}
		{
			ParameterChange & p = params_op[1];
			p.type = ParameterChange::Radius;
			p.node_ops[F] = p.node_ops[H] = vcg::Point2d();
		}
		{
			ParameterChange & p = params_op[2];
			p.type = ParameterChange::Radius;
			p.node_ops[I] = p.node_ops[K] = vcg::Point2d();
		}
		{
			ParameterChange & p = params_op[3];
			p.type = ParameterChange::Radius;
			p.node_ops[E] = p.node_ops[G] = vcg::Point2d();
		}
		{
			ParameterChange & p = params_op[4];
			p.type = ParameterChange::Radius;
			p.node_ops[A] = p.node_ops[B] = p.node_ops[C] = p.node_ops[D] = vcg::Point2d();
		}
		{
			ParameterChange & p = params_op[5];
			p.type = ParameterChange::Translation;
			vcg::Point2d JLDispl(0,1);
			p.node_ops[J] =  JLDispl;
			p.node_ops[L] = -JLDispl;
		}
		{
			ParameterChange & p = params_op[6];
			p.type = ParameterChange::Translation;
			vcg::Point2d IKDispl(1,0);
			p.node_ops[I] = -IKDispl;
			p.node_ops[K] =  IKDispl;
		}
		{
			ParameterChange & p = params_op[7];
			p.type = ParameterChange::Translation;
			p.node_ops[A] = vcg::Point2d(-M_SQRT1_2,  M_SQRT1_2);
			p.node_ops[B] = vcg::Point2d( M_SQRT1_2,  M_SQRT1_2);
			p.node_ops[C] = vcg::Point2d( M_SQRT1_2, -M_SQRT1_2);
			p.node_ops[D] = vcg::Point2d(-M_SQRT1_2, -M_SQRT1_2);
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
};

#endif // OCTACELLPATTERN_H
