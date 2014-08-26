#ifndef WIREINFLATOR2D_H
#define WIREINFLATOR2D_H

#include "InflatorParameters.h"
#include "OctaCellPattern.h"
#include "TriMeshType.h"
#include "OutMesh.h"
#include <vcg/complex/complex.h>

class WireInflator2D
{
public:
	typedef OctaCellPattern<TMesh>                               PatternGen;
	typedef typename PatternGen::PatternParameters               PatternParameters;
	typedef OutMesh<PatternParameters::NumberOfParameters, 2, 3> OutMeshType;

	static void generatePattern(const PatternParameters & inP, const TessellationParameters & inT, OutMeshType & out)
	{
		TMesh m;

		PatternGen pattern;
		pattern.params() = { inP, inT };

		pattern.generate();
		pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);

		// compute edge velocity field
		pattern.computeVelocityField(m, out.edge_fields);
	}

	static void generateTiledPattern(const Array2D<PatternParameters *> & inP, const TessellationParameters & inT, OutMeshType & out)
	{
		TMesh m;

		PatternGen pattern;
		PatternParameters p;
		pattern.params() = { p, inT };

		pattern.tile(inP);
		pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);
	}

private:
	// convert vcg mesh to outmesh structure
	static void vcgMeshToOutMesh(OutMeshType & out, TMesh & m)
	{
		typedef typename OutMeshType::IndexType  Index;

		// resize vectors
		out.nodes.resize(m.vert.size());
		out.elements.resize(m.face.size());
		out.edge_fields.clear();


		// fill in nodes
		for (size_t i=0; i<m.vert.size(); ++i)
		{
			// node position
			const TMesh::VertexType & v = m.vert[i];
			out.nodes[i] = { {double(v.cP()[0]), double(v.cP()[1])} };
		}

		// fill in elements
		for (size_t i=0; i<m.face.size(); ++i)
		{
			const TMesh::FaceType & f = m.face[i];
			out.elements[i] = {
			    { Index(vcg::tri::Index(m, f.cV(0))),
			      Index(vcg::tri::Index(m, f.cV(1))),
			      Index(vcg::tri::Index(m, f.cV(2))) }
			};
		}
	}
};

#endif // WIREINFLATOR2D_H
