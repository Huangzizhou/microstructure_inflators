#ifndef WIREINFLATOR2D_H
#define WIREINFLATOR2D_H

#include "EdgeMeshPattern.h"
#include "WireMesh2D.h"
#include "TriMeshType.h"
#include "OutMesh.h"

#include <vcg/complex/complex.h>
#include <string>

class WireInflator2D
{
public:
	typedef EdgeMeshPattern<TMesh,WireMesh2D> PatternGen;
	typedef OutMesh<2, 3>                     OutMeshType;

	WireInflator2D(const std::string & edgeMeshPath)
	    : m_pattern(edgeMeshPath)
	{
		;
	}

	void generatePattern(const CellParameters & inP,
	                     const TessellationParameters & inT,
	                     OutMeshType & out,
	                     bool genVelocityField = true)
	{
		TMesh m;

		m_pattern.params() = { inP, inT };

		m_pattern.generate();
		m_pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);

		// compute edge velocity field
		if (genVelocityField)
			m_pattern.computeVelocityField(m, out.edge_fields);
	}

	void generateTiledPattern(const Array2D<CellParameters *> & inP, const TessellationParameters & inT, OutMeshType & out)
	{
		TMesh m;

		m_pattern.params().tessellationParams = inT;

		m_pattern.tile(inP);
		m_pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);
	}

	CellParameters createParameters(void)
	{
		return CellParameters(m_pattern.numberOfParameters());
	}

	const PatternGen & patternGenerator(void)
	{
		return m_pattern;
	}

private:
	PatternGen m_pattern;

	// convert vcg mesh to outmesh structure
	static void vcgMeshToOutMesh(OutMeshType & out, TMesh & m)
	{
		typedef typename OutMeshType::IndexType Index;

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
