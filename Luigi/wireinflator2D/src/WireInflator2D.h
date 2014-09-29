#ifndef WIREINFLATOR2D_H
#define WIREINFLATOR2D_H

#include "EdgeMeshPattern.h"
#include "EdgeMeshType.h"
#include "PolyMeshType.h"
#include "PolyMeshUtils.h"
#include "TriMeshType.h"
#include "mshLoader.h"
#include "WireMesh2D.h"
#include "OutMesh.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <string>
#include <strings.h>

class WireInflator2D
{
public:
	typedef EdgeMeshPattern<TMesh, EMesh> PatternGen;
	typedef OutMesh<2, 3>                 OutMeshType;

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

	/*!
	 * \brief generateQuadsPattern generates a wiremesh pattern over a quad mesh.
	 * \param quadMeshPath file path string (can be either an OBJ or MSH file).
	 * \param quadParameters the vector of cell parameters, one per quad.
	 * \param inT the tessellation parameters
	 * \param out the output mesh structure
	 */
	void generateQuadsPattern(const std::string & quadMeshPath,
	                          const std::vector<CellParameters> & quadParameters,
	                          const TessellationParameters & inT,
	                          OutMeshType & out)
	{
		typedef PolyMeshUtils<PolyMesh> PMU;
		PolyMesh pmesh;
		bool ok = false;
		if (checkFileExt(quadMeshPath, ".OBJ"))
		{
			ok = PMU::importFromOBJ(quadMeshPath, pmesh);
			if (ok)
				generateQuadsPattern(pmesh, quadParameters, inT, out);
		}
		else if (checkFileExt(quadMeshPath, ".MSH"))
		{
			VectorF nodes;
			VectorI elements;

			ok = loadQuadMsh(quadMeshPath, nodes, elements);
			if (ok)
				generateQuadsPattern(nodes, elements, quadParameters, inT, out);
		}

		if (!ok)
		{
			std::cout << "Unable to process."<< std::endl
			          << "Input mesh is invalid." << std::endl;
		}
	}

	/*!
	 * \brief generateQuadsPattern generates a wiremesh pattern over a quad mesh.
	 * \param nodes vector of nodes (coordinates)
	 * \param elements vector of elements (node indices)
	 * \param quadParameters the vector of cell parameters, one per quad.
	 * \param inT the tessellation parameters
	 * \param out the output mesh structure
	 */
	void generateQuadsPattern(const VectorF & nodes,
	                          const VectorI & elements,
	                          const std::vector<CellParameters> & quadParameters,
	                          const TessellationParameters & inT,
	                          OutMeshType & out)
	{
		PolyMesh pmesh;
		bool ok = vectorsToPolyMesh(nodes, elements, pmesh);
		if (!ok)
		{
			std::cout << "Unable to process."<< std::endl
			          << "Input mesh is invalid." << std::endl;
			return;
		}

		generateQuadsPattern(pmesh, quadParameters, inT, out);
	}

	void generateQuadsPattern(PolyMesh & pmesh,
	                          const std::vector<CellParameters> & quadParameters,
	                          const TessellationParameters & inT,
	                          OutMeshType & out)
	{
		typedef PolyMeshUtils<PolyMesh> PMU;

		if (!PMU::isQuadMesh(pmesh))
		{
			std::cout << "Unable to process."<< std::endl
			          << "Input mesh is not a quad mesh." << std::endl;
			return;
		}

		if ((pmesh.FN() != int(quadParameters.size())))
		{
			std::cout << "Unable to process."<< std::endl
			          << "Cells and parameters number mismatch." << std::endl;
			return;
		}

		for (size_t i=0; i<quadParameters.size(); ++i)
		{
			if (quadParameters[i].numberOfParameters() != m_pattern.numberOfParameters())
			{
				std::cout << "Unable to process."<< std::endl
				          << "Some quad has an invalid number of parameters." << std::endl;
				return;
			}
		}

		// copy parameters into per face attributes
		auto faceParams = vcg::tri::Allocator<PolyMesh>::GetPerFaceAttribute<CellParameters>(pmesh, m_pattern.CellParametersAttributeName());
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			faceParams[i] = quadParameters[i];
		}

		// generate the triangulated mesh
		TMesh m;

		m_pattern.params().tessellationParams = inT;

		m_pattern.generateFromQuads(pmesh);
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

	static bool checkFileExt(const std::string & filePath, const std::string & ext)
	{
		if (filePath.length() < ext.length())
			return false;

		int a = (filePath.length() - ext.length());
		const char * cfile = filePath.c_str() + a;
		const char * cext  = ext.c_str();
		return (strncasecmp(cfile , cext, ext.length()) == 0);
	}

	static bool loadQuadMsh(const std::string & file_name, VectorF & nodes, VectorI & elements)
	{
		try
		{
			MshLoader loader(file_name);
			if (loader.get_nodes_per_element() != 4)
				return false;

			nodes    = loader.get_nodes();
			elements = loader.get_elements();

		} catch (MshLoader::ErrorCode e) { (void)e; return false; }

		return true;
	}

	static bool vectorsToPolyMesh(const VectorF & nodes, const VectorI & elements, PolyMesh & pmesh)
	{
		pmesh.Clear();

		static const int NodeDim   = 3;
		static const int FaceArity = 4;

		if ((nodes.size()    % NodeDim   != 0) ||
		    (elements.size() % FaceArity != 0))
			return false;

		// fill nodes
		vcg::tri::Allocator<PolyMesh>::AddVertices(pmesh, nodes.size() / NodeDim);
		for (size_t i=0; i<pmesh.vert.size(); ++i)
		{
			PolyMesh::CoordType & p = pmesh.vert[i].P();
			for (int k=0; k<NodeDim; k++)
			{
				p[k] = nodes[i*NodeDim + k];
			}
		}

		// fill faces
		vcg::tri::Allocator<PolyMesh>::AddFaces(pmesh, elements.size() / FaceArity);
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			PolyMesh::FaceType & f = pmesh.face[i];
			f.Alloc(FaceArity);
			for (int k=0; k<FaceArity; k++)
			{
				int nodeIdx = elements[i*FaceArity + k];
				assert(nodeIdx >= 0 && nodeIdx < pmesh.VN());
				f.V(k) = &pmesh.vert[nodeIdx];
			}
		}

		vcg::tri::UpdateBounding<PolyMesh>::Box(pmesh);

		return true;
	}
};

#endif // WIREINFLATOR2D_H
