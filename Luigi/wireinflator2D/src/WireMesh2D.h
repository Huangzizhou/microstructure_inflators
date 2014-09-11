#ifndef WIREMESH2D_H
#define WIREMESH2D_H

#include "EdgeMeshType.h"
#include "EdgeMeshUtils.h"
#include "InflatorParameters.h"

#include <string>
#include <cmath>
#include <unordered_map>

#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>

class WireMesh2D
{
public:
	WireMesh2D(void) {;}


	WireMesh2D(EMesh & em)
	{
		this->setMesh(em);
	}

	WireMesh2D(const std::string & edgeMeshPath)
	{
		this->setMesh(edgeMeshPath);
	}

	WireMesh2D(WireMesh2D & wm)
	    : m_operations(wm.m_operations)
	{
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em, wm.m_em, false, true);
	}

	void setMesh(const std::string & edgeMeshPath)
	{
		m_operations.clear();
		bool done = EdgeMeshUtils<EMesh>::importObj(m_em, edgeMeshPath);

		if (!done)
		{
			m_em.Clear();
			return;
		}

		this->setup();
	}

	void setMesh(EMesh & em)
	{
		m_operations.clear();
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em, em, false, true);
		this->setup();
	}

	bool isValid(void) const
	{
		return (m_em.VN() > 0) && (m_em.EN() > 0);
	}

	size_t numberOfParameters(void) const
	{
		return m_operations.size();
	}

	CellParameters createCellParameters(void) const
	{
		CellParameters params(this->numberOfParameters());
		for (size_t i=0; i<params.numberOfParameters(); ++i)
		{
			std::pair<double,double> range = this->getParameterRange(i);
			params.parameter(i) = (range.first + range.second) / 2;
		}

		return params;
	}

	/*!
	 * \brief parametersValid checks wheter the current set of parameters is valid for the current wire mesh.
	 * \return true if the parameters are valid, false otherwise.
	 */
	bool parametersValid(const CellParameters & params) const
	{
		if (params.numberOfParameters() != m_operations.size())
			return false;

		std::pair<double,double> range_radius = this->getRadiusParameterRange();
		std::pair<double,double> range_transl = this->getTranslationParameterRange();
		for (size_t i=0; i<m_operations.size(); ++i)
		{
			switch (m_operations.at(i).type)
			{
			case ParameterOperation::Radius:
				if (params.cParameter(i) < range_radius.first ||
				    params.cParameter(i) > range_radius.second)
					return false;
				break;
			case ParameterOperation::Translation:
				if (params.cParameter(i) < range_transl.first ||
				    params.cParameter(i) > range_transl.second)
					return false;
				break;
			default:
				assert(0);
			}
		}
		return true;
	}

	std::pair<double,double> getParameterRange(int i) const
	{
		assert(i>=0 && i<int(m_operations.size()));

		switch (m_operations.at(i).type)
		{
		case ParameterOperation::Radius:
			return this->getRadiusParameterRange();
		case ParameterOperation::Translation:
			return this->getTranslationParameterRange();
		default:
			assert(0);
		}
		return std::make_pair(0,0);
	}

	void getUnmodifiedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_em, false, true);
	}

	void getEdgeMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	const std::vector<ParameterOperation> & getParameterOperations(void) const
	{
		return m_operations;
	}

protected:
	typedef EMesh::VertexType      VertexType;
	typedef VertexType::ScalarType ScalarType;
	typedef VertexType::CoordType  CoordType;

	EMesh m_em;
	std::vector<ParameterOperation> m_operations;

	static inline ScalarType Epsilon(void) { return 0.00001; }

	void setup(void)
	{
		if (!this->isValid())
		{
			m_em.Clear();
			return;
		}

		this->normalizeMesh();

		vcg::tri::UpdateTopology<EMesh>::VertexEdge(m_em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(m_em);

		std::vector<std::vector<int> > symmetryOrbits;
		this->computeSymmetryOrbits(symmetryOrbits);
		this->generateParameterOperations(symmetryOrbits);
	}

	virtual std::pair<double,double> getRadiusParameterRange(void) const
	{
		return std::make_pair( 0.01, 0.1 );
	}

	virtual std::pair<double,double> getTranslationParameterRange(void) const
	{
		return std::make_pair(-0.15, 0.15);
	}

	// force the mesh to be 1x1 squared
	void normalizeMesh(void)
	{
		if (!this->isValid())
			return;

		vcg::tri::UpdateBounding<EMesh>::Box(m_em);

		auto & bbox = m_em.bbox;
		vcg::tri::UpdatePosition<EMesh>::Translate(m_em, bbox.min);

		vcg::tri::UpdateBounding<EMesh>::Box(m_em);
		CoordType max = bbox.max;
		if (max[0] == 0) max[0] = 1;
		if (max[1] == 0) max[1] = 1;
		for (size_t i=0; i<m_em.vert.size(); ++i)
		{
			CoordType & p = m_em.vert[i].P();
			p[0] /= max[0];
			p[1] /= max[1];
			p[2] = 0;
		}

		vcg::tri::UpdateBounding<EMesh>::Box(m_em);

		assert((bbox.min[0] == 0 || bbox.min[1] == 0) &&
		       (bbox.max[0] == 1 || bbox.max[1] == 1));
	}

	// compute symmetry orbits clusters
	void computeSymmetryOrbits(std::vector<std::vector<int> > & symmetryOrbits)
	{
		// assumes the mesh is normalized
		if (!isValid())
			return;

		static const std::string AttrName = "tag";

		// tag each vertex
		auto tag = vcg::tri::Allocator<EMesh>::AddPerVertexAttribute<int>(m_em, AttrName);
		for (int i=0; i<int(m_em.vert.size()); ++i)
		{
			tag[i] = i;
		}

		// symmetrize wrt to x and y and find matching vertices
		std::unordered_map<int, std::vector<int> > orbits;
		for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
		{
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				CoordType flipped = m_em.vert[i].cP();
				flipped[ax] = (1 - flipped[ax]);

				for (int j=0; j<int(m_em.vert.size()); ++j)
				{
					if ( (i == j) || (tag[i] == tag[j]) )
						continue;

					CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

					if (delta[0] <= Epsilon() &&
					    delta[1] <= Epsilon() &&
					    delta[2] <= Epsilon())
					{
						tag[j] = tag[i];
					}
				}
			}
		}
		for (int i=0; i<int(m_em.vert.size()); ++i)
		{
			orbits[tag[i]].push_back(i);
		}

		symmetryOrbits.clear();
		for (auto & o : orbits)
		{
			symmetryOrbits.push_back(o.second);
		}

		vcg::tri::Allocator<EMesh>::DeletePerVertexAttribute(m_em, AttrName);
	}

	void generateParameterOperations(const std::vector<std::vector<int> > & symmetryOrbits)
	{
		m_operations.clear();

		// generate radius change for each symmetry orbit
		{
			ParameterOperation op;
			op.type = ParameterOperation::Radius;
			for (auto & s : symmetryOrbits)
			{
				op.nodes = s;
				m_operations.push_back(op);
			}
		}

		// generate displacement change
		{
			ParameterOperation op;
			op.type = ParameterOperation::Translation;
			const auto & meshBbox = m_em.bbox;
			for (auto & s : symmetryOrbits)
			{
				vcg::Box3<ScalarType> bbox;
				for (int i : s)
				{
					bbox.Add(m_em.vert[i].P());
				}

				// add x and y displacement operations
				for (int ax=0; ax<2; ++ax)
				{
					if (bbox.Dim()[ax] < meshBbox.Dim()[ax] &&
					    bbox.Dim()[ax] > Epsilon())
					{
						vcg::Point2d displ(0,0);
						for (int i : s)
						{
							displ[ax] = double(this->sign(m_em.vert[i].P()[ax] - meshBbox.Center()[ax]));
							op.nodes_displ[i] = displ;
						}
						m_operations.push_back(op);
						op.nodes_displ.clear();
					}
				}
			}
		}
	}

	void applyParameterOperations(EMesh & em, const CellParameters & params)
	{
		assert(m_operations.size() == params.numberOfParameters());

		for (size_t i=0; i<m_operations.size(); ++i)
		{
			const ParameterOperation & p = m_operations[i];
			switch (p.type)
			{
			// store vertex radius into quality
			case ParameterOperation::Radius:
				for (int node : p.nodes)
				{
					em.vert[node].Q() = params.cParameter(i);
				}
				break;
			// displace the vertex
			case ParameterOperation::Translation:
				for (auto & displ : p.nodes_displ)
				{
					em.vert[displ.first].P() +=
					        CoordType(displ.second[0], displ.second[1], 0) * params.cParameter(i);
				}
				break;
			default:
				assert(0);
			}
		}
	}

	template <typename T>
	static int sign(T x)
	{
		return (T(0) < x) - (x < T(0));
	}
};

#endif // WIREMESH2D_H
