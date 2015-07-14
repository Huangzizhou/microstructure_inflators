#ifndef WIREMESH2D_H
#define WIREMESH2D_H

#include "EdgeMeshUtils.h"
#include "InflatorParameters.h"

#include <string>
#include <cmath>
#include <unordered_map>
#include <cassert>

#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>

template <class EMesh>
class WireMesh2D
{
public:
	typedef typename EMesh::VertexType      VertexType;
	typedef typename EMesh::EdgeType        EdgeType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType  CoordType;

	WireMesh2D(void) {;}

	WireMesh2D(EMesh & em)
	{
		this->setMesh(em);
	}

	WireMesh2D(const std::string & edgeMeshPath)
	{
		this->setMesh(edgeMeshPath);
	}

	WireMesh2D & operator = (WireMesh2D & wm)
	{
		m_operations      = wm.m_operations;
		m_symmetry_orbits = wm.m_symmetry_orbits;
		m_orbits_idx      = wm.m_orbits_idx;
		m_radius_range    = wm.m_radius_range;
		m_transl_range    = wm.m_transl_range;
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;

		return *this;
	}

	WireMesh2D(const WireMesh2D & wm)
	    : m_operations(wm.m_operations)
	    , m_symmetry_orbits(wm.m_symmetry_orbits)
	    , m_orbits_idx(wm.m_orbits_idx)
	    , m_radius_range(wm.m_radius_range)
	    , m_transl_range(wm.m_transl_range)
	{
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;
	}

	void setMesh(const std::string & edgeMeshPath)
	{
		m_operations.clear();
		bool done = EdgeMeshUtils<EMesh>::importObj(m_em, edgeMeshPath);

		if (!done)
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
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

	void getUnmodifiedNormalizedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_normalized_em, false, true);
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

	void getNormalizedMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedNormalizedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	vcg::Point2<ScalarType> getOriginalScale(void) const
	{
		return vcg::Point2<ScalarType>(m_bbox.DimX(), m_bbox.DimY());
	}

	const std::vector<ParameterOperation> & getParameterOperations(void) const
	{
		// all symmetry orbit radius parameters come first
		return m_operations;
	}

	int getOrbitIndexForNode(int node)
	{
		if (!this->isValid() || node < 0 || node >= m_em.VN())
			return -1;

		auto tag = vcg::tri::Allocator<EMesh>::template FindPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName());
		if (!vcg::tri::Allocator<EMesh>::IsValidHandle(m_em, tag))
			return -1;

		return m_orbits_idx[tag[node]];
	}

protected:
	EMesh                           m_em;
	EMesh                           m_normalized_em;
	vcg::Box3<ScalarType>           m_bbox; // the original bounding box
	std::vector<ParameterOperation> m_operations;
	std::vector<std::vector<int> >  m_symmetry_orbits;
	std::unordered_map<int,int>     m_orbits_idx;
	std::pair<double,double>        m_radius_range;
	std::pair<double,double>        m_transl_range;

	static inline ScalarType Epsilon(void) { return 0.00001; }

	void setup(void)
	{
		this->m_radius_range = std::make_pair( 0.01, 5.0 );
		this->m_transl_range = std::make_pair(-0.18, 0.18);

		if (!this->isValid())
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
			return;
		}

		vcg::tri::Allocator<EMesh>::CompactEveryVector(m_em);
		vcg::tri::UpdateTopology<EMesh>::VertexEdge(m_em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(m_em);

		this->normalizeMesh();

		this->computeSymmetryOrbits();
		this->generateParameterOperations(m_symmetry_orbits);
	}

	virtual std::pair<double,double> getRadiusParameterRange(void) const
	{
		return m_radius_range;
	}

	virtual std::pair<double,double> getTranslationParameterRange(void) const
	{
		return m_transl_range;
	}

	static std::string SymmetryOrbitAttributeName(void)
	{
		return "symmetry_orbit";
	}

	void normalizeMesh(void)
	{
		if (!this->isValid())
			return;

		// store bounding box for later use
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);
		m_bbox = m_em.bbox;

		vcg::tri::UpdatePosition<EMesh>::Translate(m_em, -m_em.bbox.min);
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);

		assert(m_em.bbox.Dim()[0] != 0 && m_em.bbox.Dim()[1] != 0);

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, m_em, false, true);

		// real normalization
		for (size_t i=0; i<m_normalized_em.vert.size(); ++i)
		{
			CoordType & p = m_normalized_em.vert[i].P();
			p[0] /= m_normalized_em.bbox.max[0];
			p[1] /= m_normalized_em.bbox.max[1];
			p[2] = 0;
		}
		vcg::tri::UpdateBounding<EMesh>::Box(m_normalized_em);
	}

	// compute symmetry orbits clusters
	void computeSymmetryOrbits()
	{
		if (!isValid())
			return;

		// tag each vertex
		auto tag = vcg::tri::Allocator<EMesh>::template AddPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName());
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
				flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

				for (int j=0; j<int(m_em.vert.size()); ++j)
				{
					if ( (i == j) || (tag[i] == tag[j]) )
						continue;

					CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

					if (delta[0] <= Epsilon() &&
					    delta[1] <= Epsilon() &&
					    delta[2] <= Epsilon())
					{
						tag[j] = tag[i] = std::min(tag[i], tag[j]);
					}
				}
			}
		}
		for (int i=0; i<int(m_em.vert.size()); ++i)
		{
			orbits[tag[i]].push_back(i);
		}

		m_symmetry_orbits.clear();
		m_orbits_idx.clear();
		int index = 0;
		for (auto & o : orbits)
		{
			m_symmetry_orbits.push_back(o.second);
			m_orbits_idx[o.first] = index++;
		}
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
					        CoordType(displ.second[0] * em.bbox.Dim()[0], displ.second[1] * em.bbox.Dim()[1], 0) * params.cParameter(i);
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

// Morteza: use this to override symmetry orbits.
template <class EMesh>
class WireMesh2DMorteza : public WireMesh2D<EMesh> 
{

public:
	typedef typename EMesh::VertexType      VertexType;
	typedef typename EMesh::EdgeType        EdgeType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType  CoordType;
	
	WireMesh2DMorteza(void) {;}

	WireMesh2DMorteza(EMesh & em)
	{
		this->setMesh(em);
	}

	WireMesh2DMorteza(const std::string & edgeMeshPath)
	{
		this->setMesh(edgeMeshPath);
	}

	// MHS on JUL14, 2015:
	// A new constructor with inSymmetryMode as the input ... 
	WireMesh2DMorteza(const std::string & edgeMeshPath, const int inSymmetryMode)
	{
		this->symmetryMode = inSymmetryMode;
		this->setMesh(edgeMeshPath);
	}


	WireMesh2DMorteza & operator = (WireMesh2DMorteza & wm)
	{
		m_operations      = wm.m_operations;
		m_symmetry_orbits = wm.m_symmetry_orbits;
		m_orbits_idx      = wm.m_orbits_idx;
		m_radius_range    = wm.m_radius_range;
		m_transl_range    = wm.m_transl_range;
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;

		return *this;
	}

	WireMesh2DMorteza(const WireMesh2DMorteza & wm)
	    : m_operations(wm.m_operations)
	    , m_symmetry_orbits(wm.m_symmetry_orbits)
	    , m_orbits_idx(wm.m_orbits_idx)
	    , m_radius_range(wm.m_radius_range)
	    , m_transl_range(wm.m_transl_range)
	{
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;
	}

	void setMesh(const std::string & edgeMeshPath)
	{
		m_operations.clear();
		bool done = EdgeMeshUtils<EMesh>::importObj(m_em, edgeMeshPath);

		if (!done)
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
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

	void getUnmodifiedNormalizedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_normalized_em, false, true);
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

	void getNormalizedMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedNormalizedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	vcg::Point2<ScalarType> getOriginalScale(void) const
	{
		return vcg::Point2<ScalarType>(m_bbox.DimX(), m_bbox.DimY());
	}

	const std::vector<ParameterOperation> & getParameterOperations(void) const
	{
		// all symmetry orbit radius parameters come first
		return m_operations;
	}

	int getOrbitIndexForNode(int node)
	{
		if (!this->isValid() || node < 0 || node >= m_em.VN())
			return -1;

		auto tag = vcg::tri::Allocator<EMesh>::template FindPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName());
		if (!vcg::tri::Allocator<EMesh>::IsValidHandle(m_em, tag))
			return -1;

		return m_orbits_idx[tag[node]];
	}

protected:
	EMesh                           m_em;
	EMesh                           m_normalized_em;
	vcg::Box3<ScalarType>           m_bbox; // the original bounding box
	std::vector<ParameterOperation> m_operations;
	std::vector<std::vector<int> >  m_symmetry_orbits;
	std::unordered_map<int,int>     m_orbits_idx;
	std::pair<double,double>        m_radius_range;
	std::pair<double,double>        m_transl_range;
	int symmetryMode;

	static inline ScalarType Epsilon(void) { return 0.00001; }
	static inline ScalarType OneForth(void) { return 0.25; }
	static inline ScalarType ThreeForth(void) { return 0.75; }

	void setup(void)
	{
		this->m_radius_range = std::make_pair( 0.01, 5.0 );
		this->m_transl_range = std::make_pair(-0.18, 0.18);

		if (!this->isValid())
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
			return;
		}

		vcg::tri::Allocator<EMesh>::CompactEveryVector(m_em);
		vcg::tri::UpdateTopology<EMesh>::VertexEdge(m_em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(m_em);

		this->normalizeMesh();

		this->computeSymmetryOrbits();
		this->generateParameterOperations(m_symmetry_orbits);
	}

	virtual std::pair<double,double> getRadiusParameterRange(void) const
	{
		return m_radius_range;
	}

	virtual std::pair<double,double> getTranslationParameterRange(void) const
	{
		return m_transl_range;
	}

	static std::string SymmetryOrbitAttributeName(void)
	{
		return "symmetry_orbit";
	}

	void normalizeMesh(void)
	{
		if (!this->isValid())
			return;

		// store bounding box for later use
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);
		m_bbox = m_em.bbox;

		vcg::tri::UpdatePosition<EMesh>::Translate(m_em, -m_em.bbox.min);
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);

		assert(m_em.bbox.Dim()[0] != 0 && m_em.bbox.Dim()[1] != 0);

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, m_em, false, true);

		// real normalization
		for (size_t i=0; i<m_normalized_em.vert.size(); ++i)
		{
			CoordType & p = m_normalized_em.vert[i].P();
			p[0] /= m_normalized_em.bbox.max[0];
			p[1] /= m_normalized_em.bbox.max[1];
			p[2] = 0;
		}
		vcg::tri::UpdateBounding<EMesh>::Box(m_normalized_em);
	}

	// compute symmetry orbits clusters
	void computeSymmetryOrbits()
	{
		if (!isValid())
			return;

		// tag each vertex
		auto tag = vcg::tri::Allocator<EMesh>::template AddPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName());

		// MHS JUL14, 2015:
		// TODO: Clean up the cases and get rid of unnecessary ones
		// TODO: Explain what each case means ...
		std::unordered_map<int, std::vector<int> > orbits;
		int tempTag = 1; // this is needed for case 2
		switch (symmetryMode){
			case 0:
				// symmetrize wrt to x and y and find matching vertices
				for (int i=0; i<int(m_em.vert.size()); ++i)
				{
					tag[i] = i;
				}


				for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
				{
					for (int i=0; i<int(m_em.vert.size()); ++i)
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
							tag[j] = tag[i] = std::min(tag[i], tag[j]);
						}
					}
				}
			}
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				orbits[tag[i]].push_back(i);
			}
			break;
		case 1:
			// Modified version of Luigi's symmetry orbits
			// all boundary vertices --> one orbit,
			// all vertices connected to boundary nodes --> one orbit, and
			// the remaining vertices each have thier own orbit
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				tag[i] = i;
			}

			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				double xCoord = m_em.vert[i].cP()[0];
				double yCoord = m_em.vert[i].cP()[1];
				bool midBoundaryNode = vcg::math::Abs(xCoord) <= Epsilon() ||
										vcg::math::Abs(yCoord) <= Epsilon() ||
										vcg::math::Abs(xCoord - 1.0) <= Epsilon() ||
										vcg::math::Abs(yCoord - 1.0) <= Epsilon();
				bool cornerInternalNode = vcg::math::Abs(xCoord - 0.5) > Epsilon() &&
											vcg::math::Abs(yCoord - 0.5) > Epsilon();

				bool midInternalNode = !midBoundaryNode && !cornerInternalNode;

				if (midInternalNode) 
				{
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}
				}
				else if (midBoundaryNode)
				{	
					// use this if you want the boundary verices behave the same as 
					// Luigi's code ... 
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}

					/* uncomment this if you want to have all boundary vertices be the same
					for (int j=0; j<int(m_em.vert.size()); ++j)
					{
						if ( (i == j) || (tag[i] == tag[j]) )
							continue;
						if (vcg::math::Abs(m_em.vert[j].cP()[0]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[0] - 1.0) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1] - 1.0) <= Epsilon())
						{
							tag[j] = tag[i] = std::min(tag[i], tag[j]);
						}
					}
					*/
				}
			}
			
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				orbits[tag[i]].push_back(i);
			}
			break;
		case 2:
			// Symmetry Orbits for distorded cells (this may result in non-periodic unit cells):
			// orbit 0 contains all the boundary vertices
			// orbits 1-8 each contain only one of the internal vertices 
			for (int i = 0; i < int(m_em.vert.size()); ++i)
			{
				if (m_em.vert[i].cP()[0] < OneForth() ||
					m_em.vert[i].cP()[0] > ThreeForth() ||
					m_em.vert[i].cP()[1] < OneForth() ||
					m_em.vert[i].cP()[1] > ThreeForth())
				{
					tag[i] = 0;
				}
				else
				{
					tag[i] = tempTag;
					tempTag++;
				}
			}
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				orbits[tag[i]].push_back(i);
			}
			break;
		case 3:
			// like case 1 except that the diagonal vertices are linked together
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				tag[i] = i;
			}

			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				double xCoord = m_em.vert[i].cP()[0];
				double yCoord = m_em.vert[i].cP()[1];
				bool midBoundaryNode = vcg::math::Abs(xCoord) <= Epsilon() ||
										vcg::math::Abs(yCoord) <= Epsilon() ||
										vcg::math::Abs(xCoord - 1.0) <= Epsilon() ||
										vcg::math::Abs(yCoord - 1.0) <= Epsilon();
				bool cornerInternalNode = vcg::math::Abs(xCoord - 0.5) > Epsilon() &&
											vcg::math::Abs(yCoord - 0.5) > Epsilon();

				bool midInternalNode = !midBoundaryNode && !cornerInternalNode;

				if (midInternalNode) 
				{
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}
				}
				else if (midBoundaryNode)
				{	
					// use this if you want the boundary verices behave the same as 
					// Luigi's code ... 
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}

					/* uncomment this if you want to have all boundary vertices be the same
					for (int j=0; j<int(m_em.vert.size()); ++j)
					{
						if ( (i == j) || (tag[i] == tag[j]) )
							continue;
						if (vcg::math::Abs(m_em.vert[j].cP()[0]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[0] - 1.0) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1] - 1.0) <= Epsilon())
						{
							tag[j] = tag[i] = std::min(tag[i], tag[j]);
						}
					}
					*/
				}
				else if (cornerInternalNode)
				{
				//	for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
				//	{
						CoordType flipped = m_em.vert[i].cP();
						flipped[0] = (m_em.bbox.Dim()[0] - flipped[0]);
						flipped[1] = (m_em.bbox.Dim()[1] - flipped[1]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
				//	}

				}
			}
			
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				orbits[tag[i]].push_back(i);
			}
			break;
		case 4:
			// like case 3 except that the diagonal vertices are not linked together
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				tag[i] = i;
			}

			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				double xCoord = m_em.vert[i].cP()[0];
				double yCoord = m_em.vert[i].cP()[1];
				bool midBoundaryNode = vcg::math::Abs(xCoord) <= Epsilon() ||
										vcg::math::Abs(yCoord) <= Epsilon() ||
										vcg::math::Abs(xCoord - 1.0) <= Epsilon() ||
										vcg::math::Abs(yCoord - 1.0) <= Epsilon();
				bool cornerInternalNode = vcg::math::Abs(xCoord - 0.5) > Epsilon() &&
											vcg::math::Abs(yCoord - 0.5) > Epsilon();

				bool midInternalNode = !midBoundaryNode && !cornerInternalNode;

				if (midInternalNode) 
				{
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}
				}
				else if (midBoundaryNode)
				{	
					// use this if you want the boundary verices behave the same as 
					// Luigi's code ... 
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}

					/* uncomment this if you want to have all boundary vertices be the same
					for (int j=0; j<int(m_em.vert.size()); ++j)
					{
						if ( (i == j) || (tag[i] == tag[j]) )
							continue;
						if (vcg::math::Abs(m_em.vert[j].cP()[0]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[0] - 1.0) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1] - 1.0) <= Epsilon())
						{
							tag[j] = tag[i] = std::min(tag[i], tag[j]);
						}
					}
					*/
				}
				else if (cornerInternalNode)
				{
					/*
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[0] = (m_em.bbox.Dim()[0] - flipped[0]);
						flipped[1] = (m_em.bbox.Dim()[1] - flipped[1]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}
					*/

				}
			}
			
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				orbits[tag[i]].push_back(i);
			}
			break;
		case 5:
			// like case 1 except that the diagonal vertices are linked together
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				tag[i] = i;
			}

			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				double xCoord = m_em.vert[i].cP()[0];
				double yCoord = m_em.vert[i].cP()[1];
				bool midBoundaryNode = vcg::math::Abs(xCoord) <= Epsilon() ||
										vcg::math::Abs(yCoord) <= Epsilon() ||
										vcg::math::Abs(xCoord - 1.0) <= Epsilon() ||
										vcg::math::Abs(yCoord - 1.0) <= Epsilon();
				bool cornerInternalNode = vcg::math::Abs(xCoord - 0.5) > Epsilon() &&
											vcg::math::Abs(yCoord - 0.5) > Epsilon();

				bool midInternalNode = !midBoundaryNode && !cornerInternalNode;

				if (midInternalNode) 
				{
					for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
					{
						CoordType flipped = m_em.vert[i].cP();
						flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

						for (int j=0; j<int(m_em.vert.size()); ++j)
						{
							if ( (i == j) || (tag[i] == tag[j]) )
								continue;

							CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

							if (delta[0] <= Epsilon() &&
								delta[1] <= Epsilon() &&
								delta[2] <= Epsilon())
							{
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
							}
						}
					}
				}
				else if (midBoundaryNode)
				{	
					for (int j=0; j<int(m_em.vert.size()); ++j)
					{
						if ( (i == j) || (tag[i] == tag[j]) )
							continue;

						if (vcg::math::Abs(m_em.vert[j].cP()[0]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1]) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[0] - 1.0) <= Epsilon() ||
							vcg::math::Abs(m_em.vert[j].cP()[1] - 1.0) <= Epsilon())
								tag[j] = tag[i] = std::min(tag[i], tag[j]);
					}
				}
				else if (cornerInternalNode)
				{
				}
			}
			
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				orbits[tag[i]].push_back(i);
			}
			break;

		} // end of switch(seymmetryMode)

		m_symmetry_orbits.clear();
		m_orbits_idx.clear();
		int index = 0;
		for (auto & o : orbits)
		{
			m_symmetry_orbits.push_back(o.second);
			m_orbits_idx[o.first] = index++;
		}


		// MHS JUL14, 2015: 
		// printout the orbits
		std::cout << "----------" << std::endl << "printing out the orbits" << std::endl << "----------" << std::endl;
		std::cout << std::endl << "symmetryMode is: " <<symmetryMode << std::endl << std::endl;
		for (auto & s : m_symmetry_orbits)
		{

			std::cout<< "size of the orbit is: " << s.size() << " and the points are: " << std::endl;
			for (int i :s)
			{

				std::cout<< "p" << i << ": " << "(" <<  m_em.vert[i].cP()[0] << "," <<m_em.vert[i].cP()[1] << ")" << "     ";
			}
			std::cout<<endl;
		}
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


			// MHS JUL14, 2015:
			// TODO: Clean up the cases and get rid of unnecessary ones
			// TODO: Explain what each case means ...
			switch(symmetryMode){
				case 0:
					// This is the original Luigi's code
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
					break;
				// this is supposed to fixed the issue with case 2 below
				case 1:
					for (auto & s : symmetryOrbits)
					{
						if (s.size() == 1){
							vcg::Point2d displx(0,0);
							vcg::Point2d disply(0,0);
							displx[0] = double(1.0);
							disply[1] = double(1.0);
							for (int i : s){
								// add x displs
								op.nodes_displ[i] = displx;
								m_operations.push_back(op);
								op.nodes_displ.clear();
								// add y displs
								op.nodes_displ[i] = disply;
								m_operations.push_back(op);
								op.nodes_displ.clear();
							}
						}
						else if (s.size() == 2)
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
					break;
				case 2:
					// This is the modified Morteza's code to account for nonsymmetric cells (this version results in non-periodic cells
					// in this version every x and y displacements and thicknesses of all the internal vertices can change independently (i.e., they are all pattern parameters)
					for (auto & s : symmetryOrbits)
					{
						if (s.size() == 1){
							vcg::Point2d displx(0,0);
							vcg::Point2d disply(0,0);
							displx[0] = double(1.0);
							disply[1] = double(1.0);
							for (int i : s){
								// add x displs
								op.nodes_displ[i] = displx;
								m_operations.push_back(op);
								op.nodes_displ.clear();
								// add y displs
								op.nodes_displ[i] = disply;
								m_operations.push_back(op);
								op.nodes_displ.clear();
							}
						}
					}
					break;
				case 3:
					for (auto & s : symmetryOrbits)
					{
						if (s.size() == 1){
							vcg::Point2d displx(0,0);
							vcg::Point2d disply(0,0);
							displx[0] = double(1.0);
							disply[1] = double(1.0);
							for (int i : s){
								// add x displs
								op.nodes_displ[i] = displx;
								m_operations.push_back(op);
								op.nodes_displ.clear();
								// add y displs
								op.nodes_displ[i] = disply;
								m_operations.push_back(op);
								op.nodes_displ.clear();
							}
						}
						else if (s.size() == 2)
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
					break;
				case 4:
					for (auto & s : symmetryOrbits)
					{
						if (s.size() == 1){
							vcg::Point2d displx(0,0);
							vcg::Point2d disply(0,0);
							displx[0] = double(1.0);
							disply[1] = double(1.0);
							for (int i : s){
								// add x displs
								op.nodes_displ[i] = displx;
								m_operations.push_back(op);
								op.nodes_displ.clear();
								// add y displs
								op.nodes_displ[i] = disply;
								m_operations.push_back(op);
								op.nodes_displ.clear();
							}
						}
						else if (s.size() == 2)
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
					break;
				case 5:
					for (auto & s : symmetryOrbits)
					{
						if (s.size() == 1){
							vcg::Point2d displx(0,0);
							vcg::Point2d disply(0,0);
							displx[0] = double(1.0);
							disply[1] = double(1.0);
							for (int i : s){
								// add x displs
								op.nodes_displ[i] = displx;
								m_operations.push_back(op);
								op.nodes_displ.clear();
								// add y displs
								op.nodes_displ[i] = disply;
								m_operations.push_back(op);
								op.nodes_displ.clear();
							}
						}
						else if (s.size() == 2)
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
					break;

			} // end of switch(symmetryMode)

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
					        CoordType(displ.second[0] * em.bbox.Dim()[0], displ.second[1] * em.bbox.Dim()[1], 0) * params.cParameter(i);
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
