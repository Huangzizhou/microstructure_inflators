#ifndef EDGEMESHPATTERN_H
#define EDGEMESHPATTERN_H

#include "Pattern2D.h"
#include "PolyMeshUtils.h"
#include "WireMesh2D.h"
#include "WireMeshEmbedding.h"
#include "InflatorParameters.h"

#include "tessellator2d.h"
#include <stdlib.h>
#include <cmath>
#include <utility>
#include <array>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/space/index/grid_static_ptr2d.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>

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

template <typename T>
class Array2D
{
public:
	Array2D(size_t w, size_t h)
	    : m_w(w)
	    , m_h(h)
	    , m_data(NULL)
	{
		if (w > 0 && h > 0)
			m_data = new T[w * h];
		else
			assert(0);
	}

	~Array2D(void)
	{
		delete[] m_data;
	}

	size_t width() const  { return m_w; }
	size_t height() const { return m_h; }

	const T & operator () (size_t x, size_t y) const
	{
		return m_data[ y*m_w + x ];
	}

	T & operator () (size_t x, size_t y)
	{
		return m_data[ y*m_w + x ];
	}

private:
	size_t m_w;
	size_t m_h;
	T *    m_data;

	Array2D(const Array2D<T> &);
	Array2D & operator = (const Array2D<T> &);
};


template <class TriMesh, class EMesh>
class EdgeMeshPattern : public Pattern2D<TriMesh>
{
public:
	typedef WireMesh2D<EMesh>                  WireMesh;
	typedef EdgeMeshPattern<TriMesh, EMesh>    ThisType;
	typedef Pattern2D<TriMesh>                 BaseType;
	typedef typename BaseType::ScalarType      ScalarType;
	typedef typename BaseType::Coord2Type      Coord2Type;

	typedef typename EMesh::CoordType           ECoordType;
	typedef typename EMesh::EdgeType            EEdgeType;
	typedef typename EMesh::VertexType          EVertexType;

	typedef Tessellator2DSettings              TessellatorSettings;

	EdgeMeshPattern(void)
	{
		;
	}

	EdgeMeshPattern(EdgeMeshPattern & emp)
	    : m_tri_settings(emp.m_tri_settings)
	    , m_wire(emp.m_wire)
	    , m_params(emp.m_params)
	    , BaseType::m_paths(emp.m_paths)
	{
		;
	}

	EdgeMeshPattern & operator = (EdgeMeshPattern & emp)
	{
		m_tri_settings = emp.m_tri_settings;
		m_wire = emp.m_wire;
		m_params = emp.m_params;
		this->m_paths = emp.m_paths;

		return *this;
	}

	EdgeMeshPattern(WireMesh & wm)
	    : m_wire(wm)
	{
		assert(m_wire.isValid());
		this->setup();
	}

	EdgeMeshPattern(EMesh & em)
	{
		m_wire.setMesh(em);
		assert(m_wire.isValid());
		this->setup();
	}

	EdgeMeshPattern(const std::string & edgeMeshPath)
	{
		m_wire.setMesh(edgeMeshPath);
		assert(m_wire.isValid());
		this->setup();
	}

	TessellatorSettings & getTessellationSettings(void)
	{
		return m_tri_settings;
	}

	const TessellatorSettings & getConstTessellationSettings(void) const
	{
		return m_tri_settings;
	}

	struct InputParameters
	{
		CellParameters         patternParams;
		TessellationParameters tessellationParams;
	};

	InputParameters & params(void)
	{
		return m_params;
	}

	size_t numberOfParameters(void) const
	{
		return m_wire.numberOfParameters();
	}

	bool parametersValid(const CellParameters & params) const
	{
		return m_wire.parametersValid(params);
	}

	std::pair<double,double> getParameterRange(int i) const
	{
		return m_wire.getParameterRange(i);
	}

	const std::vector<ParameterOperation> & getParameterOperations(void) const
	{
		return m_wire.getParameterOperations();
	}

	bool generate(void)
	{
		if (!m_wire.isValid())
			return false;

		this->setTessellationParameters();

		this->m_paths.clear();

		generateOneElement();

		// get bounds
		ClipperLib::IntRect bounds = getBounds(this->m_base_paths);
		this->matchPeriodicVertices(this->m_base_paths, bounds);

		// optimize path removing small edges
		simplifyContour(this->m_base_paths, bounds);

		ClipperLib::SimplifyPolygons(this->m_base_paths);
		splitAllCountourEdges(this->m_base_paths);

		// remove weird occurrencies of tiny boundary edges
		simplifyContour(this->m_base_paths, bounds, true);

		this->m_paths = this->m_base_paths;

		return true;
	}

	// works only when generating a single axis-aligned cell
	void computeVelocityField(TriMesh & mesh,
	                          std::unordered_map<std::pair<int, int>, std::vector<ScalarType> > & fields)
	{                                                              // contains one scalar for each pattern parameter
		typedef typename ParameterOperation::NodeID NodeID;

		if (!m_wire.isValid())
			return;

		// get the generating edge mesh
		EMesh em;
		m_wire.getEdgeMesh(em, this->m_params.patternParams);

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
						ECoordType n = ECoordType::Construct((v1->cP() - v0->cP()).Normalize());
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
				const EEdgeType & e = em.edge[i];
				if (e.IsD())
					continue;

				// get their vertices
				const EVertexType & v0 = *e.cV(0);
				const EVertexType & v1 = *e.cV(1);

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
				const EVertexType & v = em.vert[i];
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
		const std::vector<ParameterOperation> & params_op =
			m_wire.getParameterOperations();
		size_t numberOfParameters = params_op.size();

		fields.clear();
		typedef typename TriMesh::VertexType::NormalType NormalType;
		for (size_t p=0; p<numberOfParameters; ++p)
		{
			const ParameterOperation & par = params_op[p];

			// for each parameter compute the displacement per vertex
			if (par.type == ParameterOperation::Radius)
			{
				for (size_t i=0; i<vtxToSegs.size(); ++i)
				{
					typename TriMesh::VertexPointer v = vtx[i];
					const SegmentType * s             = vtxToSegs[i].first;
					const Coord2Type  & closest       = vtxToSegs[i].second;

					switch (s->belonging) {
					case SegmentType::EdgeSegment :
					{
						NormalType n0 = (std::find(par.nodes.begin(), par.nodes.end(), NodeID(s->index0)) == par.nodes.end()) ?
						                    NormalType(0,0,0) :
						                    NormalType::Construct(ECoordType(s->P0()[0], s->P0()[1], 0) - em.vert[s->index0].cP()).Normalize();

						NormalType n1 = (std::find(par.nodes.begin(), par.nodes.end(), NodeID(s->index1)) == par.nodes.end()) ?
						                    NormalType(0,0,0) :
						                    NormalType::Construct(ECoordType(s->P1()[0], s->P1()[1], 0) - em.vert[s->index1].cP()).Normalize();

						ScalarType t =  s->interpolationParameter(closest);
						v->N() = n0 * (1-t) + n1 * t;

						break;
					}
					case SegmentType::VertexSegment :
					{
						v->N() = (std::find(par.nodes.begin(), par.nodes.end(), NodeID(s->index0)) == par.nodes.end()) ?
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
			else if (par.type == ParameterOperation::Translation)
			{
				// for each parameter compute the displacement per vertex
				for (size_t i=0; i<vtxToSegs.size(); ++i)
				{
					typename TriMesh::VertexPointer v = vtx[i];
					const SegmentType * s             = vtxToSegs[i].first;
					const Coord2Type  & closest       = vtxToSegs[i].second;

					switch (s->belonging) {
					case SegmentType::EdgeSegment :
					{
						NormalType n0 = (par.nodes_displ.count(NodeID(s->index0)) == 0) ?
						                    NormalType(0,0,0) :
						                    NormalType(par.nodes_displ.at(NodeID(s->index0))[0], par.nodes_displ.at(NodeID(s->index0))[1], 0);

						NormalType n1 = (par.nodes_displ.count(NodeID(s->index1)) == 0) ?
						                    NormalType(0,0,0) :
						                    NormalType(par.nodes_displ.at(NodeID(s->index1))[0], par.nodes_displ.at(NodeID(s->index1))[1], 0);

						ScalarType t =  s->interpolationParameter(closest);
						v->N() = n0 * (1-t) + n1 * t;
						break;
					}
					case SegmentType::VertexSegment :
					{
						v->N() = (par.nodes_displ.count(NodeID(s->index0)) == 0) ?
						             NormalType(0,0,0) :
						             NormalType(par.nodes_displ.at(NodeID(s->index1))[0], par.nodes_displ.at(NodeID(s->index1))[1], 0);
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
				std::vector<ScalarType> & f = fields[e.first];
				if (f.size() <= 0)
				{
					f.resize(numberOfParameters);
				}
				f[p] = val;
			}
		}
	}

	bool tile(const Array2D<CellParameters *> & grid)
	{
		if (!m_wire.isValid())
			return false;

		this->setTessellationParameters();

		this->m_paths.clear();

		ClipperLib::Clipper clipper;

		for (unsigned int x=0; x<grid.width(); ++x)
		{
			for (unsigned int y=0; y<grid.height(); ++y)
			{
				if (grid(x,y) == NULL)
					continue;

				generateOneElement(*grid(x,y));
				ClipperLib::SimplifyPolygons(this->m_base_paths); // check

				// optimize path removing small edges
				ClipperLib::IntRect bounds = getBounds(this->m_base_paths);
				simplifyContour(this->m_base_paths, bounds);


				ClipperLib::IntPoint t((bounds.right - bounds.left) * x, (bounds.top - bounds.bottom) * y);
				clipper.AddPaths((this->m_base_paths + t), ClipperLib::ptSubject, true);
			}
		}

		clipper.Execute(ClipperLib::ctUnion, this->m_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		ClipperLib::SimplifyPolygons(this->m_paths);
		splitAllCountourEdges(this->m_paths);

		return true;
	}

	static std::string CellParametersAttributeName(void)
	{
		return "cell_parameters";
	}

	template <class PolyMesh>
	bool generateFromQuads(PolyMesh & pmesh)
	{
		typedef WireMeshEmbedding<EMesh, PolyMesh> WireEmbedding;

		if (!m_wire.isValid())
			return false;

		if (!PolyMeshUtils<PolyMesh>::isQuadMesh(pmesh))
			return false;

		this->setTessellationParameters();
		this->m_paths.clear();

		WireEmbedding::preprocessQuadMesh(pmesh);

		ClipperLib::Clipper clipper;
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			// generate one element per face (quad)
			typename PolyMesh::FaceType & f = pmesh.face[i];

			// retrieve cell parameters
			CellParameters params;
			auto faceParams = vcg::tri::Allocator<PolyMesh>::template FindPerFaceAttribute<CellParameters>(pmesh, CellParametersAttributeName());
			if (vcg::tri::Allocator<PolyMesh>::template IsValidHandle<CellParameters>(pmesh, faceParams))
			{
				params = faceParams[&f];
			}
			else
			{
				params = m_wire.createCellParameters();
			}

			this->generateOneElement(params, pmesh, f);

			// get the clipping quad
			ClipperLib::Path clip_poly;
			for (char i=0; i<f.VN(); i++)
			{
				clip_poly.push_back(ThisType::convertToIntPoint(Coord2Type(f.cP(i)[0], f.cP(i)[1])));
			}
			// optimize path removing small edges
			simplifyContour(this->m_base_paths, clip_poly);

			ClipperLib::SimplifyPolygons(this->m_base_paths);

			clipper.AddPaths(this->m_base_paths, ClipperLib::ptSubject, true);
		}

		clipper.Execute(ClipperLib::ctUnion, this->m_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		splitAllCountourEdges(this->m_paths);

		return true;
	}

	bool tessellate(TriMesh & mesh)
	{
		if (!m_wire.isValid() || this->m_paths.size() == 0)
			return false;

		bool ret = Tessellator2D<TriMesh>::execute(m_tri_settings, this->m_paths, ThisType::ScaleFactor, mesh);
		if (ret)
		{
			this->decimateMesh(mesh);
			vcg::tri::Allocator<TriMesh>::CompactEveryVector(mesh);
			vcg::tri::UpdateBounding<TriMesh>::Box(mesh);
		}
		return ret;
	}

protected:
	/// parameters
	TessellatorSettings m_tri_settings;
	WireMesh        m_wire;
	InputParameters m_params;

	void setup()
	{
		params().patternParams = m_wire.createCellParameters();
	}

	virtual void generateOneElement(void)
	{
		generateOneElement(this->m_params.patternParams);
	}

	template <class PolyMesh>
	void generateOneElement(const CellParameters & pars, PolyMesh & pmesh, const typename PolyMesh::FaceType & f)
	{
		typedef WireMeshEmbedding<EMesh, PolyMesh> WireEmbedding;

		if (!m_wire.isValid())
			return;

		this->m_base_paths.clear();

		// get the edge mesh embedded within the quad
		EMesh em;
		m_wire.getNormalizedMesh(em, pars);

		WireEmbedding::embedWireMesh(em, f, pmesh);

		vcg::tri::RequireVEAdjacency(em);
		vcg::tri::RequireEEAdjacency(em);

		ClipperLib::Paths profile;
		inflateEdgeMesh(em, profile);

		// get the clipping quad
		ClipperLib::Path clip_poly;
		for (char i=0; i<4; i++)
		{
			clip_poly.push_back(ThisType::convertToIntPoint(Coord2Type(f.cP(i)[0], f.cP(i)[1])));
		}

		// clip with the quad
		ClipperLib::Clipper clip;
		clip.AddPaths(profile, ClipperLib::ptSubject, true);
		clip.AddPath(clip_poly, ClipperLib::ptClip, true);
		clip.Execute(ClipperLib::ctIntersection, this->m_base_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

		// get the border segments adjacent to another quad
		clip_poly.push_back(clip_poly.back());
		ClipperLib::Paths tmp;
		for (int i=0; i<f.VN(); ++i)
		{
			if (!f.IsB(i))
			{
				ClipperLib::Path edge;
				edge.push_back(clip_poly[i]);
				edge.push_back(clip_poly[(i+1)%f.VN()]);
				tmp.push_back(edge);
			}
		}
		ClipperLib::PolyTree pt;
		clip.Clear();
		clip.AddPaths(tmp, ClipperLib::ptSubject, false);
		clip.AddPaths(profile, ClipperLib::ptClip, true);
		clip.Execute(ClipperLib::ctIntersection, pt, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

		ClipperLib::ClipperOffset offset;
		for (int i=0; i<pt.ChildCount(); i++)
		{
			this->m_base_paths.push_back(ClipperLib::Path());
			ClipperLib::PolyNode & pn = *pt.Childs.at(i);
			if (!pn.IsOpen())
			{
				assert(0);
				continue;
			}
			offset.AddPath(pn.Contour, ClipperLib::jtSquare, ClipperLib::etOpenButt);
		}
		offset.Execute(tmp, 0.0001 * ThisType::ScaleFactor);

		clip.Clear();
		clip.AddPaths(tmp, ClipperLib::ptSubject, true);
		clip.AddPaths(this->m_base_paths, ClipperLib::ptSubject, true);
		clip.Execute(ClipperLib::ctUnion, this->m_base_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	}

	void generateOneElement(const CellParameters & pars)
	{
		if (!m_wire.isValid())
			return;

		this->m_base_paths.clear();

		// get the edge mesh with displaced vertices and radius set
		EMesh em;
		m_wire.getEdgeMesh(em, pars);

		vcg::tri::RequireVEAdjacency(em);
		vcg::tri::RequireEEAdjacency(em);

		ClipperLib::Paths profile;
		inflateEdgeMesh(em, profile);

		// get the bbox
		ClipperLib::IntRect intBbox = getClipperBbox(em.bbox);
		ClipperLib::Path clip_bbox = convertToPath(intBbox);

		// clip to bounding box
		ClipperLib::Clipper clip;
		clip.AddPaths(profile, ClipperLib::ptSubject, true);
		clip.AddPath(clip_bbox, ClipperLib::ptClip, true);
		clip.Execute(ClipperLib::ctIntersection, this->m_base_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	}

	void inflateEdgeMesh(EMesh & em, ClipperLib::Paths & profile)
	{
		ClipperLib::Clipper clip;

		// generate the isosceles trapezoids corresponding to the edges
		for (size_t i=0; i<em.edge.size(); ++i)
		{
			// get the edge
			const EEdgeType & e = em.edge[i];
			if (e.IsD())
				continue;

			// get their vertices
			const EVertexType & v0 = *e.cV(0);
			const EVertexType & v1 = *e.cV(1);

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

			ClipperLib::Path edge_path;
			edge_path.push_back(ThisType::convertToIntPoint(p0));
			edge_path.push_back(ThisType::convertToIntPoint(p0 + (leftRot * r0)));
			edge_path.push_back(ThisType::convertToIntPoint(p1 + (leftRot * r1)));
			edge_path.push_back(ThisType::convertToIntPoint(p1));
			edge_path.push_back(ThisType::convertToIntPoint(p1 - (rightRot * r1)));
			edge_path.push_back(ThisType::convertToIntPoint(p0 - (rightRot * r0)));

			clip.AddPath(edge_path, ClipperLib::ptSubject, true);
		}

		// generate vertex circles
		for (size_t i=0; i<em.vert.size(); ++i)
		{
			const EVertexType & v = em.vert[i];
			if (v.IsD())
				continue;

			Coord2Type p(ScalarType(v.cP()[0]), ScalarType(v.cP()[1]));

			ClipperLib::Path circle;
			this->generateCircle(circle, p, v.cQ(), true);

			clip.AddPath(circle, ClipperLib::ptSubject, true);
		}

		// union all vertices circles and edges trapezoids
		profile.clear();
		clip.Execute(ClipperLib::ctUnion, profile, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	}

	double maxEdgeLength(void) const
	{
		const TessellatorSettings & ts = this->getConstTessellationSettings();
		if (!ts.area_constrained)
			return 1;

		return vcg::math::Sqrt(ts.max_area*2);
	}

	void setTessellationParameters(void)
	{
		// set tesselation settings from tessellation parameters
		this->m_tri_settings.area_constrained = true;
		this->m_tri_settings.max_area = m_params.tessellationParams.max_area; //  change to account for real dimension
		this->m_tri_settings.angle_constrained = true;
		this->m_tri_settings.min_angle = m_params.tessellationParams.min_angle;

		// prevent steiner points to be added on each boundary segment
		this->m_tri_settings.steiner_points = false;
		this->m_tri_settings.quiet = true;
	}

	// in order to preserve periodic boundary conditions we manually split the polygon edges in a consistent way
	void splitAllCountourEdges(ClipperLib::Paths & profile)
	{
		double maxLength = this->maxEdgeLength() * ThisType::ScaleFactor;

		ClipperLib::Paths newProfile;
		for (const ClipperLib::Path & p : profile)
		{
			ClipperLib::Path newPath;
			for (size_t i=0; i<p.size(); i++)
			{
				const ClipperLib::IntPoint & p0 = p[i];
				const ClipperLib::IntPoint & p1 = p[(i+1)%p.size()];
				newPath.push_back(p0);

				// split the edge in <refine_count> edges
				int refine_count = int(ceil(distance(p0, p1) / maxLength));

				if (refine_count > 1)
				{
					if (p1 < p0)
					{
						for (int i=(refine_count-1); i>0; i--)
						{
							double t = double(i)/refine_count;
							ClipperLib::IntPoint split_p = interpolate(p1, p0, t);
							newPath.push_back(split_p);
						}
					}
					else
					{
						for (int i=1; i<refine_count; i++)
						{
							double t = double(i)/refine_count;
							ClipperLib::IntPoint split_p = interpolate(p0, p1, t);
							newPath.push_back(split_p);
						}
					}
				}
			}

			newProfile.push_back(newPath);
		}

		profile = newProfile;
	}

	// needed to avoid vertices too close to each other, leading to artifacts in the triangulation
	struct EdgeCollapse
	{
		typedef enum
		{
			First  = 0,
			Second,
			MidPoint
		} CollapseType;

		int              index;
		ClipperLib::cInt length;
		CollapseType     type;

		bool operator < (const EdgeCollapse & rh) const
		{
			return this->length < rh.length;
		}

		ClipperLib::IntPoint collapsedPoint(const ClipperLib::Path & path) const
		{
			assert (this->index >= 0 && this->index < int(path.size()));
			switch (this->type)
			{
			case First    : return path[this->index];
			case Second   : return path[(this->index + 1) % path.size()];
			case MidPoint : return (path[this->index] + path[(this->index + 1) % path.size()]) / 2;
			}
			assert(0);
			return ClipperLib::IntPoint();
		}
	};

	template <typename BoundsType>
	void simplifyContour(ClipperLib::Paths & paths, const BoundsType & bounds, bool boundariesOnly = false) const
	{
		ClipperLib::cInt epsilon = ClipperLib::cInt(this->maxEdgeLength() * 0.4 * ThisType::ScaleFactor);

		for (size_t p=0; p<paths.size(); ++p)
		{
			while (true)
			{
				ClipperLib::Path & path = paths[p];
				std::vector<EdgeCollapse> collapseOps;

				// iterate all over the edges and insert collapse ops in an ordered set
				for (size_t i=0; i<path.size(); ++i)
				{
					const ClipperLib::IntPoint & p0 = path[i];
					const ClipperLib::IntPoint & p1 = path[(i+1)%path.size()];

					// create the edge-collapse operations preserving the points on the boundary
					bool boundary0 = isBoundary(p0, bounds);
					bool boundary1 = isBoundary(p1, bounds);

					if (!boundariesOnly)
					{
						if ((!(boundary0 && boundary1) && ClipperLib::cInt(distance(p0, p1)) < epsilon))
						{
							// generate collapse
							EdgeCollapse ec;
							ec.index  = i;
							ec.length = ClipperLib::cInt(distance(p0, p1));

							if (boundary0 == boundary1)
							{
								ec.type = EdgeCollapse::MidPoint;
							}
							else if (boundary0)
							{
								ec.type = EdgeCollapse::First;
							}
							else if (boundary1)
							{
								ec.type = EdgeCollapse::Second;
							}

							collapseOps.push_back(ec);
						}
					}
					else
					{
						// collapse boundaries
						if ((boundary0 && boundary1) && ClipperLib::cInt(distance(p0, p1)) < (epsilon/2)
						    /*&& !isBoundary((p0+p1)/2, bounds)*/)
						{
							EdgeCollapse ec;
							ec.index  = i;
							ec.length = ClipperLib::cInt(distance(p0, p1));
							ec.type   = EdgeCollapse::MidPoint;

							collapseOps.push_back(ec);
						}
					}
				}

				if (collapseOps.size() == 0)
					break;

				std::sort(collapseOps.begin(), collapseOps.end());
				// sorted by length, they represent all the candidate collapses

				// extract, in order, the collapse operations (border or not), and mark the adjacent ones as denied
				std::unordered_set<int>                       deniedIdxs;
				std::unordered_map<int, const EdgeCollapse *> allowedOps;

				for (const EdgeCollapse & ec : collapseOps)
				{
					if (deniedIdxs.count(ec.index) == 0)
					{
						// avoid reducing contours to edge or point
						if ((path.size() - allowedOps.size()) == 3)
							break;
						deniedIdxs.insert((ec.index + 1)%path.size());
						deniedIdxs.insert((ec.index - 1 + path.size())%path.size());
						allowedOps[ec.index] = &ec;
					}
				}

				// perform all the allowed collapses
				ClipperLib::Path newPath;
				for (int i=0; i<int(path.size()); ++i)
				{
					if (allowedOps.count(i) != 0)
					{
						// insert the new point
						const EdgeCollapse & ec = *allowedOps[i];
						newPath.push_back(ec.collapsedPoint(path));
					}
					else if (allowedOps.count((i-1+path.size())%path.size()) == 0)
					{
						newPath.push_back(path[i]);
					}
				}

				path = newPath;

				if (path.size() <= 3)
					break;
			}
		}
	}

	static void pathToEdgeMesh(const ClipperLib::Path & path, EMesh & em)
	{
		em.Clear();

		// fill nodes
		vcg::tri::Allocator<EMesh>::AddVertices(em, path.size());
		for (size_t i=0; i<em.vert.size(); ++i)
		{
			typename EMesh::CoordType &           p = em.vert[i].P();
			const ClipperLib::IntPoint & pp = path[i];
#ifdef use_xyz
			vcg::Point2d pd = vcg::Point2d(pp.X, pp.Y, pp.Z) / ThisType::ScaleFactor;
#else
			vcg::Point2d pd = vcg::Point2d(pp.X, pp.Y) / ThisType::ScaleFactor;
#endif
			p.X() = pd.X();
			p.Y() = pd.Y();
#ifdef use_xyz
			mp.Z() = pd.Z();
#else
			p.Z() = 0;
#endif
		}

		// fill edges
		vcg::tri::Allocator<EMesh>::AddEdges(em, path.size());
		for (size_t i=0; i<em.edge.size(); ++i)
		{
			typename EMesh::EdgeType & e = em.edge[i];

			e.V(0) = &em.vert[i];
			e.V(1) = &em.vert[(i+1) % path.size()];
		}
	}

	static bool compareClipperIntPoint (ClipperLib::IntPoint * i, ClipperLib::IntPoint * j)
	{
		return (*i < *j);
	}

	// correct small clipper errors to have identical periodical boundaries
	static void matchPeriodicVertices(ClipperLib::Paths & paths, const ClipperLib::IntRect & bounds)
	{
		static const ClipperLib::cInt Tolerance = 2000;

		bool fail = false;
		for (ClipperLib::Path & path : paths)
		{
			if (fail)
				break;

			for (int coord=0; coord<2; ++coord)
			{
				int b_pos0   = 0;                     // coordinate of border
				int b_pos1   = coord == 0 ? bounds.right : bounds.top; //ThisType::ScaleFactor; // coordinate of opposite border

				std::vector<ClipperLib::IntPoint *> bPoints0;
				std::vector<ClipperLib::IntPoint *> bPoints1;
				for (ClipperLib::IntPoint & p : path)
				{
					ClipperLib::cInt c = clipperPointCCoord(p,coord);
					if (c == b_pos0)
						bPoints0.push_back(&p);
					else if (c == b_pos1)
						bPoints1.push_back(&p);
				}

				if (bPoints0.size() != bPoints1.size())
				{
					// Fail
					std::cout << "The geometry is not 'well shaped' ! It can't be periodic!" << std::endl;

					assert(0);
					fail = true;
					break;
				}

				if (bPoints0.size() == 0)
					continue;

				// sort the points and set them coherently
				std::sort(bPoints0.begin(), bPoints0.end(), compareClipperIntPoint);
				std::sort(bPoints1.begin(), bPoints1.end(), compareClipperIntPoint);

				for (size_t i=0; i<bPoints0.size(); ++i)
				{
					ClipperLib::cInt diff = clipperPointCCoord(*bPoints1[i], 1-coord) - clipperPointCCoord(*bPoints0[i], 1-coord);
					if (abs(diff) > Tolerance)
					{
						// Fail
						std::cout << "Cannot find matching periodic counterpart for node ("
						          << bPoints0[i]->X << "," << bPoints0[i]->Y << ")" << std::endl
						          << "Candidate is (" << bPoints1[i]->X << "," << bPoints1[i]->Y << ")"
						          << "The geometry will not be periodic!" << std::endl;
						fail = true;
						break;
					}
					clipperPointCoord(*bPoints1[i], 1-coord) = clipperPointCCoord(*bPoints0[i], 1-coord);
				}
			}
		}
	}

	// decimations utility
	class SimplifyParameter : public vcg::BaseParameterClass
	{
	public:
		double maxLength;
	};

	typedef BasicVertexPair<typename TriMesh::VertexType> VertexPair;
	class EdgeCollapse2D : public vcg::tri::TriEdgeCollapse<TriMesh, VertexPair, EdgeCollapse2D>
	{
	public:
		typedef vcg::tri::TriEdgeCollapse<TriMesh, VertexPair, EdgeCollapse2D> BaseType;
		typedef typename BaseType::ScalarType                                  ScalarType;
		typedef vcg::BaseParameterClass                                        BaseParameterClass;
		inline EdgeCollapse2D(const VertexPair &p, int mark, BaseParameterClass * pp) : BaseType(p, mark, pp) {}

		inline ScalarType ComputePriority(BaseParameterClass *)
		{
			auto dist = vcg::SquaredDistance(this->pos.V(0)->cP(), this->pos.V(1)->cP());
			this->_priority = dist;
			return this->_priority;
		}

		inline bool IsFeasible(BaseParameterClass * p)
		{
			SimplifyParameter * sp = (SimplifyParameter *)p;
			// check distance
			if (vcg::Distance(this->pos.V(0)->cP(), this->pos.V(1)->cP()) >= sp->maxLength)
				return false;

			// check link conditions
			if (vcg::tri::EdgeCollapser<TriMesh,VertexPair>::LinkConditions(this->pos))
				return true;

			// if both vertices are on the border
			if (this->pos.V(0)->IsB() && this->pos.V(1)->IsB())
			{
				// check if the edge is border and belongs to a face with exaclty 2 border edges
				int countF[2] = {0, 0};
				for (int i=0; i<2; ++i)
				{
					for(vcg::face::VFIterator<typename TriMesh::FaceType> vfi(this->pos.V(i)); !vfi.End(); ++vfi)
					{
						countF[i]++;
					}
				}
				return (countF[0] == 1 && countF[1] > 1) || (countF[1] == 1 && countF[0] > 1);
			}

			return false;
		}

		inline void Execute(TriMesh & m, BaseParameterClass *)
		{
			// for inner vertices or border edge use midpoint
			// snap to border otherwise
			int bb = 0;
			if (this->pos.V(0)->IsB()) bb += 1;
			if (this->pos.V(1)->IsB()) bb += 2;
			switch(bb)
			{
			case 0 :
			case 3 :
				vcg::tri::EdgeCollapser<TriMesh,VertexPair>::Do(m, this->pos, (this->pos.V(0)->cP() + this->pos.V(1)->cP()) / 2.0);
				break;
			case 1 :
			case 2 :
				vcg::tri::EdgeCollapser<TriMesh,VertexPair>::Do(m, this->pos, this->pos.V(bb-1)->cP());
				break;
			}
		}
	};

	void decimateMesh(TriMesh & mesh)
	{
		vcg::tri::UpdateTopology<TriMesh>::VertexFace(mesh);
		vcg::tri::UpdateFlags<TriMesh>::VertexBorderFromNone(mesh);

		SimplifyParameter sp;
		sp.maxLength = this->maxEdgeLength() * 0.25;
		vcg::LocalOptimization<TriMesh> decimator(mesh, &sp);

		decimator.template Init<EdgeCollapse2D>();
		decimator.DoOptimization();
	}

	// clipper related functions
	typedef vcg::Box3<typename EMesh::ScalarType> EMeshBoxType;

	static ClipperLib::IntRect getClipperBbox(const EMeshBoxType & bbox)
	{
		ClipperLib::IntRect rect;
		ClipperLib::IntPoint min = ThisType::convertToIntPoint(Coord2Type(bbox.min[0], bbox.min[1]));
		ClipperLib::IntPoint max = ThisType::convertToIntPoint(Coord2Type(bbox.max[0], bbox.max[1]));
		rect.left    = min.X;
		rect.bottom  = min.Y;
		rect.right   = max.X;
		rect.top     = max.Y;
		return rect;
	}

	// generating circle functions
	void generateCircle(std::vector<Coord2Type> & circle, const Coord2Type & center, ScalarType radius, bool ccw = true)
	{
		static const size_t MinCircleResolution = 36;

		// snap to multiple of 4 for symmetry reasons
		size_t circle_res = std::max(MinCircleResolution,
		                             size_t(4*ceil((M_PI_2 * radius) / this->maxEdgeLength())) );

		circle.clear();
		for (size_t i=0; i<circle_res; ++i)
		{
			double angle = (M_PI * 2.0 * i) / circle_res;
			if (!ccw) angle = -angle;

			vcg::Point2d p(cos(angle), sin(angle));
			p *= radius;
			p += vcg::Point2d::Construct(center);

			circle.push_back(Coord2Type::Construct(p));
		}
	}

	void generateCircle(ClipperLib::Path & circle, const Coord2Type & center, ScalarType radius, bool ccw = true)
	{
		std::vector<Coord2Type> c;
		this->generateCircle(c, center, radius, ccw);

		circle.clear();
		for (auto it=c.begin(); it!=c.end(); ++it)
			circle.push_back(scaleToIntPoint(*it, ScalarType(ThisType::ScaleFactor)));
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

	// marker type for segments
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

	// query closest segment
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

#endif // EDGEMESHPATTERN_H
