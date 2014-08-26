#ifndef EDGEMESHPATTERN_H
#define EDGEMESHPATTERN_H

#include "Pattern2D.h"
#include "EdgeMeshType.h"
#include "InflatorParameters.h"
#include "tessellator2d.h"
#include <stdlib.h>
#include <cmath>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>

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


template <class TriMesh, typename Parameters>
class EdgeMeshPattern : public Pattern2D<TriMesh>
{
public:
	typedef EdgeMeshPattern<TriMesh, Parameters> ThisType;
	typedef Pattern2D<TriMesh>                   BaseType;
	typedef Parameters                           PatternParameters;
	typedef typename BaseType::ScalarType        ScalarType;
	typedef typename BaseType::Coord2Type        Coord2Type;

	typedef Tessellator2DSettings TessellatorSettings;

	EdgeMeshPattern(void)
	    : m_repeat_x(1)
	    , m_repeat_y(1)
	    , m_scale(1)
	{
		;
	}

	TessellatorSettings & getTessellationSettings(void)
	{
		return m_tri_settings;
	}

	size_t & repeatX(void)
	{
		return m_repeat_x;
	}

	size_t & repeatY(void)
	{
		return m_repeat_x;
	}

	double & scale(void)
	{
		return m_scale;
	}

	struct InputParameters
	{
		PatternParameters      patternParams;
		TessellationParameters tessellationParams;
	};

	InputParameters & params(void)
	{
		return m_params;
	}

	bool generate(void)
	{
		// set tesselation settings from settings
		this->m_tri_settings.area_constrained = true;
		this->m_tri_settings.max_area = m_params.tessellationParams.max_area;
		this->m_tri_settings.angle_constrained = true;
		this->m_tri_settings.min_angle = m_params.tessellationParams.min_angle;

		// prevent steiner points to be added on each boundary segment
		this->m_tri_settings.steiner_points = false;
		this->m_tri_settings.quiet = true;

		this->m_paths.clear();

		generateOneElement();

		assert(m_scale > 0);

		ClipperLib::Clipper clipper;
		ClipperLib::IntRect bounds = { 0, ThisType::ScaleFactor, ThisType::ScaleFactor, 0 };

		for (unsigned int x=0; x<m_repeat_x; ++x)
		{
			for (unsigned int y=0; y<m_repeat_y; ++y)
			{
				ClipperLib::IntPoint t((bounds.right - bounds.left)*x, (bounds.top - bounds.bottom)*y);
				clipper.AddPaths((this->m_base_paths + t), ClipperLib::ptSubject, true);
			}
		}

		clipper.Execute(ClipperLib::ctUnion, this->m_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

		splitAllCountourEdges(this->m_paths);
		this->m_paths *= m_scale;

		return true;
	}

	bool tile(const Array2D<PatternParameters *> & grid)
	{
		// set tesselation settings from settings
		this->m_tri_settings.area_constrained = true;
		this->m_tri_settings.max_area = m_params.tessellationParams.max_area;
		this->m_tri_settings.angle_constrained = true;
		this->m_tri_settings.min_angle = m_params.tessellationParams.min_angle;

		// prevent steiner points to be added on each boundary segment
		this->m_tri_settings.steiner_points = false;
		this->m_tri_settings.quiet = true;

		this->m_paths.clear();

		assert(m_scale > 0);

		ClipperLib::Clipper clipper;

		ClipperLib::IntRect bounds = { 0, ThisType::ScaleFactor, ThisType::ScaleFactor, 0 };

		for (unsigned int x=0; x<grid.width(); ++x)
		{
			for (unsigned int y=0; y<grid.height(); ++y)
			{
				if (grid(x,y) == NULL)
					continue;

				generateOneElement(*grid(x,y));

				ClipperLib::IntPoint t((bounds.right - bounds.left)*x, (bounds.top - bounds.bottom)*y);
				clipper.AddPaths((this->m_base_paths + t), ClipperLib::ptSubject, true);
			}
		}

		clipper.Execute(ClipperLib::ctUnion, this->m_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

		splitAllCountourEdges(this->m_paths);
		this->m_paths *= m_scale;

		return true;
	}

	bool tessellate(TriMesh & mesh)
	{
		return Tessellator2D<TriMesh>::execute(m_tri_settings, this->m_paths, this->ScaleFactor, mesh);
	}

protected:
	/// parameters
	TessellatorSettings m_tri_settings;
	size_t m_repeat_x;
	size_t m_repeat_y;
	double m_scale;

	InputParameters m_params;

	struct ParameterChange
	{
		typedef enum
		{
			Translation,
			Radius
		} OperationType;

		OperationType             type;
		std::map<int, Coord2Type> node_ops;
	};

	virtual void getParameterOperations(std::vector<ParameterChange> & params_op) const = 0;

	virtual void generateUnmodifiedEdgeMesh(EMesh & em) = 0;

	virtual void getEdgeMesh(EMesh & em)
	{
		getEdgeMesh(em, m_params.patternParams);
	}

	virtual void generateOneElement(void)
	{
		generateOneElement(this->m_params.patternParams);
	}

	void generateOneElement(const PatternParameters & pars)
	{
		this->m_base_paths.clear();

		EMesh em;
		getEdgeMesh(em, pars);

		vcg::tri::RequireVEAdjacency(em);
		vcg::tri::RequireEEAdjacency(em);

		ClipperLib::Clipper clip;

		// generate the isosceles trapezoids corresponding to the edges
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

			ClipperLib::Path edge_path;
			edge_path.push_back(ThisType::convertToIntPoint(p0));
			edge_path.push_back(ThisType::convertToIntPoint(p0 + (leftRot * r0)));
			edge_path.push_back(ThisType::convertToIntPoint(p1 + (leftRot * r1)));
			edge_path.push_back(ThisType::convertToIntPoint(p1));
			edge_path.push_back(ThisType::convertToIntPoint(p1 - (rightRot * r1)));
			edge_path.push_back(ThisType::convertToIntPoint(p0 - (rightRot * r0)));

			clip.AddPath(edge_path, ClipperLib::ptSubject, true);
		}

		// get the bbox
		ClipperLib::IntRect bbox = { 0, ThisType::ScaleFactor, ThisType::ScaleFactor, 0 };
		ClipperLib::Path clip_bbox = convertToPath(bbox);

		// generate vertex circles
		for (size_t i=0; i<em.vert.size(); ++i)
		{
			const EMesh::VertexType & v = em.vert[i];
			if (v.IsD())
				continue;

			Coord2Type p(ScalarType(v.cP()[0]), ScalarType(v.cP()[1]));

			ClipperLib::Path circle;
			this->generateCircle(circle, p, v.cQ(), true);

			clip.AddPath(circle, ClipperLib::ptSubject, true);
		}

		// union all vertices circles and edges trapezoids
		ClipperLib::Paths tmp;
		clip.Execute(ClipperLib::ctUnion, tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		clip.Clear();

		// clip to bounding box
		clip.AddPaths(tmp, ClipperLib::ptSubject, true);
		clip.AddPath(clip_bbox, ClipperLib::ptClip, true);
		tmp.clear();
		clip.Execute(ClipperLib::ctIntersection, this->m_base_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

		// optimize path
		optimizeContours(this->m_base_paths);
	}

	void getEdgeMesh(EMesh & em, const Parameters & pars)
	{
		typedef typename EMesh::CoordType CoordType;

		this->generateUnmodifiedEdgeMesh(em);

		std::vector<ParameterChange> params_op;
		this->getParameterOperations(params_op);
		assert(params_op.size() == Parameters::NumberOfParameters);

		for (size_t i=0; i<params_op.size(); ++i)
		{
			const ParameterChange & p = params_op[i];
			switch (p.type)
			{
			case ParameterChange::Radius:
				for (auto & op : p.node_ops)
				{
					em.vert[op.first].Q() = pars.cParameter(i);
				}
				break;
			case ParameterChange::Translation:
				for (auto & op : p.node_ops)
				{
					em.vert[op.first].P() +=
					        CoordType(op.second[0], op.second[1], 0) * pars.cParameter(i);
				}
				break;
			default:
				assert(0);
			}
		}

		vcg::tri::UpdateTopology<EMesh>::VertexEdge(em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(em);
	}

	double maxEdgeLength(void)
	{
		const TessellatorSettings & ts = this->getTessellationSettings();
		if (!ts.area_constrained)
			return 1;

		return vcg::math::Sqrt(ts.max_area*2);
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
				const ClipperLib::IntPoint * p0 = &p[i];
				const ClipperLib::IntPoint * p1 = &p[(i+1)%p.size()];
				newPath.push_back(*p0);

				bool swapped = (*p1 < *p0);
				if (swapped)
					std::swap(p1, p0);

				// split the edge in <refine_count> edges
				int refine_count = int(distance(*p0, *p1) / maxLength/* + 0.5*/) + 1;

				if (swapped)
				{
					for (int i=(refine_count-1); i>0; i--)
					{
						double t = double(i)/refine_count;
						ClipperLib::IntPoint split_p = interpolate(*p0, *p1, t);
						newPath.push_back(split_p);
					}
				}
				else
				{
					for (int i=1; i<refine_count; i++)
					{
						double t = double(i)/refine_count;
						ClipperLib::IntPoint split_p = interpolate(*p0, *p1, t);
						newPath.push_back(split_p);
					}
				}
			}
			newProfile.push_back(newPath);
		}

		profile = newProfile;
	}

	// needed to avoid vertices too close to each other, leading to artifacts in the triangulation
	void optimizeContours(ClipperLib::Paths & profile)
	{
		// TODO fix: multiple null segments are generated
		ClipperLib::cInt epsilon = ClipperLib::cInt(double(0.00075) * ThisType::ScaleFactor);

		// get bounds
		ClipperLib::IntRect bounds;
		{
			ClipperLib::Clipper c;
			c.AddPaths(profile, ClipperLib::ptSubject, true);
			bounds = c.GetBounds();
		}

		ClipperLib::Paths newProfile;
		for (const ClipperLib::Path & p : profile)
		{
			ClipperLib::Path newPath;
			for (size_t i=0; i<p.size(); i++)
			{
				const ClipperLib::IntPoint p0 = p[i];
				const ClipperLib::IntPoint p1 = p[(i+1)%p.size()];

				if (ClipperLib::cInt(distance(p0, p1)) < epsilon)
				{
					// merge
					i++;
					// choose to preserve the point on the boundary
					if (p1.X == bounds.left   || p1.X == bounds.right ||
					    p1.Y == bounds.bottom || p1.Y == bounds.left)
					{
						newPath.push_back(p1);
					}
					else
					{
						newPath.push_back(p0);
					}
					// if the index wraps (p1 is p[0]) we don't mind, we add an additional redundant point
				}
				else
				{
					newPath.push_back(p0);
				}
			}
			newProfile.push_back(newPath);
		}

		profile = newProfile;
	}

	void generateCircle(std::vector<Coord2Type> & circle, const Coord2Type & center, ScalarType radius, bool ccw = true)
	{
		static const size_t MinCircleResolution = 36;

		size_t circle_res = std::max(MinCircleResolution,
		                             size_t(ceil((M_PI * 2.0 * radius) / this->maxEdgeLength())));

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

	void generateCircle(ClipperLib::Path & circle, const Coord2Type & center, ScalarType radius, bool ccw = true, int scaleFactor = ThisType::ScaleFactor)
	{
		std::vector<Coord2Type> c;
		this->generateCircle(c, center, radius, ccw);

		circle.clear();
		for (auto it=c.begin(); it!=c.end(); ++it)
			circle.push_back(scaleToIntPoint(*it, ScalarType(scaleFactor)));
	}

};

#endif // EDGEMESHPATTERN_H
