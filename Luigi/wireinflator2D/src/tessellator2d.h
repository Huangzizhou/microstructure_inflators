#ifndef TESSELLATOR2D_H
#define TESSELLATOR2D_H

#include "clipperHelper.h"
#include "triangle.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/clean.h>

#include <iostream>
#include <iomanip>

struct Tessellator2DSettings
{
	bool plsg;
	bool area_constrained;
	double max_area;
	bool angle_constrained;
	double min_angle;
	bool steiner_points;
	bool quiet;

	Tessellator2DSettings(void)
	    : plsg             (true)
	    , area_constrained (true)
	    , max_area         (0.001)
	    , angle_constrained(true)
	    , min_angle        (30)
	    , steiner_points   (true)
	    , quiet            (false)
	{
		;
	}

	std::string getTriangleSwitches(void) const
	{
		// z  -> uses 0 based index
		// p  -> triangulates a Planar Straight Line Graph
		// a  -> set maximum triangle area
		// q  -> quality tessellation (min angle constraint)
		// Y  -> no steiner points (AVOID)
		// YY -> steiner point added except on boundary segments
		// Q  -> quiet

		std::ostringstream ss;
		ss << 'z';
		if (plsg) ss << 'p';
		if (quiet) ss << 'Q';
		if (!steiner_points) ss << "YY";
		if (angle_constrained)
		{
			ss << 'q';
			if (min_angle > 0)
				ss << std::fixed << std::setprecision(9) << min_angle;
		}
		if (area_constrained)
		{
			ss << 'a';
			if (max_area > 0)
				ss << std::fixed << std::setprecision(9) << max_area;
		}
		return ss.str();
	}

	bool isValid() const
	{
		return (!area_constrained  || max_area > 0)
		    && (!angle_constrained || (min_angle > 0 && min_angle <= 33));
	}
};

template <class VcgMesh>
class Tessellator2D
{
public:
	typedef typename VcgMesh::ScalarType ScalarType;
	typedef typename VcgMesh::CoordType  CoordType;

	static bool execute(const Tessellator2DSettings & settings, const ClipperLib::Paths & polygon, const ScalarType clipperScaleFactor, VcgMesh & out_mesh)
	{
		// Do triangulation using only boundary segments (a Planar Straight Line Graph)
		if (polygon.size() == 0) return false;

		// create in and out structs for triangle
		triangulateio in, out;
		memset(&in , 0, sizeof(triangulateio));
		memset(&out, 0, sizeof(triangulateio));

		// count points
		for (size_t i=0; i<polygon.size(); ++i)
		{
			const ClipperLib::Path & p = polygon[i];
			in.numberofpoints   += p.size();
		}
		in.numberofsegments  = in.numberofpoints;

		// initialize lists
		in.pointlist         = (REAL *) malloc(in.numberofpoints   * 2 * sizeof(REAL));
		in.segmentlist       = (int *)  malloc(in.numberofsegments * 2 * sizeof(int));
		in.segmentmarkerlist = (int *)  malloc(in.numberofsegments * sizeof(int));

		// fill triangle input structure with points and segments
		int startIdx = 0;
		for (size_t i=0; i<polygon.size(); ++i)
		{
			const ClipperLib::Path & p = polygon[i];
			for (size_t j=0; j<p.size(); ++j)
			{
				int idx = startIdx + j;
				in.pointlist[idx*2]     = REAL(p[j].X) / clipperScaleFactor;
				in.pointlist[idx*2 + 1] = REAL(p[j].Y) / clipperScaleFactor;

				in.segmentlist[idx*2]     = startIdx + j;
				in.segmentlist[idx*2 + 1] = startIdx + ((j+1)%p.size());
				in.segmentmarkerlist[idx]   = 1; // mark each segment as boundary
			}

			startIdx += p.size();
		}

		//find holes
		ClipperLib::Path holePoints = getPointsInHoles(polygon);
		if (holePoints.size() > 0)
		{
			in.numberofholes = holePoints.size();
			in.holelist = (REAL *) malloc(in.numberofholes * 2 * sizeof(REAL));
			for (int i=0; i<in.numberofholes; ++i)
			{
				in.holelist[i*2]     = REAL(holePoints[i].X) / clipperScaleFactor;
				in.holelist[i*2 + 1] = REAL(holePoints[i].Y) / clipperScaleFactor;
			}
		}

		// launch the triangulation with parameters from settings
		std::string cmd = settings.getTriangleSwitches();
		triangulate((char *)cmd.c_str(), &in, &out, (triangulateio *) NULL);

		// output check
		if (out.numberoftriangles <= 0)
		{
			std::cout << "unable to perform trianglulation!" << std::endl << std::flush;
			return false;
		}

		if (!settings.quiet)
			std::cout << std::flush;

		// create the vcg mesh
		out_mesh.Clear();
		vcg::tri::Allocator<VcgMesh>::AddVertices(out_mesh, out.numberofpoints);
		vcg::tri::Allocator<VcgMesh>::AddFaces(out_mesh, out.numberoftriangles);
		for (int i=0; i<out.numberofpoints; i++)
		{
			// vertices
			out_mesh.vert[i].P() = CoordType(ScalarType(out.pointlist[i*2]), ScalarType(out.pointlist[(i*2)+1]), 0);
		}
		for (int i=0; i<out.numberoftriangles; i++)
		{
			// faces
			int index0 = out.trianglelist[i*3];
			int index1 = out.trianglelist[(i*3)+1];
			int index2 = out.trianglelist[(i*3)+2];
			out_mesh.face[i].V(0) = &out_mesh.vert[index0];
			out_mesh.face[i].V(1) = &out_mesh.vert[index1];
			out_mesh.face[i].V(2) = &out_mesh.vert[index2];
		}
		deallocateTriangulation(out);

		// Clean steps
		vcg::tri::Clean<VcgMesh>::RemoveUnreferencedVertex(out_mesh);
		vcg::tri::Clean<VcgMesh>::RemoveDuplicateVertex(out_mesh);
		vcg::tri::Allocator<VcgMesh>::CompactEveryVector(out_mesh);

		// Update normals and bbox
		vcg::tri::UpdateNormal<VcgMesh>::PerFaceNormalized(out_mesh);
		vcg::tri::UpdateNormal<VcgMesh>::PerVertexFromCurrentFaceNormal(out_mesh);
		vcg::tri::UpdateBounding<VcgMesh>::Box(out_mesh);

		return true;
	}

	static void deallocateTriangulation(triangulateio & triangleio)
	{
		// deallocate the triangle library output
		if (triangleio.edgelist)              trifree((VOID *)triangleio.edgelist);
		if (triangleio.edgemarkerlist)        trifree((VOID *)triangleio.edgemarkerlist);
		if (triangleio.holelist)              trifree((VOID *)triangleio.holelist);
		if (triangleio.neighborlist)          trifree((VOID *)triangleio.neighborlist);
		if (triangleio.normlist)              trifree((VOID *)triangleio.normlist);
		if (triangleio.pointattributelist)    trifree((VOID *)triangleio.pointattributelist);
		if (triangleio.pointlist)             trifree((VOID *)triangleio.pointlist);
		if (triangleio.pointmarkerlist)       trifree((VOID *)triangleio.pointmarkerlist);
		if (triangleio.regionlist)            trifree((VOID *)triangleio.regionlist);
		if (triangleio.segmentlist)           trifree((VOID *)triangleio.segmentlist);
		if (triangleio.segmentmarkerlist)     trifree((VOID *)triangleio.segmentmarkerlist);
		if (triangleio.trianglearealist)      trifree((VOID *)triangleio.trianglearealist);
		if (triangleio.triangleattributelist) trifree((VOID *)triangleio.triangleattributelist);
		if (triangleio.trianglelist)          trifree((VOID *)triangleio.trianglelist);
	}
};

#endif // TESSELLATOR2D_H
