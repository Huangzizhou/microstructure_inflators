#ifndef EDGEMESHTYPE_H
#define EDGEMESHTYPE_H

#include <vcg/complex/complex.h>

class EdgeType;
class VertexType;

struct EUsedTypes : public vcg::UsedTypes< vcg::Use<VertexType>::AsVertexType, vcg::Use<EdgeType>::AsEdgeType > {};

class VertexType : public vcg::Vertex< EUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::Qualityd, vcg::vertex::BitFlags, vcg::vertex::VEAdj> {};
class EdgeType   : public vcg::Edge< EUsedTypes, vcg::edge::VertexRef, vcg::edge::Qualityd, vcg::edge::BitFlags, vcg::edge::EEAdj, vcg::edge::VEAdj> {};
class EMesh      : public vcg::tri::TriMesh< std::vector<VertexType>, std::vector<EdgeType> > {};


#include <vcg/complex/algorithms/update/bounding.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/import.h>

template <class MESH>
class EdgeMeshUtils
{
public:
	static void cleanMesh(MESH & mesh)
	{
		vcg::tri::Clean<MESH>::RemoveUnreferencedVertex(mesh);
		vcg::tri::Clean<MESH>::RemoveDuplicateVertex(mesh);
	}

	static void compactMesh(MESH & mesh)
	{
		vcg::tri::Allocator<MESH>::CompactEveryVector(mesh);
	}

	///
	/// \brief exportObj
	/// \param mesh the mesh to be exported.
	/// \param filePath the file path to which export to.
	/// \param export_attrib if false it does not export per vertex attributes (i.e. exports geometry only).
	/// \return true if the export was successful, false otherwise.
	///
	static bool exportObj(MESH & mesh, const std::string & filePath, bool export_attrib = true)
	{
		typedef typename vcg::tri::io::ExporterOBJ<MESH> Exporter;

		int cap = Exporter::GetExportMaskCapability();
		int mask = export_attrib ? createExportMask(cap, mesh) : 0;

		int res = Exporter::Save(mesh, filePath.c_str(), mask);

		return res == Exporter::E_NOERROR;
	}

	static bool importObj(MESH & mesh, const std::string & filePath)
	{
		typedef typename vcg::tri::io::ImporterOBJ<MESH> Importer;

		int load_mask;

		if (!Importer::LoadMask(filePath.c_str(), load_mask))
			return false;

		int result = Importer::Open(mesh, filePath.c_str(), load_mask);
		if (Importer::ErrorCritical(result))
			return false;

		vcg::tri::UpdateBounding<MESH>::Box(mesh);

		return true;
	}

	static int createExportMask(const int exportCapabilities, const MESH & mesh)
	{
		int mask = vcg::tri::io::Mask::IOM_NONE;

		if (vcg::tri::HasPerVertexNormal   (mesh)) mask |= (vcg::tri::io::Mask::IOM_VERTNORMAL & exportCapabilities);

		if (vcg::tri::HasPerVertexColor    (mesh)) mask |= (vcg::tri::io::Mask::IOM_VERTCOLOR & exportCapabilities);

		if (vcg::tri::HasPerVertexTexCoord (mesh)) mask |= (vcg::tri::io::Mask::IOM_VERTTEXCOORD & exportCapabilities);

		return mask;
	}
};

#endif // EDGEMESHTYPE_H
