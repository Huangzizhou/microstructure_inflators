#ifndef WIREMESHEMBEDDING_H
#define WIREMESHEMBEDDING_H

#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>

#include <cstdlib>

template <class EMesh, class PolyMesh>
class WireMeshEmbedding
{
public:
	typedef WireMeshEmbedding<EMesh, PolyMesh> ThisType;

	typedef typename EMesh::CoordType          ECoordType;

	typedef typename PolyMesh::CoordType       PCoordType;
	typedef typename PolyMesh::FaceType        PFaceType;
	typedef typename PolyMesh::FacePointer     PFacePointer;\


	struct QuadParametrization
	{
		// these represents the indices of the vertices
		// defining the positive X-axis (V(index0) ->V(index1))
		char index0;

		PCoordType interpolate(const vcg::Point2d & t, const PFaceType & face) const
		{
			char index1 = (index0+1)%4;
			char index2 = (index0+2)%4;
			char index3 = (index0+3)%4;
			      PCoordType   p0 = face.cP(index0);
			      PCoordType   p1 = face.cP(index1);
			const PCoordType & p2 = face.cP(index2);
			const PCoordType & p3 = face.cP(index3);

			// bilinear interpolation
			p0 = p0 * (1 - t.X()) + p1 * (t.X());
			p1 = p3 * (1 - t.X()) + p2 * (t.X());
			return p0 * (1 - t.Y()) + p1 * (t.Y());
		}
	};

	typedef typename PolyMesh::template PerFaceAttributeHandle<QuadParametrization> QuadParamHandle;

	static QuadParamHandle getQuadParametrizationHandle(PolyMesh & pmesh)
	{
		QuadParamHandle h = vcg::tri::Allocator<PolyMesh>::template FindPerFaceAttribute<QuadParametrization>(pmesh, ThisType::ParametrizationAttributeName());
		assert(vcg::tri::Allocator<PolyMesh>::IsValidHandle(pmesh, h));
		return h;
	}

	// used to planarize and obtain a coherent mapping of edge meshes among quads
	static void preprocessQuadMesh(PolyMesh & pmesh)
	{
		if (pmesh.VN() == 0)
			return;

		// planarize
		planarizeQuadMesh(pmesh);

		// compute adjacency
		vcg::tri::UpdateTopology<PolyMesh>::FaceFace(pmesh);
		// computer border
		vcg::tri::UpdateFlags<PolyMesh>::FaceBorderFromFF(pmesh);

		// create a coherent (if possibile) parametrization for the quad mesh
		createParametrization(pmesh);
	}

	// assume em is a normalized edge mesh
	static void embedWireMesh(EMesh & em, const PFaceType & f, PolyMesh & pmesh)
	{
		QuadParametrization qpar = getQuadParametrizationHandle(pmesh)[&f];
		for (size_t i=0; i<em.vert.size(); ++i)
		{
			ECoordType & p = em.vert[i].P();
			vcg::Point2d t(p[0], p[1]);
			p = ECoordType::Construct(qpar.interpolate(t, f));
		}
	}

private:
	static std::string ParametrizationAttributeName(void)
	{
		return "parametrization";
	}

	static void planarizeQuadMesh(PolyMesh & pmesh)
	{
		for (auto vi=pmesh.vert.begin(); vi!=pmesh.vert.end(); vi++)
		{
			vi->P()[2] = 0;
		}
	}

	static void createParametrization(PolyMesh & pmesh)
	{
		vcg::tri::RequireFFAdjacency(pmesh);
		vcg::tri::RequirePerFaceFlags(pmesh);

		// clear visited flag
		vcg::tri::UpdateFlags<PolyMesh>::FaceClearV(pmesh);

		QuadParamHandle handle =
		        vcg::tri::Allocator<PolyMesh>::template GetPerFaceAttribute<QuadParametrization>(pmesh, ThisType::ParametrizationAttributeName());

		// propagate if we encounter faces not visited
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			PFacePointer fp = &pmesh.face[i];
			if (!fp->IsV())
			{
				// set the parametrization for the first face
				setOptimalParametrization(handle, fp);

				// propagate to all other quads in the same connected component
				std::vector<PFacePointer> toPropagate;
				toPropagate.push_back(fp);

				for (size_t k=0; k<toPropagate.size(); ++k)
				{
					PFacePointer f = toPropagate[k];
					for (char e=0; e<4; ++e)
					{
						PFacePointer fAdj = setAdjacentParametrization(handle, f, e);
						if (fAdj)
							toPropagate.push_back(fAdj);
					}
				}
			}
		}
	}

	static void setOptimalParametrization(QuadParamHandle handle, PFacePointer f)
	{
		// TODO CHANGE TO OPTIMAL PARAMETRIZATION
		// right now, it gets the first quad and assumes a configuration like this:
		//
		//        Y
		//        ^
		//        |
		//  3-----------2
		//  |     |     |
		//  |     |     |
		//  |     +---- | --> X
		//  |           |
		//  |           |
		//  0-----------1

		assert(f);

		QuadParametrization & qp = handle[f];
		qp.index0 = 0;

		f->SetV();
	}

	static PFacePointer setAdjacentParametrization(QuadParamHandle & handle, const PFacePointer face, char adjIdx)
	{
		assert(adjIdx>=0 && adjIdx<4);

		PFaceType * adjFace = face->FFp(adjIdx);
		if (adjFace != NULL && !adjFace->IsV())
		{
			char edgeAdj = face->cFFi(adjIdx);
			handle[adjFace].index0 = (handle[face].index0 - adjIdx + (edgeAdj+2) + 4)%4;

			adjFace->SetV();
			return adjFace;
		}
		return NULL;
	}
};

#endif // WIREMESHEMBEDDING_H
