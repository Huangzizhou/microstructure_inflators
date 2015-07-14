#ifndef WIREMESHEMBEDDING_H
#define WIREMESHEMBEDDING_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cstdlib>
#include <math.h>

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
		if (pmesh.VN() == 0 || pmesh.FN() == 0)
			return;

		// planarize
		planarizeQuadMesh(pmesh);

		// compute adjacencies
		vcg::tri::UpdateTopology<PolyMesh>::FaceFace(pmesh);
		// vf adj does not exist for polygonal meshes
		auto vf = vcg::tri::Allocator<PolyMesh>::template GetPerVertexAttribute<std::vector<PFacePointer> >(pmesh, VFAttributeName());
		for (auto fi=pmesh.face.begin(); fi!=pmesh.face.end(); ++fi)
		{
			for (int i=0; i<fi->VN(); ++i)
				vf[fi->V(i)].push_back(&(*fi));
		}

		// compute border
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

	static std::vector<PFacePointer> adjacentFaceEdge(PolyMesh & pmesh, const PFaceType & f, unsigned char edgeIdx)
	{
		std::vector<PFacePointer> ret;
		if (edgeIdx >= 0 && edgeIdx < f.VN())
		{
			QuadParamHandle qpar = getQuadParametrizationHandle(pmesh);
			unsigned char actualIndex = (qpar[&f].index0 + edgeIdx) % 4;
			PFacePointer af = f.FFp(actualIndex);
			if (af != &f)
				ret.push_back(af);
		}
		return ret;
	}

	static std::vector<PFacePointer> adjacentFaceVertex(PolyMesh & pmesh, const PFaceType & f, unsigned char vertIdx)
	{
		std::vector<PFacePointer> ret;
		if (vertIdx >= 0 && vertIdx < f.VN())
		{
			QuadParamHandle qpar = getQuadParametrizationHandle(pmesh);
			unsigned char actualIndex = (qpar[&f].index0 + vertIdx) % 4;

			auto vf = vcg::tri::Allocator<PolyMesh>::template GetPerVertexAttribute<std::vector<PFacePointer> >(pmesh, VFAttributeName());
			std::vector<PFacePointer> & adjFaces = vf[f.V(actualIndex)];
			for (PFacePointer af : adjFaces)
			{
				if (af != &f)
					ret.push_back(af);
			}
		}
		return ret;
	}

private:
	static std::string ParametrizationAttributeName(void)
	{
		return "parametrization";
	}

	static std::string VFAttributeName(void)
	{
		return "vf_adj";
	}

	static void planarizeQuadMesh(PolyMesh & pmesh)
	{
		for (auto vi=pmesh.vert.begin(); vi!=pmesh.vert.end(); vi++)
		{
			vi->P()[2] = 0;
		}
	}

	// MHS on JUL14, 2015:
	// This method returns the deformation corresponding to the equivalent parallelogram for each
	// bilinear quad in the quad mesh.
	// TODO: maybe you should make this public so that you can call it from outside ...
	// TODO: there are a lot of cout statements for debug purposes that should be deleted later ...
	static void getDefs(PolyMesh & pmesh)
	{
		for (auto fc = pmesh.face.begin(); fc != pmesh.face.end(); fc++)
		{

			/* PCoordType myV1 = fc->cP(1) - fc->cP(0); */
			/* PCoordType myV2 = fc->cP(2) - fc->cP(1); */

			/* double crossProduct = myV1[0] * myV2[1] - myV1[1] * myV2[0]; */

			/* cout << endl; */
			/* if (crossProduct > 0) */
			/* 	cout << "CCW" << endl; */
			/* else */ 
			/* 	cout << "CW" << endl; */

			/* cout << fc->cP(0)[0] << ", " << fc->cP(0)[1] << endl; */
			/* cout << fc->cP(1)[0] << ", " << fc->cP(1)[1] << endl; */
			/* cout << fc->cP(2)[0] << ", " << fc->cP(2)[1] << endl; */
			/* cout << fc->cP(3)[0] << ", " << fc->cP(3)[1] << endl; */


			
			/* char lowerLeftVertIndx = 0; */
			/* for (int i = 1; i < fc->VN(); i++) */
			/* 	if (fc->cP(i)[0] <= fc->cP(lowerLeftVertIndx)[0] && */ 
			/* 		fc->cP(i)[1] <= fc->cP(lowerLeftVertIndx)[1]) */
			/* 		lowerLeftVertIndx = i; */

			/* cout << endl << "lower left vertex id is ..." << lowerLeftVertIndx << endl; */


			/* PCoordType p1 = fc->cP((lowerLeftVertIndx + 2)%4); */
			/* PCoordType p2 = fc->cP((lowerLeftVertIndx + 3)%4); */
			/* PCoordType p3 = fc->cP((lowerLeftVertIndx + 0)%4); */
			/* PCoordType p4 = fc->cP((lowerLeftVertIndx + 1)%4); */

			PCoordType p1 = fc->cP(0);
			PCoordType p2 = fc->cP(1);
			PCoordType p3 = fc->cP(2);
			PCoordType p4 = fc->cP(3);

			// Explain how the equivalent parallelograms are computed:
			// refer to the corresponding paper
			// assume face vertices are defined in a clock-wise manner ... 
			// note that in the paper it is counter-clock-wise ...
			
			PCoordType p1_equivalent = (p1 + p1 + p1) / 4.0 + (p2 - p3 + p4) / 4.0; 
			PCoordType p2_equivalent = (p2 + p2 + p2) / 4.0 + (p3 - p4 + p1) / 4.0;
			PCoordType p3_equivalent = (p3 + p3 + p3) / 4.0 + (p4 - p1 + p2) / 4.0;
			PCoordType p4_equivalent = (p4 + p4 + p4) / 4.0 + (p1 - p2 + p3) / 4.0;

			PCoordType alpha = (p2_equivalent - p3_equivalent); //  / 2.0 by construction but this is droped because the undefored cell in the above papers calculations has size 2;
			PCoordType beta  = (p4_equivalent - p3_equivalent); //  / 2.0;
			

			/* cout << endl << endl; */
			/* //cout << "Face Vertices of the Original Quads:------" << endl; */
			/* cout << p1[0] << ", " << p1[1] << endl; */
			/* cout << p2[0] << ", " << p2[1] << endl; */
			/* cout << p3[0] << ", " << p3[1] << endl; */
			/* cout << p4[0] << ", " << p4[1] << endl; */

			/* //cout << "Face Vertices of the Equivalent Quads:------" << endl; */
			/* cout << p1_equivalent[0] << ", " << p1_equivalent[1] << endl; */
			/* cout << p2_equivalent[0] << ", " << p2_equivalent[1] << endl; */
			/* cout << p3_equivalent[0] << ", " << p3_equivalent[1] << endl; */
			/* cout << p4_equivalent[0] << ", " << p4_equivalent[1] << endl; */

			/* cout << endl << endl; */

			Eigen::Matrix<double, 2, 2> jacobian;
			jacobian(0, 0) = alpha[0];
			jacobian(0, 1) = beta[0];
			jacobian(1, 0) = alpha[1];
			jacobian(1, 1) = beta[1];

			/* double maxElement = jacobian.array().abs().matrix().maxCoeff(); */
			/* jacobian = jacobian / maxElement; */

			jacobian = jacobian / sqrt(jacobian.determinant());

			cout  << jacobian(0, 0) << "\t" << jacobian(0, 1) << "\t" << jacobian(1, 0) << "\t" << jacobian(1, 1) << endl; 

			/*
			cout << endl << " F is ... " << endl << jacobian << endl;
			cout << endl << " F^T F is ... " << endl << jacobian.transpose() * jacobian << endl;
			cout << endl << " U is ... " << endl << stretch << endl;
			cout << endl << " U^2 is ... " << endl << stretch * stretch << endl;

			cout << endl << "sanity check FTF - U2 is ... " << endl << (jacobian.transpose() * jacobian - stretch * stretch) << endl;

			cout << "----------" << endl << endl;
			*/

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
