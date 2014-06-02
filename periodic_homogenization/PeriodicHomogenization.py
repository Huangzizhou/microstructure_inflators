import numpy as np
import scipy.sparse
from scipy.sparse.linalg import spsolve

from MatrixUtils import *

class PeriodicHomogenization:
    ''' Implements a linear tet-based periodic homogenization.
    '''

    def __init__(self, mesh, assembler):
        self.__mesh = mesh;
        self.__assembler = assembler;
        dim = self.__mesh.dim;

    def solveCellProblems(self, identifiedVertices):
        ''' Solve the periodic homogenization cell problems:
                -div(E : (strain(w_ij) + e_ij)) = 0 in base cell
                w_ij periodic over base cell (w on identified vertices match)
                int w_ij = 0
            for each (flattened) e_ij = (1, 0, 0, 0, 0, 0), (0, 1, 0, ...) ...

            @param[in] identifiedVertices   list of pairs of boundary vertices
                                            to identify.
            @return    [w_00, w_11, w_22, w_12, w_02, w_01]
        '''
        dim = self.__mesh.dim;
        assert(dim == 3);
        n = self.__mesh.num_vertices;

        K = self.__assembler.stiffness_matrix;
        assert(K.shape[0] == dim * n);
        assert(K.shape[1] == dim * n);

        C = self.__assembler.rigid_motion_matrix;
        # Only use the translation part of the no rigid motion constraint.
        C = scipy.sparse.vstack([C.getrow(0), C.getrow(1), C.getrow(2)]);
        
        i, j, v = [], [], []
        numPeriodicConstraints = 0;
        for (v1, v2) in identifiedVertices:
            for d in range(dim):
                i.append(numPeriodicConstraints);
                j.append(dim * v1 + d);
                v.append(1.0);
                i.append(numPeriodicConstraints);
                j.append(dim * v2 + d);
                v.append(-1.0);

                numPeriodicConstraints = numPeriodicConstraints + 1;
        PC = scipy.sparse.coo_matrix((v, (i, j)),
                                     shape=(numPeriodicConstraints, n * dim));
        C = scipy.sparse.vstack([C, PC]);
        C_rhs = np.zeros(3 + numPeriodicConstraints);

        assert(C.shape[0] == C_rhs.shape[0]);

        # Assemble lagrange multiplier system.
        nc = C.shape[0]
        A = stack([[K, C], [C.T, Z(nc, nc)]]);
        # save_matrix(A, 'A');

        B = self.__assembler.displacement_strain_matrix;
        D = self.__assembler.strain_stress_matrix;
        v = self.__mesh.get_attribute('voxel_volume').ravel();
        V = scipy.sparse.dia_matrix((np.repeat(v, 6), 0), shape=D.shape);

        # save_matrix(D, 'D');

        w_ij = []
        # Figure out how to do only back/forward sub for subsequent solves.
        for e_ij in np.eye(6):
            # Double shear components so B.T * V * D * e_ij computes nodal load
            e_ij = np.multiply(e_ij, np.array([1, 1, 1, 2, 2, 2]));
            # Constant macroscopic strain e_ij is repeated for each tet
            # (i.e. pull out "column" ij of each tet's elasticity tensor)
            source_term = -B.T * (V * (D *
                np.tile(e_ij, self.__mesh.num_voxels)));
            b = np.hstack([source_term, C_rhs]);
            assert(A.shape[0] == b.shape[0]);
            soln = spsolve(A, b);
            w_ij.append(soln[0:n*dim]);

        return w_ij

    def periodicHomogenize(self):
        ''' Run periodic homogenization, returning the effective elasticity
            tensor.
        '''
        dim = self.__mesh.dim;
        assert(dim == 3);
        import PeriodicBoundaryCondition as pbc

        numElems = self.__mesh.num_voxels

        # Get pairs of vertices on opposite periodic cell boundary faces
        bmin,bmax = self.__mesh.bbox;
        vtx = self.__mesh.vertices;
        identifiedPairs = []
        for d in range(self.__mesh.dim):
            minFaceVertices, = np.where(np.abs(vtx[:,d] - bmin[d]) < 1e-5);
            maxFaceVertices, = np.where(np.abs(vtx[:,d] - bmax[d]) < 1e-5);
            # print self.__mesh.vertices[minFaceVertices]
            # print self.__mesh.vertices[maxFaceVertices]
            assert(minFaceVertices.shape == maxFaceVertices.shape);
            identifiedPairs = identifiedPairs + pbc.match_boundaries(
                    self.__mesh, minFaceVertices, maxFaceVertices);

        # print identifiedPairs;
        identifiedPairs = pbc.breakIdentificationCycles(identifiedPairs)
        # print identifiedPairs;
        # print len(identifiedPairs);
        w_ijs = self.solveCellProblems(identifiedPairs);
        assert(len(w_ijs) == 6);

        B = self.__assembler.displacement_strain_matrix;
        D = self.__assembler.strain_stress_matrix;
        v = self.__mesh.get_attribute('voxel_volume').ravel();
        V = scipy.sparse.dia_matrix((np.repeat(v, 6), 0), shape=D.shape);

        # shearDoubler = scipy.sparse.dia_matrix(
        #         (np.tile([1, 1, 1, 2, 2, 2], numElems), 0),
        #         shape=(6 * numElems, 6 * numElems));

        e = np.tile(np.eye(6), (numElems, 1));
        # Integrate base elasticity tensor over the microstructure
        elasticity = np.dot(e.T, (V * (D * e)));

        # Add in fluctuation contribution to each row of the elasticity tensor
        for row in range(6):
            elasticity[row, :] = elasticity[row, :] + \
                     np.dot(e.T, V * (D * (B * w_ijs[row])));

        cell_volume = np.prod(bmax - bmin);
        return elasticity / cell_volume;
