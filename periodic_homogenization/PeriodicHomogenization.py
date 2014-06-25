import numpy as np
import scipy.sparse
import PyMesh
from scipy.sparse.linalg import spsolve
from numpy.linalg import norm
import itertools

from MatrixUtils import *

class PeriodicHomogenization:
    ''' Implements a linear tet-based periodic homogenization.
    '''

    def __init__(self, mesh, assembler):
        self.__mesh = mesh;
        self.__assembler = assembler;
        dim = self.__mesh.dim;

    def solveCellProblems(self):
        ''' Solve the periodic homogenization cell problems:
                -div(E : (strain(w_ij) + e_ij)) = 0 in base cell
                w_ij periodic over base cell (w on identified vertices match)
                int w_ij = 0
            for each (flattened) e_ij = (1, 0, 0, 0, 0, 0), (0, 1, 0, ...) ...

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

        vtx = self.__mesh.vertices;
        # Find identified vertices (sets of vertices on the periodic boundary)
        import PeriodicBoundaryCondition as pbc
        pbcPairs = []
        bmin,bmax = self.__mesh.bbox;
        for d in range(self.__mesh.dim):
            minFaceVertices, = np.where(np.abs(vtx[:,d] - bmin[d]) < 1e-8);
            maxFaceVertices, = np.where(np.abs(vtx[:,d] - bmax[d]) < 1e-8);
            # print self.__mesh.vertices[minFaceVertices]
            # print self.__mesh.vertices[maxFaceVertices]
            assert(minFaceVertices.shape == maxFaceVertices.shape);
            pbcPairs += [pbc.match_boundaries(self.__mesh, minFaceVertices,
                maxFaceVertices)];

        identifiedPairs = list(itertools.chain(*pbcPairs))
        identifiedPairs = pbc.breakIdentificationCycles(identifiedPairs)
        
        i, j, v = [], [], []
        numPeriodicConstraints = 0;
        for (v1, v2) in identifiedPairs:
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
        # TODO: Figure out how to do only back/forward sub for subsequent solves.
        # Loop over probing strains:
        #     [1, 0, 0, 0, 0, 0], ..., [0, 0, 0, 0, 0, 0.5]
        for e_ij in np.eye(6):
            e_ij[3:] *= .5;
            # B.T * V * S * D * e_ij computes nodal load (S is the shear
            # doubling matrix).
            # Constant macroscopic strain e_ij is repeated for each tet
            # (i.e. pull out "column" ij of each tet's elasticity tensor)

            numElems = self.__mesh.num_voxels
            shearDoubler = scipy.sparse.dia_matrix(
                    (np.tile([1, 1, 1, 2, 2, 2], numElems), 0),
                    shape=(6 * numElems, 6 * numElems));

            sourceTerm = -B.T * (V * (shearDoubler * (D *
                np.tile(e_ij, numElems))));

            ij = len(w_ij)
            # # Check for force balance on periodic boundary for x loading
            # # scenario
            # if (ij == 0):
            #     for (u, v) in pbcPairs[0]:
            #         f1 = np.array(sourceTerm[3 * u : 3 * u + 3])
            #         f2 = np.array(sourceTerm[3 * v : 3 * v + 3])
            #         print "total force on periodic vertices (" + str(u) + ", " + \
            #                 str(v) + "): " + str(f1[0]) + " + " + str(f2[0]) + \
            #                 " = " + str(norm(f1 + f2))

            name = "rhs " + str(ij);
            self.__mesh.add_attribute(name);
            self.__mesh.set_attribute(name, sourceTerm);
            b = np.hstack([sourceTerm, C_rhs]);
            assert(A.shape[0] == b.shape[0]);
            soln = spsolve(A, b);

            # # Check how accurately periodic boundary conditions are enforced:
            # maxNorm = 0;
            # for (u, v) in identifiedPairs:
            #     maxNorm = max(maxNorm, abs(soln[3 * u + 0] - soln[3 * v + 0]))
            #     maxNorm = max(maxNorm, abs(soln[3 * u + 1] - soln[3 * v + 1]))
            #     maxNorm = max(maxNorm, abs(soln[3 * u + 2] - soln[3 * v + 2]))
            # print "max pbc violation for w_ij " + str(ij) + \
            #       ": " + str(maxNorm)

            # Print displacements on the x-periodic boundary for x loading
            # Also output strains
            if (ij == 0):
                # for (u, v) in pbcPairs[0]:
                #     d1 = np.array(soln[3 * u : 3 * u + 3])
                #     d2 = np.array(soln[3 * v : 3 * v + 3])
                #     print "displacement on periodic vertices (" + str(u) + ", " + \
                #             str(v) + "): " + str(d1) + ", " + str(d2)
                strain = B * soln[0:n * dim];
                for i in range(6):
                   strain_i = strain[6 * np.arange(numElems) + i]
                   name = "strain_0 component " + str(i)
                   self.__mesh.add_attribute(name)
                   self.__mesh.set_attribute(name, strain_i)

            w_ij.append(soln[0:n*dim]);

        return w_ij

    def homogenizedElasticityTensor(self, w_ij):
        '''
        Compute the homogenized elasticity tensor given solutions to the cell
        problems.
        @param[in]  w_ij    6 cell problem solutions (fluctuation displacements)
        @return     homogenized elasticity tensor
        '''
        numElems = self.__mesh.num_voxels

        B = self.__assembler.displacement_strain_matrix;
        D = self.__assembler.strain_stress_matrix;
        v = self.__mesh.get_attribute('voxel_volume').ravel();
        V = scipy.sparse.dia_matrix((np.repeat(v, 6), 0), shape=D.shape);

        # we need both shear-doubled and true e_ij to pull out symmetric
        # flattened elasticity tensor:
        # Eijkl = e_ij : E : e_kl = flat(e_ij) * S * D * flat(e_kl)
        # Since D maps true strain to true stress
        eSDouble = np.tile(np.eye(6), (numElems, 1));
        e = np.eye(6);
        e[3:, 3:] *= .5;
        e = np.tile(e, (numElems, 1))
        # Integrate base elasticity tensor over the microstructure
        elasticity = np.dot(eSDouble.T, (V * (D * e)));

        # Add in fluctuation contribution to each row of the elasticity tensor
        for row in range(6):
            elasticity[row, :] = elasticity[row, :] + \
                     np.dot(eSDouble.T, V * (D * (B * w_ij[row])));

        bmin,bmax = self.__mesh.bbox;
        cell_volume = np.prod(bmax - bmin);
        return elasticity / cell_volume;


    def periodicHomogenize(self, out_file):
        ''' Run periodic homogenization, returning the effective elasticity
            tensor.
        '''
        dim = self.__mesh.dim;
        assert(dim == 3);

        w_ij = self.solveCellProblems()
        if (out_file != None):
            mesh = self.__mesh
            writer = PyMesh.MeshWriter.create_writer(out_file)
            for ij in range(6):
                name = "w_ij " + str(ij)
                mesh.add_attribute(name)
                mesh.set_attribute(name, w_ij[ij])
                writer.with_attribute(name)
                rhs_name = "rhs " + str(ij) # Already in attributes...
                writer.with_attribute(rhs_name)
            for ij in range(6):
                name = "strain_0 component " + str(ij)
                writer.with_attribute(name)

            writer.write_mesh(mesh.raw_mesh)
        return self.homogenizedElasticityTensor(w_ij)


    def homogenizedElasticityCoefficientShapeDerivative(self, Etarget, w_ij):
        ''' Compute the normal velocity that gives a steepest descent of the
            objective:
            ||Etarget - Eh||_fro^2
        '''
        return