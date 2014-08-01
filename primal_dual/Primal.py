import LinearElasticitySettings

import numpy as np

import scipy.sparse
from scipy.sparse.linalg import spsolve

class Primal(object):
    def __init__(self, mesh, assembler):
        self.mesh = mesh;
        self.assembler = assembler;
        self.source_term = np.zeros(self.mesh.num_vertices * self.mesh.dim);

    def solve(self):
        num_dof = self.mesh.num_vertices * self.mesh.dim;
        K = self.assembler.stiffness;
        C = self.C;
        A = scipy.sparse.bmat([
            [K, C.T], [C, None] ]).tocsc();
        rhs = np.vstack((self.source_term.reshape((-1, 1)), self.Cb));
        sol = spsolve(A, rhs);

        self.displacement = sol[:num_dof];

    def add_dirichlet_bc(self, bc):
        dim = self.mesh.dim;
        rows = [];
        cols = [];
        data = [];
        rhs = [];
        count = 0;
        for idx, pos, area in zip(*bc):
            col_idx = np.arange(dim);
            if np.any(pos.mask):
                col_idx = col_idx[~pos.mask];
            for i in col_idx:
                rows.append(count);
                cols.append(idx*dim + i);
                data.append(1.0);
                rhs.append(pos[i]);
                count+=1;

        self.C = scipy.sparse.coo_matrix((data, (rows, cols)),
                shape=(count, dim*self.mesh.num_vertices)).tocsc();
        self.Cb = np.reshape(rhs, (count, 1));

    def add_neumann_bc(self, bc):
        dim = self.mesh.dim;
        for idx, val, area in zip(*bc[:3]):
            self.source_term[dim*idx:dim*(idx+1)] += \
                    np.array(val) * area;

    @property
    def strain(self):
        B = self.assembler.displacement_strain;
        return B * self.displacement;
