import LinearElasticitySettings

import numpy as np

import scipy.sparse
from scipy.sparse.linalg import spsolve

class Dual:
    def __init__(self, mesh, assembler):
        self.mesh = mesh;
        self.assembler = assembler;
        self.source_term = np.zeros(self.mesh.num_vertices * self.mesh.dim);

    def solve(self):
        num_dof = self.mesh.num_vertices * self.mesh.dim;
        K = self.assembler.stiffness;
        R = self.assembler.rigid_motion;
        A = scipy.sparse.bmat([
            [K, R.T], [R, None] ]).tocsc();
        rhs = np.hstack((self.source_term, self.Rb));
        sol = spsolve(A, rhs);

        self.displacement = sol[:num_dof];

    def add_neumann_bc(self, bc):
        dim = self.mesh.dim;
        for idx, val, area in zip(*bc[:3]):
            self.source_term[dim*idx:dim*(idx+1)] += \
                    np.array(val) * area;

    def set_rigid_motion_as(self, target_displacement):
        R = self.assembler.rigid_motion;
        self.Rb = R * target_displacement;

    @property
    def stress(self):
        B = self.assembler.displacement_strain;
        E = self.assembler.elasticity;
        return E * B * self.displacement;
