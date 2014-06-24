import numpy as np

class MaterialOptimizer(object):
    def __init__(self, mesh):
        self.mesh = mesh;
        num_dof = self.mesh.num_vertices * self.mesh.dim;
        self.source_term = np.zeros(num_dof);
        self.target_displacement = np.zeros(num_dof);
        self.bd_index_map = np.ones(self.mesh.num_vertices, dtype=int) * -1;
        self.bd_index_map[self.mesh.boundary_nodes] = np.arange(self.mesh.num_boundary_nodes);

    def add_neumann_bc(self, neumann_bc):
        self.applied_idx, self.applied_traction, self.applied_node_area\
                = neumann_bc;
        self.applied_traction = np.array(self.applied_traction);
        assert(self.mesh.dim == self.applied_traction.shape[1]);

        neumann_term = np.zeros((self.mesh.num_vertices, self.mesh.dim));
        neumann_term[self.applied_idx] = self.applied_traction;
        self.source_term += neumann_term.ravel(order="C");

    def add_dirichlet_bc(self, dirichlet_bc):
        self.fixed_idx, self.fixed_position, self.fixed_node_area\
                = dirichlet_bc;

        dirichlet_term = np.zeros((self.mesh.num_vertices, self.mesh.dim));
        dirichlet_term[self.fixed_idx] = self.fixed_position;
        self.target_displacement += dirichlet_term.ravel(order="C");
        self.target_areas = np.zeros(self.mesh.num_vertices);
        self.target_areas[self.fixed_idx] = self.fixed_node_area;

    def optimize(self, max_iterations):
        raise NotImplementedError("This method is abstract");

