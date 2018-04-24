import LinearElasticitySettings
import PyAssembler
from MatrixUtils import format

class Assembler(object):
    def __init__(self, mesh, material_state):
        self.mesh = mesh;
        self.material_state = material_state;
        self.assembler = PyAssembler.FEAssembler.create(
                self.mesh.raw_mesh, self.material_state.material);
        self.__init_matrices();

    def __init_matrices(self):
        self.rigid_motion = format(self.assembler.assemble("rigid_motion"));
        self.displacement_strain = format(self.assembler.assemble("displacement_strain"));

    def update(self):
        self.stiffness = format(self.assembler.assemble("stiffness"));
        self.elasticity = format(self.assembler.assemble("elasticity_tensor"));
