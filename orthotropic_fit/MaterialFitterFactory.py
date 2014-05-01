from OrthotropicMaterialFitter2D import OrthotropicMaterialFitter2D
from OrthotropicMaterialFitter3D import OrthotropicMaterialFitter3D

class MaterialFitterFactory:
    def __init__(self, mesh, material_file):
        self.mesh = mesh;
        self.material_file = material_file;

    def create(self, fitter_type="orthotropic"):
        if fitter_type == "orthotropic":
            if self.mesh.dim == 2:
                return OrthotropicMaterialFitter2D(self.mesh, self.material_file);
            elif self.mesh.dim == 3:
                return OrthotropicMaterialFitter3D(self.mesh, self.material_file);
            else:
                raise RuntimeError("Only 2D and 3D meshes are supported");
        else:
            raise NotImplementedError("Unknow fitter type: {}".format(
                fitter_type));
