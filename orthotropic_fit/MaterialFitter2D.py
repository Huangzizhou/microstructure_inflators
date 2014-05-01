import LinearElasticitySettings
from MaterialFitter import MaterialFitter
import PredefinedBoundaryConditions as PredefinedBC

class MaterialFitter2D(MaterialFitter):
    def __init__(self, mesh, material_file=None):
        super(MaterialFitter2D, self).__init__(mesh, material_file);

    @property
    def bc_configs(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        return [PredefinedBC.compress_x(bbox_min, bbox_max, eps),
                PredefinedBC.compress_y(bbox_min, bbox_max, eps),
                PredefinedBC.compress_xy(bbox_min, bbox_max, eps)];
