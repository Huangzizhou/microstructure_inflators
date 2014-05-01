from numpy.linalg import norm

import LinearElasticitySettings
from MaterialFitter import MaterialFitter
import PredefinedBoundaryConditions as PredefinedBC

class MaterialFitter3D(MaterialFitter):
    def __init__(self, mesh, material_file=None):
        super(MaterialFitter3D, self).__init__(mesh, material_file);

    @property
    def bc_configs(self):
        bbox_min, bbox_max = self.mesh.bbox;
        eps = min(2e-1, norm(bbox_max - bbox_min) * 0.01);
        return [PredefinedBC.compress_x(bbox_min, bbox_max, eps),
                PredefinedBC.compress_y(bbox_min, bbox_max, eps),
                PredefinedBC.compress_z(bbox_min, bbox_max, eps),
                PredefinedBC.compress_xy(bbox_min, bbox_max, eps),
                PredefinedBC.compress_yz(bbox_min, bbox_max, eps),
                PredefinedBC.compress_zx(bbox_min, bbox_max, eps) ];

