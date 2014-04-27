import numpy as np
import PyMeshSetting
import PyMesh

class UniquePointExtractor():
    @classmethod
    def extract(cls, pts, eps=1e-8):
        if len(pts) == 0:
            return np.array(pts);

        unique_pts = [];
        grid = PyMesh.HashGrid.create(eps);
        for i,p in enumerate(pts):
            nearby_pts = grid.get_items_near_point(p);
            if len(nearby_pts) == 0:
                grid.insert(i, p);
                unique_pts.append(p);
        return np.vstack(unique_pts);

