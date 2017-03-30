import os

def inflator(dim = 3):
#    if (dim == 3): return os.environ['MICRO_DIR'] + "/pattern_optimization/Inflator_cli"
#    return os.environ['MICRO_DIR'] + "/pattern_optimization/Inflator2D_cli"
    return os.environ['MICRO_DIR'] + "/isosurface_inflator/isosurface_cli"

def homogenizer(dim = 3):
    return os.environ['MeshFEM'] + "/PeriodicHomogenization_cli"

def pattern(pat, dim = 3):
    if (dim == 3):
        if (pat == 'truncated_octahedron'): return os.environ['MICRO_DIR'] + "/patterns/3D/truncated_octahedron.wire"
        return os.environ['MICRO_DIR'] + ("/patterns/3D/reference_wires/pattern%04i.wire" % pat)
    if (dim == 2):
        if (pat > 3000):
            return os.environ['MICRO_DIR'] + ("/patterns/2D/topologies/twice/%04i.obj" % pat)
        if (pat > 1000):
            return os.environ['MICRO_DIR'] + ("/patterns/2D/topologies/splitted/%04i.obj" % pat)
        if (pat > 0):
            return os.environ['MICRO_DIR'] + ("/patterns/2D/topologies/%04i.obj" % pat)
        if (pat == -1): 
            return os.environ['MICRO_DIR'] + "/pattern_optimization/experiments/horizontal_cover/octa_cell_splitted.obj"
        if (pat == -2): 
            return os.environ['MICRO_DIR'] + "/pattern_optimization/experiments/horizontal_cover/octa_cell_double_splitted.obj"
        return os.environ['MICRO_DIR'] + "/Luigi/wireinflator2D/meshes/octa_cell.obj"

def material(name):
    return os.environ['MICRO_DIR'] + ("/materials/%s.material" % name)
        
def optimizer(dim = 3):
    return os.environ['MICRO_DIR'] + "/worst_case_stress/WCSOptimization_cli"

def jobGenerator(dim = 3):
    return os.environ['MICRO_DIR'] + "/pattern_optimization/GenIsosurfaceJob"
