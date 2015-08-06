import inflation, paths
import tempfile, os, subprocess, re

def homogenize(pattern, params, material, dim = 3, deg = 2, constraints=[]):
    mesh_fd, mesh_path = tempfile.mkstemp(suffix='.msh')
    os.close(mesh_fd)
    try:
        inflation.inflate(pattern, params, mesh_path, dim=dim, constraints=constraints)
        homog_out = subprocess.check_output([paths.homogenizer(dim),
            '-m', material, '-d', str(deg), mesh_path])
        young, shear, poisson, anisotropy = [], [], [], 0
        for l in homog_out.strip().split('\n'):
            m = re.match('^Approximate Young moduli:(.+)$', l)
            if (m): young = map(float, m.group(1).strip().split('\t'))

            if (dim == 3): m = re.match('^Approximate shear moduli:(.+)$', l)
            else:          m = re.match('^Approximate shear modulus:(.+)$', l)
            if (m): shear = map(float, m.group(1).strip().split('\t'))

            if (dim == 3): m = re.match('^v_yx, v_zx, v_zy:(.+)$', l)
            else:          m = re.match('^v_yx, v_xy:(.+)$', l)
            if (m): poisson = map(float, m.group(1).strip().split('\t'))

            m = re.match('^Anisotropy:(.+)$', l)
            if (m): anisotropy = float(m.group(1).strip().split('\t')[0])
        return (young, shear, poisson, anisotropy)
    finally: os.unlink(mesh_path)

def homogenize_cubic(pattern, params, material, dim = 3, deg = 2, constraints=[]):
    young, shear, poisson, anisotropy = homogenize(pattern, params, material, dim, deg, constraints)
    return (young[0], poisson[0], anisotropy)
