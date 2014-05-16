import numpy as np

def compress_x(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 0, 0, 1, eps),
                generate_force(bbox_min, bbox_max, 0, 0,-1, eps),
                ]
            };
    return bc_config;

def compress_y(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 1, 1, 1, eps),
                generate_force(bbox_min, bbox_max, 1, 1,-1, eps),
                ]
            };
    return bc_config;

def compress_z(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 2, 2, 1, eps),
                generate_force(bbox_min, bbox_max, 2, 2,-1, eps),
                ]
            };
    return bc_config;

def compress_xy(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 0, 0, 1, eps),
                generate_force(bbox_min, bbox_max, 0, 0,-1, eps),
                generate_force(bbox_min, bbox_max, 1, 1, 1, eps),
                generate_force(bbox_min, bbox_max, 1, 1,-1, eps),
                ]
            };
    return bc_config;

def compress_yz(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 1, 1, 1, eps),
                generate_force(bbox_min, bbox_max, 1, 1,-1, eps),
                generate_force(bbox_min, bbox_max, 2, 2, 1, eps),
                generate_force(bbox_min, bbox_max, 2, 2,-1, eps),
                ]
            };
    return bc_config;

def compress_zx(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 2, 2, 1, eps),
                generate_force(bbox_min, bbox_max, 2, 2,-1, eps),
                generate_force(bbox_min, bbox_max, 0, 0, 1, eps),
                generate_force(bbox_min, bbox_max, 0, 0,-1, eps),
                ]
            };
    return bc_config;

def shear_xy(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 0, 1, 1, eps),
                generate_force(bbox_min, bbox_max, 0, 1,-1, eps),
                generate_force(bbox_min, bbox_max, 1, 0, 1, eps),
                generate_force(bbox_min, bbox_max, 1, 0,-1, eps),
                ]
            };
    return bc_config;

def shear_yz(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 1, 2, 1, eps),
                generate_force(bbox_min, bbox_max, 1, 2,-1, eps),
                generate_force(bbox_min, bbox_max, 2, 1, 1, eps),
                generate_force(bbox_min, bbox_max, 2, 1,-1, eps),
                ]
            };
    return bc_config;

def shear_zx(bbox_min, bbox_max, eps=1e-3):
    bc_config = {
            "neumann": [
                generate_force(bbox_min, bbox_max, 2, 0, 1, eps),
                generate_force(bbox_min, bbox_max, 2, 0,-1, eps),
                generate_force(bbox_min, bbox_max, 0, 2, 1, eps),
                generate_force(bbox_min, bbox_max, 0, 2,-1, eps),
                ]
            };
    return bc_config;

def generate_force(bbox_min, bbox_max, normal_idx, force_idx, sign, eps):
    bbox = np.vstack((bbox_min, bbox_max));
    if sign > 0:
        bound = bbox_max;
    else:
        bound = bbox_min;
    bbox[0, normal_idx] = bound[normal_idx] - eps;
    bbox[1, normal_idx] = bound[normal_idx] + eps;
    force = [0.0, 0.0, 0.0];
    force[force_idx] = -sign;
    spec = { "box": bbox.ravel(order="F").tolist(),
            "force": force };
    return spec;

