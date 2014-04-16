import numpy as np
from datetime import datetime
import os
import os.path
from subprocess import check_call

import LinearElasticitySettings
from mesh_io import form_mesh, save_mesh, load_mesh
from scipy.spatial import Delaunay, ConvexHull
from timethis import timethis

@timethis
def generate_box_mesh(box_min, box_max, num_samples):
    dim = len(box_min);
    if dim == 2:
        return generate_2D_box_mesh(box_min, box_max, num_samples);
    elif dim == 3:
        return generate_3D_box_mesh(box_min, box_max, num_samples);

def generate_2D_box_mesh(box_min, box_max, num_samples):
    vertices = [];
    faces = [];
    for i in range(num_samples+1):
        for j in range(num_samples+1):
            x_ratio = float(i) / float(num_samples);
            y_ratio = float(j) / float(num_samples);
            x = x_ratio * box_max[0] + (1.0 - x_ratio) * box_min[0];
            y = y_ratio * box_max[1] + (1.0 - y_ratio) * box_min[1];
            vertices.append([x,y]);

    row_size = num_samples+1;
    for i in range(num_samples):
        for j in range(num_samples):
            idx = [i*row_size+j,
                    i*row_size+j+1,
                    (i+1)*row_size+j,
                    (i+1)*row_size+j+1 ];
            faces.append([idx[0], idx[2], idx[3]]);
            faces.append([idx[0], idx[3], idx[1]]);

    vertices = np.array(vertices, dtype=float);
    faces = np.array(faces, dtype=int);
    mesh = form_mesh(vertices, faces);
    return mesh;

def reorientate_triangles(vertices, faces):
    """ This only works for convex shapes
    """
    centroid = np.mean(vertices, axis=0);
    face_centers = np.mean(vertices[faces], axis=1);
    out_dir = face_centers - centroid;
    edge_1= vertices[faces[:,1]] - vertices[faces[:,0]]
    edge_2= vertices[faces[:,2]] - vertices[faces[:,0]]
    normals = np.cross(edge_1, edge_2);

    pointing_inwards = np.sum(normals * out_dir, axis=1) < 0.0;
    faces[pointing_inwards, :] = faces[pointing_inwards][:,[0,2,1]];
    return faces;

def reorientate_tets(vertices, tets):
    edge_1 = vertices[tets[:,1]] - vertices[tets[:,0]];
    edge_2 = vertices[tets[:,2]] - vertices[tets[:,0]];
    edge_3 = vertices[tets[:,3]] - vertices[tets[:,0]];
    volume = np.sum(np.cross(edge_1, edge_2) * edge_3, axis=1);
    bad = np.absolute(volume) < 1e-6;
    print(vertices[tets[bad]]);
    inverted = volume < 0.0;
    tets[inverted] = tets[inverted][:, [0,1,3,2]];
    tets = tets[np.logical_not(bad)];
    return tets;

def generate_3D_box_mesh_old_4(box_min, box_max, num_samples):
    tmp_dir = "/tmp/"
    x_coord = np.linspace(box_min[0], box_max[0], num_samples+1);
    y_coord = np.linspace(box_min[1], box_max[1], num_samples+1);
    z_coord = np.linspace(box_min[2], box_max[2], num_samples+1);

    coord = np.meshgrid(x_coord, y_coord, z_coord);
    X = coord[0].ravel();
    Y = coord[1].ravel();
    Z = coord[2].ravel();
    vertices = np.array([X, Y, Z]).T;

    hull = ConvexHull(vertices);
    faces = hull.simplices;

    stamp = datetime.isoformat(datetime.now());
    tmp_obj = os.path.join(tmp_dir, stamp+".obj");
    tmp_msh = os.path.join(tmp_dir, stamp+".msh");

    mesh = form_mesh(vertices, np.array([]));
    save_mesh(tmp_obj, mesh);

    command = "tetgen.py {} {}".format(tmp_obj, tmp_msh);
    check_call(command.split());

    mesh = load_mesh(tmp_msh);
    os.remove(tmp_obj);
    os.remove(tmp_msh);
    return mesh;

def generate_3D_box_mesh_old(box_min, box_max, num_samples):
    x_coord = np.linspace(box_min[0], box_max[0], num_samples+1);
    y_coord = np.linspace(box_min[1], box_max[1], num_samples+1);
    z_coord = np.linspace(box_min[2], box_max[2], num_samples+1);

    coord = np.meshgrid(x_coord, y_coord, z_coord);
    X = coord[0].ravel();
    Y = coord[1].ravel();
    Z = coord[2].ravel();
    vertices = np.array([X, Y, Z]).T;

    faces = np.zeros(0);

    generator = Delaunay(vertices);
    tets = generator.simplices.copy();
    tets = reorientate_tets(vertices, tets);

    mesh = form_mesh(vertices, faces, tets);
    return mesh;

def generate_3D_box_mesh_b(box_min, box_max, num_samples):
    tmp_dir = "/tmp/"
    x_coord = np.linspace(box_min[0], box_max[0], num_samples+1);
    y_coord = np.linspace(box_min[1], box_max[1], num_samples+1);
    z_coord = np.linspace(box_min[2], box_max[2], num_samples+1);

    vertices = np.array([[x,y,z]
        for x in x_coord
        for y in y_coord
        for z in z_coord]);

    tets = [];
    faces = [];
    step_size = [(num_samples+1)**2, num_samples+1, 1];
    for i in range(num_samples):
        for j in range(num_samples):
            for k in range(num_samples):
                idx_0 = np.dot([i  ,j  ,k+1], step_size);
                idx_1 = np.dot([i  ,j  ,k  ], step_size);
                idx_2 = np.dot([i+1,j  ,k  ], step_size);
                idx_3 = np.dot([i+1,j  ,k+1], step_size);
                idx_4 = np.dot([i  ,j+1,k+1], step_size);
                idx_5 = np.dot([i  ,j+1,k  ], step_size);
                idx_6 = np.dot([i+1,j+1,k  ], step_size);
                idx_7 = np.dot([i+1,j+1,k+1], step_size);

                if i == 0:
                    faces.append([idx_1, idx_4, idx_5]);
                    faces.append([idx_1, idx_0, idx_4]);
                if i == num_samples-1:
                    faces.append([idx_2, idx_7, idx_3]);
                    faces.append([idx_2, idx_6, idx_7]);

                if j == 0:
                    faces.append([idx_0, idx_1, idx_2]);
                    faces.append([idx_0, idx_2, idx_3]);
                if j == num_samples-1:
                    faces.append([idx_4, idx_6, idx_5]);
                    faces.append([idx_4, idx_7, idx_6]);

                if k == 0:
                    faces.append([idx_1, idx_6, idx_2]);
                    faces.append([idx_1, idx_5, idx_6]);
                if k == num_samples-1:
                    faces.append([idx_0, idx_3, idx_7]);
                    faces.append([idx_0, idx_7, idx_4]);

    faces = np.array(faces, dtype=int);

    stamp = datetime.isoformat(datetime.now());
    tmp_obj = os.path.join(tmp_dir, stamp+".obj");
    tmp_msh = os.path.join(tmp_dir, stamp+".msh");

    mesh = form_mesh(vertices, faces);
    save_mesh(tmp_obj, mesh);

    command = "tetgen.py {} {}".format(tmp_obj, tmp_msh);
    check_call(command.split());

    mesh = load_mesh(tmp_msh);
    os.remove(tmp_obj);
    os.remove(tmp_msh);
    print(mesh.voxels);
    return mesh;


def generate_3D_box_mesh(box_min, box_max, num_samples):
    x_coord = np.linspace(box_min[0], box_max[0], num_samples+1);
    y_coord = np.linspace(box_min[1], box_max[1], num_samples+1);
    z_coord = np.linspace(box_min[2], box_max[2], num_samples+1);

    vertices = np.array([[x,y,z]
        for x in x_coord
        for y in y_coord
        for z in z_coord]);

    tets = [];
    faces = [];
    step_size = [(num_samples+1)**2, num_samples+1, 1];
    for i in range(num_samples):
        for j in range(num_samples):
            for k in range(num_samples):
                idx_0 = np.dot([i  ,j  ,k  ], step_size);
                idx_1 = np.dot([i  ,j  ,k+1], step_size);
                idx_2 = np.dot([i  ,j+1,k  ], step_size);
                idx_3 = np.dot([i  ,j+1,k+1], step_size);
                idx_4 = np.dot([i+1,j  ,k  ], step_size);
                idx_5 = np.dot([i+1,j  ,k+1], step_size);
                idx_6 = np.dot([i+1,j+1,k  ], step_size);
                idx_7 = np.dot([i+1,j+1,k+1], step_size);

                tets.append([idx_0, idx_2, idx_6, idx_7]);
                tets.append([idx_0, idx_2, idx_7, idx_3]);
                tets.append([idx_0, idx_3, idx_7, idx_1]);
                tets.append([idx_5, idx_0, idx_7, idx_1]);
                tets.append([idx_5, idx_0, idx_4, idx_7]);
                tets.append([idx_7, idx_0, idx_4, idx_6]);

                #if i == 0:
                #    faces.append([idx_1, idx_4, idx_5]);
                #    faces.append([idx_1, idx_0, idx_4]);
                #elif i == num_samples-1:
                #    faces.append([idx_2, idx_7, idx_3]);
                #    faces.append([idx_2, idx_6, idx_7]);

                #if j == 0:
                #    faces.append([idx_0, idx_1, idx_2]);
                #    faces.append([idx_0, idx_2, idx_3]);
                #elif j == num_samples-1:
                #    faces.append([idx_4, idx_6, idx_5]);
                #    faces.append([idx_4, idx_7, idx_6]);

                #if k == 0:
                #    faces.append([idx_1, idx_6, idx_2]);
                #    faces.append([idx_1, idx_5, idx_6]);
                #elif k == num_samples-1:
                #    faces.append([idx_0, idx_3, idx_7]);
                #    faces.append([idx_0, idx_7, idx_4]);

    faces = np.array(faces, dtype=int);
    tets = np.array(tets, dtype=int);
    mesh = form_mesh(vertices, faces, tets);
    return mesh;


