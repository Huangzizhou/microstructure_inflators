import numpy as np
from datetime import datetime
import os
import os.path
from scipy.spatial import Delaunay, ConvexHull
from subprocess import check_call

import LinearElasticitySettings
from mesh_io import form_mesh, save_mesh, load_mesh
import PyMeshUtils
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

def generate_3D_box_mesh(bbox_min, bbox_max, num_samples, keep_symmetry=False):
    step_size = np.divide((bbox_max - bbox_min),
            [num_samples, num_samples, num_samples]);

    num_vertices = 0;
    vertices = [];
    tets = [];
    for i in range(num_samples):
        for j in range(num_samples):
            for k in range(num_samples):
                p = np.multiply([i,j,k], step_size);
                corners = [
                        [p[0]             , p[1]             , p[2]             ],
                        [p[0]+step_size[0], p[1]             , p[2]             ],
                        [p[0]+step_size[0], p[1]+step_size[1], p[2]             ],
                        [p[0]             , p[1]+step_size[1], p[2]             ],
                        [p[0]             , p[1]             , p[2]+step_size[2]],
                        [p[0]+step_size[0], p[1]             , p[2]+step_size[2]],
                        [p[0]+step_size[0], p[1]+step_size[1], p[2]+step_size[2]],
                        [p[0]             , p[1]+step_size[1], p[2]+step_size[2]],
                        ];
                if keep_symmetry:
                    cell_vertices, cell_tets =\
                            split_hex_into_tets_symmetrically(corners);
                else:
                    cell_vertices, cell_tets =\
                            split_hex_into_tets(corners);
                vertices.append(cell_vertices);
                tets.append(cell_tets + num_vertices);
                num_vertices += len(cell_vertices);

    vertices = np.vstack(vertices);
    tets = np.vstack(tets);

    vertices, tets = remove_duplicated_vertices(vertices, tets);
    vertices, tets = remove_isolated_vertices(vertices, tets);

    faces = np.array([], dtype=int);
    mesh = form_mesh(vertices, faces, tets);
    return mesh;

def split_hex_into_tets(corners):
    """ Convert a hex to 24 tets.
    Algorithm: Each hex is split into 6 tets.  The resulting tet mesh does not
    symmetric connectivities.

    Corner ordering:
         7 _______ 6
          /:     /|        z
       4 /______/ |        |
        |  :   5| |        |  y
        |  :... |.|        | /
        | . 3   | /2       |/
        |_______|/         /-------x
        0        1

    """
    vertices = np.array(corners);
    tets = np.array([
        [0, 3, 7, 6],
        [0, 3, 6, 2],
        [0, 2, 6, 1],
        [5, 0, 6, 1],
        [5, 0, 4, 6],
        [6, 0, 4, 7],
        ]);
    return vertices, tets;

def split_hex_into_tets_symmetrically(corners):
    """ Convert a hex to 24 tets.
    Algorithm: Form 6 pyramids by using cell centroid and each face.  Break each
    pyramid into 4 tets by spliting the base (the original cell face) into 4 triangles.

    Corner ordering:
         7 _______ 6
          /:     /|        z
       4 /______/ |        |
        |  :   5| |        |  y
        |  :... |.|        | /
        | . 3   | /2       |/
        |_______|/         /-------x
        0        1

    """
    assert(len(corners) == 8);
    corners = np.array(corners);
    centroid = np.mean(corners, axis=0);
    faces = np.array([
        [0, 3, 2, 1], # bottom
        [4, 5, 6, 7], # top
        [1, 2, 6, 5], # right
        [0, 4, 7, 3], # left
        [0, 1, 5, 4], # front
        [7, 6, 2, 3], # back
        ]);

    face_centers = np.array([np.mean(corners[face], axis=0)
        for face in faces]);

    vertices = np.vstack([corners, face_centers, centroid]);
    tets = np.array([
        # tets based on bottom face
        [0, 1, 8, 14],
        [1, 2, 8, 14],
        [2, 3, 8, 14],
        [3, 0, 8, 14],
        # tets based on top faces
        [5, 4, 9, 14],
        [4, 7, 9, 14],
        [7, 6, 9, 14],
        [6, 5, 9, 14],
        # tets based on right face
        [2, 1, 10, 14],
        [1, 5, 10, 14],
        [5, 6, 10, 14],
        [6, 2, 10, 14],
        # tets based on left face
        [0, 3, 11, 14],
        [3, 7, 11, 14],
        [7, 4, 11, 14],
        [4, 0, 11, 14],
        # tets based on front face
        [0, 4, 12, 14],
        [4, 5, 12, 14],
        [5, 1, 12, 14],
        [1, 0, 12, 14],
        # tets based on back face
        [3, 2, 13, 14],
        [2, 6, 13, 14],
        [6, 7, 13, 14],
        [7, 3, 13, 14],
        ]);

    return vertices, tets;

def remove_duplicated_vertices(vertices, tets):
    remover = PyMeshUtils.DuplicatedVertexRemoval(vertices, tets);
    remover.run(1e-6);
    return remover.get_vertices(), remover.get_faces();

def remove_isolated_vertices(vertices, tets):
    remover = PyMeshUtils.IsolatedVertexRemoval(vertices, tets);
    remover.run();
    return remover.get_vertices(), remover.get_faces();


