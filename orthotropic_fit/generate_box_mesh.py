import numpy as np
from mesh_io import form_mesh
from timethis import timethis

@timethis
def generate_box_mesh(box_min, box_max, num_samples):
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

