import PyTetgen

def tetgen(flags, vertices, faces):
    tetgen = PyTetgen.TetgenWrapper(vertices, faces);
    tetgen.run(flags);

    tet_vertices = tetgen.get_vertices();
    tet_faces = tetgen.get_faces();
    tet_voxels = tetgen.get_voxels();

    return tet_vertices, tet_faces, tet_voxels;
