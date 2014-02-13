import os.path
import sys
import numpy as np

# Update path to import PyMesh
py_mesh_path = os.environ.get("PYMESH_PATH");
if py_mesh_path == None:
    raise ImportError("Please set PYMESH_PATH to the correct lib path.");
sys.path.append(os.path.join(py_mesh_path, "lib"));
sys.path.append(os.path.join(py_mesh_path, "swig"));
import PyMesh

class Mesh(object):
    """
    A very thin wrapper around the mesh structure defined by PyMesh module.
    """
    def __init__(self, raw_mesh):
        self.__mesh = raw_mesh;

    def add_attribute(self, name):
        self.__mesh.add_attribute(name);

    def has_attribute(self, name):
        return self.__mesh.has_attribute(name);

    def get_attribute(self, name):
        return self.__mesh.get_attribute(name);

    def get_vertex_attribute(self, name):
        return self.__mesh.get_attribute(name).reshape(
                (self.num_vertices, -1), order="");

    def get_face_attribute(self, name):
        return self.__mesh.get_attribute(name).reshape(
                (self.num_faces, -1), order="C");

    def get_voxel_attribute(self, name):
        return self.__mesh.get_attribute(name).reshape(
                (self.num_voxels, -1), order="C");

    def set_attribute(self, name, val):
        self.__mesh.set_attribute(name, val);

    def get_attribute_names(self):
        return self.__mesh.get_attribute_names();

    @property
    def vertices(self):
        return self.__mesh.get_vertices().reshape(
                (-1,self.dim), order="C");

    @property
    def faces(self):
        return self.__mesh.get_faces().reshape(
                (-1, self.vertex_per_face), order="C");

    @property
    def voxels(self):
        if self.num_voxels == 0:
            return self.__mesh.get_voxels();
        else:
            return self.__mesh.get_voxels().reshape(
                    (-1, self.vertex_per_voxel), order="C");

    @property
    def num_vertices(self):
        return self.__mesh.get_num_vertices();

    @property
    def num_faces(self):
        return self.__mesh.get_num_faces();

    @property
    def num_voxels(self):
        return self.__mesh.get_num_voxels();

    @property
    def dim(self):
        return self.__mesh.get_dim();

    @property
    def vertex_per_face(self):
        return self.__mesh.get_vertex_per_face();

    @property
    def vertex_per_voxel(self):
        return self.__mesh.get_vertex_per_voxel();

    @property
    def bbox(self):
        vts = self.vertices.reshape((-1,self.dim), order='C');
        bmin = np.amin(vts, axis=0);
        bmax = np.amax(vts, axis=0);
        return bmin, bmax;

    @property
    def raw_mesh(self):
        return self.__mesh;

