import PyMeshSetting
import PyMeshUtils

class Subdivision:
    def __init__(self, algorithm="simple"):
        self.__subdivision = PyMeshUtils.Subdivision.create(algorithm);

    def subdivide(self, vertices, faces, num_iterations=1):
        self.__subdivision.subdivide(vertices, faces, num_iterations);
        vertices = self.__subdivision.get_vertices();
        faces = self.__subdivision.get_faces();
        face_indices = self.__subdivision.get_face_indices();
        return vertices, faces, face_indices;
