import os
import sys

# Update path to import PyMesh
py_mesh_path = os.environ.get("PYMESH_PATH");
if py_mesh_path == None:
    raise ImportError("Please set PYMESH_PATH to the correct lib path.");
sys.path.append(os.path.join(py_mesh_path, "lib"));
sys.path.append(os.path.join(py_mesh_path, "swig"));

