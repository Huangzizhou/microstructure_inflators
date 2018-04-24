import os
import sys

microstructure_path = os.environ.get("MICROSTRUCTURES_PATH");
if microstructure_path is None:
    raise ImportError("Please set MICROSTRUCTURES_PATH to the correct lib path");
sys.path.append(os.path.join(microstructure_path, "Luigi/wireinflator2D/lib"));
sys.path.append(os.path.join(microstructure_path, "Luigi/wireinflator2D/swig"));

