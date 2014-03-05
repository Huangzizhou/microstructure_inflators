import os
import sys

linear_elasticity_path = os.environ.get("LINEAR_ELASTICITY_PATH");
if linear_elasticity_path is None:
    raise ImportError("Please set LINEAR_ELASTICITY_PATH the correct path");
if linear_elasticity_path not in sys.path:
    sys.path.append(linear_elasticity_path);
    sys.path.append(os.path.join(linear_elasticity_path, "PyUtils"));
