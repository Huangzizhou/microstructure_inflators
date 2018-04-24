import json
from timethis import timethis
from mesh_io import save_mesh

@timethis
def add_attribute_to_mesh(mesh, attr_name, attr_value):
    if mesh.has_attribute(attr_name):
        raise RuntimeError("Attribute {} already exists.".format(attr_name));
    mesh.add_attribute(attr_name);
    mesh.set_attribute(attr_name, attr_value);

@timethis
def save_mesh_fields(mesh_file, mesh, pressures=None, displacements=None,
        stress_traces=None):
    num_fields = len(displacements);
    attributes_to_save = [];
    for i in range(num_fields):
        if pressures is not None:
            pressure_attr_name = "pressure_{}".format(i);
            add_attribute_to_mesh(mesh, pressure_attr_name, pressures[i]);
            attributes_to_save.append(pressure_attr_name);

        if displacements is not None:
            disp_attr_name = "displacement_{}".format(i);
            add_attribute_to_mesh(mesh, disp_attr_name, displacements[i]);
            attributes_to_save.append(disp_attr_name);

        if stress_traces is not None:
            strs_attr_name = "stress_trace_{}".format(i);
            add_attribute_to_mesh(mesh, strs_attr_name, stress_traces[i]);
            attributes_to_save.append(strs_attr_name);
    save_mesh(mesh_file, mesh, *attributes_to_save);

@timethis
def save_parameters(parameter_file, fitter):
    parameter_dict = {
            "elasticity_tensor": fitter.elasticity_tensor.tolist(),
            "residual_error": fitter.residual_error,
            "condition_num": fitter.condition_num };
    if hasattr(fitter, "youngs_modulus"):
        parameter_dict["youngs_modulus"] = fitter.youngs_modulus.tolist();
    if hasattr(fitter, "poisson_ratio"):
        parameter_dict["poisson_ratio"] = fitter.poisson_ratio.tolist();
    if hasattr(fitter, "shear_modulus"):
        parameter_dict["shear_modulus"] = fitter.shear_modulus.tolist();

    with open(parameter_file, 'w') as fout:
        json.dump(parameter_dict, fout, indent=4);

