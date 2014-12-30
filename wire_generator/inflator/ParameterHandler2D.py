import numpy as np

import core.PyWireInflator2DSetting
import PyWireInflator2D

class ParameterHandler2D(object):
    def __init__(self, wire_network, inflator):
        self.wire_network = wire_network;
        self.inflator = inflator;

    def convert_to_flattened_parameters(self, parameters, **kwargs):
        if (parameters.num_dofs != self.inflator.get_num_parameters()):
            raise RuntimeError(
                    "Pattern parameter is incompatible with 2D inflator\n" +
                    "Expect {} dofs, but get {} dofs".format(
                        self.inflator.get_num_parameters(),
                        parameters.num_dofs));

        raw_parameters = parameters.raw_parameters;
        thickness_map = raw_parameters.get_thickness_dof_map().ravel();
        offset_map = raw_parameters.get_offset_dof_map();

        if len(kwargs) > 0:
            source_parameters = self.__evaluate_formula(parameters, **kwargs);
        else:
            source_parameters = parameters.dofs;

        target_parameters = np.zeros(self.inflator.get_num_parameters());

        num_parameters = self.inflator.get_num_parameters();
        for i in range(num_parameters):
            param_type = self.inflator.get_parameter_type(i);
            affected_vertices = self.inflator.get_affected_vertex_orbit(i).ravel();

            if param_type == PyWireInflator2D.WireInflatorFacade.VERTEX_OFFSET:
                offset_dir = self.inflator.get_offset_direction(i);
                mask = np.amax(np.absolute(offset_dir), axis=0) > 0.0;
                dof_indices = offset_map[affected_vertices, mask];
                assert(np.all(dof_indices == dof_indices[0]));
                target_parameters[i] = source_parameters[dof_indices[0]];

            elif param_type == PyWireInflator2D.WireInflatorFacade.THICKNESS:
                dof_indices = thickness_map[affected_vertices];
                assert(np.all(dof_indices == dof_indices[0]));
                target_parameters[i] = source_parameters[dof_indices[0]];

            else:
                raise NotImplementedError(
                        "Unknown parameter type: {}".format(param_type));
        return target_parameters;

    def __evaluate_formula(self, parameters, **kwargs):
        formulas = parameters.raw_parameters.get_formulas();
        values = parameters.dofs;
        dofs = [];
        for value, formula in zip(values, formulas):
            if formula == "":
                dofs.append(value);
            else:
                dofs.append(float(kwargs[formula]));

        return np.array(dofs);

