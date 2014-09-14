#pragma once
#include <iostream>
#include <memory>
#include <sstream>

#include "InflatorParameters.h"
#include "WireInflator2D.h"

#include "EigenTypedef.h"
#include "Exception.h"

class WireInflatorFacade {
    public:
        enum ParameterType {
            THICKNESS,
            VERTEX_OFFSET
        };

    public:
        WireInflatorFacade(std::string wire_file) : m_inflator(wire_file) {
            m_t_params.max_area = 0.001;
        }

        void set_dimension(size_t rows, size_t cols) {
            m_rows = rows;
            m_cols = cols;
            m_p_params = ParameterGridPtr(new ParameterGrid(m_cols, m_rows));
            for (size_t row = 0; row < m_rows; row++) {
                for (size_t col=0; col < m_cols; col++) {
                    (*m_p_params)(col, row) = NULL;
                }
            }
        }

        size_t get_num_parameters() {
            const WireInflator2D::PatternGen& pattern_gen = m_inflator.patternGenerator();
            return pattern_gen.numberOfParameters();
        }

        ParameterType get_parameter_type(size_t i) {
            const ParameterOperation& param_op = m_inflator.patternGenerator().getParameterOperations()[i];
            switch (param_op.type) {
                case ParameterOperation::Radius:
                    return THICKNESS;
                case ParameterOperation::Translation:
                    return VERTEX_OFFSET;
                default:
                    throw NotImplementedError("Unknown parameter type detected.");
            }
        }

        VectorI get_affected_vertex_orbit(size_t i) {
            VectorI vertex_orbit;
            const ParameterOperation& param_op = m_inflator.patternGenerator().getParameterOperations()[i];
            size_t count = 0;
            switch (param_op.type) {
                case ParameterOperation::Radius:
                    vertex_orbit.resize(param_op.nodes.size());
                    std::copy(param_op.nodes.begin(), param_op.nodes.end(), vertex_orbit.data());
                    break;
                case ParameterOperation::Translation:
                    vertex_orbit.resize(param_op.nodes_displ.size());
                    for (auto itr : param_op.nodes_displ) {
                        vertex_orbit[count] = itr.first;
                        count++;
                    }
                    break;
                default:
                    throw NotImplementedError("Unknown parameter type detected.");
            }
            assert(vertex_orbit.size() > 0);
            return vertex_orbit;
        }

        MatrixF get_offset_direction(size_t i) {
            const ParameterOperation& param_op = m_inflator.patternGenerator().getParameterOperations()[i];
            if (param_op.type == ParameterOperation::Translation) {
                const size_t num_nodes = param_op.nodes_displ.size();
                MatrixF offset_dir(num_nodes, 2);
                size_t count = 0;
                for (auto itr : param_op.nodes_displ) {
                    offset_dir.coeffRef(count, 0) = itr.second[0];
                    offset_dir.coeffRef(count, 1) = itr.second[1];
                    count++;
                }
                return offset_dir;
            } else {
                std::stringstream err_msg;
                err_msg << "Parameter " << i << " is not offset parameter";
                throw RuntimeError(err_msg.str());
            }
        }

        void set_parameter(size_t row, size_t col, const VectorF& param) {
            assert(row < m_rows && col < m_cols);
            if (row >= m_rows) {
                std::stringstream err_msg;
                err_msg << "Out of bound: " << row << " >= " << m_rows;
                throw RuntimeError(err_msg.str());
            }
            if (col >= m_cols) {
                std::stringstream err_msg;
                err_msg << "Out of bound: " << col << " >= " << m_cols;
                throw RuntimeError(err_msg.str());
            }
            const WireInflator2D::PatternGen& pattern_gen = m_inflator.patternGenerator();
            CellParameters* p = new CellParameters(pattern_gen.numberOfParameters());

            const size_t num_param = p->numberOfParameters();
            for (size_t i=0; i<num_param; i++) {
                std::pair<double, double> range = pattern_gen.getParameterRange(i);
                if (param[i] >= range.first && param[i] <= range.second) {
                    p->parameter(i) = param[i];
                } else {
                    std::stringstream err_msg;
                    err_msg << "param " << i << "=" << param[i]
                        << " is out of range ["
                        << range.first << ", " << range.second << "]";
                    throw RuntimeError(err_msg.str());
                }
            }

            (*m_p_params)(col, row) = p;
        }

        void set_max_triangle_area(Float max_area) {
            m_t_params.max_area = max_area;
        }

        void generate_periodic_pattern() {
            assert((*m_p_params)(0, 0) != NULL);
            m_inflator.generatePattern(*(*m_p_params)(0, 0), m_t_params, m_mesh);
        }

        void generate_tiled_pattern() {
            std::cout << "generating "
                << m_p_params->height()
                << "x"
                << m_p_params->width()
                << " tiled pattern" << std::endl;
            m_inflator.generateTiledPattern(*m_p_params, m_t_params, m_mesh);
        }

        VectorF get_vertices() {
            const size_t dim = 2;
            typedef WireInflator2D::OutMeshType MeshType;
            MeshType::NodeVector nodes = m_mesh.nodes;
            VectorF vertices(nodes.size() * dim);
            for (size_t i=0; i<nodes.size(); i++) {
                for (size_t j=0; j<dim; j++) {
                    vertices[i*dim + j] = nodes[i][j];
                }
            }
            return vertices;
        }

        VectorI get_triangles() {
            const size_t dim = 2;
            const size_t nodes_per_elem = 3;
            typedef WireInflator2D::OutMeshType MeshType;
            MeshType::ElementVector elems = m_mesh.elements;

            VectorI faces(elems.size() * nodes_per_elem);
            for (size_t i=0; i<elems.size(); i++) {
                for (size_t j=0; j<nodes_per_elem; j++) {
                    faces[i*nodes_per_elem + j] = elems[i][j];
                }
            }
            return faces;
        }

        MatrixFr get_boundary_velocity() {
            const size_t dim = 2;
            typedef WireInflator2D::OutMeshType MeshType;

            MatrixFr bd_velocity = MatrixF::Zero(m_mesh.nodes.size(), get_num_parameters());
            size_t count = 0;
            for (auto itr = m_mesh.edge_fields.begin();
                    itr != m_mesh.edge_fields.end(); itr++) {
                MeshType::EdgeType edge = itr->first;
                MeshType::Fields edge_velocity = itr->second;
                for (size_t i=0; i<edge_velocity.size(); i++) {
                    std::cout << edge_velocity[i] << " ";
                }
                std::cout << std::endl;
                std::copy(edge_velocity.begin(), edge_velocity.end(),
                        bd_velocity.row(edge.first).data());
                std::copy(edge_velocity.begin(), edge_velocity.end(),
                        bd_velocity.row(edge.second).data());
                count++;
            }

            return bd_velocity;
        }

    private:
        typedef Array2D<CellParameters*> ParameterGrid;
        typedef std::shared_ptr<ParameterGrid> ParameterGridPtr;
        size_t m_rows;
        size_t m_cols;
        TessellationParameters m_t_params;
        ParameterGridPtr  m_p_params;
        WireInflator2D::OutMeshType m_mesh;
        WireInflator2D m_inflator;
};
