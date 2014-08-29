#pragma once
#include <sstream>
#include <iostream>

#include "WireInflator2D.h"
#include "EigenTypedef.h"
#include "Exception.h"

class WireInflatorFacade {
    public:
        WireInflatorFacade(size_t rows, size_t cols) : 
            m_rows(rows), m_cols(cols), m_p_params(cols, rows) {
                m_t_params.max_area = 0.001;
                for (size_t row = 0; row < m_rows; row++) {
                    for (size_t col=0; col < m_cols; col++) {
                        m_p_params(col, row) = NULL;
                    }
                }
            }

        size_t get_num_parameters() const {
            WireInflator2D::PatternParameters p;
            return p.numberOfParameters();
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
            WireInflator2D::PatternParameters* p = new WireInflator2D::PatternParameters();

            const size_t num_param = p->numberOfParameters();
            for (size_t i=0; i<num_param; i++) {
                std::pair<double, double> range = p->parameterRange(i);
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

            m_p_params(col, row) = p;
        }

        void set_max_triangle_area(Float max_area) {
            m_t_params.max_area = max_area;
        }

        void generate_periodic_pattern() {
            assert(m_p_params(0, 0) != NULL);
            WireInflator2D::generatePattern(*m_p_params(0, 0), m_t_params, m_mesh);
        }

        void generate_tiled_pattern() {
            std::cout << "generating "
                << m_p_params.height()
                << "x"
                << m_p_params.width()
                << " tiled pattern" << std::endl;
            WireInflator2D::generateTiledPattern(m_p_params, m_t_params, m_mesh);
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
        size_t m_rows;
        size_t m_cols;
        TessellationParameters m_t_params;
        Array2D<WireInflator2D::PatternParameters*> m_p_params;
        WireInflator2D::OutMeshType m_mesh;
};
