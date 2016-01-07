////////////////////////////////////////////////////////////////////////////////
// InflatorBase.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      CRTP "abstract base class" for inflators.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/10/2015 22:16:39
////////////////////////////////////////////////////////////////////////////////
#ifndef INFLATORBASE_HH
#define INFLATORBASE_HH

#include <vector>
#include <MeshIO.hh>
#include <Functions.hh> 

enum class ParameterType { Thickness, Offset };

// Normal shape velocities are (discontinuous) piecwise linear scalar fields.
template<size_t N>
using NormalShapeVelocity = std::vector<Interpolant<Real, N - 1, 1>>;

template<class _Derived>
class InflatorBase {
public:
    void clear() {
        m_elements.clear();
        m_vertices.clear();
    }

    const std::vector<MeshIO::IOElement> &elements() const { return m_elements; }
    const std::vector<MeshIO::IOVertex>  &vertices() const { return m_vertices; }

protected:
    std::vector<MeshIO::IOElement> m_elements;
    std::vector<MeshIO::IOVertex>  m_vertices;
};

#endif /* end of include guard: INFLATORBASE_HH */
