////////////////////////////////////////////////////////////////////////////////
// ParameterConstraints.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Parses and holds constraints relating the pattern parameters.
//      Constraints should be LINEAR and of the form:
//      p0  = p1 + c
//      p2 <= p1 + c
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/11/2015 14:02:54
////////////////////////////////////////////////////////////////////////////////
#ifndef PARAMETERCONSTRAINTS_HH
#define PARAMETERCONSTRAINTS_HH
#include <string>
#include <Types.hh>
#include <util.h>
#include <ExpressionVector.hh>

class ParameterConstraint {
public:
    enum class Type { EQ, LEQ, LE, GEQ, GE, UNKNOWN };
    ParameterConstraint(size_t numParams, const std::string &constraint)
        : m_numParams(numParams), m_type(Type::UNKNOWN)
    {
        std::runtime_error invalid("Invalid constraint.");

        // NOTE: Must be arranged longest to shortest.
        static const std::vector<std::string> tokens =
            { "<=", ">=", "<", ">", "=" };
        static const std::vector<Type> types =
            { Type::LEQ, Type::GEQ, Type::LE, Type::GE, Type::EQ };

        std::string lhs, rhs;
        for (size_t i = 0; i < tokens.size(); ++i) {
            const auto &token = tokens[i];
            auto pos = constraint.find(token);
            if (pos == std::string::npos) continue;
            m_type = types[i];
            lhs = trim(constraint.substr(0, pos));
            rhs = trim(constraint.substr(pos + token.size(), std::string::npos));
            break;
        }
        if (m_type == Type::UNKNOWN) throw invalid;
        if ((lhs.size() == 0) || (rhs.size() == 0)) throw invalid;

        // Make sure there weren't any other (in)equality operators
        for (const auto &token : tokens) {
            if ((lhs.find(token) != std::string::npos) ||
                (rhs.find(token) != std::string::npos)) {
                throw invalid;
            }
        }

        Expression lhsExpr(lhs), rhsExpr(rhs);
        ExpressionEnvironment env;

        env.setVectorValue("p", std::vector<Real>(m_numParams, 0.0));
        // Collect all the constants onto the RHS
        Real lhsConst = lhsExpr.eval(env);
        Real rhsConst = rhsExpr.eval(env);
        m_rhsConst = rhsConst - lhsConst;

        // Probe for the parameter coefficients, moving them to the LHS
        m_lhsCoeffs.assign(m_numParams, 0.0);
        for (size_t p = 0; p < m_numParams; ++p) {
            env.setValue("p" + std::to_string(p), 1.0);
            m_lhsCoeffs[p] = (lhsExpr.eval(env) - lhsConst) -
                             (rhsExpr.eval(env) - rhsConst);
            env.setValue("p" + std::to_string(p), 0.0);
        }
    }

    // Get this constraint's row in an augmented system of the form (LHS | RHS)
    std::vector<Real> augmentedRow() const {
        std::vector<Real> result(m_lhsCoeffs);
        result.push_back(m_rhsConst);
        return result;
    }

    Type type() const { return m_type; }
    bool isEqualityConstraint() const { return m_type == Type::EQ; }

    friend std::ostream & operator<<(std::ostream &os, const ParameterConstraint &c) {
        for (size_t p = 0; p < c.m_lhsCoeffs.size(); ++p) {
            os << (p ? " + " : "") << c.m_lhsCoeffs[p] << " * p" << p;
        }

        switch (c.m_type) {
            case Type::EQ:  os <<  " = "; break;
            case Type::LE:  os <<  " < "; break;
            case Type::GE:  os <<  " > "; break;
            case Type::LEQ: os << " <= "; break;
            case Type::GEQ: os << " >= "; break;
            default: throw std::runtime_error("Invalid constraint.");
        }

        os << c.m_rhsConst;
        return os;
    }

private:
    size_t m_numParams;
    Type m_type;
    Real m_rhsConst;
    std::vector<Real> m_lhsCoeffs;
};

#endif /* end of include guard: PARAMETERCONSTRAINTS_HH */
