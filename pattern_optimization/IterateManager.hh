////////////////////////////////////////////////////////////////////////////////
// IterateManager.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Wraps the creation and re-use of iterates for particular parameter
//      values. (Repeated calls to get() with the same parameters only create
//      one iterate).
//
//      IterateManagerBase provides a non-templated interface to the iterate
//      manager.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/28/2016 17:30:15
////////////////////////////////////////////////////////////////////////////////
#ifndef ITERATEMANAGER_HH
#define ITERATEMANAGER_HH

#include <memory>
#include "PatternOptimizationIterate.hh"
#include "IterateManagerBase.hh"

namespace PatternOptimization {

// Implementation of IterateManagerBase for a particular IterateFactory.
template<class _ItFactory>
struct IterateManager : public IterateManagerBase {
    using Iterate  = typename _ItFactory::Iterate;

    IterateManager(std::shared_ptr<_ItFactory> itFactory)
        : m_iterFactory(itFactory) { }

    Iterate &get(size_t nParams, const double * const params) {
        m_currIterate = m_iterFactory->getIterate(std::move(m_currIterate), nParams, &params[0]);
        return *m_currIterate;
    }

          Iterate &get()       { assert(m_currIterate); return *m_currIterate; }
    const Iterate &get() const { assert(m_currIterate); return *m_currIterate; }

          Iterate *getPtr()       { return m_currIterate.get(); }
    const Iterate *getPtr() const { return m_currIterate.get(); }

    virtual size_t numParameters() const { return m_iterFactory->numParameters(); }

    virtual ~IterateManager() { }
private:
    std::shared_ptr<_ItFactory> m_iterFactory;
    std::unique_ptr<Iterate>    m_currIterate;
};

template<class _IF>
std::shared_ptr<IterateManager<_IF>> make_iterate_manager(std::shared_ptr<_IF> itFactory) {
    return std::make_shared<IterateManager<_IF>>(itFactory);
}

}

#endif /* end of include guard: ITERATEMANAGER_HH */
