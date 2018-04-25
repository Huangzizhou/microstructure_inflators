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

#include "PatternOptimizationIterate.hh"
#include <optimizers/IterateManagerBase.hh>
#include <memory>

namespace PatternOptimization {

// Implementation of IterateManagerBase for a particular IterateFactory.
template<class _ItFactory>
struct IterateManager : public IterateManagerBase {
    using Iterate  = typename _ItFactory::Iterate;

    IterateManager(std::unique_ptr<_ItFactory> itFactory)
        : m_iterFactory(std::move(itFactory)) { }

    Iterate &get(size_t nParams, const double * const params) override {
        m_currIterate = m_iterFactory->getIterate(std::move(m_currIterate), nParams, &params[0]);
        return *m_currIterate;
    }

    virtual       Iterate &get()          override { assert(m_currIterate); return *m_currIterate; }
    virtual const Iterate &get()    const override { assert(m_currIterate); return *m_currIterate; }

    virtual       Iterate *getPtr()       override { return m_currIterate.get(); }
    virtual const Iterate *getPtr() const override { return m_currIterate.get(); }

    virtual size_t numParameters()  const override { return m_iterFactory->numParameters(); }
    virtual bool    isParametric()  const override { return m_iterFactory->isParametric(); }

    virtual ~IterateManager() { }
private:
    std::unique_ptr<_ItFactory> m_iterFactory;
    std::unique_ptr<Iterate>    m_currIterate;
};

template<class _IF>
std::shared_ptr<IterateManager<_IF>> make_iterate_manager(std::unique_ptr<_IF> itFactory) {
    return std::make_shared<IterateManager<_IF>>(std::move(itFactory));
}

}

#endif /* end of include guard: ITERATEMANAGER_HH */