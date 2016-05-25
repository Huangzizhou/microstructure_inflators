#include "MakeInflator.hh"

#include <Future.hh>
#include <stdexcept>
#include <Utilities/ci_string.hh>

#include "Inflator.hh"
#include "inflators/IsoinflatorWrapper.hh"
#include "inflators/JamesInflatorWrapper.hh"
#include "inflators/LuigiInflatorWrapper.hh"
#include "inflators/LpHoleInflator.hh"
#include "inflators/EqualityConstrainedInflator.hh"

using namespace std;
namespace po = boost::program_options;

template<typename T>
T extract_required(po::variables_map &opts, const string &key) {
    if (opts.count(key)) {
        T val(opts[key].as<T>());
        opts.erase(key);
        return val;
    }
    throw runtime_error("Options passed to makeInflator missing required option '" + key + "'");
}

template<typename T>
T extract_defaulted(po::variables_map &opts, const string &key, T default_val) {
    if (opts.count(key)) {
        T val = opts[key].as<T>();
        opts.erase(key);
        return val;
    }
    return default_val;
}

bool extract_flag(po::variables_map &opts, const string &key) {
    if (opts.count(key)) { opts.erase(key); return true; }
    return false;
}

template<typename T>
void extract_notify(po::variables_map &opts, const string &key, const T &func) {
    static_assert(std::is_same<typename function_traits<T>::result_type, void>::value, "Notify function must return void.");
    static_assert(function_traits<T>::arity == 1, "Notify function must take a single argument.");
    using OptionType = typename function_traits<T>::template arg<0>::type;
    if (opts.count(key)) {
        func(opts[key].as<OptionType>());
        opts.erase(key);
    }
}

unique_ptr<InflatorBase> make_inflator(const string &name, po::variables_map opts,
                                       const vector<string> &constraints) {
    unique_ptr<InflatorBase> infl;
    if (ci_string("Isosurface2D") == name.c_str()) {
        infl = Future::make_unique<IsoinflatorWrapper<2>>(
                extract_required<string>(opts, "pattern"),

                extract_flag(opts, "isotropicParameters"),
                extract_flag(opts, "vertexThickness"));
    }
    else if (ci_string("Isosurface3D") == name.c_str()) {
        infl = Future::make_unique<IsoinflatorWrapper<3>>(
                extract_required<string>(opts, "pattern"),
                extract_flag(opts, "isotropicParameters"),
                extract_flag(opts, "vertexThickness"));
    }
    else if (ci_string("James") == name.c_str()) {
        infl = Future::make_unique<JamesInflatorWrapper>(
                extract_required<string>(opts, "pattern"),
                extract_defaulted<double>(opts, "cell_size", 5.0),
                0.5 * sqrt(2),
                extract_flag(opts, "isotropicParameters"),
                extract_flag(opts, "vertexThickness"));
        infl->configureSubdivision(extract_defaulted<string>(opts, "sub_algorithm", "simple"),
                                   extract_defaulted<size_t>(opts,     "subdivide",        0));
    }
    else if (ci_string("Luigi") == name.c_str()) {
        infl = Future::make_unique<LuigiInflatorWrapper>(extract_required<string>(opts, "pattern"));
    }
    else if (ci_string("LpHole") == name.c_str()) {
        auto lphole_infl = Future::make_unique<LpHoleInflator>();
        extract_notify(opts, "hole_segments", [&](size_t hs) { lphole_infl->setNumSubdiv(hs); });
        infl = std::move(lphole_infl);
    }
    else throw runtime_error("Invalid inflator: " + name);

    infl->setReflectiveInflator(!extract_flag(opts, "fullCellInflator"));
    extract_notify(opts, "meshingOptions", [&](const string &v) { infl->loadMeshingOptions(v); });
    extract_notify(opts, "max_volume",     [&](double        v) { infl->setMaxElementVolume(v); });

    // Report any ignored options.
    for (const auto &entry : opts)
        cerr << "WARNING: inflator " << name << " ignored option " << entry.first << endl;

    // Wrap the inflator in an EqualityConstrainedInflator if there are
    // constraints.
    if (constraints.size()) {
        InflatorBase *rawPtr = infl.release();
        if (auto infl2D = dynamic_cast<Inflator<2> *>(rawPtr))
            infl = Future::make_unique<EqualityConstrainedInflator<2>>(unique_ptr<Inflator<2>>(infl2D), constraints);
        else if (auto infl3D = dynamic_cast<Inflator<3> *>(rawPtr))
            infl = Future::make_unique<EqualityConstrainedInflator<3>>(unique_ptr<Inflator<3>>(infl3D), constraints);
        else throw runtime_error("Invalid inflator dimension.");
    }

    return infl;
}

// Extract the options that are meant to be passed to the inflator constructor:
po_vm filterInflatorOptions(const po_vm &opts) {
    auto keys = {
        "pattern",
        "isotropicParameters",
        "vertexThickness",
        "cell_size",
        "hole_segments",
        "max_volume",
        "meshingOptions",
        "subdivide",
        "sub_algorithm"
    };
    po_vm filtered;
    for (const string &key : keys)
        if (opts.count(key)) filtered.emplace(key, opts[key]);

    return filtered;
}