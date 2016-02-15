#ifndef INFLATORTRAITS_INL
#define INFLATORTRAITS_INL

template<size_t N>
struct InflatorTraits {
    template<class Sim>
    using Iterate = WCStressOptimization::Iterate<Sim>;

    template<class type>
    static SField initParams(shared_ptr<type> iptr,
                         const po::variables_map &/* args */,
                         const PatternOptimization::Job<N> *job) {
        if (job->numParams() != iptr->numParameters()) {
            for (size_t i = 0; i < iptr->numParameters(); ++i) {
                cout << "param " << i << " role: " <<
                    (iptr->parameterType(i) == ParameterType::Offset ? "Offset" : "Thickness")
                    << endl;
            }
            throw runtime_error("Invalid number of parameters.");
        }

        SField params(job->initialParams);
        for (const auto &boundEntry : job->varLowerBounds) {
            if (boundEntry.first > params.domainSize())
                cerr << "WARNING: bound on nonexistent variable" << endl;
        }

        for (size_t p = 0; p < params.domainSize(); ++p) {
            if (job->varLowerBounds.count(p)) {
                 if ((params[p] < job->varLowerBounds.at(p)) ||
                     (params[p] > job->varUpperBounds.at(p))) {
                    throw std::runtime_error("Initial point infeasible");
                 }
            }
        }
        return params;
    }

    template<class type>
    static void finalize(shared_ptr<type> /* iptr */, const std::vector<Real> &result,
                         const po::variables_map &/* args */,
                         const PatternOptimization::Job<N> * /* job */) {
        std::cout << "Final p:";
        for (size_t i = 0; i < result.size(); ++i)
            cout << "\t" << result[i];
        cout << endl;
    }
};

template<size_t N>
struct InflatorTraitsConstrainedInflator;

template<size_t N>
struct InflatorTraitsLpHole;

template<>
struct InflatorTraitsConstrainedInflator<2> : public InflatorTraits<2> {
    using type = ConstrainedInflator<2>;
    static shared_ptr<type> construct(const po::variables_map &args, const PatternOptimization::Job<2> *job) {
        auto inflator_ptr = make_shared<type>(
                job->parameterConstraints,
                args["pattern"].as<string>(),
                args["cell_size"].as<double>(),
                0.5 * sqrt(2),
                args.count("isotropicParameters"),
                args.count("vertexThickness"));
        if (args.count("max_volume"))
            inflator_ptr->setMaxElementVolume(args["max_volume"].as<double>());
        return inflator_ptr;
    }
    static void finalize(shared_ptr<type> iptr, const std::vector<Real> &result,
                         const po::variables_map &args,
                         const PatternOptimization::Job<2> * job) {
        InflatorTraits<2>::finalize(iptr, result, args, job);
    }
};

template<>
struct InflatorTraitsConstrainedInflator<3> : public InflatorTraits<3> {
    using type = ConstrainedInflator<3>;
    static shared_ptr<type> construct(const po::variables_map &args, const PatternOptimization::Job<3> *job) {
        auto inflator_ptr = make_shared<type>(
                job->parameterConstraints,
                args["pattern"].as<string>(),
                args["cell_size"].as<double>(),
                0.5 * sqrt(2),
                args.count("isotropicParameters"),
                args.count("vertexThickness"));
        inflator_ptr->configureSubdivision(args["sub_algorithm"].as<string>(),
                                           args["subdivide"].as<size_t>());
        inflator_ptr->setReflectiveInflator(args.count("fullCellInflator") == 0);
        if (args.count("dofOut"))
            inflator_ptr->setDoFOutputPrefix(args["dofOut"].as<string>());
        if (args.count("max_volume"))
            inflator_ptr->setMaxElementVolume(args["max_volume"].as<double>());

        return inflator_ptr;
    }

    static void finalize(shared_ptr<type> iptr, const std::vector<Real> &result,
                         const po::variables_map &args,
                         const PatternOptimization::Job<3> * job) {
        InflatorTraits<3>::finalize(iptr, result, args, job);
        if (args.count("dofOut"))
            iptr->writePatternDoFs(args["dofOut"].as<string>() + ".final.dof", result);
    }
};

template<>
struct InflatorTraitsLpHole<2> : public InflatorTraits<2> {
    using type = LpHoleInflator;
    static shared_ptr<type> construct(const po::variables_map &args, const PatternOptimization::Job<2> * /* job */) {
        if (args.count("pattern"))
            std::cerr << "WARNING: pattern argument not expected for LpHoleInflator" << std::endl;
        auto inflator_ptr = make_shared<type>();
        inflator_ptr->setNumSubdiv(args["hole_segments"].as<size_t>());
        if (args.count("max_volume"))
            inflator_ptr->setMaxElementVolume(args["max_volume"].as<double>());
        return inflator_ptr;
    }

    static void finalize(shared_ptr<type> iptr, const std::vector<Real> &result,
                         const po::variables_map &args,
                         const PatternOptimization::Job<2> * job) {
        InflatorTraits<2>::finalize(iptr, result, args, job);
    }
};

template<size_t N>
struct InflatorTraitsBoundaryPerturbation : public InflatorTraits<N> {
    using type = BoundaryPerturbationInflator<N>;

    // Special iterate for BoundaryPerturbationInflator
    template<class Sim>
    using Iterate = WCStressOptimization::BoundaryPerturbationIterate<Sim>;

    static shared_ptr<type> construct(const po::variables_map &args, const PatternOptimization::Job<N> * /* job */) {
        std::vector<MeshIO::IOVertex>  vertices;
        std::vector<MeshIO::IOElement> elements;
        MeshIO::load(args["pattern"].as<string>(), vertices, elements);

        auto inflator_ptr = make_shared<type>(vertices, elements);
        return inflator_ptr;
    }

    // For boundary pertrubation inflator, ignore the initial params set in job
    // file unless the sizes match (it's hard to specify all params).
    static SField initParams(shared_ptr<type> iptr,
                         const po::variables_map &args,
                         const PatternOptimization::Job<N> *job) {
        if (job->numParams() != iptr->numParameters()) {
            SField params(iptr->numParameters());
            params.clear();
            return params;
        }
        return InflatorTraits<N>::initParams(iptr, args, job);
    }

    static void finalize(shared_ptr<type> /* iptr */, const std::vector<Real> &/* result */,
                         const po::variables_map &/* args */,
                         const PatternOptimization::Job<N> * /* job */) {
    }
};

#endif /* end of include guard: INFLATORTRAITS_INL */
