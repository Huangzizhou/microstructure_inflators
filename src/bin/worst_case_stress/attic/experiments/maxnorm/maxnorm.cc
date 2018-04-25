#include "quartic_form.hh"
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/Types.hh>
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>

#include <algorithm>
#include <vector>

#include <boost/algorithm/string.hpp>

template<size_t N>
void getSampleMesh(std::vector<MeshIO::IOVertex> &vertices,
                   std::vector<MeshIO::IOElement> &elements,
                   std::vector<Real> &theta1,
                   std::vector<Real> &theta2) {
    if (N == 3) { MeshIO::load("sphere.off", vertices, elements); }

    if (N == 2) {
        constexpr size_t nsubdivs = 1024;
        theta1.clear(), theta2.clear();
        theta1.reserve(nsubdivs);
        for (size_t i = 0; i < nsubdivs; ++i) {
            Real theta = M_PI * Real(i) / nsubdivs; // half circle
            vertices.emplace_back(cos(theta), sin(theta));
            theta1.push_back(theta);
        }
        for (size_t i = 0; i < nsubdivs - 1; ++i)
            elements.emplace_back(i, i + 1);
        // elements.emplace_back(vertices.size() - 1, 0);
    }
}

template<size_t N>
void execute(const std::vector<std::vector<std::vector<Real>>> &tensors,
             const std::string &geometryPath) {
    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    std::vector<Real> theta1, theta2;
    getSampleMesh<N>(vertices, elements, theta1, theta2);

    std::vector<VectorND<N>> pts;
    for (auto &v : vertices)
        pts.emplace_back(truncateFrom3D<VectorND<N>>(v.point));

    ScalarField<Real> maxima(tensors.size());
    for (size_t t = 0; t < tensors.size(); ++t) {
        const auto &components = tensors[t];
        ElasticityTensor<Real, N> T;
        for (size_t i = 0; i < flatLen(N); ++i) {
            for (size_t j = i;  j < flatLen(N); ++j) {
                if (std::abs(components[i][j] - components[j][i]) > 1e-8)
                    throw std::runtime_error("Only support major-symmetric tensors.");
                T.D(i, j) = components[i][j];
            }
        }

        QuarticForm<Real, N> q(T);

        ScalarField<Real> qval(vertices.size());
        Real maxq = 0.0;
        for (size_t i = 0; i < vertices.size(); ++i) {
            qval[i] = q(pts[i]);
            maxq = std::max(maxq, qval[i]);
        }

        maxima[t] = maxq;
        

        // MSHFieldWriter writer("result.msh", vertices, elements);
        // writer.addField("qval", qval, DomainType::PER_NODE);

        // std::ofstream qvalOut("qvalues.txt");
        // if (!qvalOut) throw std::runtime_error("Couldn't open 'qvalues.txt'");
        // for (size_t i = 0; i < vertices.size(); ++i) {
        //     if (theta1.size() == vertices.size()) qvalOut << theta1[i] << '\t';
        //     if (theta2.size() == vertices.size()) qvalOut << theta2[i] << '\t';
        //     qvalOut << qval[i] << std::endl;
        // }
    }

    if (geometryPath.size()) {
        std::vector<MeshIO::IOVertex > outVertices;
        std::vector<MeshIO::IOElement> outElements;
        MeshIO::load(geometryPath, outVertices, outElements);
        if (outElements.size() != tensors.size())
            throw std::runtime_error("Expected one tensor per element");
        MSHFieldWriter writer("maxima.msh", outVertices, outElements);
        writer.addField("maxq", maxima, DomainType::PER_ELEMENT);

        ScalarField<Real> maxStress(maxima.domainSize());
        for (size_t i = 0; i < maxima.domainSize(); ++i)
            maxStress[i] = std::sqrt(maxima[i]);

        writer.addField("max stress", maxStress, DomainType::PER_ELEMENT);
    }
}

int main(int argc, const char *argv[]) {
    if ((argc != 2) && (argc != 3)) {
        std::cerr << "usage: ./maxnorm tensors.txt [geometry.msh]" << std::endl;
        exit(1);
    }

    std::string tensorFilePath = argv[1];
    std::string geometryPath = (argc > 2) ? argv[2] : "";

    std::ifstream tensorFile(tensorFilePath);
    if (!tensorFile.is_open()) { throw std::runtime_error("Couldn't open tensor"); }

    // Read flattened 2- or 3D 4th order tensors
    std::vector<std::vector<std::vector<Real>>> tensors;
    std::runtime_error parseError("Invalid tensor; failed to parse 2- or 3D 4th order tensor.");
    std::string line;
    size_t tensorSize = 0;
    tensors.emplace_back();
    while (std::getline(tensorFile, line)) {
        boost::trim(line);
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of("\t "),
                     boost::token_compress_on);
        if (tensorSize == 0) {
            tensorSize = tokens.size();
            if (!((tensorSize == 3) || (tensorSize == 6)))
                throw parseError;
        }
        if ((tokens.size() != tensorSize)) throw parseError;
        if (tensors.back().size() == tensorSize) tensors.emplace_back();
        auto &tensor = tensors.back();
        tensor.push_back(std::vector<Real>(tensorSize));
        std::transform(tokens.begin(), tokens.end(),
                       tensor.back().begin(),
                       [](const std::string &s) { return std::stod(s); });
    }
    if (tensors.back().size() != tensorSize) throw parseError;

    if (tensorSize == 3) { execute<2>(tensors, geometryPath); }
    else                 { execute<3>(tensors, geometryPath); }

    return 0;
}
