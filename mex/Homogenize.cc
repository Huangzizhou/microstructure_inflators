#include "Homogenize.hh"
#include "LinearElasticity.hh"
#include "Materials.hh"
#include "PeriodicHomogenization.hh"
#include <vector>
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "MeshIO.hh"

using namespace std;
using namespace PeriodicHomogenization;

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N, size_t _FEMDegree>
void execute(const vector<MeshIO::IOVertex> &inVertices, 
             const vector<MeshIO::IOElement> &inElements,
             const char *materialPath,
             Eigen::MatrixXd jacobian,
             std::vector<double> &moduli,
             Eigen::MatrixXd &elasticityTensor)
{
    auto &mat = HMG<_N>::material;
    mat.setFromFile(materialPath);
    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;
    Simulator sim(inElements, inVertices);

    if ((jacobian.rows() == 0) && (jacobian.cols() == 0)) {
        jacobian = Eigen::MatrixXd::Identity(_N, _N);
    }
    if ((jacobian.rows() != 3) || (jacobian.cols() != 3)) {
        auto nstr = std::to_string(_N);
        throw std::runtime_error("Jacobian should be " + nstr + "x" + nstr);
    }

    // Morteza's transformation formulas
    mat.setTensor(mat.getTensor().transform(jacobian.inverse()));
    vector<typename Simulator::VField> w_ij;
    solveCellProblems(w_ij, sim);
    auto EhDefo = homogenizedElasticityTensor(w_ij, sim).transform(jacobian);
    EhDefo.getOrthotropicParameters(moduli);

    elasticityTensor.resize(flatLen(_N), flatLen(_N));
    for (size_t i = 0; i < flatLen(_N); ++i)
        for (size_t j = 0; j < flatLen(_N); ++j)
            elasticityTensor(i, j) = EhDefo.D(i, j);
}

void homogenize(const char *meshPath, const char *materialPath,
                const Eigen::MatrixXd &jacobian,
                std::vector<double> &moduli,
                Eigen::MatrixXd &elasticityTensor)
{
    vector<MeshIO::IOVertex>  inVertices;
    vector<MeshIO::IOElement> inElements;
    auto type = load(meshPath, inVertices, inElements,
            MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

    // Infer dimension from mesh type.
    size_t dim;
    if      (type == MeshIO::MESH_TET) dim = 3;
    else if (type == MeshIO::MESH_TRI) dim = 2;
    else    throw std::runtime_error("Mesh must be triangle or tet.");

    // Look up and run appropriate instantiation.
    // We always run degree 2 FEM for now.
    auto exec = (dim == 3) ? execute<3, 2> : execute<2, 2>;
    exec(inVertices, inElements, materialPath, jacobian, moduli,
         elasticityTensor);
}
