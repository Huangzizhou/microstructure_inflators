#ifndef HOMOGENIZE_HH
#define HOMOGENIZE_HH

#include <cstdlib> // size_t
#include <vector>
#include <Eigen/Dense>

void homogenize(const char *meshPath, const char *materialPath,
                const Eigen::MatrixXd &jacobian,
                std::vector<double> &moduli,
                Eigen::MatrixXd &elasticityTensor);

#endif /* end of include guard: HOMOGENIZE_HH */
