////////////////////////////////////////////////////////////////////////////////
// mex_homogenize.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Top-level mex code for homogenization.
//      Due to difficulties with c++11 in mex, this is simiply a wrapper
//      that calls code from Homogenize.cc, which is compiled using a c++11
//      compliant compiler.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  05/20/2015 14:11:23
////////////////////////////////////////////////////////////////////////////////
#include "mex.h"
#include <stdexcept>
#include <string>
#include "Homogenize.hh"
#include <Eigen/Dense>

// #if !((DIMENSION == 2) || (DIMENSION == 3))
// #error "must set preprocessor constant DIMENSION"
// #endif

typedef double real;
using namespace std;

void mexFunction(
    int             nlhs,
    mxArray      *plhs[],
    int             nrhs,
    const mxArray *prhs[])
{
    // Parsing input
    if ((nrhs != 2) && (nrhs != 3)) {
        mexErrMsgTxt("usage: [moduli, etensor] = homogenize(meshPath, materialPath[, defoJacobian])");
    }

    const char *meshPath = mxArrayToString(prhs[0]);
    const char *materialPath = mxArrayToString(prhs[1]);

    Eigen::MatrixXd jacobian(0, 0);
    if (nrhs == 3) {
        // Size checking will be done in homogenizer
        jacobian.resize(mxGetM(prhs[2]), mxGetN(prhs[2]));
        // Column major by default, like MATLAB
        jacobian = Eigen::Map<Eigen::MatrixXd>(mxGetPr(prhs[2]),
                jacobian.rows(), jacobian.cols());
    }

    if (nlhs < 1 || nlhs > 2)
        mexErrMsgTxt("Expects one or two output arguments.");

    try {
        std::vector<double> moduli;
        Eigen::MatrixXd etensor;
        homogenize(meshPath, materialPath, jacobian, moduli, etensor);

        plhs[0] = mxCreateDoubleMatrix(moduli.size(), 1, mxREAL); // moduli
        double *moduli_out = mxGetPr(plhs[0]);
        std::copy(moduli.begin(), moduli.end(), moduli_out);

        if (nlhs == 2) {
            // Output full elasticity tensor if requested
            plhs[1] = mxCreateDoubleMatrix(etensor.rows(), etensor.cols(), mxREAL); // elasticity tensor
            double *etensor_out = mxGetPr(plhs[1]);
            size_t size = etensor.rows() * etensor.cols();
            std::copy(etensor.data(), etensor.data() + size, etensor_out);
        }
    }
    catch(exception &e) {
        mexErrMsgTxt((string("Exception: ") + e.what()).c_str());
    }
    catch(...) {
        mexErrMsgTxt("Unknown exception");
    }

    return;
}
