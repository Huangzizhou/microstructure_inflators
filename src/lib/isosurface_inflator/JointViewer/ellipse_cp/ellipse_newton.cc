////////////////////////////////////////////////////////////////////////////////
// ellipse_newton.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Solves for the closest point on an axis-aligned ellipse with major axis
//  a and minor axis b <= a:
//  e(t) = [a cos(t), b sin(t)]
//  (Note: the major axis is aligned with the first coordinate).
//  Applies Newton's method with a carefully chosen starting point that
//  guarantees convergence.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/31/2016 14:48:46
////////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <fstream>

#include <iostream>
#include "ellipse_newton.hh"

#ifndef OUTPUT_FILE
#define OUTPUT_FILE 1
#endif


int main(int argc, const char *argv[]) {
    if (!((argc == 3) || (argc == 5))) {
        std::cerr << "usage: ./ellipse_newton a b [x y]" << std::endl;
        exit(-1);
    }

    double a = std::stod(argv[1]),
           b = std::stod(argv[2]);

    double aspect = std::min(a / b, 3.0);
    size_t gridSizeX = 1024,
           gridSizeY = ceil(gridSizeX / aspect);
    if (gridSizeY < 2) throw std::runtime_error("Invalid aspect ratio");
    double gridSpacing = 1.5 * a / (gridSizeX - 1);
    // For now, always sample square domain--makes creating animations easier.
    gridSizeY = gridSizeX;

    std::cout << std::setprecision(19);
    std::cerr << std::setprecision(19);

#if OUTPUT_FILE
    std::ofstream iterateFile("iterates.txt");
    iterateFile << std::setprecision(19);
#endif

    if (argc > 3) {
        double x = std::stod(argv[3]),
               y = std::stod(argv[4]);
        auto cp = ellipseClosestPoint<double, true>(a, b, x, y);
        std::cout << std::get<0>(cp) << "\t" << std::get<1>(cp) << std::endl;
    }
    else {
        for (size_t xi = 0; xi < gridSizeX; ++xi) {
            for (size_t yi = 0; yi < gridSizeY; ++yi) {
                double x = xi * gridSpacing,
                       y = yi * gridSpacing;
                auto cp = ellipseClosestPoint<double>(a, b, x, y);
#if OUTPUT_FILE
                iterateFile << x
                    << '\t' << y
                    << '\t' << std::get<0>(cp) // cp_x
                    << '\t' << std::get<1>(cp) // cp_y
                    << '\t' << std::get<2>(cp) // dist
                    << '\t' << std::get<3>(cp) // iterations
                    << '\t' << std::get<4>(cp) // backtracking iterations
                    << '\n';
#endif
            }
#if OUTPUT_FILE
            iterateFile  << '\n'; // gnuplot wants blank lines between x values
#endif
        }
    }
    return 0;
}
