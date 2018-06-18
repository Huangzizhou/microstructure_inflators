#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cassert>

using Real = double;
#ifndef OUTPUT_FILE
#define OUTPUT_FILE 1
#endif

int maxIterations = 0;

Real RobustLength(Real a, Real b) {
    if (std::abs(b) > std::abs(a))
        std::swap(a, b);
    // |a| > |b|
    b /= a;
    return std::abs(a) * sqrt(1 + b * b);
}

Real Sqr(Real a) { return a * a; }

Real GetRoot(Real r0 , Real z0 , Real z1 , Real g)
{
    Real n0 = r0*z0;
    Real s0 = z1 - 1, s1 = (g < 0 ? 0 : RobustLength(n0, z1) - 1); Real s = 0;
    for (int i = 0; i < maxIterations; ++i)
    {
        s = (s0 + s1) / 2.0;
        if (s == s0 || s == s1) { break; }
        Real ratio0 = n0/(s + r0), ratio1 = z1/(s + 1.0);
        g = Sqr(ratio0) + Sqr(ratio1) - 1.0;
        if (g > 0) { s0 = s; } else if (g < 0) { s1 = s; } else { break; }
    }
    return s ;
}

Real DistancePointEllipse( Real e0 , Real e1 , Real y0 , Real y1 , Real &x0 , Real &x1) {
    Real distance;
    assert(e0 >= e1);
    if (y1 > 0){
        if (y0 > 0){
            Real z0 = y0 / e0; 
            Real z1 = y1 / e1; 
            Real g = z0*z0+z1*z1 - 1.0;
            if (g != 0){
                Real r0 = (e0/e1)*(e0/e1);
                Real sbar = GetRoot(r0 , z0 , z1 , g);
                x0 = r0 * y0 /( sbar + r0 );
                x1 = y1 /( sbar + 1 );
                distance = sqrt( (x0-y0)*(x0-y0) + (x1-y1)*(x1-y1) );
            }else{
                x0 = y0; 
                x1 = y1;
                distance = 0;
            }
        }
        else {// y0 == 0
            x0 = 0; x1 = e1; distance = std::abs( y1 - e1 );
        }
    }else{ // y1 == 0
        Real numer0 = e0*y0 , denom0 = e0*e0 - e1*e1;
        if ( numer0 < denom0 ){
            Real xde0 = numer0/denom0;
            x0 = e0*xde0 ; x1 = e1*sqrt(1 - xde0*xde0 );
            distance = sqrt( (x0-y0)*(x0-y0) + x1*x1 );
        }else{
            x0 = e0; 
            x1 = 0; 
            distance = std::abs( y0 - e0 );
        }
    }
    return distance;
}

int main(int argc, const char *argv[]) {
    if (!((argc == 4) || (argc == 6))) {
        std::cerr << "usage: ./bisection a b max_iterations [x y]" << std::endl;
        exit(-1);
    }

    double a = std::stod(argv[1]),
           b = std::stod(argv[2]);
    maxIterations = std::stoi(argv[3]);

    double aspect = std::min<Real>(a / b, 3.0);
    size_t gridSizeX = 1024,
           gridSizeY = ceil(gridSizeX / aspect);
    if (gridSizeY < 2) throw std::runtime_error("Invalid aspect ratio");
    double gridSpacing = 1.5 * a / (gridSizeX - 1);
    // For now, always sample square domain--makes creating animations easier.
    gridSizeY = gridSizeX;

    std::cout << std::setprecision(19);
    std::cerr << std::setprecision(19);

    if (argc > 4) {
        double x = std::stod(argv[4]),
               y = std::stod(argv[5]);
        Real p_x, p_y;
        Real dist = DistancePointEllipse(a, b, x, y, p_x, p_y);
        std::cout << x
                << '\t' << y
                << '\t' << p_x
                << '\t' << p_y
                << '\t' << dist
                << '\n';
        return 0;
    }

#if OUTPUT_FILE
    std::ofstream iterateFile("bisection_iterates.txt");
    iterateFile << std::setprecision(19);
#endif

    for (size_t xi = 0; xi < gridSizeX; ++xi) {
        for (size_t yi = 0; yi < gridSizeY; ++yi) {
            double x = xi * gridSpacing,
                   y = yi * gridSpacing;
            Real p_x, p_y;
            volatile Real dist = DistancePointEllipse(a, b, x, y, p_x, p_y);
#if OUTPUT_FILE
            iterateFile << x
                << '\t' << y
                << '\t' << p_x
                << '\t' << p_y
                << '\t' << dist
                << '\n';
#endif
        }
#if OUTPUT_FILE
        iterateFile  << '\n'; // gnuplot wants blank lines between x values
#endif
    }
    return 0;
}
