////////////////////////////////////////////////////////////////////////////////
// AutomaticDifferentiation.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Includes and functions needed for automaic differentiation (using Adept)
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 15:26:43
////////////////////////////////////////////////////////////////////////////////
#ifndef AUTOMATICDIFFERENTIATION_HH
#define AUTOMATICDIFFERENTIATION_HH

#include <adept.h>
#include <algorithm>

// A note on Eigen's norm() vs squaredNorm():
// Adept's sqrt overload is not visible to Eigen, so we must use
// sqrt(x.squaredNorm()) (or include Adept before Eigen).

#endif /* end of include guard: AUTOMATICDIFFERENTIATION_HH */