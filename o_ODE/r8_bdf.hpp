
// r8_bdf.hpp
#ifndef R8_BDF_HPP
#define R8_BDF_HPP

#include "ode.hpp"  // extern void F(double t, double *V, double *VP);
// at top of r8_bdf.hpp

double integrate_bdf1(
    double t0,
    double t1,
    double *V,
    int neqn,
    double h,
    double tol     = 1e-5,
    int    max_iter= 5
    );

#endif // R8_BDF_HPP
