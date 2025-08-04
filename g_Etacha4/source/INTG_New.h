// INTG_New.h
#pragma once

#include <functional>
#include <vector>

namespace quad_new {

/// Adaptive Gauss–Kronrod integrator ported from INTG.FOR
/// @param F       integrand function
/// @param BOUND   integration bound or mapping parameter
/// @param INF     0=fixed [0,BOUND], 1=semi-infinite, 2=full infinite
/// @param EPSABS  absolute error tolerance
/// @param EPSREL  relative error tolerance
/// @param RESULT  computed integral (output)
/// @param ABSERR  estimated absolute error (output)
/// @param NEVAL   number of function evaluations (output)
/// @param IER     error code (0=success, >0=failure)
/// @param LIMIT   maximum number of subintervals (default 1000)
void INTG_New(
    const std::function<double(double)>& F,
    double BOUND,
    int INF,
    double EPSABS,
    double EPSREL,
    double& RESULT,
    double& ABSERR,
    int& NEVAL,
    int& IER,
    int LIMIT = 1000
    );

} // namespace quad_new
