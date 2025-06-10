
// r8_bdf.cpp
#include "r8_bdf.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

extern void F(double t, double *V, double *VP);

double integrate_bdf1(
    double t0,
    double t1,
    double *V,
    int neqn,
    double h,
    double tol,
    int max_iter
    ) {
    int n = neqn;
    double t = t0;

    // adaptivity parameters
    const double safety  = 0.1;
    const double fac_min = 0.1;
    const double fac_max = 2.0;
    const int    order   = 1;   // BDF-1 is first order
    const double eps_fd  = 1e-6;

    std::vector<double> Vold(n), Vnew(n), VP(n), VP2(n),
        J(n), Rold(n), Rnew(n), dV(n);

    // initialize
    for (int i = 0; i < n; ++i)
        Vold[i] = Vnew[i] = V[i];

    // time‐stepping
    while (t + 1e-14 < t1) {
        double h_step = std::min(h, t1 - t);
        bool accepted = false;

        while (!accepted) {
            double tnext = t + h_step;
            double dt    = h_step;

            // reset Newton guess
            for (int i = 0; i < n; ++i)
                Vnew[i] = Vold[i];

            // build FD Jacobian diagonal & initial residual
            F(tnext, Vnew.data(), VP.data());
            for (int i = 0; i < n; ++i) {
                Rold[i] = Vold[i] + dt * VP[i] - Vnew[i];
                double tmp = Vnew[i];
                Vnew[i] = tmp + eps_fd;
                F(tnext, Vnew.data(), VP2.data());
                J[i]    = (VP2[i] - VP[i]) / eps_fd;
                Vnew[i] = tmp;
            }

            // Newton iterations
            for (int iter = 0; iter < max_iter; ++iter) {
                double maxR = 0.0;
                for (int i = 0; i < n; ++i) {
                    double denom = 1.0 - dt * J[i];
                    if (std::fabs(denom) < 1e-12)
                        denom = (denom < 0 ? -1e-12 : 1e-12);
                    dV[i]    = Rold[i] / denom;
                    Vnew[i] += dV[i];
                    maxR = std::max(maxR, std::fabs(Rold[i]));
                }
                if (maxR < tol) break;

                F(tnext, Vnew.data(), VP.data());
                for (int i = 0; i < n; ++i) {
                    Rnew[i]  = Vold[i] + dt * VP[i] - Vnew[i];
                    Rold[i]  = Rnew[i];
                }
            }

            // compute BE residual + scale it
            F(tnext, Vnew.data(), VP.data());
            double err = 0.0;
            for (int i = 0; i < n; ++i) {
                double resi  = (Vnew[i] - Vold[i]) / dt - VP[i];
                double scale = tol + tol * std::fabs(Vnew[i]);  // atol=tol, rtol=tol
                err = std::max(err, std::fabs(resi) / scale);
            }

            if (err <= 1.0) {
                // accept
                for (int i = 0; i < n; ++i) {
                    Vold[i] = Vnew[i];
                    V[i]    = Vnew[i];
                }
                t = tnext;
                // step‐size update (exponent = 1/(order+1) = 1/2)
                double fac = safety * std::pow(1.0 / (err + 1e-16), 1.0 / (order + 1.0));
                fac = std::min(std::max(fac, fac_min), fac_max);
                h   = h_step * fac;
                accepted = true;
            } else {
                // reject & shrink trial
                double fac = safety * std::pow(1.0 / (err + 1e-16), 1.0 / (order + 1.0));
                fac = std::min(std::max(fac, fac_min), 1.0);
                h_step *= fac;
            }
        }
    }
}
