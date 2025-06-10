#include "o_ODE/r8_Euler.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

extern int ODEsteps;
extern void F(double t, double *y, double *yp);

// Tuning constants:
static constexpr double SAFETY  = 0.3;
static constexpr double FAC_MIN = 0.1;
static constexpr double FAC_MAX = 1.0;
static constexpr int    MAX_STEPS = 500000;

double integrate_euler(
    double t0,
    double t1,
    double *y,
    int     neqn,
    double  h_init,
    double  relTol,
    double  absTol)
{
    std::vector<double> yp(neqn);
    std::vector<double> yTent(neqn);
    std::vector<double> yMid(neqn);
    std::vector<double> yHalf(neqn);

    double t = t0;
    double h = (t1 > t0 ? std::fabs(h_init) : -std::fabs(h_init));
    int stepCount = 0;
    ODEsteps = 0;

    while ((t1 - t) * h > 0.0) {
        if (++stepCount > MAX_STEPS) {
            throw std::runtime_error("integrate_euler: too many steps");
        }
        ODEsteps = stepCount;

        if (std::fabs(h) > std::fabs(t1 - t)) {
            h = t1 - t;
        }

        F(t, y, yp.data());

        // Tentative full-step Euler
        for (int i = 0; i < neqn; ++i) {
            yTent[i] = y[i] + h * yp[i];
        }

        // Midpoint (half-step Euler)
        double h2 = 0.5 * h;
        for (int i = 0; i < neqn; ++i) {
            yMid[i] = y[i] + h2 * yp[i];
        }

        F(t + h2, yMid.data(), yp.data());

        for (int i = 0; i < neqn; ++i) {
            yHalf[i] = yMid[i] + h2 * yp[i];
        }

        // Error estimation
        double errSum = 0.0;
        for (int i = 0; i < neqn; ++i) {
            double sc = absTol + relTol * std::max(std::fabs(yTent[i]), std::fabs(yHalf[i]));
            double delta = yTent[i] - yHalf[i];
            errSum += (delta * delta) / (sc * sc);
        }
        double err = std::sqrt(errSum / static_cast<double>(neqn));

        if (err <= 1.0) {
            for (int i = 0; i < neqn; ++i) {
                y[i] = yHalf[i];
            }
            t += h;
        }

        double fac;
        if (err == 0.0) {
            fac = FAC_MAX;
        } else {
            fac = SAFETY * std::exp(-0.5 * std::log(err));
            fac = std::clamp(fac, FAC_MIN, FAC_MAX);
        }

        h *= fac;
    }

    return t;
}
