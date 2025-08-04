#include "INTG_New.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace quad_new {

// Machine constants (from Fortran R2MACH)
static double r2mach(int i) {
    static const double RMACH[6] = {0.0, 1.18e-38, 3.40e+38, 0.595e-07, 1.19e-07, 0.30102999566};
    return (i >= 1 && i <= 5) ? RMACH[i] : 0.0;
}

// Error handler
static void xerreur(const char* msg, int ier) {
    std::cerr << "INTG_New error (IER=" << ier << "): " << msg << std::endl;
}

// 15-point Gauss–Kronrod rule (QK15I)
static void qk15i(
    const std::function<double(double)>& F,
    double BOUND, int INF,
    double A, double B,
    double& RESULT, double& ABSERR,
    double& RESABS, double& RESASC
    ) {
    static const double XGK[8] = { 0.0000000000000000,
                                  0.2077849550078985, 0.4058451513773972,
                                  0.5860872354676911, 0.7415311855993944,
                                  0.8648644233597691, 0.9491079123427585,
                                  0.9914553711208126 };
    static const double WGK[8] = { 0.2094821410847278,
                                  0.2044329400752989, 0.1903505780647854,
                                  0.1690047266392679, 0.1406532597155259,
                                  0.1047900103222502, 0.0630920926299786,
                                  0.0229353220105292 };
    static const double WG[4] = { 0.4179591836734694,
                                 0.3818300505051189, 0.2797053914892767,
                                 0.1294849661688697 };

    double epmach = r2mach(4), uflow = r2mach(1);
    double half = 0.5 * (B - A), centre = 0.5 * (A + B);
    RESULT = 0.0; RESABS = 0.0; RESASC = 0.0;

    auto evalF = [&](double t) {
        double x = centre + half * t;
        if (INF != 0) {
            double dinf = (INF == 2 ? 1.0 : 1.0);
            double mapped = BOUND + dinf * (1.0 - x) / x;
            x = mapped;
        }
        return F(x);
    };

    // central point
    double fv0 = evalF(0.0);
    double resg = WG[0] * fv0;
    double resk = WGK[0] * fv0;
    RESABS += WGK[0] * std::fabs(fv0);

    // other Kronrod/Gauss points
    for (int j = 1; j < 8; ++j) {
        double f1 = evalF(XGK[j]);
        double f2 = evalF(-XGK[j]);
        double sum = f1 + f2;
        resk += WGK[j] * sum;
        if (j % 2 == 1) {
            resg += WG[(j - 1) / 2] * sum;
        }
        RESABS += WGK[j] * (std::fabs(f1) + std::fabs(f2));
    }

    RESULT = resk * half;
    double diff = (resk - resg) * half;
    ABSERR = std::fabs(diff);

    // compute RESASC
    double reskh = resk * 0.5;
    auto accum = [&](double t) {
        double fv = evalF(t);
        return std::fabs(fv - reskh);
    };
    RESASC += WGK[0] * accum(0.0);
    for (int j = 1; j < 8; ++j) {
        RESASC += WGK[j] * (accum(XGK[j]) + accum(-XGK[j]));
    }
    RESASC *= half;

    // refine error
    if (RESASC > 0.0 && ABSERR > 0.0) {
        double ratio = std::pow((200.0 * ABSERR / RESASC), 1.5);
        ABSERR = std::max(ABSERR, RESASC * std::min(1.0, ratio));
    }
    ABSERR = std::max(ABSERR, epmach * 50.0 * RESABS);
}

// Extrapolation (EA)
static void ea(bool &NEWFLG, double SVALUE, int LIMEXP,
               double &RESULT, double &ABSERR,
               std::vector<double> &EPSTAB, int &IERR) {
    if (LIMEXP < 3) {
        IERR = 1;
        xerreur("EA: LIMEXP<3", 1);
        return;
    }
    IERR = 0;
    static double RES3LA[3];
    if (NEWFLG) {
        EPSTAB[0] = SVALUE;
        RESULT = SVALUE;
        ABSERR = std::fabs(SVALUE);
        NEWFLG = false;
        return;
    }
    int N = static_cast<int>(EPSTAB[LIMEXP + 2]);
    int NRES = static_cast<int>(EPSTAB[LIMEXP + 3]);
    EPSTAB[N] = SVALUE;
    // simple 3-term epsilon algorithm
    if (N >= 2) {
        double e0 = EPSTAB[N-2];
        double e1 = EPSTAB[N-1];
        double e2 = EPSTAB[N];
        double delta1 = e1 - e0;
        double delta2 = e2 - e1;
        double ss = 1.0 / delta1 + 1.0 / delta2;
        double res = e1 + 1.0/ss;
        RESULT = res;
        ABSERR = std::fabs(res - e2) + std::fabs(e2 - e1);
    }
    EPSTAB[LIMEXP + 2] = N + 1;
}

// Priority-queue sort (QPSRT)
static void qpsrt(int LIMIT, int LAST, int &MAXERR, double &ERMAX,
                  std::vector<double> &ELIST, std::vector<int> &IORD, int &NRMAX) {
    // simple bubble to maintain IORD sorted by ELIST desc
    for (int i = 0; i < LAST; ++i) IORD[i] = i;
    std::sort(IORD.begin(), IORD.begin() + LAST,
              [&](int a, int b){ return ELIST[a] > ELIST[b]; });
    MAXERR = IORD[0];
    ERMAX = ELIST[MAXERR];
    NRMAX = LAST;
}

// Adaptive engine (INTGE)
static void intge(
    const std::function<double(double)>& F,
    double BOUND, int INF,
    double EPSABS, double EPSREL,
    int LIMIT,
    double &RESULT, double &ABSERR,
    int &NEVAL, int &IER,
    std::vector<double> &ALIST,
    std::vector<double> &BLIST,
    std::vector<double> &RLIST,
    std::vector<double> &ELIST,
    std::vector<int>    &IORD,
    int &LAST
    ) {
    // initial subinterval [0,1]
    double resabs, resasc;
    qk15i(F, BOUND, INF, 0.0, 1.0, RESULT, ABSERR, resabs, resasc);
    NEVAL = 15;
    LAST = 1;
    ALIST[0] = 0.0;  BLIST[0] = 1.0;
    RLIST[0] = RESULT;  ELIST[0] = ABSERR;
    IORD[0]  = 0;

    double area = RESULT;
    double errsum = ABSERR;
    double errbnd = std::max(EPSABS, EPSREL * std::fabs(area));
    if (errsum <= errbnd) { IER = 0; return; }
    if (LIMIT == 1) { IER = 1; return; }

    const int LIMEXP = 50;
    std::vector<double> epstab(LIMEXP+5);
    bool newflg = true;
    int ierr_ea = 0;

    for (LAST = 2; LAST <= LIMIT; ++LAST) {
        // find subinterval with max error
        int maxidx = 0;
        for (int i = 1; i < LAST; ++i) {
            if (ELIST[i] > ELIST[maxidx]) maxidx = i;
        }
        double a = ALIST[maxidx];
        double b = BLIST[maxidx];
        double mid = 0.5*(a + b);
        double r1,e1,r2,e2,dumy;
        qk15i(F, BOUND, INF, a, mid, r1, e1, dumy, dumy);
        qk15i(F, BOUND, INF, mid, b, r2, e2, dumy, dumy);
        RLIST[maxidx] = r1; RLIST[LAST-1] = r2;
        ELIST[maxidx] = e1; ELIST[LAST-1] = e2;
        ALIST[LAST-1] = mid; BLIST[LAST-1] = b;
        BLIST[maxidx] = mid;

        area  = area - (r1 + r2) + r1 + r2;
        errsum= errsum - ELIST[maxidx] + (e1 + e2);
        errbnd = std::max(EPSABS, EPSREL * std::fabs(area));
        NEVAL += 30;

        if (errsum <= errbnd) {
            RESULT = area;
            ABSERR = errsum;
            IER = 0;
            return;
        }
        if (LAST == LIMIT) { IER = 1; return; }

        // optional extrapolation
        if (LAST <= LIMEXP && ierr_ea == 0) {
            double ea_res, ea_err;
            ea(newflg, area, LAST, ea_res, ea_err, epstab, ierr_ea);
            if (ea_err < errsum) {
                area = ea_res;
                errsum = ea_err;
            }
        }
    }
    // if here, failed to converge
    RESULT = area;
    ABSERR = errsum;
    IER = 1;
}

// Public entry
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
    int LIMIT
    ) {
    std::vector<double> ALIST(LIMIT), BLIST(LIMIT), RLIST(LIMIT), ELIST(LIMIT);
    std::vector<int> IORD(LIMIT);
    int LAST = 0;
    intge(F, BOUND, INF, EPSABS, EPSREL,
          LIMIT, RESULT, ABSERR, NEVAL, IER,
          ALIST, BLIST, RLIST, ELIST, IORD, LAST);
    if (IER != 0) {
        xerreur("INTG_New did not converge", IER);
    }
}

} // namespace quad_new
