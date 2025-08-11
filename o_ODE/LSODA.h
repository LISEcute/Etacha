#ifndef LSODE_H
#define LSODE_H

#include <array>
#include <cmath>
#include <memory>
#include <vector>

using real = double;
using namespace std;

/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Type definition of LSODA ode system. See the file test_LSODA.cpp
 * for an example.
 *
 * @Param time, double
 * @Param y, array of double.
 * @Param dydt, array of double
 * @Param data, void*
 *
 * @Returns void
 */
/* ----------------------------------------------------------------------------*/
typedef void (*LSODA_ODE_SYSTEM_TYPE)(real t, real* y, real* dydt, void*);

class LSODA {

public:
    LSODA();
    ~LSODA();

    size_t idamax1(
        const vector<real>& dx, const size_t n, const size_t offset);

    void dscal1(const real da, vector<real>& dx, const size_t n,
        const size_t offset);

    real ddot1(const vector<real>& a, const vector<real>& b,
        const size_t n, const size_t offsetA, const size_t offsetB);

    void daxpy1(const real da, const vector<real>& dx, vector<real>& dy,
        const size_t n, const size_t offsetX, const size_t offsetY);

    void dgesl(const vector<vector<real>>& a, const size_t n,
               vector<int>& ipvt, vector<real>& b, const size_t job);


    void dgefa(vector<vector<real>>& a, const size_t n, vector<int>& ipvt,
        size_t* const info);
    void dgbfa(vector<vector<real>>& a, const size_t lda, const size_t n,
               const size_t ml, const size_t mu, vector<int>& ipvt, size_t* info);
    void dgbsl(vector<vector<real>>& a, const size_t lda, const size_t n,
               const size_t ml, const size_t mu, vector<int>& ipvt, vector<real>& b,
               const size_t job);

    void prja(const size_t neq, vector<real>& y, LSODA_ODE_SYSTEM_TYPE f,
        void* _data);

    void lsoda(LSODA_ODE_SYSTEM_TYPE f, const size_t neq, vector<real>& y,
                real* t, real tout, int itask, int* istate, int iopt, int jt,
                array<int, 7>& iworks, array<real, 4>& rworks, void* _data);


    void correction(const size_t neq, vector<real>& y,
                    LSODA_ODE_SYSTEM_TYPE f, size_t* corflag, real pnorm, real* del,
                    real* delp, real* told, size_t* ncf, real* rh, size_t* m,
        void* _data);

    void stoda(const size_t neq, vector<real>& y, LSODA_ODE_SYSTEM_TYPE f,
        void* _data);

    // We call this function in VoxelPools::
    void lsoda_update(
        LSODA_ODE_SYSTEM_TYPE f, const size_t neq, vector<real>& y,
        std::vector<real>& yout, real* t, const real tout, int* istate,
        void* const _data, real rtol = 1e-6f, real atol = 1e-6f // Tolerance
    );

    void terminate(int* istate);
    void terminate2(vector<real>& y, real* t);
    void successreturn(vector<real>& y, real* t, int itask, int ihit,
                       real tcrit, int* istate);
    void _freevectors(void);
    void ewset(const vector<real>& ycur);
    void resetcoeff(void);
    void solsy(vector<real>& y);
    void endstoda(void);
    void orderswitch(
        real* rhup, real dsm, real* pdh, real* rh, size_t* orderflag);
    void intdy(real t, int k, vector<real>& dky, int* iflag);
    void corfailure(real* told, real* rh, size_t* ncf, size_t* corflag);
    void methodswitch(real dsm, real pnorm, real* pdh, real* rh);
    void cfode(int meth_);
    void scaleh(real* rh, real* pdh);
    real fnorm(
        int n, const vector<vector<real>>& a, const vector<real>& w);
    real vmnorm(
        const size_t n, const vector<real>& v, const vector<real>& w);
    static bool abs_compare(real a, real b);

private:
    size_t ml, mu, imxer;
    real sqrteta;

    // NOTE: initialize in default constructor. Older compiler e.g. 4.8.4 would
    // produce error if these are initialized here. With newer compiler,
    // initialization can be done here.
    array<size_t, 3> mord;
    array<real, 13> sm1;

    array<real, 14> el; // = {0};
    array<real, 13> cm1; // = {0};
    array<real, 6> cm2; // = {0};


    array<array<real, 14>, 13> elco;
    array<array<real, 4>, 13> tesco;

    size_t illin, init, ierpj, iersl, jcur, l, miter, maxord, maxcor, msbp,
        mxncf;

    int kflag, jstart;

    size_t ixpr = 0, jtyp, mused, mxordn = 12, mxords = 5;
    size_t meth_;

    size_t n, nq, nst, nfe, nje, nqu;
    size_t mxstep, mxhnil;
    size_t nslast, nhnil, ntrep, nyh;

    real ccmax, el0, h_ = .0;
    real hmin, hmxi, hu, rc, tn_ = 0.0;
    real tsw, pdnorm;
    real conit, crate, hold, rmax;

    size_t ialth, ipup, lmax;
    size_t nslp;
    real pdest, pdlast, ratio;
    int icount, irflag;

    vector<real> ewt;
    vector<real> savf;
    vector<real> acor;
    vector<vector<real>> yh_;
    vector<vector<real>> wm_;

    vector<int> ipvt;

private:
    int itol_ = 2;
    std::vector<real> rtol_;
    std::vector<real> atol_;

public:
    void* param = nullptr;
};

#endif /* end of include guard: LSODE_H */
