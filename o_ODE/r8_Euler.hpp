#ifndef R8_EULER_HPP
#define R8_EULER_HPP
double integrate_euler(
    double t0,
    double t1,
    double *y,
    int     neqn,
    double  h_init,
    double  relTol,
    double  absTol
    );
#endif // R8_EULER_HPP
