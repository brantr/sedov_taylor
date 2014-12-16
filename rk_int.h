/*! \file rk_int.h
 *  \brief Function declarations for numerical integration routines.
 */
#ifndef BRANT_RK_INT
#define BRANT_RK_INT

/*! \fn double integrate(double (*FUNC)(double,void*), void *fp ,int np,double a,double b,double dxinit, double eps)
 *  \brief Numerical integration routine using 5th-order Runge-Kutta.
 */
double integrate(double (*FUNC)(double,void*), void *fp ,int np,double a,double b,double dxinit, double eps);
/*! \fn void RUNGE5VAR(double *y,double dydx,double *x,double htry,double eps,double yscale,double *hnext,double (*DERIVS)(double,void*), void *fp)
 *  \brief Runge-Kutta step in numerical integration routine
 */
void RUNGE5VAR(double *y,double dydx,double *x,double htry,double eps,double yscale,double *hnext,double (*DERIVS)(double,void*), void *fp);
/*! \fn void RUNGE(double y,double dydx,double x,double h,double *yout,double *yerr,double (*DERIVS)(double,void*),void *fp)
 *  \brief Function to advance the RK solution in the numerical integration.
 */
void RUNGE(double y,double dydx,double x,double h,double *yout,double *yerr,double (*DERIVS)(double,void*),void *fp);



/*! \fn double midpoint_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level)
 *  \brief Midpoint rule integration
 */
double midpoint_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);

/*! \fn double trapezoid_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);
 *  \brief Trapezoid rule integration
 */
double trapezoid_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);

/*! \fn double romberg_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int m, int k);
 *  \brief Romberg integration
 */
double romberg_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int m, int k);

/*! \fn double gauss_quadrature_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int degree);
 *  \brief Gauss quadrature integration
 */
double gauss_quadrature_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int degree);


//jacobi polynomials

/*! \fn double double jacobi_poly(double x, double alpha, double beta, int degree)
 *  \brief Function to calculate the Jacobi polynomials.
 */
double jacobi_poly(double x, double alpha, double beta, int degree);

/*! \fn double jacobi_poly_deriv(double x, double alpha, double beta, int degree)
 *  \brief Function to calculate the derivative of Jacobi polynomials.
 */
double jacobi_poly_deriv(double x, double alpha, double beta, int degree);

/*! \fn void jacobi_zeros(double *z, double alpha, double beta, int degree)
 *  \brief Function to find zeros of Jacobi polynomials. 
 */
void   jacobi_zeros(double *z, double alpha, double beta, int degree);

/*! \fn void jacobi_zeros_and_weights(double *z, double *w, double alpha, double beta, int degree)
 *  \brief Function to find zeros of the Jacobi polynomials, and weights for Gauss quadrature int.
 */
void   jacobi_zeros_and_weights(double *z, double *w, double alpha, double beta, int degree);
#endif
