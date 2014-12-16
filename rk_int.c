/*! \file rk_int.c
 *  \brief Function definitions for numerical integration routines.
 */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_sf_gamma.h>
#include"rk_int.h"
#include"routines.h"
/********************************************************
*
*     Numerical Integration Subroutines.
*
*********************************************************/


/*! \fn double integrate(double (*FUNC)(double,void*), void *fp ,int np,double a,double b,double dxinit, double eps)
 *  \brief Numerical integration routine using 5th-order Runge-Kutta.
 *
 *   Quadrature using fifth order Runge-Kutta with adaptive step size.
 *   Based on Press et al, Numerical Recipes in C, 2nd ed, pp 719-722.
 * 
 *   Runge-Kutta driver with adaptive stepsize control.  Integrate starting
 *   value y from a to b with accuracy eps, storing intermediate results in
 *   global variables.  dxinit should be set as a guessed first stepsize.
 * 
 *   Pass a second parameter to FUNC in fparm.
 *  
 *   Original fortan routine by M.A.K. Gross, C implementation by Brant Robertson
 *  
 *   func is the function to be integrated.
 *   parameters are passed in fp(np) array.
 *   func(x,xp,np).
 */
double integrate(double (*FUNC)(double,void*), void *fp ,int np,double a,double b,double dxinit, double eps)
{
	int maxsteps=10000000;
      
	double x, dx, dxnext, y, dydx, yscale;
      
	int Nstep;

      
	x     = a;
	dx    = dxinit;
	y     = 0.0;
	Nstep = 0;

      
	do 
	{
        
		Nstep = Nstep + 1;
		dydx = FUNC(x,fp);

		//yscale is the scaling used to monitor accuracy.  This general-purpose
		//choice can be modified if need be.
		yscale = fmax(fabs(y) + fabs(dx*dydx), 1.e-12);

		if ((x+dx-b)*(x+dx-a)>0.0)  //! If stepsize overshoots, decrease it.
			dx = b - x;
        
		RUNGE5VAR(&y,dydx,&x,dx,eps,yscale,&dxnext,FUNC,fp);
        
		dx = dxnext;
	}while (((x-b)*(b-a)<0.0) && (Nstep<maxsteps));

      
	if (Nstep>=maxsteps)
	{
		printf("Failed to converge in integral!\n");
		exit(-1);
	}
	return y;
}


/*! \fn void RUNGE5VAR(double *y,double dydx,double *x,double htry,double eps,double yscale,double *hnext,double (*DERIVS)(double,void*), void *fp)
 *  \brief Runge-Kutta step in numerical integration routine
 *
 *
 *   Fifth-order Runge-Kutta step with monitoring of local truncation error
 *   to ensure accuracy and adjust stepsize.  Input are the dependent
 *   variable y and its derivative dydx at the starting value of the
 *   independent variable x.  Also input are the stepsize to be attempted
 *   htry, the required accuracy eps, and the value yscale, against which the
 *   error is scaled.  On output, y and x are replaced by their new values.
 *   hdid is the stepsize that was actually accomplished, and hnext is the
 *   estimated next stepsize.  DERIVS is the user-supplied routine that
 *   computes right-hand-side derivatives.  The argument fparm is for an
 *   optional second argument to DERIVS (NOT integrated over).
 * 
 * 
 *   Original fortran by M.A.K. Gross, c implementation by Brant Robertson
 */
void RUNGE5VAR(double *y,double dydx,double *x,double htry,double eps,double yscale,double *hnext,double (*DERIVS)(double,void*), void *fp)
{
      
	//external DERIVS
	double errmax,h,hold,htemp,xnew,yerr,ytemp;
	double safety=0.9;
	double pgrow=-0.2;
	double pshrink=-0.25;
	double errcon=1.89e-4;
      
	yerr = 0.0;
      
	h = htry;  //! Set stepsize to initial accuracy.
      
	errmax = 10.0;
      
	do 
	{
		RUNGE(*y,dydx,*x,h,&ytemp,&yerr,DERIVS,fp);

		errmax = fabs(yerr/yscale)/eps;// ! Scale relative to required accuracy.
         
		if (errmax>1.0)
		{
			//! Truncation error too large; reduce h
          
			htemp = safety*h*pow(errmax,pshrink);
			hold = h;
			//h = sign(fmax(fabs(htemp),0.1*fabs(h)),h);  //! No more than factor of 10
			if(h<0.0)
			{
			
				//! No more than factor of 10
				h = -1.0*fmax(fabs(htemp),0.1*fabs(h));
			}else{
				//! No more than factor of 10
				h = 1.0*fmax(fabs(htemp),0.1*fabs(h));
			}
			xnew = *x + h;
          
			if (xnew == *x)
			{
				h = hold;
				errmax = 0.0;
			}
		}
	}while(errmax>1.0);
      
	if (errmax>errcon)
	{
        
		*hnext = safety*h*pow(errmax,pgrow);
	}else{
        
		*hnext = 5.0 * h;//! No more than factor of 5 increase.
	}
      
	*x = *x + h;
      
	*y = ytemp;
}

/*! \fn void RUNGE(double y,double dydx,double x,double h,double *yout,double *yerr,double (*DERIVS)(double,void*),void *fp)
 *  \brief Function to advance the RK solution in the numerical integration.
 *
 *  Given values for a variable y and its derivative dydx known at x, use
 *  the fifth-order Cash-Karp Runge-Kutta method to advance the solution
 *  over an interval h and return the incremented variables as yout.  Also
 *  return an estimate of the local truncation error in yout using the
 *  embedded fourth order method.  The user supplies the routine
 *  DERIVS(x,y,dydx), which returns derivatives dydx at x.
 * 
 *  Original fortran by M.A.K. Gross, c implementation by Brant Robertson.
 */
void RUNGE(double y,double dydx,double x,double h,double *yout,double *yerr,double (*DERIVS)(double,void*),void *fp)
{
      
	double ak3, ak4, ak5 ,ak6;
	double a2,a3,a4,a5,a6;
	double c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6;

      
	a2  =    0.2;
	a3  =    0.3;
	a4  =    0.6;
	a5  =    1.0;
	a6  =    0.875;
	c1  =   37.0/378.0;
	c3  =  250.0/621.0;
	c4  =  125.0/594.0;
	c6  =  512.0/1771.0;
	dc1 = c1 -  2825.0/27648.0;
	dc3 = c3 - 18575.0/48384.0;
	dc4 = c4 - 13525.0/55296.0;
	dc5 = -277.0/14336.0;
	dc6 = c6 -     0.25;

      
	ak3 = DERIVS(x+a3*h,fp);
	ak4 = DERIVS(x+a4*h,fp);
	ak5 = DERIVS(x+a5*h,fp);
	ak6 = DERIVS(x+a6*h,fp);

	//Estimate the fifth order value.
      
	*yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6);

	//Estimate error as difference between fourth and fifth order
      
	*yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6);
}
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/

/*! \fn double midpoint_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level)
 *  \brief Midpoint rule integration
 *
 *   mid point rule integration, see karniadakis and kirby, s4.2
 */
double midpoint_rule_integration(double(*func)(double, void*), void *fp, int np, double a, double b, int level)
{
	int nsteps = (int) pow(2.0,level)-1;
	double h   = (b-a)/pow(2.0,level);
	double sum = 0.0;

	for(int i=0;i<=nsteps;i++)
		sum += func( a + (i+0.5)*h, fp);
	sum*=h;

	return sum;
}

/*! \fn double trapezoid_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);
 *  \brief Trapezoid rule integration
 *
 *  trapezoid rule integration, see karniadakis and kirby, s4.2
 */
double trapezoid_rule_integration(double(*func)(double, void*), void *fp, int np, double a, double b, int level)
{
	int nsteps = (int) pow(2.0,level)-1;
	double h   = (b-a)/pow(2.0,level);
	double sum = 0.0;

	for(int i=1;i<=nsteps;i++)
		sum += func( a + i*h, fp);
	sum*=2;

	//add the first and the last point to the sum

	sum += func(a,fp) + func(b,fp);

	sum*= 0.5*h;

	return sum;
}

/*! \fn double romberg_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int m, int k);
 *  \brief Romberg integration
 *
 *  Romberg integration, see karniadakis and kirby, s4.2 
 */
double romberg_integration(double(*func)(double, void*), void *fp, int np, double a, double b, int m, int k)
{
	double RI, I1, I2;

	double coeff = pow(4.0,m);

	if(k<m)
	{
		printf("in romberg integration, k must be >=m, but k=%d, m=%d; setting k=m\n",k,m);	
		fflush(stdout);
		k=m;
	}

	if(m==0)
	{
		RI = trapezoid_rule_integration(func, fp, np, a, b, k);
	}else{
		I1 = romberg_integration(func, fp, np, a, b, m-1, k); 
		I2 = romberg_integration(func, fp, np, a, b, m-1, k-1); 
		RI = (coeff*I1 - I2)/(coeff-1.0);
	}
	return RI;
}

/*! \fn double double jacobi_poly(double x, double alpha, double beta, int degree)
 *  \brief Function to calculate the Jacobi polynomials.
 */
double jacobi_poly(double x, double alpha, double beta, int degree)
{
	double value;
	double tmp, degm1;

	double a1=0.,a2=0.,a3=0.,a4=0.;

	switch(degree)
	{
		case 0:
			value=1.0;
			break;

		case 1:
			value = 0.5*(alpha-beta+(alpha+beta+2.0)*x);
			break;

		default:
			degm1 = degree-1.0;
			tmp   = 2.0*degm1 + alpha + beta;
			a1    = 2.0*(degm1+1)*(degm1+alpha+beta+1)*tmp;
			a2    = (tmp+1)*(alpha*alpha - beta*beta);
			a3    = tmp*(tmp+1.0)*(tmp+2.0);
			a4    = 2.0*(degm1+alpha)*(degm1+beta)*(tmp+2.0);

			value = ((a2+a3*x)*jacobi_poly(x,alpha,beta,degree-1) - a4*jacobi_poly(x,alpha,beta,degree-2))/a1;
	}

	return value;
}

/*! \fn double jacobi_poly_deriv(double x, double alpha, double beta, int degree)
 *  \brief Function to calculate the derivative of Jacobi polynomials.
 */
double jacobi_poly_deriv(double x, double alpha, double beta, int degree)
{
	double value;

	double tmp;

	double b1,b2,b3;

	switch(degree)
	{
		case 0:
			value = 0.0;
			break;

		default:
			tmp = 2.0*degree + alpha + beta;

			b1  = tmp*(1.0-x*x);
	
			b2  = degree*(alpha-beta-tmp*x);

			b3  = 2.0*(degree+alpha)*(degree+beta);

			value = (b2*jacobi_poly(x,alpha,beta,degree) + b3*jacobi_poly(x,alpha,beta,degree-1))/b1;
	}

	return value;
}

/*! \fn void jacobi_zeros(double *z, double alpha, double beta, int degree)
 *  \brief Function to find zeros of Jacobi polynomials. 
 */
void jacobi_zeros(double *z, double alpha, double beta, int degree)
{
	int i,j,k;

	const int maxit=30;

	const double EPS = 1.0e-14;
	
	double dth = M_PI/(2.0*degree);

	double poly, pder, rlast=0.0;

	double sum, delr, r;

	double one = 1.0;
	double two = 2.0;


	//if the poly is degree zero or less, then no roots

	if(degree<=0)
		return;

	for(k=0;k<degree;k++)
	{
		r = -cos((two*k + one)*dth);

		if(k)
			r = 0.5*(r+rlast);


		for(j=1;j<maxit;++j)
		{

			poly = jacobi_poly(r,alpha,beta,degree);
			pder = jacobi_poly_deriv(r,alpha,beta,degree);

			sum =0.0;

			for(i=0;i<k;i++)
				sum+= one/(r-z[i]);

			delr = -poly / (pder - sum*poly);

			r += delr;

			if( fabs(delr) < EPS) 
				break;
		}

		z[k] = r;

		rlast = r;

	}

	return;
}

/*! \fn void jacobi_zeros_and_weights(double *z, double *w, double alpha, double beta, int degree)
 *  \brief Function to find zeros of the Jacobi polynomials, and weights for Gauss quadrature int.
 */
void jacobi_zeros_and_weights(double *z, double *w, double alpha, double beta, int degree)
{

	int i;

	double fac, one=1.0, two=2.0, apb = alpha+beta;

	jacobi_zeros(z, alpha, beta, degree);

	for(i=0;i<degree;i++)
		w[i] = jacobi_poly_deriv(z[i],alpha,beta,degree);

	fac = pow(two,apb+one)*gsl_sf_gamma(alpha+degree+one)*gsl_sf_gamma(beta+degree+one);
	fac /= gsl_sf_gamma(degree+one)*gsl_sf_gamma(apb+degree+one);

	for(i=0;i<degree;i++)
		w[i] = fac/(w[i]*w[i]*(one-z[i]*z[i]));

	return;	
}

/*! \fn double gauss_quadrature_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int degree);
 *  \brief Gauss quadrature integration
 */
double gauss_quadrature_integration(double(*func)(double, void*), void *fp, int np, double a, double b, int degree)
{
	double sum=0.;
	double *x, *y, *w;

	x = calloc_double_array(degree);
	y = calloc_double_array(degree);
	w = calloc_double_array(degree);

	jacobi_zeros_and_weights(x,w,0,0,degree);


	for(int k=0;k<degree;k++)
	{
		//scale x coordinates

		x[k] = 0.5*(b-a) + 0.5*(b-a)*x[k] + a;

		//scale weights

		w[k] *= 0.5*(b-a);

		//calculate y[k]'s

		y[k] = func(x[k],fp);

		//sum partial integrand

		sum+=y[k]*w[k];
	}

	free(x);
	free(y);
	free(w);

	return sum;

}
