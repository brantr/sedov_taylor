#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include "routines.h"
Interpolation::Interpolation(void)
{
	flag_init = 0;
}
Interpolation::~Interpolation(void)
{
	FreeMemory();
}
double Interpolation::y(double x)
{
	size_t ix;
	if(!flag_init)	
	{
		printf("ERROR! must initialize!\n");
	}
	//printf("itnerp x %e xt[0] %e yt[n-1] %e\n",x,xt[0],xt[n-1]);
	if(x<=xt[0])
		return yt[0];
	if(x>xt[n-1])
		return yt[n-1];

	ix = gsl_interp_bsearch(xt,x,0,n-2);
	return (yt[ix+1]-yt[ix])*(x-xt[ix])/(xt[ix+1]-xt[ix]) + yt[ix];
	//return gsl_spline_eval(spline_x_y,x,acc_x_y);
	//return gsl_interp_eval(linear_x_y,xt,yt,x,acc_l_x_y);
}
void Interpolation::ShowTable(void)
{
	for(int i=0;i<n;i++)
		printf("x %e y %e\n",xt[i],yt[i]);
}
void Interpolation::Initialize(double *x_in, double *y_in, int n_in)
{

	n = n_in;
	xt = calloc_double_array(n);
	yt = calloc_double_array(n);

	for(int i=0;i<n;i++)
	{
		xt[i] = x_in[i];
		yt[i] = y_in[i];
		//printf("xt %e yt %e\n",xt[i],yt[i]);
	}


	acc_x_y    = gsl_interp_accel_alloc();
	spline_x_y = gsl_spline_alloc(gsl_interp_akima,n);
	gsl_spline_init(spline_x_y,xt,yt,n);

	//acc_l_x_y = gsl_interp_accel_alloc();
	//linear_x_y = gsl_interp_alloc(gsl_interp_linear,n);
	//gsl_interp_init(linear_x_y,xt,yt,n);	

	flag_init = 1;
}
double *Interpolation::MakeLog10Table(int n_in, double x_min, double x_max)
{
	double *x = calloc_double_array(n_in);

	for(int i=0;i<n_in;i++)
		x[i] = pow(10.0, (log10(x_max)-log10(x_min))*((double) i)/((double) (n_in-1)) + log10(x_min));

	return x;
}
double *Interpolation::MakeLinearTable(int n_in, double x_min, double x_max)
{
	double *x = calloc_double_array(n_in);

	for(int i=0;i<n_in;i++)
		x[i] = (x_max-x_min)*((double) i)/((double) (n_in-1)) + x_min;

	return x;
}
void Interpolation::FreeMemory(void)
{
	if(flag_init)
	{
		free(xt);
		free(yt);

		//gsl_spline_free(spline_x_y);
		//gsl_interp_accel_free(acc_x_y);
		gsl_interp_free(linear_x_y);
		gsl_interp_accel_free(acc_l_x_y);
	}
}
