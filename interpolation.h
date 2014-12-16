#ifndef INTERPOLATION
#define INTERPOLATION
#include<gsl/gsl_spline.h>
class Interpolation
{
	public:

		int n;

		double *xt; 
		double *yt; 

		int flag_init;

		Interpolation(void);
		~Interpolation(void);

		void Initialize(double *x_in, double *y_in, int n_in);
		void ShowTable(void);
		void FreeMemory(void);

		gsl_spline *spline_x_y;
		gsl_interp_accel *acc_x_y;
		gsl_interp *linear_x_y;
		gsl_interp_accel *acc_l_x_y;

		double *MakeLinearTable(int n_in, double x_min, double x_max);
		double *MakeLog10Table(int n_in, double x_min, double x_max);


		//return the interpolated value
		double y(double x);

	
};
#endif //INTERPOLATION
