/*! \file routines.c
 *  \brief Function definitions for utility routines.
 */
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_spline.h>
#include<time.h>
#include"routines.h"


/*! \fn double double_log10_index(int i, int n, double xmin, double xmax)
 *  \brief Provides the i^th out of n log10 incremented value between xmin and xmax.
 *	   Useful for creating a ordinate array for an interpolation.
 */
double double_log10_index(int i, int n, double xmin, double xmax)
{
	return pow(10.0, (log10(xmax) - log10(xmin))*((double) i)/((double) (n-1)) + log10(xmin));
}
/*! \fn double double double_linear_index(int i, int n, double xmin, double xmax)
 *  \brief Provides the i^th out of n linear incremented value between xmin and xmax.
 *	   Useful for creating a ordinate array for an interpolation.
 */
double double_linear_index(int i, int n, double xmin, double xmax)
{
	return (xmax - xmin)*((double) i)/((double) (n-1)) + xmin;
}
/*! \fn void check_args(int argc, char **argv, int num_args)
 *  \brief Ensures that the number of command line arguements equals num_args.
 */
void check_args(int argc, char **argv, int num_args)
{
	if(argc!=num_args)
	{
		printf("Error: argc == %d, expected num_args == %d.\n",argc,num_args);
		fflush(stdout);
		exit(-1);
	}
}

/*! \fn FILE *fopen_brant(char fname[], const char *mode)
 *  \brief Safe method for opening a FILE pointer.
 */
FILE *fopen_brant(char fname[], const char *mode)
{
	FILE *fp;

	if(!(fp = fopen(fname,mode)))
	{
		printf("Error opening %s.\n",fname);
		fflush(stdout);
		exit(-1);
	}

	return fp;
}
/*! \fn double *calloc_double_array(int n)
 *  \brief Safe method for callocing a double array
 */
double *calloc_double_array(int n)
{
	double *f;

	if(!(f = (double *) calloc(n,sizeof(double))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(double)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}


/*! \fn float *calloc_float_array(int n)
 *  \brief Safe method for callocing a float array
 */
float *calloc_float_array(int n)
{
	float *f;

	if(!(f = (float *) calloc(n,sizeof(float))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(float)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}


/*! \fn int *calloc_int_array(int n)
 *  \brief Safe method for callocing an int array
 */
int *calloc_int_array(int n)
{
	int *f;

	if(!(f = (int *) calloc(n,sizeof(int))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(int)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}

/*! \fn size_t *calloc_size_t_array(int n)
 *  \brief Safe method for callocing a size_t array
 */
size_t *calloc_size_t_array(int n)
{
	size_t *f;

	if(!(f = (size_t *) calloc(n,sizeof(size_t))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(size_t)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}
/*! \fn double max_three(double a, double b, double c)
 *  \brief Returns the max of 3 numbers */
double max_three(double a, double b, double c)
{
	return GSL_MAX_DBL(GSL_MAX_DBL(a,b),c);
}

/*! \fn double **two_dimensional_array(int n, int l)
 *  \brief Allocate a two dimensional (n x l) array
 */
double **two_dimensional_array(int n, int l)
{
	double **x;

	x = new double *[n];
	for(int i=0;i<n;i++)
		x[i] = new double[l];

	return x;
}
/*! \fn void deallocate_two_dimensional_array(double **x, int n, int l)
 *  \brief De-allocate a two dimensional (n x l) array
 */
void deallocate_two_dimensional_array(double **x, int n, int l)
{
	for(int i=0;i<n;i++)
		delete[] x[i];
	delete x;
}
/*! \fn double ***three_dimensional_array(int n, int l, int m)
 *  \brief Allocate a three dimensional (n x l x m) array 
 */
double ***three_dimensional_array(int n, int l, int m)
{
	double ***x;

	x = new double **[n];
	for(int i=0;i<n;i++)
	{
		x[i] = new double *[l];
		for(int j=0;j<l;j++)
		{
			x[i][j] = new double [m];
		}
	}

	return x;
}
/*! \fn void deallocate_three_dimensional_array(double ***x, int n, int l, int m)
 *  \brief De-allocate a three dimensional (n x l x m) array
 */
void deallocate_three_dimensional_array(double ***x, int n, int l, int m)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<l;j++)
			delete[] x[i][j];
		delete[] x[i];
	}
	delete x;
}

/*! \fn int ***three_dimensional_int_array(int n, int l, int m)
 *  \brief Allocate a three dimensional (n x l x m) int array
 */
int ***three_dimensional_int_array(int n, int l, int m)
{
	int ***x;

	x = new int **[n];
	for(int i=0;i<n;i++)
	{
		x[i] = new int *[l];
		for(int j=0;j<l;j++)
		{
			x[i][j] = new int [m];
		}
	}

	return x;
}
/*! \fn void deallocate_three_int_dimensional_array(int ***x, int n, int l, int m)
 *  \brief De-allocate a three dimensional (n x l x m) int array.
 */
void deallocate_three_int_dimensional_array(int ***x, int n, int l, int m)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<l;j++)
			delete[] x[i][j];
		delete[] x[i];
	}
	delete x;
}


/*! \fn int compare_doubles(const void *a, const void *b)
 *  \brief Function to compare doubles, for use with qsort.
 */
int compare_doubles(const void *a, const void *b)
{
	//compare doubles, from GNU C Library Ref. Manual
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	return (*da > *db) - (*da < *db);
}

/*! \fn double time_in_seconds(clock_t A, clock_t B)
 *  \brief Returns time difference B-A in seconds for two clock_t
 */
double time_in_seconds(clock_t A, clock_t B)
{
	//stupid work around for icpc
	long bl = (long) B;
	long al = (long) A;
	double BA = (double) (bl-al);
	double time_check = BA;
	time_check = BA/CLOCKS_PER_SEC;
	return time_check;

	//return (double) (B-A)/CLOCKS_PER_SEC;
}





/*! \fn void create_log10_spline(double (*func)(double, void *), double *log10x, double *&log10y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
 *  \brief Routine to make a spline, interpolating in log10.
 */
void create_log10_spline(double (*func)(double, void *), double *log10x, double *&log10y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
{
	//log10x contains the locations at which
	//wish to evaluate the function func()

	//log10y is input unallocated
	//we will allocate it and it will
	//contain log10(func(x))

	//n is the size of the array of log10x and log10y
	
	//params is the array of parameters to pass to
	//func(x,params)

	//spline is the unallocated gsl_spline

	//acc is the unallocated gsl_interp_accel

	double x;
	double y;

	//allocate log10y
	log10y = calloc_double_array(n);

	//evaluate func at x locations
	for(int i=0;i<n;i++)
	{
		x = pow(10.0,log10x[i]);

		y = func(x,params);

		if(y<=0)
		{
			printf("func(%e) <= 0 (%e), cannot use log10 spline here.\n",x,y);
			fflush(stdout);
			exit(-1);
		}
		log10y[i] = log10(y);
	}


	//allocate interpolants
	spline = gsl_spline_alloc(gsl_interp_cspline,n);

	acc    = gsl_interp_accel_alloc();

	//create interpoolation
	gsl_spline_init(spline, log10x, log10y, n);
}

/*! \fn void create_linear_spline(double (*func)(double, void *), double *x, double *&y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
 *  \brief Routine to create a spline, interpolated linear.
 */
void create_linear_spline(double (*func)(double, void *), double *x, double *&y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
{
	//x contains the locations at which
	//wish to evaluate the function func()

	//y is input unallocated
	//we will allocate it and it will
	//contain func(x)

	//n is the size of the array of x and y
	
	//params is the array of parameters to pass to
	//func(x,params)

	//spline is the unallocated gsl_spline

	//acc is the unallocated gsl_interp_accel

	//allocate y
	y = calloc_double_array(n);

	//evaluate func at x locations
	for(int i=0;i<n;i++)
		y[i] = func(x[i],params);


	//allocate interpolants
	spline = gsl_spline_alloc(gsl_interp_cspline,n);

	acc    = gsl_interp_accel_alloc();

	//create interpoolation
	gsl_spline_init(spline, x, y, n);
}

double array_max(double *x, int n)
{
	double m = x[0];
	for(int i=1;i<n;i++)
	{
		m = GSL_MAX_DBL(m,x[i]);
	}
	return m;
}
double array_min(double *x, int n)
{
	double m = x[0];
	for(int i=1;i<n;i++)
	{
		m = GSL_MIN_DBL(m,x[i]);
	}
	return m;
}
