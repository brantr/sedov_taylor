/*! \file routines.h
 *  \brief Declarations for a variety of utility routines.
 */
#ifndef  BRANT_ROUTINES
#define  BRANT_ROUTINES
#include<stdio.h>
#include<gsl/gsl_spline.h>
#include<time.h>

/*! \fn void check_args(int argc, char **argv, int num_args)
 *  \brief Ensures that the number of command line arguements equals num_args.
 */
void      check_args(int argc, char **argv, int num_args);

/*! \fn FILE *fopen_brant(char fname[], const char *mode)
 *  \brief Safe method for opening a FILE pointer.
 */
FILE     *fopen_brant(char fname[], const char *mode);
/*! \fn double *calloc_double_array(int n)
 *  \brief Safe method for callocing a double array
 */
double   *calloc_double_array(int n);
/*! \fn float *calloc_float_array(int n)
 *  \brief Safe method for callocing a float array
 */
float    *calloc_float_array(int n);
/*! \fn int *calloc_int_array(int n)
 *  \brief Safe method for callocing an int array
 */
int      *calloc_int_array(int n);
/*! \fn size_t *calloc_size_t_array(int n)
 *  \brief Safe method for callocing a size_t array
 */
size_t   *calloc_size_t_array(int n);
/*! \fn double **two_dimensional_array(int n, int l)
 *  \brief Allocate a two dimensional (n x l) array
 */
double  **two_dimensional_array(int n, int l);
/*! \fn void deallocate_two_dimensional_array(double **x, int n, int l)
 *  \brief De-allocate a two dimensional (n x l) array
 */
void      deallocate_two_dimensional_array(double **x, int n, int l);
/*! \fn double ***three_dimensional_array(int n, int l, int m)
 *  \brief Allocate a three dimensional (n x l x m) array 
 */
double ***three_dimensional_array(int n, int l, int m);
/*! \fn void deallocate_three_dimensional_array(double ***x, int n, int l, int m)
 *  \brief De-allocate a three dimensional (n x l x m) array
 */
void      deallocate_three_dimensional_array(double ***x, int n, int l, int m);
/*! \fn int ***three_dimensional_int_array(int n, int l, int m)
 *  \brief Allocate a three dimensional (n x l x m) int array
 */
int    ***three_dimensional_int_array(int n, int l, int m);
/*! \fn void deallocate_three_int_dimensional_array(int ***x, int n, int l, int m)
 *  \brief De-allocate a three dimensional (n x l x m) int array.
 */
void      deallocate_three_dimensional_int_array(int ***x, int n, int l, int m);
/*! \fn double max_three(double a, double b, double c)
 *  \brief Returns the max of 3 numbers */
double    max_three(double a, double b, double c);
/*! \fn int compare_doubles(const void *a, const void *b)
 *  \brief Function to compare doubles, for use with qsort.
 */
int       compare_doubles(const void *a, const void *b);
/*! \fn double time_in_seconds(clock_t A, clock_t B)
 *  \brief Returns time difference B-A in seconds for two clock_t
 */
double    time_in_seconds(clock_t A, clock_t B);
/*! \fn double double double_linear_index(int i, int n, double xmin, double xmax)
 *  \brief Provides the i^th out of n linear incremented value between xmin and xmax.
 *	   Useful for creating a ordinate array for an interpolation.
 */
double double_linear_index(int i, int n, double xmin, double xmax);
/*! \fn double double_log10_index(int i, int n, double xmin, double xmax)
 *  \brief Provides the i^th out of n log10 incremented value between xmin and xmax.
 *	   Useful for creating a ordinate array for an interpolation.
 */
double double_log10_index(int i, int n, double xmin, double xmax);
/*! \fn void create_linear_spline(double (*func)(double, void *), double *x, double *&y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
 *  \brief Routine to create a spline, interpolated linear.
 */
void     create_linear_spline(double (*func)(double, void *), double *x, double *&y, int n, double *params, gsl_spline *spline, gsl_interp_accel *acc);
/*! \fn void create_log10_spline(double (*func)(double, void *), double *log10x, double *&log10y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
 *  \brief Routine to make a spline, interpolating in log10.
 */
void     create_log10_spline(double (*func)(double, void *), double *log10x, double *&log10y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc);


/*! \fn double array_max(double *x, int n)
 *  \brief Find the maximum of array x */
double array_max(double *x, int n);

/*! \fn double array_min(double *x, int n)
 *  \brief Find the minimum of array x */
double array_min(double *x, int n);

#endif //BRANT_ROUTINES
