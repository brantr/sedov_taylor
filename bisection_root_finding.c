#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bisection_root_finding.h"


double bisection_root_finding(double (*f)(double x, void *params), void *params, double x_min_start, double x_max_start)
{
  /* function that uses bisection search to find a root */

  double x_min = x_min_start;  /* minimum x in bracket */
  double x_mid;		      /* 1/2 between x_min and x_max */
  double x_max = x_max_start;  /* maximum x in bracket */
  double y_min;		      /* f at x_min */
  double y_mid;		      /* f at x_mid */
  double y_max;		      /* f at x_max */

  int flag    = 0;	      /* a flag for triggering */
  int iter    = 0;	      /* current number of iterations */
  int imax    = 100000;	      /* maximum number of iterations */

  /* check our initial bracket to make sure it's valid */
  flag = check_initial_values(f, params, x_min, x_max);
  if(flag)
  {
    /* we made a lucky guess and we're done */
    if(flag==1)	    /*x_min is the root*/
      return x_min;
    if(flag==2)     /*x_max is the root*/
      return x_max;
  }

  /*instead, if flag==0, we continue on our way to find the root*/

  /* loop until we find the root or we reach the maximum number of iterations*/
  while(!flag)
  {
    /* find value midway between bracket values */
    x_mid = 0.5*(x_min + x_max);

    /* determine the function value at this point */
    y_mid = f(x_mid,params);

    /* if the absolute value of f(x) is < tolerance, then exit */
    if(fabs(y_mid) < TOLERANCE || fabs(x_max-x_min) < 1.0e-14)
    {
      /* trigger an exit */
      flag = 1;

    }else{

      /*x_mid is not a root*/
  
      y_min = f(x_min,params);
      y_max = f(x_max,params);

      /* print guesses */
      //print_guess(x_min,y_min);
      //print_guess(x_mid,y_mid);
      //print_guess(x_max,y_max);

      //printf("xmin %15.14e ymin %15.14e xmid %15.14e ymid %15.14e xmax %15.14e ymax %15.14e\n",x_min,y_min,x_mid,y_mid,x_max,y_max);

      /* if the product of the function at the midpoint and at one of the endpoints */
      /* is greater than zero, replace this end point */
      //if( f(x_min,params)*f(x_mid,params) > 0)
      if( y_min * y_mid > 0)
      {
	/*replace x_min with x_mid */
	x_min = x_mid;
      }else{
	/*replace x_max with x_mid */
	x_max = x_mid;
      }
      

      /* count the iteration */
      iter++;

      /* if we have exceeded the max number of iterations, exit*/
      if(iter>=imax)
      {
	printf("Exceeded %d iterations, exiting (x=%f, f=%f)...\n",imax,x_mid,y_mid);
	fflush(stdout);	
	exit(-1);
      }
    }
  }

  /* return x_mid == root */
  return x_mid;
  
}


int check_initial_values(double (*f)(double x, void *params), void *params, double x_min, double x_max)
{

  /*check to be sure our initial guesses are valid*/

  double y_min = f(x_min,params);
  double y_max = f(x_max,params);

 // printf("CIV x_min %15.14e y_min %15.14e x_max %15.14e y_max %15.14e\n",x_min,y_min,x_max,y_max);

  /*check for zero crossing*/
  if( y_min*y_max >= 0.0 )
  {
    printf("No zero crossing found in range [%15.14e, %15.14e].\n",x_min,x_max);
    printf("f(%15.14e) = %15.14e, f(%15.14e) = %15.14e\n",x_min,y_min,x_max,y_max);
    printf("Change initial values.\n");
    fflush(stdout);
    exit(-1);
  }

  /*If x_min is a root, then return a flag == 1*/
  if( fabs(y_min) < TOLERANCE)
    return 1;

  /*If x_max is a root, then return a flag == 2*/
  if( fabs(y_max) < TOLERANCE)
    return 2;

  /*if we reach this point, the bracket is valid.*/
  return 0;
}
