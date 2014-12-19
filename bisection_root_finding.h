#ifndef BISECTION_H
#define BISECTION_H
#define TOLERANCE 1.e-9

/*check to be sure our initial guesses are valid*/
int   check_initial_values(double (*f)(double x, void *params), void *params, double x_min, double x_max);

/* function that uses bisection search to find a root */
double bisection_root_finding(double (*f)(double x, void *params), void *params, double x_min_start, double x_max_start);

#endif /*BISECTION_H*/
