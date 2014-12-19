#include <stdio.h>
#include "sedov_taylor.h"

int main(int argc, char **argv)
{
	SedovTaylor ST;

	double r_min = 1.0e-5;
	double r_max = 0.6;
	double r;
	int n = 1000;
	double gamma = 5./3.;
	double xi;

	ST.t = 0.06;
	ST.E = 1.0;

	ST.Initialize(gamma);


	r_max = 0.999999999*ST.R(ST.t);
	r_min = 1.0e-2;//*ST.R(ST.t);
	//r_max = 1.1*ST.R(ST.t);
	//printf("rho_1 %e\n",ST.rho_1);
	printf("%d\n",n);
	for(int i=0;i<n;i++)
	{
		r = (r_max-r_min)*((double) i)/((double) (n-1)) + r_min;
    //printf("rhere %e r_min %e r_max %e\n",r,r_min,r_max);
    fflush(stdout);
		printf("%e\t%e\t%e\t%e\n",r,ST.rho(r),ST.p(r),ST.v(r));
		xi = r/ST.R(ST.t);
		//printf("%e\t%e\t%e\t%e\n",xi,ST.G(xi),ST.Z(xi),ST.V(xi));
	}
	return 0;
}
