#include<math.h>
#include"routines.h"
#include"rk_int.h"
#include"sedov_taylor.h"
Interpolation V_xi_Interpolant;
SedovTaylor::SedovTaylor(void)
{
	flag_init = 0;

	n_V_table = 1000;

	//dummy values
	rho_1 = 1;
	E   = 1;
	t   = 1;
}
SedovTaylor::~SedovTaylor(void)
{
}
void SedovTaylor::Initialize(double gamma_in)
{
	double *V_in;
	double *xi_in;

	flag_init = 1;
	
	gamma = gamma_in;

	V_min     = 1./gamma;
	V_max     = 2./(gamma+1);

	v1 = -(13*gamma*gamma - 7*gamma + 12)/((3*gamma-1)*(2*gamma+1));
	v2 = 5*(gamma-1)/(2*gamma+1);
	v3 = 3/(2*gamma+1);
	v4 = -v1/(2-gamma);
	v5 = -2/(2-gamma);


	//V_in  = V_xi_Interpolant.MakeLog10Table(n_V_table,V_min,V_max);	
	V_in =  calloc_double_array(n_V_table);

	for(int i=0;i<n_V_table/2;i++)
		V_in[i] = pow(10., (log10(1.0001/gamma) - log10(1.0/gamma + 1.0e-6))*((double) i)/((double) n_V_table/2 -1) + log10(1.0/gamma+1.0e-6));

		//V_in[i] = (1.0001/gamma - 1./gamma)*((double) i)/((double) n_V_table/2 -1) + 1./gamma;

	for(int i=n_V_table/2;i<n_V_table;i++)
		V_in[i] = (2./(gamma+1) - 1.00011/gamma)*((double) i-n_V_table/2)/((double) n_V_table/2 -1) + 1.00011/gamma;
	

	xi_in = calloc_double_array(n_V_table);

	//printf("n_V_table %d\n",n_V_table);
	for(int i=0;i<n_V_table;i++)
	{
		xi_in[i] = xi_V(V_in[i]);
		//printf("i %d xi %e V %e\n",i,xi_in[i],V_in[i]);
	}

	//printf("HERE\n");
	//fflush(stdout);

	V_xi_Interpolant.Initialize(xi_in,V_in,n_V_table);
	//V_xi_Interpolant.ShowTable();
	//printf("Done!\n");
	//exit(0);
	//printf("V(1) %e check %e\n",V(1.),2./(gamma+1));
	//printf("Z(1) %e check %e\n",Z(1.),2*gamma*(gamma-1)/pow(gamma+1,2));
	//printf("G(1) %e check %e\n",G(1.),(gamma+1)/(gamma-1));
	
	//fflush(stdout);

	DetermineBeta();

	free(xi_in);
	free(V_in);
}
double SedovTaylor::xi_V(double V)
{
	double fp[6];
	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;
	return xi_V_SedovTaylor(V,fp);
}
double xi_V_SedovTaylor(double V, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    = fp[3];
	double v4    = fp[4];
	double v5    = fp[5];
	double A = pow( 0.5*(gamma+1)*V, -2.);
	double B = pow( (gamma+1)/(7-gamma)*(5-(3*gamma-1)*V),v1);
	double C = (gamma+1)/(gamma-1)*(gamma*V-1);

	//printf("V %e A %e B %e C %e v2 %e\n",V,A,B,C,v2);
	return pow(A*B*pow(C,v2),0.2);
}
double SedovTaylor::V(double xi)
{
	double fp[6];
	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;
	return V_SedovTaylor(xi,fp);
}
double SedovTaylor::Z(double xi)
{
	double fp[6];
	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;
	return Z_SedovTaylor(xi,fp);
}
double SedovTaylor::G(double xi)
{
	double fp[6];

	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;

	return G_SedovTaylor(xi,fp);
}
void SedovTaylor::DetermineBeta(void)
{
	
	double fp[6];
	double b_int;
	double xi_min = 0.00;
	double xi_max = 1;
	double dinit = 0.1*(xi_max - xi_min);
	double eps = 1.0e-6;
	double K = 16.*M_PI/25;

	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;

	b_int = integrate(Beta_Integrand,fp,6,xi_min,xi_max,dinit,eps);
	
	beta = pow(K*b_int,-0.2);
	//printf("beta = %e b_int %e K %e\n",beta,b_int,K);
}
double SedovTaylor::R(double t)
{
	//shock radius
	return beta*pow(E*t*t/rho_1,0.2);
}
double SedovTaylor::dRdt(double t)
{
	//shock radius time derivative
	return 0.4*beta*pow(E/rho_1,0.2)*pow(t,-0.6);
}
double SedovTaylor::c_squared(double r)
{
	double xi = r/R(t);
	return 4*r*r*Z(xi)/(25.*t*t);
}
double SedovTaylor::v(double r)
{
	double xi = r/R(t);
	if(xi<=1.0)
	{
		return 2.0*r*V(xi)/(5.*t);
	}else{
		return 0;
	}
}
double SedovTaylor::p(double r)
{
	double xi = r/R(t);
	if(xi<=1.0)
	{
		return rho(r)*c_squared(r)/gamma;
	}else{
		return 0;
	}
}
double SedovTaylor::rho(double r)
{
	double xi = r/R(t);
	if(xi<1.0)
	{
		return rho_1*G(xi);
	}else{
		return rho_1;
	}
}
double V_SedovTaylor(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	//printf("VST xi %e\n",xi);
	//fflush(stdout);
	double V = V_xi_Interpolant.y(xi);
	//printf("x %e V %e\n",xi,V);
	//fflush(stdout);
	if(V<1.000001/gamma)
	{
		double xi_s = xi_V_SedovTaylor(1.000001/gamma,params);
		//double G_s  = G_SedovTaylor(xi_s,params);

		//printf("xi %e xis %e V %e\n",xi,xi_s,V);
		V = 1.0/gamma + 0.000001*pow(xi/xi_s,5./v2) + 1.0e-12;
	}
	return V;
}
double Z_SedovTaylor(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double V = V_SedovTaylor(xi,fp);
	double Z = gamma*(gamma-1)*(1-V)*V*V/(2*(gamma*V-1));


	double xi_s = xi_V_SedovTaylor(1./gamma+1.0e-6,params);

	//printf("xi %e xi_s %e V %e Z %e\n",xi,xi_s,V,Z);

	if(xi<xi_s)
	{
		double G_s  = G_SedovTaylor(xi_s,params);
		///printf("G_s %e\n",G_s);
		//fflush(stdout);
		double G    = G_SedovTaylor(xi,params);
		double Z_s  = Z_SedovTaylor(xi_s,params);
		return Z_s*(G_s/G)*pow((xi+1.0e-15)/xi_s,-2.0);
	}
	
	
	//printf("ZZ V %e 1./gamma %e Z %e\n",V,1./gamma,Z);
	return Z;
}
double G_SedovTaylor(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    = fp[3];
	double v4    = fp[4];
	double v5    = fp[5];
	//printf("xi %e\n",xi);
	//fflush(stdout);
	double V = V_SedovTaylor(xi,fp);
	//printf("V %e\n",V);
	//fflush(stdout);
	//double V = V_xi_Interpolant.y(xi);
	//if(V<1.000001/gamma)
	//	V=1.000001/gamma;
		
	double A = (gamma+1)/(gamma-1);
	double B = pow(A*(gamma*V-1),v3);
	double C = pow((gamma+1)/(7-gamma)*(5-(3*gamma-1)*V),v4);
	double D = pow(A*(1-V),v5);
	//printf("xi %e V %e check %e\n",xi,V,2./(gamma+1));
	//printf("GG V %e A %E B %e C %e D %e\n",V,A,B,C,D);
	//if(xi<=0.25)
	//	return 0;

	return A*B*C*D;
}
double Beta_Integrand(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	//printf("xi %e gamma %e\n",xi,gamma);
	//fflush(stdout);
	double G = G_SedovTaylor(xi,fp);
	//printf("HERE G %e\n",G);
	//fflush(stdout);
	double V = V_SedovTaylor(xi,fp);
	//	printf("HERE V %e\n",V);
	//fflush(stdout);
	double Z = Z_SedovTaylor(xi,fp);
	//	printf("HERE Z %e\n",Z);
	//fflush(stdout);
	double bi = G*(0.5*V*V + Z/(gamma*(gamma-1)))*xi*xi*xi*xi;
	bi = fabs(bi);
	//printf("xi %e G %e V %e Z %e 1./gamma %e bi %e\n",xi,G,V,Z,1./gamma,bi);
	//fflush(stdout);
	return bi;
}
