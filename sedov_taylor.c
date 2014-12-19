#include <math.h>
#include "routines.h"
#include "rk_int.h"
#include "sedov_taylor.h"
#include "bisection_root_finding.h"
SedovTaylor::SedovTaylor(void)
{
	flag_init = 0;

	n_V_table = 20000;

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
	V_max     = 2./(gamma+1.);

	v1 = -(13.*gamma*gamma - 7.*gamma + 12.)/((3.*gamma-1.)*(2.*gamma+1.));
	v2 = 5.*(gamma-1.)/(2.*gamma+1.);
	v3 = 3./(2.*gamma+1.);
	v4 = -v1/(2.-gamma);
	v5 = -2./(2.-gamma);


	DetermineBeta();

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
	double xi_V_r = xi_V_SedovTaylor(V,fp);
	return xi_V_r;
}
double SedovTaylor::V_xi(double xi)
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
double xi_V_SedovTaylor(double V, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    = fp[3];
	double v4    = fp[4];
	double v5    = fp[5];
	double A = pow( 0.5*(gamma+1.)*V, -2.);
	double B = pow( (gamma+1.)/(7.-gamma)*(5.-(3.*gamma-1.)*V),v1);
	double C = (gamma+1.)/(gamma-1.)*(gamma*V-1.);
	double xi_V_r = pow(A*B*pow(C,v2),0.2);
	return xi_V_r;
}

double f_V_xi(double dV, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    = fp[3];
	double v4    = fp[4];
	double v5    = fp[5];
	double xi    = fp[6];
	double xiV = xi_V_SedovTaylor(1./gamma + dV,params);
	double dxi = xi-xiV;// 
	return dxi/xi;
}
double df_V_xi_dV(double dV, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    = fp[3];
	double v4    = fp[4];
	double v5    = fp[5];
	double xi    = fp[6];
	double V     = 1./gamma;
	double A = pow( 0.5*(gamma+1.)*(V+dV), -2.);
	double dAdV = -2.0*(0.5*gamma+1.)*pow( 0.5*(gamma+1.)*(V+dV), -3.);
	double B = pow( (gamma+1)/(7-gamma)*(5-(3*gamma-1)*(V+dV)),v1);
	double dBdV = v1*( (gamma+1)/(7-gamma)*(3*gamma-1) )*pow( (gamma+1)/(7-gamma)*(5-(3*gamma-1)*(V+dV)),v1-1.);
	double C = (gamma+1)/(gamma-1)*(gamma*(V+dV)-1);
	double dCdV = (gamma+1)/(gamma-1)*gamma;
	double dfdV = 0.2*( dAdV*B*pow(C,v2) + A*dBdV*pow(C,v2) + v2*A*B*pow(C,v2-1.)*dCdV )*pow(A*B*pow(C,v2),0.2);
	return dfdV;
}
/*

[ A * B * C^v2 ]^0.2 - xi = 0
dAdV
0.2*[ dAdV * B * C^v2 + A*dBdV * C^v2 + v2*A*B*C^(v2-1)*dCdV ]

*/


double SedovTaylor::V(double xi)
{
	double fp[6];
	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;
	return V_xi(xi);
//	return V_SedovTaylor(xi,fp);
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
	double V = V_xi(xi);
	double Z = gamma*(gamma-1.)*(1.-V)*V*V/(2*(gamma*V-1.));
	return Z;
	//return Z_SedovTaylor(xi,fp);
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
	double Vv = V_xi(xi);
	double A = (gamma+1.)/(gamma-1.);
	double B = pow(A*(gamma*Vv-1.),v3);
	double C = pow((gamma+1.)/(7.-gamma)*(5.-(3.*gamma-1.)*Vv),v4);
	double D = pow(A*(1.-Vv),v5);

	return A*B*C*D;
}
void SedovTaylor::DetermineBeta(void)
{
	
	double fp[6];
	double b_int;
	double xi_min = 1.0e-12;
	double xi_max = 1;
	double dinit = 0.1*(xi_max - xi_min);

	double eps = 1.0e-7;
	double K = 16.*M_PI/25.;

	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;



	b_int = integrate(Beta_Integrand,fp,6,xi_min,xi_max,dinit,eps);
	
	beta = pow(K*b_int,-0.2);

}
double SedovTaylor::R(double t)
{
	//shock radius
	double Rs;
	Rs = beta*pow(E*t*t/rho_1,0.2);
	return Rs;
}
double SedovTaylor::dRdt(double t)
{
	//shock radius time derivative
	return 0.4*beta*pow(E/rho_1,0.2)*pow(t,-0.6);
}
double SedovTaylor::c_squared(double r)
{
	double xi = r/R(t);
	return 4.*r*r*Z(xi)/(25.*t*t);
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

	double Gg;
	if(xi<1.0)
	{
		Gg = G(xi);
	}else{
		Gg=1;
	}

	return rho_1*Gg;
}
double V_SedovTaylor(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
	double gp[7];
	for(int i=0;i<6;i++)
		gp[i] = fp[i];
	gp[6] = xi;


	double V_return;
	double V2;
	double xi_min = 1.0e-16;
	if(xi < 1.0e-2)
	{
		if(xi<xi_min)
			return 1./gamma;
		gp[6] = 1.0e-2;
		V2 = bisection_root_finding(f_V_xi,gp,1.0e-14,(1.0-1/gamma)-1.e-14);
		V_return = 1./gamma + V2*( log(xi) - log(xi_min) )/(log(1.0e-2) - log(xi_min));

	}else{
		V_return = 1./gamma + bisection_root_finding(f_V_xi,gp,1.0e-14,(1.0-1/gamma)-1.e-14);
	}

	return V_return;
}
double Z_SedovTaylor(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];

	double V = V_SedovTaylor(xi,fp);
	double Z = gamma*(gamma-1.)*(1.-V)*V*V/(2*(gamma*V-1.));

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

	double V = V_SedovTaylor(xi,fp);

		
	double A = (gamma+1.)/(gamma-1.);
	double B = pow(A*(gamma*V-1.),v3);
	double C = pow((gamma+1.)/(7.-gamma)*(5.-(3.*gamma-1.)*V),v4);
	double D = pow(A*(1.-V),v5);

	return A*B*C*D;
}
//double Beta_Integrand(double xi, void *params)
double Beta_Integrand(double lnxi, void *params)
{
	//double xi = exp(lnxi);
	double xi = lnxi;
	double *fp = (double *) params;
	double gamma = fp[0];
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    = fp[3];
	double v4    = fp[4];
	double v5    = fp[5];

	double gp[7];
	for(int i=0;i<6;i++)
		gp[i] = fp[i];
	gp[6] = xi;


	double G = G_SedovTaylor(xi,gp);

	double V = V_SedovTaylor(xi,gp);

	double Z = Z_SedovTaylor(xi,gp);
;
	double bi = G*(0.5*V*V + Z/(gamma*(gamma-1.)))*xi*xi*xi*xi;
	bi = fabs(bi);

	return bi;
}
