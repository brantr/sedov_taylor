#include <math.h>
#include "routines.h"
#include "rk_int.h"
#include "sedov_taylor.h"
#include "bisection_root_finding.h"
Interpolation V_xi_Interpolant;
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


	//V_in  = V_xi_Interpolant.MakeLog10Table(n_V_table,V_min,V_max);	
	V_in =  calloc_double_array(n_V_table);


	double xi_1 = 1.0e-9;
	double xi_2 = 1.1*xi_1;
	double xi_3 = 1.0e-6;
	double xi_4 = 1.1*xi_3;


	for(int i=0;i<n_V_table/4;i++)
		V_in[i] = 1. + xi_1*((double) i)/((double) (n_V_table/4-1)) + 2.5e-14;

	for(int i=0;i<3*n_V_table/4;i++)
		V_in[i + n_V_table/4] = ( 2.*gamma/(gamma+1.) - (1. + xi_2) )*((double) i)/((double) (3*n_V_table/4-1)) + (1. + xi_2);

	//for(int i=0;i<n_V_table/4;i++)
	//	V_in[i + n_V_table/4] = pow(10., ( log10(1.+xi_3)- log10(1.+xi_2))*((double) i)/((double) (n_V_table/4 -1)) + log10(1. + xi_2));
		//V_in[i + n_V_table/4] = (1.+xi_3- (1.+xi_2))*((double) i)/((double) (n_V_table/4 -1)) + (1. + xi_2);

		//V_in[i] = pow(10., (log10(1.0001/gamma) - log10(1.0/gamma + 1.0e-6))*((double) i)/((double) n_V_table/2 -1) + log10(1.0/gamma+1.0e-6));

		//V_in[i] = (1.0001/gamma - 1./gamma)*((double) i)/((double) n_V_table/2 -1) + 1./gamma;

	//for(int i=n_V_table/2;i<n_V_table;i++)
	//	V_in[i] = (2.*gamma/(gamma+1) - (1. + xi_4))*((double) i-n_V_table/2)/((double) n_V_table/2 -1) + (1. + xi_4);

	xi_in = calloc_double_array(n_V_table);

	//printf("0.6 xi %e\n",xi_V(0.6+1.0e-9));
	//exit(0);
	//printf("n_V_table %d\n",n_V_table);
	for(int i=0;i<n_V_table;i++)
	{	
		V_in[i] /= gamma;
		//printf("i %d V_in %15.14e\n",i,V_in[i]-1./gamma);
		//fflush(stdout);
		xi_in[i] = xi_V(V_in[i]);

		//if(xi_in[i]<2.0e-2)
			//printf("BUILD i %d Vi %15.14e xi %15.14e Vr %15.14e\n",i,V_in[i],xi_in[i],V_xi(xi_in[i]));
		//fflush(stdout);
		//printf("i %d xi %15.14e r %15.14e V %15.14e Vmax %e\n",i,xi_in[i],xi_in[i]*R(t),V_in[i],2./(gamma+1));
	}
	//fflush(stdout);
	//exit(0);

	//printf("HERE\n");
	//fflush(stdout);

	//V_xi_Interpolant.Initialize(xi_in,V_in,n_V_table);
	//V_xi_Interpolant.ShowTable();
	//printf("Done!\n");
	//exit(0);
	//printf("V(1) %e check %e\n",V(1.),2./(gamma+1));
	//printf("Z(1) %e check %e\n",Z(1.),2*gamma*(gamma-1)/pow(gamma+1,2));
	//printf("G(1) %e check %e\n",G(1.),(gamma+1)/(gamma-1));
	
	//fflush(stdout);

	DetermineBeta();

/*
	for(int i=0;i<n_V_table;i++)
	{	
		printf("i %d xi %15.14e r %15.14e V %15.14e Vmax %e\n",i,xi_in[i],xi_in[i]*R(t),V_in[i],2./(gamma+1));
	}
	fflush(stdout);
	exit(0);
*/
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
	double xi_V_r = xi_V_SedovTaylor(V,fp);
	//printf("%15.14e\n",xi_V_r);
	//printf("V %15.14e %15.4e V %15.4e\n",V,xi_V_r,V);
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
//	fp[6] = xi;
	//printf("V_xi  xi %15.14e\n",xi);
	//return 1./gamma + newton_raphson(f_V_xi,df_V_xi_dV,fp,1.0e-10);
	//return 1./gamma + bisection_root_finding(f_V_xi,fp,1.0e-12,(1.0-1/gamma)-1.e-12);
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

	//printf("xi_V_SedovTaylor V %15.14e A %15.14e B %15.14e C %15.14e v2 %15.14e xi_V_r %15.14e\n",V,A,B,C,v2,xi_V_r);
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
	double dxi = xi-xiV;// - xi;
	//printf("f_V_xi dV %15.14e V %15.14e xi %15.14e xiV %15.14e dxi %15.14e\n",dV,1./gamma+dV,log(xi),log(xiV),dxi);
	//returns 
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
	//printf("V %e A %e B %e C %e dAdV %e dBdV %e dCdV %e v2 %e dfdV %e\n",V+dV,A,B,C,dAdV,dBdV,dCdV,v2,dfdV);
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
	//printf("xi %e V %e check %e\n",xi,V,2./(gamma+1));
	//printf("GG V %e A %E B %e C %e D %e\n",V,A,B,C,D);
	//if(xi<=0.25)
	//	return 0;

	return A*B*C*D;
	//return G_SedovTaylor(xi,fp);
}
void SedovTaylor::DetermineBeta(void)
{
	
	double fp[6];
	double b_int;
	double xi_min = 1.0e-12;
	double xi_max = 1;
	double dinit = 0.1*(xi_max - xi_min);
  //double dinit = 0.1*(log(xi_max) - log(xi_min));

	double eps = 1.0e-7;
	double K = 16.*M_PI/25.;

	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;



	b_int = integrate(Beta_Integrand,fp,6,xi_min,xi_max,dinit,eps);

	//b_int = integrate(Beta_Integrand,fp,6,log(xi_min),log(xi_max),dinit,eps);

	//printf("gamma %e v1 %e v2 %e v3 %e v4 %e v5 %e b_int %e K %e \n",gamma,v1,v2,v3,v4,v5,b_int,K);

	
	beta = pow(K*b_int,-0.2);
	//printf("beta = %15.14e b_int %15.14e K %15.14e\n",beta,b_int,K);
	//beta = 1.15;
	//printf("R = %e\n",R(t));
	//exit(0);
}
double SedovTaylor::R(double t)
{
	//shock radius
	double Rs;
	Rs = beta*pow(E*t*t/rho_1,0.2);
	//printf("beta %e E %e t %e t %e rho_1 %e Rs %e\n",beta,E,t,t,rho_1,Rs);
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

	//printf("in rho xi %e r %e rho_1 %e G %e rho_1*G %e\n",xi,r,rho_1,Gg,rho_1*Gg);
	//fflush(stdout);
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
	/*
	double v1    = fp[1];
	double v2    = fp[2];
	double v3    =
		double fp[7];
	fp[0] = gamma;
	fp[1] = v1;
	fp[2] = v2;
	fp[3] = v3;
	fp[4] = v4;
	fp[5] = v5;
	fp[6] = xi;*/
	//#error HERE

	//return 1./gamma + newton_raphson(f_V_xi,df_V_xi_dV,fp,1.0e-10);

	//printf("gp gamma %15.14e xi %15.14e\n",gamma,xi);

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
	
	/*
	printf("V_ST  xi %15.14e V_ret %15.14e\n",xi,V_return);
	printf("**********************************\n");
	printf("**********************************\n");
	printf("**********************************\n");
	*/

	fflush(stdout);
	return V_return;
	//printf("VST xi %e\n",xi);
	//fflush(stdout);
	double V = 0;//V_xi_Interpolant.y(xi);
	double V_alt = 0.0;
	//printf("x %e V %e\n",xi,V);
	//fflush(stdout);


	
	if(V<1.0/gamma + 0.5e-11)//HERE
	{
		double xi_s = xi_V_SedovTaylor(1./gamma + 1.0e-11,params);
		//double G_s  = G_SedovTaylor(xi_s,params);
		double V_s = V_SedovTaylor(xi_s,params);

		V_alt = (V_s - 1.0/gamma)*(xi/xi_s) + 1./gamma;// + 1.0e-11 *xi/xi_s;//pow(xi/xi_s,5./v2);
		//printf("xi %e xis %e V %15.14e V_s %15.14e V_alt %15.14e\n",xi,xi_s,V,V_s,V_alt);
		V = V_alt;
	}
	
	return V;
}
double Z_SedovTaylor(double xi, void *params)
{
	double *fp = (double *) params;
	double gamma = fp[0];
/*		double gp[7];
	for(int i=0;i<6;i++)
		gp[i] = fp[i];
	gp[6] = xi;*/
	double V = V_SedovTaylor(xi,fp);
	double Z = gamma*(gamma-1.)*(1.-V)*V*V/(2*(gamma*V-1.));

	return Z;

	double xi_s = xi_V_SedovTaylor(1./gamma+1.0e-11,params);
	double G_s = 0;
	double G   = 0;
	double Z_s = 0;
	double V_s = 0;


	//printf("Z_SedovTaylor xi %e xi_s %e V %e Z %e\n",xi,xi_s,V,Z);

	if(xi<xi_s)
	{
		//G_s  = G_SedovTaylor(xi_s,params);
		//G    = G_SedovTaylor(xi,params);
		//printf("G_s %e G %e\n",G_s,G);
		//fflush(stdout);
		V_s = V_SedovTaylor(xi_s,fp);
		//printf("V_s %e\n",V_s);
		Z_s  = Z_SedovTaylor(xi_s,params);
		Z = Z_s + (Z-Z_s)*(xi-xi_s)/(xi_s);
		//Z = Z_s*(G_s/G)*pow((xi+1.0e-15)/xi_s,-2.0);
	}
	
	
	//printf("ZZ V %e 1./gamma %e Z %e g*(g-1) %e (1-V)*V*V %e 2*(g*V-1) %e xi %e G %e G_s %e Z_s %e\n",V,1./gamma,Z,gamma*(gamma-1.),(1.0-V)*V*V,2*(gamma*V-1.),xi,G,G_s,Z_s);
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
/*			double gp[7];
	for(int i=0;i<6;i++)
		gp[i] = fp[i];
	gp[6] = xi;*/
	//double xi    = fp[6];
	//printf("G+ST xi %15.14e\n",xi);
	//fflush(stdout);
	double V = V_SedovTaylor(xi,fp);
	//printf("V %e\n",V);
	//fflush(stdout);
	//double V = V_xi_Interpolant.y(xi);
	//if(V<1.000001/gamma)
	//	V=1.000001/gamma;
		
	double A = (gamma+1.)/(gamma-1.);
	double B = pow(A*(gamma*V-1.),v3);
	double C = pow((gamma+1.)/(7.-gamma)*(5.-(3.*gamma-1.)*V),v4);
	double D = pow(A*(1.-V),v5);
	//printf("xi %e V %e check %e\n",xi,V,2./(gamma+1));
	//printf("GG V %e A %E B %e C %e D %e\n",V,A,B,C,D);
	//if(xi<=0.25)
	//	return 0;

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

	/*if(xi<1.0e-2)
	{
		return Beta_Integrand(log(1.0e-2),params);
	}*/

	//printf("xi %e gamma %e\n",xi,gamma);
	//fflush(stdout);
	double G = G_SedovTaylor(xi,gp);
	//printf("HERE G %e\n",G);
	//fflush(stdout);
	double V = V_SedovTaylor(xi,gp);
	//printf("HERE V %e\n",V);
	//fflush(stdout);
	double Z = Z_SedovTaylor(xi,gp);
	//printf("HERE Z %e\n",Z);
	//fflush(stdout);
	double bi = G*(0.5*V*V + Z/(gamma*(gamma-1.)))*xi*xi*xi*xi;
	bi = fabs(bi);
	//printf("Beta_Integrand xi %15.14e G %e V %15.14e Z %e 1./gamma %e bi %e\n",xi,G,V,Z,1./gamma,bi);
	//fflush(stdout);
	return bi;
}
