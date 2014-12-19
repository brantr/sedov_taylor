#ifndef SEDOV_TAYLOR
#define SEDOV_TAYLOR
#include "interpolation.h"
extern Interpolation V_xi_Interpolant;
class SedovTaylor
{
	public:
		double v1;
		double v2;
		double v3;
		double v4;
		double v5;

		double gamma;

		double beta;

		int flag_init;

		void Initialize(double gamma_in);

		void DetermineBeta(void);
		SedovTaylor(void);
		~SedovTaylor(void);

		int n_V_table;
		double V_min;
		double V_max;
		double xi_V(double Vx);
		double V_xi(double xV);


		double Z(double xi);
		double G(double xi);
		double V(double xi);

		double t;
		double rho_1;
		double E;
		double rho(double r);
		double p(double r);
		double v(double r);
		double c_squared(double r);
		double R(double t);
		double dRdt(double t);
};
double f_V_xi(double V, void *params);
double df_V_xi_dV(double V, void *params);
double xi_V_SedovTaylor(double V, void *params);
double Z_SedovTaylor(double xi, void *params);
double G_SedovTaylor(double xi, void *params);
double V_SedovTaylor(double xi, void *params);
double Beta_Integrand(double xi, void *params);
#endif //SEDOV_TAYLOR
