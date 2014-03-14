#include "kepler.h"
#include "function.h"
#include "vector.h"
#include "misc.h"

#include <math.h>
#include <stdlib.h>

double ad_ekepler(double ecc, double E) {
	return E - ecc * sin(E);
}

double ad_gkepler(double r, double eta, double zeta, double beta, double X) {
	return r * X + eta * stumpff_G(2, beta, X) 
			+ zeta * stumpff_G(3, beta, X);
}

int ad_solve_ekepler(double tol, int maxit, double ecc, double M, double *E) {
	double f, df, ddf, x=M;
	int i=0;
	do  {
		/* use a suggestion by Seppo Mikkola */
		i++;
		f = ad_ekepler(ecc, x) - x; 
		df = 1.0 - ecc*cos(x);
		ddf = ecc*sin(x);
		x = x - f/sqrt(df*df - f*ddf);
	} while (fabs(f) > tol && i < maxit);

	*E = x;
	if (fabs(f) > tol)
		return 1;

	return 0;
}

/* use a 3rd degree Householder method */
int ad_solve_gkepler(double tol, int maxit, 
		double r, double eta, double zeta, double beta, double t,
		double X0,
		double *res) {
	/* function to minimize and derivatives w.r.t X */
	double f, df, ddf, dddf;
	double sm1, s0, s1, s2, s3;
	int i=0;
	double x = X0;

	do  {
		i++;
		/* sm2 = -2*beta(cos(beta*x*x)*x + sin(beta*x*x)); */
		sm1 = -sin(beta*x*x)*2*beta*x;
		s0 = stumpff_G(0, beta, x);
		s1 = stumpff_G(1, beta, x);
		s2 = stumpff_G(2, beta, x);
		s3 = stumpff_G(3, beta, x);
		f = r*x + eta*s2 + zeta*s3 - t;
		df = r + eta*s1 + zeta*s2;
		ddf = eta*s0 + zeta*s1;
		dddf = eta*sm1 + zeta*s0;
		x = x - (6*f*df*df - 3*f*f*ddf)
			/(6*df*df*df - 6*f*df*ddf + f*f*dddf);
	} while (fabs(f) > tol && i < maxit);

	*res = x;
	if (fabs(f) > tol)
		return 1;

	return 0;
}

void ad_kepler_propagate(double tol, int maxit, double dt,
		double *X, double *rv, double *vv) {
	double r0 = ad_vec3_norm(rv), v0 = ad_vec3_norm(vv);
	double tau = dt; /* assume mu = 1 */
	double eta = ad_vec3_iprod(rv, vv);
	double alpha = 2/r0 - v0*v0;
	double zeta = 1 - alpha*r0;
	double x0 = (X == NULL ? 0.0 : *X);
	double s1, s2, s3;
	double x, r, f, g, df, dg;
	double nrv[3];

	if (ad_solve_gkepler(tol, maxit, 
		r0, eta, zeta, alpha, tau, x0, &x) != 0)
		ad_warn("ad_kepler_propagate", 
				"ad_solve_gkepler: no convergence");
	if (X != NULL)
		*X = x;

	/* dt/dX = r */
	s1 = stumpff_G(1, alpha, x);
	s2 = stumpff_G(2, alpha, x);
	s3 = stumpff_G(3, alpha, x);
	r = r0 + eta*s1 + zeta*s2;
	f = 1 - s2/r0;
	g = tau - s3;
	df = -s1/(r*r0);
	dg = 1-s2/r;

	ad_vec3_lc2(f, rv, g, vv, nrv);
	ad_vec3_lc2(df, rv, dg, vv, vv);
	ad_vec3_cpy(rv, nrv);
}
