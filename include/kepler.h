#ifndef AD_KEPLER_H
#define AD_KEPLER_H

/* elliptic Kepler's equation. returns mean anomaly. */
double ad_ekepler(double ecc, double E);

/* general variable formulation of Kepler's equation. returns tau =
 * sqrt(mu) * time */
double ad_gkepler(double r, double eta, double zeta, double beta, double X);

/* solve elliptic Kepler's equation for eccentric anomaly, and return it
 * in *res. returns 0 on convergence, -1 on non-convergence */
int ad_solve_ekepler(double tol, int maxit, double ecc, double M, double *res);

/* solve general Kepler's equation for variable X, and return it
 * in *res. returns 0 on convergence, -1 on non-convergence */
int ad_solve_gkepler(double tol, int maxit, 
		double r, double eta, double zeta, double beta, double t,
		double X0,
		double *res);

/* propagate position and velocity by dt along a keplerian path.
 * previous value of X or NULL can be given in *X, which will also be
 * updated. assumes a coordinate system in which mu = G(m1+m2) = 1 */
void ad_kepler_propagate(double tol, int maxit, double dt,
		double *X, double *rv, double *vv);

#endif /* kepler.h */
