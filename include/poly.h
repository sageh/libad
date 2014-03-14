#ifndef AD_POLY_H
#define AD_POLY_H

#include <stdio.h>

/* maximal degree of static polynomials. */
#define POLY_MAXDEG	5

/* polynomial with real coefficients reserved from stack. */
typedef struct _polyr {
	int deg;
	double coeff[POLY_MAXDEG+1];
} polyr;

/* a convenient hook for haskell FFI */
polyr *poly_create(int deg, const double *coeff);
int poly_deg(const polyr *p);
double poly_eval(const polyr *p, double x); 
void poly_diff(const polyr *p, polyr *res); 
void poly_print(FILE *stream, const polyr *p);

/* return 0 on error (such as failure to converge), non-zero on success */
int root_laguerre(double tol, int maxit, const polyr *p, 
		double x0, double *res); 

#endif /* poly.h */
