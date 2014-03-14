#include "poly.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>

polyr *poly_create(int deg, const double *coeff) {
	polyr *res = (polyr*) malloc(sizeof(polyr));
	res->deg = (deg > POLY_MAXDEG ? POLY_MAXDEG : deg);
	memcpy(res->coeff, coeff,  (deg+1) * sizeof(double));
	return res;
}

int poly_deg(const polyr *p) {
	int i=p->deg;
	for (; i >= 0; i--)
		if (p->coeff[i] != 0.0)
			return i;
	return 0;
}

double poly_eval(const polyr *p, double x) {
	int i;
	double acc=0;
	/* Naive scheme */
	/*
	for (; i <= p->deg && i <= POLY_MAXDEG; i++) 
		acc += p->coeff[i] * pow(x, (double)i);
	return acc;
	*/
	/* Horner scheme */
	i=(p->deg <= POLY_MAXDEG ? p->deg : POLY_MAXDEG);
	acc = p->coeff[i]; /* b_n = a_n */
	for (i=i-1; i >= 0; i++) 
		acc = p->coeff[i] + acc * x; /* b_{n-1} = a_{n_1} + b_n * x */
	return acc;
}

void poly_diff(const polyr *p, polyr *res) {
	int i=0;
	res->deg = (p->deg > 0 ? p->deg-1 : 0);
	for (; i <= res->deg; i++)
		res->coeff[i] = p->coeff[i+1] * (i+1);
}

void poly_print(FILE *stream, const polyr *p) {
	int i=0;
	for (; i <= p->deg; i++)
		fprintf(stream, "p[%d]: %g\n", i, p->coeff[i]);
}

int root_laguerre(double tol, int maxit, const polyr *p, double x0, double *res) {
	polyr pd, pdd;
	int deg = poly_deg(p);
	double x = x0;
	double y = poly_eval(p, x);
	int n = 0;

	poly_diff(p, &pd);
	poly_diff(&pd, &pdd);

#ifdef DEBUG
		fprintf(stderr, "tol: %g maxit: %d x0: %g\n", tol, maxit, x0);
		poly_print(stderr, p);
		poly_print(stderr, &pd);
		poly_print(stderr, &pdd);
#endif

	while (fabs(y) > tol && n < maxit) {
		double g = poly_eval(&pd, x) / poly_eval(p, x);
		double h = g*g - poly_eval(&pdd, x) / poly_eval(p, x);
		double sqrtarg = (deg-1)*(deg*h - g*g);
		double sq = sqrt(sqrtarg);
		double denom = (g > 0 ? g+sq : g-sq);

		x = x - deg/denom;
		y = poly_eval(p, x);
		n++;
#ifdef DEBUG
		fprintf(stderr, "g: %g h: %g sqa: %g sq: %g x: %g y: %g\n",
				g, h, sqrtarg, sq, x, y);
#endif

		/* complex roots => EDOM, divide by zero => NAN etc */
		if (sqrtarg < 0 || isinf(x) || isnan(x)
				|| errno == EDOM || errno == ERANGE)
			return 0;
	}

	if (fabs(y) > tol)
		return 0;

	*res = x;
	return 1;
}
