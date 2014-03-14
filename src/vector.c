#include "vector.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* n-vectors */

void ad_vec_cpy(int n, double *to, const double *from) {
	memcpy(to, from, n*sizeof(double));
}

void ad_vec_add(int n, const double *v1, const double *v2, double *res) {
	int i=0;
	for (; i < n; i++) res[i] = v1[i]+v2[i];
}

void ad_vec_sub(int n, const double *v1, const double *v2, double *res) {
	int i=0;
	for (; i < n; i++) res[i] = v1[i]-v2[i];
}

void ad_vec_smul(int n, double *v, double s) {
	int i=0;
	for (; i < n; i++) v[i] *= s;
}

double ad_vec_iprod(int n, const double *v1, const double *v2) {
	int i=0;
	double acc=0;
	for (; i < n; i++) acc += v1[i]*v2[i];
	return acc;
}

double ad_vec_norm(int n, const double *v) {
	return sqrt(ad_vec_iprod(n, v, v));
}


/* 3-vectors */

const double ad_vec3_zero[3] = {0., 0., 0.};

void ad_vec3_cpy(double *to, const double *from) {
	ad_vec_cpy(3, to, from);
}

void ad_vec3_add(const double *v1, const double *v2, double *res) {
	ad_vec_add(3, v1, v2, res);
}

void ad_vec3_sub(const double *v1, const double *v2, double *res) {
	ad_vec_sub(3, v1, v2, res);
}

void ad_vec3_smul(double *v, double s) {
	ad_vec_smul(3, v, s);
}

double ad_vec3_iprod(const double *v1, const double *v2) {
	return ad_vec_iprod(3, v1, v2);
}

double ad_vec3_norm(const double *v) {
	return ad_vec_norm(3, v);
}

void ad_vec3_cross(const double *a, const double *b, double *ret) {
	ret[0] = a[1]*b[2] - a[2]*b[1];
	ret[1] = -a[0]*b[2] + a[2]*b[0];
	ret[2] = a[0]*b[1] - a[1]*b[0];
}

void ad_vec3_lc2(const double a, const double *x, 
		const double b, const double *y, double *res) {
	int i=0;
	for (; i < 3; i++) res[i] = a*x[i] + b*y[i];
}

void ad_vec3_add_lc2(const double a, const double *x, 
		const double b, const double *y, double *res) {
	int i=0;
	for (; i < 3; i++) res[i] += a*x[i] + b*y[i];
}

void ad_rodr_rot(const double *z, double th, const double *v, double *ret) {
	double vpar[3];
	double vtr[3];
	double vz[3];

	/* the formula is v' = cos t * v + sin t * (z x v) 
	 * 			+ (z . v)*(1 - cos t) * z
	 */

	ad_vec3_cpy(vpar, v);
	ad_vec3_smul(vpar, cos(th));

	ad_vec3_cross(z, v, vtr);
	ad_vec3_smul(vtr, sin(th));

	ad_vec3_cpy(vz, z);
	ad_vec3_smul(vz, ad_vec3_iprod(z, v)*(1 - cos(th)));

	ad_vec3_cpy(ret, ad_vec3_zero);
	ad_vec3_add(vpar, vtr, ret);
	ad_vec3_add(vz, ret, ret);
}
