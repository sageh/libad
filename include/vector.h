#ifndef AD_VECTOR_H
#define AD_VECTOR_H

/* n-vectors */

void ad_vec_cpy(int n, double *to, const double *from);
void ad_vec_add(int n, const double *v1, const double *v2, double *res);
void ad_vec_sub(int n, const double *v1, const double *v2, double *res); 
void ad_vec_smul(int n, double *v, double s); 
double ad_vec_iprod(int n, const double *v1, const double *v2); 
double ad_vec_norm(int n, const double *v); 

/* 3-vectors */

extern const double ad_vec3_zero[]; /* zero vector */
void ad_vec3_cpy(double *to, const double *from);
void ad_vec3_add(const double *v1, const double *v2, double *res);
void ad_vec3_sub(const double *v1, const double *v2, double *res); 
void ad_vec3_smul(double *v, double s); 
double ad_vec3_iprod(const double *v1, const double *v2); 
double ad_vec3_norm(const double *v); 
void ad_vec3_cross(const double *a, const double *b, double *ret);
/* linear combination of 2 vectors. */
void ad_vec3_lc2(const double a, const double *x, 
		const double b, const double *y, double *res);
void ad_vec3_add_lc2(const double a, const double *x, 
		const double b, const double *y, double *res);

/* Rodrigues' rotation formula. */
void ad_rodr_rot(const double *axis, double angle, 
		const double *vec, double *res);


#endif /* vector.h */
