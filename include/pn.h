#ifndef AD_PN_H
#define AD_PN_H

/* calculate sum of non-spin dependent PN acceleration terms of orders
 * 1, 2 & 2.5. vector sum is returned in a. */
void ad_pn(double gc, double c, double m1, double m2, 
		double *rv, double *vv, double *a);

/* calculate sum of spin dependent PN acceleration terms SO and (orders
 * 1.5 and 1). spin vector in s is to be given as a cartesian 3-vector
 * sum is returned in a. */
void ad_spin_pn(double gc, double c, double m1, double m2, 
		double chi, double q, 
		double *rv, double *vv, double *s, double *a);

#endif /* pn.h */
