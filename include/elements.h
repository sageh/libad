#ifndef AD_ELEMENTS_H
#define AD_ELEMENTS_H

/* calculate standard orbital elements: 
 * (semi-major axis, eccentricity, inclination, longitude of ascending
 * node, argument of pericenter, mean anomaly at epoch)
 * from position, velocity and standard gravitational parameter */
void ad_rv2oe_gp(double mu, const double *r, const double *v, double *res);

/* as above, with mu = 1 */
void ad_rv2oe(const double *r, const double *v, double *res);

/* calculate position and velocity from orbital elements */
void ad_oe2rv_gp(double mu, const double *elements, double *r, double *v);
/* as above, with mu = 1 */
void ad_oe2rv(const double *elements, double *r, double *v);


#endif /* elements.h */
