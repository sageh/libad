#include "elements.h"
#include "vector.h"
#include "kepler.h"
#include "misc.h"

#include <math.h>

void ad_rv2oe_gp(double mu, const double *rv, const double *vin, double *res) {
	double vv[3], ev[3], hv[3], rhat[3];
	double r = ad_vec3_norm(rv);
	double alpha, eta, zeta, h, b;

	/* vv = v/mu */
	ad_vec3_cpy(vv, vin);
	ad_vec3_smul(vv, 1/sqrt(mu));

	/* alpha = 2/r - vv . vv, reciprocal of semimajor axis */
	alpha = 2/r - ad_vec3_iprod(vv, vv);
	res[0] = 1/alpha;

	/* eta = rv . vv */
	eta = ad_vec3_iprod(rv, vv);

	/* zeta = 1 - alpha * r */
	zeta = 1 - alpha*r;

	/* hv = rv x vv, specific angular momentum */
	ad_vec3_cross(rv, vv, hv);
	h = ad_vec3_norm(hv);

	/* rhat = rv / r, unit vector in r direction */
	ad_vec3_cpy(rhat, rv);
	ad_vec3_smul(rhat, 1/r);

	/* ev = vv x hv - rhat, eccentricity vector */
	ad_vec3_cross(vv, hv, ev);
	ad_vec3_sub(ev, rhat, ev);
	res[1] = ad_vec3_norm(ev);

	/* length of the projection of angular momentum on xy-plane */
	b = sqrt(hv[0]*hv[0] + hv[1]*hv[1]);

	/* inclination */
	res[2] = atan2(b, hv[2]);

	/* longitude of ascending node */
	res[3] = atan2(hv[0], -hv[1]);

	/* argument of pericentre */
	res[4] = atan2(ev[2]*h, ev[1]*hv[0] - ev[0]*hv[1]);

	/* mean anomaly at epoch */
	res[5] = atan2(eta*sqrt(alpha), zeta) - eta*sqrt(alpha);
}

void ad_rv2oe(const double *r, const double *v, double *res) {
	ad_rv2oe_gp(1.0, r, v, res);
}

void ad_oe2rv_gp(double mu, const double *el, double *rv, double *vv) {
	double pv[3], qv[3], wv[3]; /* orbital plane coordinate system axes */
	double av[3], bv[3]; /* semi-major & minor axes */
	double E, r, vc;

	pv[0] =  cos(el[4])*cos(el[3]) - sin(el[4])*sin(el[3])*cos(el[2]);
	pv[1] =  cos(el[4])*sin(el[3]) + sin(el[4])*cos(el[3])*cos(el[2]);
	pv[2] =  sin(el[4])*sin(el[2]);

	qv[0] = -sin(el[4])*cos(el[3]) - cos(el[4])*sin(el[3])*cos(el[2]);
	qv[1] = -sin(el[4])*sin(el[3]) + cos(el[4])*cos(el[3])*cos(el[2]);
	qv[2] =  cos(el[4])*sin(el[2]);

	wv[0] =  sin(el[3])*sin(el[2]);
	wv[1] = -cos(el[3])*sin(el[2]);
	wv[2] =  cos(el[2]);

	/* TODO: make tolerance and maxit settable */
	/* solve elliptic kepler M = E - e*sin(E) */
	if (ad_solve_ekepler(1e-15, 10, el[1], el[5], &E) != 0)
		ad_warn("ad_oe2rv_gp", "ad_solve_ekepler: no convergence\n");

	/* av = a * pv */
	ad_vec3_cpy(av, pv);
	ad_vec3_smul(av, el[0]);

	/* bv = sqrt(1-e^2) * wv x av */
	ad_vec3_cross(wv, av, bv);
	ad_vec3_smul(bv, sqrt(1-el[1]*el[1]));

	/* rv = (cos E - e) * av + sin E * bv */
	ad_vec3_lc2(cos(E)-el[1], av, sin(E), bv, rv);
	r = ad_vec3_norm(rv);

	/* vv = sqrt(mu/(r*sqrt(a))) * [-sin E * av + cos E * bv ] */
	vc = sqrt(mu)/(r*sqrt(el[0]));
	ad_vec3_lc2(-vc*sin(E), av, cos(E), bv, vv);
}

void ad_oe2rv(const double *elements, double *r, double *v) {
	ad_oe2rv_gp(1.0, elements, r, v);
}

