#include "pn.h"
#include "vector.h"

#include <math.h>

void ad_pn(double gc, double c, double m1, double m2, 
		double *rv, double *vv, double *a) {
	int i;
	double m = m1+m2, eta = m1*m2/(m*m), gm = gc*m, gmperc2 = gm/(c*c);
	double r = ad_vec3_norm(rv), v = ad_vec3_norm(vv);
	double rd = ad_vec3_iprod(rv, vv)/r;
	double nbrac[3], vbrac[3];
	double a1[3], a2[3], a2_5[3];

	/* a1 */
	nbrac[0] = -2*(2+eta)*gm/r + (1+3*eta)*v*v - 3/2*eta*rd*rd;
	vbrac[0] = -2*(2-eta)*rd;

	/* a2 */
	nbrac[1] = 3/4*(12+29*eta)*gm*gm/(r*r) + eta*(3-4*eta)*pow(v,4)
		+ 15/8*eta*(1-3*eta)*pow(rd,4)
		- 3/2*eta*(3-4*eta)*v*v*rd*rd
		- 1/2*eta*(13-4*eta)*(gm/r)*v*v
		- (2+25*eta+2*eta*eta)*(gm/r)*rd*rd;
	vbrac[1] = eta*(15+4*eta)*v*v - (eta+41*eta+8*eta*eta)*(gm/r)
		- 3*eta*(3+2*eta)*rd*rd;

	/* a2.5 */
	nbrac[2] = 9*v*v + 17*gm/r;
	vbrac[2] = 3*v*v + 9*gm/r;

	/* collect contributions */
	ad_vec3_lc2(nbrac[0]/r, rv, vbrac[0], vv, a1);
	ad_vec3_smul(a1, -gmperc2/(r*r));

	ad_vec3_lc2(nbrac[1]/r, rv, -(rd/2)*vbrac[1], vv, a2);
	ad_vec3_smul(a2, -gmperc2/(c*c*r*r));

	ad_vec3_lc2(nbrac[2]*rd/r, rv, -vbrac[2], vv, a2_5);
	ad_vec3_smul(a2_5, 8/15*gmperc2*gmperc2*eta/(c*r*r));

	for (i=0; i < 3; i++) 
		a[i] = a1[i]+a2[i]+a2_5[i];
}

void ad_spin_pn(double gc, double c, double m1, double m2, 
		double chi, double q, 
		double *rv, double *vv, double *s, double *a) {
	int i;
	double m = m1+m2, eta = m1*m2/(m*m), gm = gc*m;
	double r = ad_vec3_norm(rv);
	double rd = ad_vec3_iprod(rv, vv)/r;
	double nbrac[2], nsbrac, vsbrac;
	double nv[3], ns[3], vs[3], nvvv[3], nvs;
	double SO[3], Q[3];

	/* nv = rhat */
	ad_vec3_cpy(nv, rv);
	ad_vec3_smul(nv, 1/r);

	/* ns = nv x s */
	ad_vec3_cross(nv, s, ns);

	/* vs = vv x s */
	ad_vec3_cross(vv, s, vs);

	/* nvvv = nv x vv */
	ad_vec3_cross(nv, vv, nvvv);

	/* nvs = nv.s */
	nvs = ad_vec3_iprod(nv, s);

	/* SO */
	nbrac[0] = 12 * ad_vec3_iprod(s, nvvv);
	nsbrac = (9+3*sqrt(1-4*eta))*rd;
	vsbrac = 7 + sqrt(1-4*eta);

	/* Q */
	nbrac[1] = 5 * nvs*nvs - 1;

	/* collect results */
	for (i=0; i < 3; i++) {
		SO[i] = gm*gm/(r*r*r*c*c*c)*(1+sqrt(1-4*eta))/4 * chi
			* (nbrac[0]*nv[i] + nsbrac*ns[i] - vsbrac*vs[i]);
		Q[i] = -q*chi*chi*3*gc*gc*gc*m1*m1*m/(2*pow(c,4)*pow(r,4))
			* (nbrac[1]*nv[i] - 2*nvs*s[i]);
	}
}

