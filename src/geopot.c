#include "geopot.h"
#include "misc.h"

int ad_geo_idx(int n, int m) {
	/* correct base row is found from the sequence
	 * 0,3,3+4,3+4+5,3+4+5+6, ..., 3+4+...+n = (n-2)(n+3)/2 (when n >= 2) */
	if (n < 2 || m < 0 || m > n) {
		ad_warn("ad_geo_idx", 
			"expected n >= 2, m >= 0, got n: %d m: %d\n", n, m);
		return 0; /* a safe-ish return value */
	}
	return (n-2)*(n+3)/2 + m;
}
