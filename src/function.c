#include "function.h"
#include "misc.h"


long factorial(int n) {
	if (n < 0) {
		ad_warn("factorial", 
			"negative argument %ld treated as zero\n", n);
		return 1L;
	}
	if (n < 2) return 1L;

	{
		/* TODO: use a better algorithm for larger values of n*/
		long res=1;
		int i=2;
		for (; i <= n; i++)
			res *= i;
		return res;
	}
}

double stumpff_c(int k, double x) {
	double sqx = sqrt(x);
	double ret;
	if (k < 0) {
		ad_warn("stumpff_c", "negative order %d treated as zero\n", k);
		k = 0;
	}
	/* calculate lower orders explicitly, higher orders with
	 * recursion */
	switch (k) {
		case 0: return cos(sqrt(x));
		case 1: return (x >= 0 ? sin(sqx)/sqx : sinh(sqx)/sqx);
		case 2: return (x >= 0 ? (1.0-cos(sqx))/x : 
					(-1.0 + cosh(sqx))/x);
		case 3: return (x >= 0 ? (sqx - sin(sqx))/(sqx*sqx*sqx) :
					(-sqx + sinh(sqx))/(sqx*sqx*sqx));
		default: break;
	}
	/* k > 3, use recursion formula iteratively */
	/* TODO:
	 * - use angle doubling and halving instead?
	 * - calculate all k orders at once?
	 * - VALIDATE!
	 */ 	
	{
		long fac = factorial(k-2);
		double xpow = x;
		int n=0;
		ret = (n % 2 ? 1.0 : -1.0) / (fac*xpow);
		k -= 2;
		while (k > 3) {
			xpow *= x;
			fac /= (k-2)*(k-1);
			k -= 2;
			n++;
			ret += (n % 2 ? 1.0 : -1.0) / (fac*xpow);
		}; 
		ret += ((n+1) % 2 ? 1.0 : -1.0) / xpow * stumpff_c(k,x);
	}
	return ret;
}

double stumpff_G(int k, double b, double x) {
	if (k < 0) {
		ad_warn("stumpff_G", "negative order %d treated as zero\n", k);
		k = 0;
	}
	return pow(x, k) * stumpff_c(k, b*x*x);
}
