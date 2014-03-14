#ifndef AD_STUMPFF_H
#define AD_STUMPFF_H

#include <math.h>

/* special functions. */

/* calculates factorial of n */
long factorial(int n);

/* calculates Stumpff's c-function of order k */
double stumpff_c(int k, double x);

/* calculates Stumpff's G-function of order k, with parameter beta */
double stumpff_G(int k, double beta, double x);

#endif /* stumpff.h */
