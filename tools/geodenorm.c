/* denormalizes geopotential coefficients from the geodesic standard to
 * the conventional gravitational ones. the relation between normalized
 * coefficients An_nm (A = C or S) and denormalized coefficients A_nm
 * is An_nm = sqrt( (n+m)! / ( (n-m)! * (2n+1) * k ) ) A_nm, where k = 1
 * if m = 0 and k = 2 otherwise. 
 * 
 * reads input from stdin, outputs to stdout.
 *
 * input format: n m Cn_nm Sn_nm
 * n, m: integers
 * Cn_nm, Sn_nm: valid ANSI C floating point numbers
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
	char buf[200];
	char fmtout[] = "%5d %5d %25.15e %25.15e\n";
	int read;
	int n, m, k, i;
	double c, s, denorm;
	int in=0, out=0;

	for (;;) {
		fgets(buf, 200, stdin);
		read = sscanf(buf, "%d %d %lf %lf", &n, &m, &c, &s);
		in++;

		if (read != 4)
			break;

		k = (m == 0 ? 1 : 2);
		denorm = (2.0*n+1)*k;
		for (i=n-m+1; i <= n+m; i++)
			denorm *= 1.0/i;
		denorm = sqrt(denorm);
		printf(fmtout, n, m, c*denorm, s*denorm);
		out++;
	}

	fprintf(stderr, "read %d lines, wrote %d lines\n", in, out);

	if (read != EOF) {
		fprintf(stderr, "exit due to invalid input\n");
		return EXIT_FAILURE;
	}

	fprintf(stderr, "conversion successful\n");
	return EXIT_SUCCESS;
}
