#include "misc.h"

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>


int ad_debug_printf(char *fmt, ...) {
#ifdef DEBUG
	int ret;
	va_list ap;
	va_start(ap, fmt);
	ret = vfprintf(stderr, fmt, ap);
	va_end(ap);
	return ret;
#else
	return 0;
#endif
}

int ad_warn(char *fname, char *fmt, ...) {
	int ret;
	va_list ap;

	fprintf(stderr, "%s: ", fname);
	va_start(ap, fmt);
	ret = vfprintf(stderr, fmt, ap);
	va_end(ap);
	return ret;
}

void ad_error(char *fname, char *fmt, ...) {
	ad_warn(fname, fmt);
	exit(EXIT_FAILURE);
}

