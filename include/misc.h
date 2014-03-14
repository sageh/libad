#ifndef AD_MISC_H
#define AD_MISC_H

/* prints debugging messages to stderr if DEBUG is defined. */
int ad_debug_printf(char *fmt, ...);

/* prints a formatted warning to stderr prefixed by funcname and a
 * colon. */
int ad_warn(char *funcname, char *fmt, ...);

#endif /* misc.h */
