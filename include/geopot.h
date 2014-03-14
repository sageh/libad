#ifndef AD_GEOPOT_H
#define AD_GEOPOT_H

#include "geopot_egm08.h"

/* find geopotential table index from given degree and order */
int ad_geo_idx(int n, int m);

/* same, as a more unsafe macro */
#define AD_GEO_IDX(n,m) ((((n)-2)*((n)+3)/2) + (m))

#endif /* geopot.h */
