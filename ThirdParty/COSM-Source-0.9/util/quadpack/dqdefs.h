#ifndef _DQ_DEFS_H
#define _DQ_DEFS_H

#define uflow   DBL_MIN
#define oflow   DBL_MAX
#define epmach  DBL_EPSILON
#define LIMIT   500
#define MAXP1   21
#define _USE_MATH_DEFINES
#define Pi      M_PI
#define COSINE  1
#define SINE    2
                                                                             
#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif

#endif  // _DQ_DEFS_H
