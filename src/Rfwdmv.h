/*
 *  rfwdmv.h
 *  rfwdmv
 *
 *  Created by Kjell Konis on 12/09/2006.
 *  Copyright 2006. All rights reserved.
 *
 */


#include "S.h"
#include "Rmath.h"
#include "R_ext/Utils.h"

#ifdef USING_R
  typedef double Sfloat;
  typedef int Sint;
  #define SINT_MAX INT_MAX
  #define SINT_MIN INT_MIN
#else
  typedef double Sfloat;
  typedef long Sint;
  #define SINT_MAX LONG_MAX
  #define SINT_MIN LONG_MIN
#endif

typedef enum {false, true} bool;

void FSfwdmv(Sfloat* X, Sint* pn, Sint* pp, Sint* bsb, Sint* nbsb,
             Sint* monitor, Sint* pscaled, Sfloat* Distances,
             Sfloat* Center, Sfloat* Cov, Sfloat* Determinant,
             Sint* Unit, Sfloat* Max, Sfloat* Mth, Sfloat* Min,
             Sfloat* Mpo, Sfloat* QR, Sfloat* A, Sfloat* center,
             Sfloat* distances, Sfloat* wrk1n, Sfloat* wrk2n,
             Sint* newbsb);

void FSfwdmvRandomStart(Sfloat* X, Sint* pn, Sint* pp, Sint* bsbmat,
                        Sint* pnstarts, Sint* ppbsb, Sint* pscaled,
                        Sfloat* mpomat);

/* Uncategorized */

double dsum(Sint, Sfloat*, Sint, Sfloat*);

/* The BLAS and LAPACK */

void F77_NAME(dcopy)(Sint*, Sfloat*, Sint*, Sfloat*, Sint*);
void F77_NAME(dgeqrf)(Sint*, Sint*, Sfloat*, Sint*, Sfloat*, Sfloat*, Sint*, Sint*);
void F77_CALL(dtrsm)(char*, char*, char*, char*, Sint*, Sint*, Sfloat*, Sfloat*, Sint*, Sfloat*, Sint*);
void F77_NAME(dtptrs)(char*, char*, char*, Sint*, Sint*, Sfloat*, Sfloat*, Sint*, Sint*);
double F77_NAME(ddot)(Sint*, Sfloat*, Sint*, Sfloat*, Sint*);
void F77_NAME(drotg)(Sfloat*, Sfloat*, Sfloat*, Sfloat*);
void F77_NAME(drot)(Sint*, Sfloat*, Sint*, Sfloat*, Sint*, Sfloat*, Sfloat*);
void F77_CALL(dgemm)(char*, char*, Sint*, Sint*, Sint*, Sfloat*, Sfloat*, Sint*, Sfloat*,
                     Sint*, Sfloat*, Sfloat*, Sint*);
double F77_CALL(dnrm2)(Sint*, Sfloat*, Sint*);

/*In R header file somewhere

double R_pow_di(Sfloat, Sint);
double R_pow(Sfloat, Sfloat);

*/


