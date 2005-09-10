/*
 *  FSfwdmvRandomStarts.c
 *  Rfwdmv
 *
 *  Created by Kjell Konis on 9/14/06.
 *  Copyright 2006. All rights reserved.
 *
 */

#include "Rfwdmv.h"

void FSfwdmvRandomStart(Sfloat* X, Sint* pn, Sint* pp, Sint* bsbmat, Sint* pnstarts, Sint* ppbsb,
                        Sint* pscaled, Sfloat* mpomat)
{
  Sint n = *pn, p = *pp, nstarts = *pnstarts, pbsb = *ppbsb, scaled = *pscaled;
  Sint start = 0, j = 0;
  Sint IONE = 1;
  Sint steps = n-pbsb;

  /* monitor only mpo */
  Sint monitor[] = {0, 0, 0, 0, 0, 0, 0, 0, 1};
  Sfloat Distances = 0.0;
  Sfloat Center = 0.0;
  Sfloat Cov = 0.0;
  Sfloat Determinant = 0.0;
  Sint Unit = 0;
  Sfloat Max = 0.0;
  Sfloat Mth = 0.0;
  Sfloat Min = 0.0;

  /* required working space */
  Sfloat* QR = Salloc(n*p, Sfloat);
  Sfloat* A = Salloc(n*p, Sfloat);
  Sfloat* center = Salloc(p, Sfloat);
  Sfloat* distances = Salloc(n, Sfloat);
  Sfloat* wrk1n = Salloc(n, Sfloat);
  Sfloat* wrk2n = Salloc(n, Sfloat);
  Sint* bsb = Salloc(n, Sint);
  Sint* newbsb = Salloc(n, Sint);

  for(start = 0; start < nstarts; start++) {

    /* copy row from bsbmat in bsb */
    for(j = 0; j < pbsb; j++)
      bsb[j] = bsbmat[start+j*nstarts];

    /* compute forward search and store mpo */
    FSfwdmv(X, &n, &p, bsb, &pbsb, monitor, &scaled, &Distances, &Center, &Cov,
            &Determinant, &Unit, &Max, &Mth, &Min, mpomat+start*steps, QR, A,
            center, distances, wrk1n, wrk2n, newbsb);
  }
}


