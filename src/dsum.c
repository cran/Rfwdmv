/*
 *  dsum.c
 *  Rfwdmv
 *
 *  Created by Kjell Konis on 9/13/06.
 *  Copyright 2006. All rights reserved.
 *
 */

#include "Rfwdmv.h"

double dsum(Sint n, Sfloat* x, Sint incx, Sfloat* wrkn)
{
  Sint i = 0, j = 0;

  if(n == 1)
    return x[0];

  while(i < n/2) {
    wrkn[i] = x[2*incx*i] + x[(2*i+1)*incx];
    i++;
  }

  if(2*i < n)
    wrkn[i-1] = wrkn[i-1] + x[2*incx*i];

  return dsum(i, wrkn, 1, wrkn+i);
}

