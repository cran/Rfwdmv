/*
 *  FSfwdmv.c
 *  FSfwdmv
 *
 *  Created by Kjell Konis on 12/09/2006.
 *  Copyright 2006. All rights reserved.
 *
 */

#include "Rfwdmv.h"

void FSfwdmv(Sfloat* X, Sint* pn, Sint* pp, Sint* bsb, Sint* nbsb,
             Sint* monitor, Sint* pscaled, Sfloat* Distances,
             Sfloat* Center, Sfloat* Cov, Sfloat* Determinant,
             Sint* Unit, Sfloat* Max, Sfloat* Mth, Sfloat* Min,
             Sfloat* Mpo, Sfloat* QR, Sfloat* A, Sfloat* center,
             Sfloat* distances, Sfloat* wrk1n, Sfloat* wrk2n,
             Sint* newbsb)
{
	Sint n = *pn, p = *pp, m = *nbsb, scaled = *pscaled;
  char char_N = 'N', char_U = 'U', char_R = 'R';
  Sfloat DZERO = 0.0, DONE = 1.0;
  Sint IONE = 1;
	Sint i = 0, j = 0, k = 0, ifail = 0, step = 0, nsteps = n-m+1;
  Sfloat det = 1.0;

  Sint* temp = NULL;
  Sfloat* cov = NULL;
  if(monitor[2])
    cov = Salloc(p*p, Sfloat);

	for(m = m; m <= n; m++) {

    /* copy the rows of X which are in the subset to A */
    for(i = 0; i < m; i++)
      F77_CALL(dcopy)(&p, X+bsb[i], &n, A+i, &n);

    /* compute the means of the columns of A */
    for(j = 0; j < p; j++)
      center[j] = dsum(m, A+j*n, 1, wrk1n) / (Sfloat) m;

    /* sweep means from A */
    for(j = 0; j < p; j++)
      for(i = 0; i < m; i++)
        A[i+j*n] -= center[j];

    /* compute the QR factorization of A */
    for(j = 0; j < p; j++)
      F77_CALL(dcopy)(&m, A+j*n, &IONE, QR+j*n, &IONE);
    F77_CALL(dgeqrf)(&m, &p, QR, &n, wrk1n, wrk2n, &n, &ifail);

    /* center X and store in A */
    for(j = 0; j < p; j++)
      for(i = 0; i < n; i++)
        A[i+j*n] = X[i+j*n] - center[j];

    /* monitor determinant */
    if(scaled || monitor[3]) {
      det = 1.0;
      for(i = 0; i < p; i++)
        det = det * QR[i+i*n];
      det = det < 0.0 ? -1.0*det : det;
    }
    
    if(monitor[3])
      Determinant[step] = (det*det) / R_pow_di((Sfloat) (m-1), p);

    if(scaled)
      det = R_pow(det, 1.0 / (Sfloat) p) / sqrt((Sfloat) (m-1));

    /* compute Mahalanobis distances for all units */
    F77_CALL(dtrsm)(&char_R, &char_U, &char_N, &char_N, &n, &p, &DONE, QR, &n, A, &n);
    if(scaled)
      for(i = 0; i < n; i++)
        distances[i] = det * sqrt((Sfloat)(m-1) * F77_CALL(ddot)(&p, A+i, &n, A+i, &n));
    else
      for(i = 0; i < n; i++)
        distances[i] = sqrt((Sfloat)(m-1) * F77_CALL(ddot)(&p, A+i, &n, A+i, &n));

    /* monitor distances */
    if(monitor[0])
      F77_CALL(dcopy)(&n, distances, &IONE, Distances+step*n, &IONE);

    /* monitor center */
    if(monitor[1])
      F77_CALL(dcopy)(&p, center, &IONE, Center+step, &nsteps);

    /* monitor cov */
    if(monitor[2]) {
      /* zero the lower triangle of the R factor in QR */
      for(i = 1; i < p; i++)
        for(j = 0; j < i; j++)
          QR[i+j*n] = 0.0;

      /* copy the transpose of the R factor into the top of A */
      for(i = 0; i < p; i++)
        F77_CALL(dcopy)(&p, QR+i*n, &IONE, A+i, &n);

      /* compute the covariance matrix (t(R) %*% R)/(m-1) */
      F77_CALL(dgemm)(&char_N, &char_N, &p, &p, &p, &DONE, A, &n, QR, &n, &DZERO, cov, &p);
      for(i = 0; i < p*p; i++)
        cov[i] = cov[i] / ((Sfloat) (m-1));

      /* copy cov into Cov in packed storage */
      k = 0;
      for(j = 0; j < p; j++)
        for(i = 0; i <= j; i++) {
          Cov[step+k*nsteps] = cov[i+j*p];
          k++;
        }
      }

    /* monitor unit */
    if(monitor[4]) {
      for(i = 0; i < m; i++)
        Unit[step*n+bsb[i]] = 1;
    }

    /* monitor max */
    if(monitor[5]) {
      Max[step] = 0.0;
      for(i = 0; i < m; i++)
        Max[step] = distances[bsb[i]] > Max[step] ? distances[bsb[i]] : Max[step];
    }

    /* monitor min */
    if(monitor[7] && m < n) {
      R_isort(bsb, m);
      Min[step] = DBL_MAX;
      k = 0;
      for(i = 0; i < n; i++) {
        if(k < m && bsb[k] == i)
          k++;
        else {
          Min[step] = distances[i] < Min[step] ? distances[i] : Min[step];
        }
      }
    }

    /* update subset*/
    for(i = 0; i < n; i++)
      newbsb[i] = i;
    rsort_with_index(distances, newbsb, n);

    /* monitor mth */
    if(monitor[6])
      Mth[step] = distances[m-1];

    /* monitor mpo */
    if(monitor[8] && m < n)
      Mpo[step] = distances[m];

    /* get ready for next loop */
    temp = bsb;
    bsb = newbsb;
    newbsb = temp;
    step++;
  }
}



