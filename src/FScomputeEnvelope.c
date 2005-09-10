/*
 *  FScomputeEnvelope.c
 *  Rfwdmv
 *
 *  Created by Kjell Konis on 22/09/2006.
 *  Copyright 2006. All rights reserved.
 *
 */

#include "Rfwdmv.h"

Sfloat FScdfk(Sfloat x, Sint n, Sint p, Sint k, Sfloat* nCj, Sfloat* wrk1,
              Sfloat* wrk2)
{
  Sint i = 0;
  Sfloat G = 0.0, oneMinusG  = 0.0;

  G = pf(x, (Sfloat) p, (Sfloat) (k-p-1), 1, 0);
  oneMinusG = 1.0 - G;

  for(i = k; i <= n; i++)
    wrk1[i] = nCj[i] * R_pow_di(G, i) * R_pow_di(oneMinusG, (n-i));

  return dsum(n-k+1, wrk1+k, 1, wrk2);
}


void FScomputeEnvelope(Sint* pn, Sint* pp, Sint* pm, Sfloat* probs, Sint* pnprobs,
                       Sfloat* quantiles, Sfloat* ptol)
{
  Sint k = 0;
  Sint m = *pm, n = *pn, p = *pp, nprobs = *pnprobs;
  Sint i = 0, j = 0, nsteps = 0;
  Sint IONE = 1;

  Sfloat tol = *ptol;
  Sfloat d_eps = 0.0, p_eps = 0.0, q_eps = 0.0, sigma = 0.0;
  Sfloat x = 0.0, u = 0.0, leftx = 0.0, rightx = 0.0, leftu = 0.0, rightu = 0.0;
  Sfloat delta = 0.0;

  Sfloat* nCj = Salloc(n+1, Sfloat);
  Sfloat* wrk1 = Salloc(n+1, Sfloat);
  Sfloat* wrk2 = Salloc(n+1, Sfloat);

  nsteps = n-m;

  for(i = 0; i <= n; i++)
    nCj[i] = choose((Sfloat) n, (Sfloat) i);

  for(i = m; i < n; i++) {

    /* compute asymptotic result */
    p_eps = ((Sfloat) i-(3.0/8.0)+1.0)/((Sfloat) n-(3.0/4.0)+1.0);
    q_eps = qchisq(p_eps, (Sfloat) p, 1, 0);
    d_eps = R_pow_di(dchisq(q_eps, p, 0), 2);
    sigma = sqrt(p_eps*(1.0-p_eps)/((Sfloat) n*d_eps));

    /* improve estimate using exact proceedure */

    for(j = 0; j < nprobs; j++) {

      /* use asymptotic result as initial value */
      leftx = q_eps + qnorm(probs[j], 0.0, sigma, 1, 0);
      leftu = FScdfk(leftx, n, p, i, nCj, wrk1, wrk2);
      rightx = 0.0;
      rightu = 0.0;

      /* bracket the target probability */
      if(leftu > probs[j]) {
        rightx = leftx;
        leftx = 0.0;
        rightu = leftu;
        leftu = 0.0;
      }
      else {
        rightx = 2.0*leftx;
        rightu = FScdfk(rightx, n, p, i, nCj, wrk1, wrk2);
        while(rightu < probs[j]) {
          leftx = rightx;
          leftu = rightu;
          rightx = 2.0*leftx;
          rightu = FScdfk(rightx, n, p, i, nCj, wrk1, wrk2);
        }
      }

      /* find x with a bisection search */
      while(rightu-leftu > tol) { 
        x = (leftx+rightx)/2.0;
        u = FScdfk(x, n, p, i, nCj, wrk1, wrk2);

        if(u < probs[j]) {
          leftx = x;
          leftu = u;
        }
        else {
          rightx = x;
          rightu = u;
        }
      }
      quantiles[i-m+j*nsteps] = (leftx+rightx)/2.0;
    }
  }
}
