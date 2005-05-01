/* OUTCOR Core Program - written in Visual Developer Studio (MS VC++ 32 Bit)  */
/* (C) Aldo Corbellini 1 Nov. 1998. */
/* input: xr,yr,length(xr); output: $X,$Y(spline points),$z(length of $X=$Y) */
/* Note! Watch out for <storage.mode(x) <- "double"> */
/* and <as.integer(x)> for a correct working parameter passing routines */


#include <math.h>


void b_spline(x,X,y,Y,z,n)
double *x,*X;
double *y,*Y;
long *n,*z;
{
const int N=30;
int i;
int j;
double t;
long k=1;
double xa,xb,xc,xd;
double ya,yb,yc,yd;
double a3,a2,a1,a0;
double b3,b2,b1,b0;
  for (i = 1; i < *n-1; i++)  {
  xa = x[i - 1];  xb = x[i];  xc = x[i + 1];  xd = x[i + 2];
  ya = y[i - 1];  yb = y[i];  yc = y[i + 1];  yd = y[i + 2];
  a3 = ( - xa + 3 * (xb - xc) + xd)/6;   b3 = ( - ya + 3 * (yb - yc) + yd)/6;
  a2 = (xa - 2 * xb + xc)/2;             b2 = (ya - 2 * yb + yc)/2;
  a1 = (xc - xa)/2;                      b1 = (yc - ya)/2;
  a0 = (xa + 4 * xb + xc)/6;             b0 = (ya + 4 * yb + yc)/6;
         for (j = 0; j <= N; j++)  {
                   t=(float)j / (float)N;
                   *z=k;
           /*      printf("%9.2f\n",t);
                   printf("%d\n",k);
                   printf("%d\n",N);           */
                   X[k]=(((a3 * t + a2) * t + a1) * t + a0);
                   Y[k]=(((b3 * t + b2) * t + b1) * t + b0);
                   k=k+1;

      }
   }
}


void outcor(X,Out,eps,x,y,wx,wy,mx,my,nn,ll,mm)

double *X,*Out,*eps;     /* output vector containing number of outliers */
double *x,*y;				/* points outside hull 90% */
double *wx,*wy;				/* spline list of points */
double *mx,*my;				/* median point*/
long *nn,*ll,*mm;			/* counters */
{
double xa,xb,xc,xd,xk;
double ya,yb,yc,yd,yk;
double den,r,s,dp,dp1,dp2;
long n,j,k,kk;
	k = 0;
	kk = 0;
	n = 0;
	while(n < *nn) {
	xa = mx[0];     xb = x[n];
	ya = my[0];     yb = y[n];
		j = 1;
		while(j < *ll) {
		xc = wx[j - 1];   xd = wx[j];
		yc = wy[j - 1];   yd = wy[j];
		den = (xb - xa) * (yd - yc) - (yb - ya) * (xd - xc);
		r = ((ya - yc) * (xd - xc) - (xa - xc) * (yd - yc))/den;
		s = ((ya - yc) * (xb - xa) - (xa - xc) * (yb - ya))/den;
		    if((s >= 0 && s <= 1) && (r >= 0 && r <= 1)) {
				xk = xa + r * (xb - xa);
				yk = ya + r * (yb - ya);
				dp1 = sqrt((xk-xa)*(xk-xa)+(yk-ya)*(yk-ya));
				dp2 = sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));
				dp = sqrt((dp1-dp2)*(dp1-dp2));
			    	if (dp > *eps){
						X[n] = 1;
						
						k++;
					}
					else{
					X[n] = 0;

					kk++;};
			}
			j++;
		}
		n++;
	}
	 	   Out[0]=k;
	   *mm=k;
}


