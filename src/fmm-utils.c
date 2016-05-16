#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "adap_fmm.h"

void FMMClean(void)
{
  free(g_fmmsources);
  free(g_fmmcharges);
  free(g_fmmtargets);
  free(g_fmmpotential);
  free(g_fmmfield);
  if (g_fmminnernormal)
    free(g_fmminnernormal);
  if (g_fmmouternormal)
    free(g_fmmouternormal);
}

void MakeUList(const int nexpo, const dcomplex * const expo, 
	       const int * const list, const int nlist, const int * const xoff, 
	       const int * const yoff, const dcomplex * const xs, 
	       const dcomplex * const ys, dcomplex * const mexpo)
{
  if (nlist) {
    for (int i = 0; i < nlist; i++) {
      int offset = list[i]*nexpo; 
      for (int j = 0; j < nexpo; j++) {
	dcomplex zmul = 1; 
	if (xoff[i] > 0) 
	  zmul *= xs[3*j + xoff[i] - 1];
	if (xoff[i] < 0) 
	  zmul *= conj(xs[3*j - xoff[i] - 1]);
	if (yoff[i] > 0) 
	  zmul *= ys[3*j + yoff[i] - 1];
	if (yoff[i] < 0) 
	  zmul *= conj(ys[3*j - yoff[i] - 1]);
	mexpo[j] += zmul*expo[offset + j];
      }
    }
  }
}

void MakeDList(const int nexpo, const dcomplex * const expo, 
	       const int * const list, const int nlist, const int * const xoff, 
	       const int * const yoff, const dcomplex * const xs, 
	       const dcomplex * const ys, dcomplex * const mexpo)
{
  if ( nlist ) {
    for (int i = 0; i < nlist; i++) {
      int offset = list[i]*nexpo; 
      for (int j = 0; j < nexpo; j++) {
	dcomplex zmul = 1; 
	if (xoff[i] > 0) 
	  zmul *= conj(xs[3*j + xoff[i] - 1]);
	if (xoff[i] < 0) 
	  zmul *= xs[3*j - xoff[i] - 1];
	if (yoff[i] > 0) 
	  zmul *= conj(ys[3*j + yoff[i] - 1]);
	if (yoff[i] < 0) 
	  zmul *= ys[3*j - yoff[i] - 1];
	mexpo[j] += zmul*expo[offset + j];
      }
    }
  }
}

void lgndr(const int nmax, const double x, double *y)
{
  int m, n;
  n = (nmax + 1)*(nmax + 1);
  for (m = 0; m < n; m++) 
    y[m] = 0.0;

  double u = -sqrt(1 - x*x);
  y[0] = 1;

  // m = 0 
  y[1] = x*y[0]; 
  for (n = 2; n <= nmax; n++) 
    y[n] = ((2*n - 1)*x*y[n - 1] - (n - 1)*y[n - 2])/n;

  // m = 1:nmax - 1
  int offset1 = nmax + 2;
  for (m = 1; m <= nmax - 1; m++) {
    int offset2 = m*offset1;
    y[offset2] = y[offset2 - offset1]*u*(2*m - 1);
    y[offset2 + 1] = y[offset2]*x*(2*m + 1);
    for (n = m + 2; n <= nmax; n++) {
      int offset3 = n + m*(nmax + 1);
      y[offset3] = ((2*n - 1)*x*y[offset3 - 1] - (n + m - 1)*y[offset3 - 2])/(n - m);
    }
  }

  // m = nmax 
  y[nmax + nmax*(nmax + 1)] = y[nmax - 1 + (nmax - 1)*(nmax + 1)]*u*(2*nmax - 1);
}

void rotz2y(const dcomplex *multipole, const double *rd, dcomplex *mrotate)
{
  dcomplex *mwork = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 1, sizeof(dcomplex));
  if (mwork == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  ephi[0] = 1.0;
  for (int m =1; m <= g_pterms; m++) 
    ephi[m] = -ephi[m - 1]*_Complex_I;
 
  for (int m = 0; m <= g_pterms; m++) {
    int offset = m*(g_pterms + 1);
    for (int ell = m; ell <= g_pterms; ell++) {
      int index = offset + ell;
      mwork[index] = ephi[m]*multipole[index];
    }
  }

  for (int m = 0; m <= g_pterms; m++) {
    int offset = m*(g_pterms + 1);
    for (int ell = m; ell <= g_pterms; ell++) {
      int index = ell + offset;
      mrotate[index] = mwork[ell]*rd[ell + (m + g_pterms)*g_pgsz]; 
      for (int mp = 1; mp <= ell; mp++) {
	int index1 = ell + mp*(g_pterms + 1);
	mrotate[index] += mwork[index1]*rd[ell + mp*(g_pterms + 1) + (m + g_pterms)*g_pgsz] +
	  conj(mwork[index1])*rd[ell + mp*(g_pterms + 1) + (-m + g_pterms)*g_pgsz];
      }
    }
  }
  free(ephi);
  free(mwork);
}

void roty2z(const dcomplex *multipole, const double *rd, dcomplex *mrotate)
{
  dcomplex *mwork = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(1 + g_pterms, sizeof(dcomplex));
  if (mwork == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  ephi[0] = 1.0;
  for (int m = 1; m <= g_pterms; m++) 
    ephi[m] = ephi[m - 1]*_Complex_I;
  
  for (int m = 0; m <= g_pterms; m++) {
    int offset = m*(g_pterms + 1);
    for (int ell = m; ell <= g_pterms; ell++) {
      int index = ell + offset;
      mwork[index] = multipole[ell]*rd[ell + (m + g_pterms)*g_pgsz];
      for (int mp = 1; mp <= ell; mp++) {
	int index1 = ell + mp*(g_pterms + 1);
	mwork[index] += multipole[index1]*rd[ell + mp*(g_pterms + 1) + (m + g_pterms)*g_pgsz] +
	  conj(multipole[index1] )*rd[ell + mp*(g_pterms + 1) + (g_pterms - m)*g_pgsz];
      }
    }
  }
 
  for (int m = 0; m <= g_pterms; m++) {
    int offset = m*(g_pterms + 1);
    for (int ell = m; ell <= g_pterms; ell++) {
      int index = ell + offset;
      mrotate[index] = ephi[m]*mwork[index];
    }
  }

  free(ephi);
  free(mwork);
}

void rotz2x(const dcomplex *multipole, const double *rd, dcomplex *mrotate)
{
  int offset1 = g_pterms*g_pgsz; 
  for (int m = 0; m <= g_pterms; m++) {
    int offset2 = m*(g_pterms + 1);
    int offset3 = m*g_pgsz + offset1;
    int offset4 = -m*g_pgsz + offset1; 
    for (int ell = m; ell <= g_pterms; ell++) {
      mrotate[ell + offset2] = multipole[ell]*rd[ell + offset3];
      for (int mp = 1; mp <= ell; mp++) {
	int offset5 = mp*(g_pterms + 1);
	mrotate[ell + offset2] += multipole[ell + offset5]*rd[ell + offset3 + offset5] + 
	  conj(multipole[ell + offset5])*rd[ell + offset4 + offset5];
      }
    }
  }
}

void in(const double scal, const double x, const int nb, double *b, int *ncalc)
{
  const double ensig = 1.0e-4; 
  const double enmten = 1.0e-300; 

  for (int i = 0; i <= nb; i++) 
    b[i] = 0; 

  if (x < 0) 
    exit(-1); 

  if (x <= ensig){
    double xscal = x/scal;
    double term1 = 1.0;
    double term2 = 0.5*x*x;
    b[0] = term1*(1.0 + term2/3.0);
    for (int i = 1; i <= nb; i++) {
      term1 = term1*xscal/(2*i + 1);
      term1 = ( term1 <= enmten ? 0.0 : term1 );
      b[i] = term1*(1.0 + term2/(2*i + 3));
    }
    *ncalc = nb + 1;
  } else if (x > 1.0e2) {
    for (int i = 0; i <= nb; i++)
      b[i] = 0;
    *ncalc = nb + 1;
  } else {    
    double constant = sqrt(M_PI_2/x);
    double alpha = 0.5;
    int ize = 1;
    int nb1 = nb + 1;
    ribesl_(&x, &alpha, &nb1, &ize, b, ncalc);
    for (int i = 0; i <= nb; i++) {
      b[i] = b[i]*constant;
      constant = constant/scal;
      constant = (fabs(b[i]) <= enmten ? 0 : constant);
    }
  }
}

void numthetahalf(int *numfour)
{
  if (g_nlambs == 9) {
    numfour[0] = 2;
    numfour[1] = 4;
    numfour[2] = 4;
    numfour[3] = 6;
    numfour[4] = 6;
    numfour[5] = 4;
    numfour[6] = 6;
    numfour[7] = 4;
    numfour[8] = 2;
  } else if (g_nlambs == 18) {
    numfour[0] = 4;
    numfour[1] = 6;
    numfour[2] = 6;
    numfour[3] = 8;
    numfour[4] = 8;
    numfour[5] = 8;
    numfour[6] = 10;
    numfour[7] = 10;
    numfour[8] = 10;
    numfour[9] = 10;
    numfour[10] = 12;
    numfour[11] = 12;
    numfour[12] = 12;
    numfour[13] = 12;
    numfour[14] = 12;
    numfour[15] = 12;
    numfour[16] = 8;
    numfour[17] = 2;
  }
}

void numthetafour(int *numphys)
{
  if (g_nlambs == 9) {
    numphys[0] = 4;
    numphys[1] = 8;
    numphys[2] = 12;
    numphys[3] = 16;
    numphys[4] = 20;
    numphys[5] = 20;
    numphys[6] = 24;
    numphys[7] = 8;
    numphys[8] = 2;
  } else if (g_nlambs == 18) {
    numphys[0] = 6;
    numphys[1] = 8;
    numphys[2] = 12;
    numphys[3] = 16;
    numphys[4] = 20;
    numphys[5] = 26;
    numphys[6] = 30;
    numphys[7] = 34;
    numphys[8] = 38;
    numphys[9] = 44;
    numphys[10] = 48;
    numphys[11] = 52;
    numphys[12] = 56;
    numphys[13] = 60;
    numphys[14] = 60;
    numphys[15] = 52;
    numphys[16] = 4;
    numphys[17] = 2;
  }
}

void bnlcft(double *c, const int p)
{
  for (int n = 0; n <= p; n++)
    c[n] = 1.0;

  for (int m = 1; m <= p; m++) {
    int offset = m*(p + 1);
    int offset1 = offset - p - 1;
    c[m + offset] = 1.0;
    for (int n = m + 1; n <= p; n++) 
      c[n + offset] = c[n - 1 + offset] + c[n - 1 + offset1];
  }

  for (int m = 1; m <= p; m++) {
    int offset = m*(p + 1);
    for (int n = m + 1; n <= p; n++) {
      c[n + offset] = sqrt(c[n + offset]);
    }
  }
}

