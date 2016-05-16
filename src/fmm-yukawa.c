#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "adap_fmm.h"

void YukFMMInit(const double beta)
{
  g_beta = beta;

  g_yifl[0] = 3; 
  g_yifl[1] = 4; 
  g_yifl[2] = 2;
  g_yifl[3] = 1; 
  g_yifl[4] = 3; 
  g_yifl[5] = 4; 
  g_yifl[6] = 2; 
  g_yifl[7] = 1; 

  g_yuknumfour = (int *)calloc(g_nlambs, sizeof(int)); 
  g_yuknumphys = (int *)calloc(g_nlambs*g_nslev, sizeof(int)); 
  g_yukwhts = (double *)calloc(g_nlambs, sizeof(double));
  g_yukrlams = (double *)calloc(g_nlambs, sizeof(double));
  g_yukytop = (double *)calloc(g_pgsz, sizeof(double));
  g_yukrdplus = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  g_yukrdminus = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  g_yukrdsq3 = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  g_yukrdmsq3 = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));

  if (g_yuknumfour == 0 || g_yuknumphys == 0 || g_yukwhts == 0 || 
      g_yukrlams == 0 || g_yukytop == 0 || g_yukrdplus == 0 || g_yukrdminus == 0 ||
      g_yukrdsq3 == 0 || g_yukrdmsq3 == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  yhfrmini();
  yhrotgen();
  yukvwts();
  numthetahalf(g_yuknumfour); 
  numthetafour(g_yuknumphys); 

  for (int i = 0; i < g_nlambs; i++) {
    double test1 = g_yukrlams[i]; 
    double test2 = sqrt(test1*test1 + 2.0*test1*g_beta*g_size);
    int indd = i; 
    int mmax = g_yuknumphys[i]; 
    for (int j = i; j < g_nlambs; j++) {
      if ( test2 <= g_yukrlams[j] ) {
	indd = j; 
	break; 
      } else {
	mmax = MAX(g_yuknumphys[j], mmax); 
      }
    }
    g_yuknumphys[i] = MAX(g_yuknumphys[indd], mmax);
  }

  g_yuknexptot = 0; 
  g_yuknthmax = 0; 
  g_yukmnexptotp = 0; 

  for (int i = 0; i < g_nlambs; i++) {
    g_yuknexptot += g_yuknumfour[i]; 
    if (g_yuknumfour[i] > g_yuknthmax) 
      g_yuknthmax = g_yuknumfour[i]; 
    g_yukmnexptotp += g_yuknumphys[i]; 
  }

  g_yukmnexptotp /= 2.0; 
  g_yuknexpmax = MAX(g_yuknexptot, g_yukmnexptotp) + 1;
  //g_yukdcpgsz = pow(2*g_pterms + 1, 3);
  g_yukdcpgsz = (2*g_pterms + 1)*(2*g_pterms+1)*(2*g_pterms+1);

  yuk_multipole = (dcomplex *)calloc((1 + g_nsboxes)*g_pgsz, sizeof(dcomplex));
  yuk_local = (dcomplex *)calloc((1 + g_ntboxes)*g_pgsz, sizeof(dcomplex)); 
  yuk_expu = (dcomplex *)calloc((1 + g_nsboxes)*g_yuknexpmax, sizeof(dcomplex)); 
  yuk_expd = (dcomplex *)calloc((1 + g_nsboxes)*g_yuknexpmax, sizeof(dcomplex)); 
  yuk_expn = (dcomplex *)calloc((1 + g_nsboxes)*g_yuknexpmax, sizeof(dcomplex)); 
  yuk_exps = (dcomplex *)calloc((1 + g_nsboxes)*g_yuknexpmax, sizeof(dcomplex)); 
  yuk_expe = (dcomplex *)calloc((1 + g_nsboxes)*g_yuknexpmax, sizeof(dcomplex)); 
  yuk_expw = (dcomplex *)calloc((1 + g_nsboxes)*g_yuknexpmax, sizeof(dcomplex)); 

  yuk_fexpe = (dcomplex *)calloc(15000*g_nslev, sizeof(dcomplex));
  yuk_fexpo = (dcomplex *)calloc(15000*g_nslev, sizeof(dcomplex));
  yuk_fexpback = (dcomplex *)calloc(15000*g_nslev, sizeof(dcomplex));
  g_yukvnexptot = (int *)calloc(g_nslev, sizeof(int));
  g_yukvnexptotp = (int *)calloc(1 + g_nslev, sizeof(int)); 
  g_yukvnthmax = (int *)calloc(g_nslev, sizeof(int)); 
  g_yukdcu = (double *)calloc(g_yukdcpgsz*(1 + g_nslev), sizeof(double)); 
  g_yukdcd = (double *)calloc(g_yukdcpgsz*(1 + g_nslev), sizeof(double));
  g_yuksfactor = (double *)calloc(1 + g_nslev, sizeof(double)); 
  g_yuksfactor2 = (double *)calloc(1 + g_nslev, sizeof(double)); 
  g_yukbetascale = (double *)calloc(1 + g_nslev, sizeof(double)); 
  yuk_zs = (double *)calloc(g_yuknexpmax*3*g_nslev, sizeof(double)); 
  g_yukrlsc = (double *)calloc(g_pgsz*g_nlambs*g_nslev, sizeof(double)); 
  yuk_xs = (dcomplex *)calloc(g_yuknexpmax*3*g_nslev, sizeof(dcomplex));
  yuk_ys = (dcomplex *)calloc(g_yuknexpmax*3*g_nslev, sizeof(dcomplex)); 

  if (yuk_multipole == 0 || yuk_local == 0 || yuk_expu == 0 || yuk_expd == 0 || 
      yuk_expn == 0 || yuk_exps == 0 || yuk_expe == 0 || yuk_expw == 0 || 
      yuk_fexpe == 0 || yuk_fexpo == 0 || yuk_fexpback == 0 || g_yukvnexptot == 0 || 
      g_yukvnexptotp == 0 || g_yukvnthmax == 0 || g_yukdcu == 0 || g_yukdcd == 0 || 
      g_yuksfactor == 0 || g_yuksfactor2 == 0 || g_yukbetascale == 0 || yuk_zs == 0 || 
      g_yukrlsc == 0 || yuk_xs == 0 || yuk_ys == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  g_yukscale = (g_beta*g_size > 1.0 ? 1.0/g_size : g_beta);

  for (int i = 0; i <= g_nslev; i++) {
    double size = g_size/pow(2, i); 
    double r0 = sqrt(3.0)*size/4; 
    g_yukbetascale[i] = g_beta*size/2.0; 
    g_yuksfactor[i] = g_yukscale*size; 
    g_yuksfactor2[i] = g_yuksfactor[i]*0.5; 
    ympshftcoef(r0, i); 
    ylcshftcoef(r0, i); 
  }

  for (int lev = 0; lev < g_nslev; lev++) {
    numthetafour(&g_yuknumphys[lev*g_nlambs]);
    for (int i = 0; i < g_nlambs; i++) {
      double test1 = g_yukrlams[i]; 
      double test2 = sqrt(test1*test1 + 2.0*test1*g_yuksfactor2[lev]); 
      int indd = i; 
      int mmax = g_yuknumphys[lev*g_nlambs + i]; 
      for (int j = i; j < g_nlambs; j++ ) {
	if ( test2 <= g_yukrlams[j] ) {
	  indd = j; 
	  break; 
	} else {
	  mmax = MAX(mmax, g_yuknumphys[lev*g_nlambs + j]);
	}
      }
      g_yuknumphys[lev*g_nlambs + i] = MAX(mmax, g_yuknumphys[lev*g_nlambs + indd]);
    }

    g_yukvnexptot[lev] = 0;
    g_yukvnthmax[lev] = 0; 
    g_yukvnexptotp[lev] = 0; 

    for (int i = 0; i < g_nlambs; i++) {
      g_yukvnexptot[lev] += g_yuknumfour[i]; 
      if (g_yuknumfour[i] > g_yukvnthmax[lev]) 
	g_yukvnthmax[lev] = g_yuknumfour[i]; 
      g_yukvnexptotp[lev] += g_yuknumphys[lev*g_nlambs + i]; 
    }

    g_yukvnexptotp[lev] /= 2.0; 

    ymkfexp(lev); 
    yrlscini(lev); 
    ymkexps(lev); 
  }
}

void YukFMMClean(void)
{
  free(yuk_xs);
  free(yuk_ys);
  free(yuk_zs);
  free(yuk_fexpe);
  free(yuk_fexpo);
  free(yuk_fexpback);

  free(yuk_multipole);
  free(yuk_local);
  free(yuk_expu);
  free(yuk_expd);
  free(yuk_expn);
  free(yuk_exps);
  free(yuk_expe);
  free(yuk_expw);

  free(g_yuknumphys);
  free(g_yukvnexptot);
  free(g_yukvnexptotp);
  free(g_yukvnthmax); 
  free(g_yukdcu);
  free(g_yukdcd);
  free(g_yuksfactor);
  free(g_yuksfactor2);
  free(g_yukbetascale);
  free(g_yukrlsc);

  free(g_yuknumfour);
  free(g_yukwhts);
  free(g_yukrlams);
  free(g_yukytop);
  free(g_yukrdplus);
  free(g_yukrdminus);
  free(g_yukrdsq3);
  free(g_yukrdmsq3);
}

void sYukSourceToMultipole(const int ibox)
{
  int ptr = g_sboxes[ibox].addr;
  int level = g_sboxes[ibox].level;
  dcomplex *multipole = &yuk_multipole[g_pgsz*ibox];
  double *sources = &g_fmmsources[3*ptr]; 
  double *charges = &g_fmmcharges[ptr]; 
  int nsources = g_sboxes[ibox].npts; 
  double scale = g_yuksfactor[level]; 
  double x0y0z0[3];
  double h = g_size/pow(2, level + 1);
  int ix = g_sboxes[ibox].idx;
  int iy = g_sboxes[ibox].idy;
  int iz = g_sboxes[ibox].idz;
  x0y0z0[0] = g_bbcorner[0] + (2*ix + 1)*h;
  x0y0z0[1] = g_bbcorner[1] + (2*iy + 1)*h;
  x0y0z0[2] = g_bbcorner[2] + (2*iz + 1)*h;
  double *p = (double *)calloc(g_pgsz, sizeof(double)); 
  double *bi = (double *)calloc(g_pterms + 2, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 2, sizeof(dcomplex));
  if (p == 0 || bi == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }
  const double precision=1.0e-14;
  for (int i = 0; i < nsources; i++) {
    double rx = sources[3*i] - x0y0z0[0];
    double ry = sources[3*i + 1] - x0y0z0[1];
    double rz = sources[3*i + 2] - x0y0z0[2];
    double proj = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = ( d <= precision ? 1.0 : rz/d );
    ephi[0] = ( proj <= precision*d ? 1.0 : rx/proj - _Complex_I*ry/proj);
    for (int ell = 1; ell < g_pterms; ell++) { 
      ephi[ell] = ephi[ell-1]*ephi[0];
    }

    double rk = d*g_beta;
    int ncalc; 
    in(scale, rk, g_pterms, bi, &ncalc); 
    multipole[0] += charges[i]*bi[0];
    lgndr(g_pterms, ctheta, p);
    for (int ell = 1; ell <= g_pterms; ell++) {
      dcomplex cp = charges[i]*p[ell]*bi[ell]*g_yukytop[ell];
      multipole[ell] += cp;
    }

    for (int m = 1; m <= g_pterms; m++) {
      int offset = m*(g_pterms + 1);
      for (int ell = m; ell <= g_pterms; ell++) {
	double cp = charges[i]*bi[ell]*g_yukytop[ell + offset]*p[ell + offset];
	multipole[ell + offset] += cp*ephi[m - 1];
      }
    }
  }
  free(p);
  free(bi);
  free(ephi);
}

void dYukSourceToMultipole(const int ibox)
{
  int ptr = g_sboxes[ibox].addr;
  int level = g_sboxes[ibox].level;
  dcomplex *multipole = &yuk_multipole[g_pgsz*ibox];
  double *sources = &g_fmmsources[3*ptr];
  double *charges = &g_fmmcharges[ptr];
  double *dn = &g_fmminnernormal[3*ptr];
  int nsources = g_sboxes[ibox].npts;
  double scale = g_yuksfactor[level];
  double h = g_size/pow(2, level + 1);
  int ix = g_sboxes[ibox].idx;
  int iy = g_sboxes[ibox].idy;
  int iz = g_sboxes[ibox].idz;
  double x0y0z0[3];
  x0y0z0[0] = g_bbcorner[0] + (2*ix + 1)*h; 
  x0y0z0[1] = g_bbcorner[1] + (2*iy + 1)*h; 
  x0y0z0[2] = g_bbcorner[2] + (2*iz + 1)*h;

  double *p = (double *)calloc(g_pgsz, sizeof(double));
  double *bi = (double *)calloc(g_pterms + 2, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 2, sizeof(dcomplex));
  if (p == 0 || bi == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  const double precision = 1.0e-14;
  int i, ell, m, ncalc, index;
  double rk, d, cp, proj, rx, ry, rz, rr, ctheta;
  dcomplex cdn1, cdn2, cpz, cpz2;

  for (i = 0; i < nsources; i++) {
    rx = sources[3*i] - x0y0z0[0];
    ry = sources[3*i + 1] - x0y0z0[1];
    rz = sources[3*i + 2] - x0y0z0[2];
    proj = rx*rx + ry*ry;
    rr = proj + rz*rz;
    proj = sqrt(proj);
    d = sqrt(rr);
    ctheta = (d <= precision ? 1.0: rz/d);

    ephi[0] = (proj <= precision*d ? 1.0: rx/proj - _Complex_I*ry/proj);
    for (ell = 1; ell <= g_pterms; ell++) {
      ephi[ell] = ephi[ell - 1]*ephi[0];
    }

    rk = d*g_beta;
    in(scale, rk, g_pterms, bi, &ncalc);

    cdn1 = dn[3*i] - dn[3*i + 1]*_Complex_I;
    cdn2 = conj(cdn1);

    cp = charges[i]*bi[0]*g_beta/scale;
    multipole[g_pterms + 2] += cp/2.0*cdn1;
    multipole[1] -= cp*dn[3*i + 2];

    lgndr(g_pterms, ctheta, p);
    cp = charges[i]*p[1]*bi[1]*g_yukytop[1]*g_beta/3.0;
    multipole[g_pterms + 3] += cp/2.0/scale*cdn1;
    multipole[0] -= cp*scale*dn[3*i + 2];
    multipole[2] -= cp/scale*2.0*dn[3*i + 2];

    for (ell = 2; ell <= g_pterms-1; ell++) {
      cp = charges[i]*p[ell]*bi[ell]*g_yukytop[ell]*g_beta/(2*ell + 1);
      cpz = cp/2.0*cdn1;
      cp = cp*dn[3*i + 2];
      index = ell + g_pterms + 1; 
      multipole[index - 1] -= cpz*scale;
      multipole[index + 1] += cpz/scale;
      multipole[ell - 1] -= cp*scale*ell;
      multipole[ell + 1] -= cp/scale*(ell + 1);
    }

    ell = g_pterms; 
    cp = charges[i]*p[ell]*bi[ell]*g_yukytop[ell]*g_beta/(2*ell + 1)*scale;
    multipole[ell + g_pterms] -= cp/2.0*cdn1;
    multipole[ell - 1] -= cp*ell*dn[3*i + 2];

    cp = charges[i]*bi[1]*g_yukytop[g_pterms + 2]*p[g_pterms + 2]*g_beta/3.0;
    cpz = cp*ephi[0];
    cp = 2.0*creal(cpz*cdn2);
    multipole[0] += cp*scale;
    multipole[2] -= cp/scale;
    cpz = cpz/scale;
    multipole[2*g_pterms + 4] += cpz/2.0*cdn1;
    multipole[g_pterms + 3] -= cpz*dn[3*i + 2];

    for (ell = 2; ell <= g_pterms-1; ell++) {
      for (m = 1; m <= ell-2; m++) {
	index = ell + m*(g_pterms + 1); 
	cp = charges[i]*bi[ell]*g_yukytop[index]*p[index]*g_beta/(2*ell + 1);
	cpz = cp*ephi[m - 1];
	cpz2 = cpz*cdn2;
	if (m == 1) {
	  cp = ell*(ell + 1)*creal(cpz2);
	  multipole[ell - 1] += cp*scale;
	  multipole[ell + 1] -= cp/scale;
	} else {
	  cpz2 = cpz2/2.0;
	  index = ell + (m - 1)*(g_pterms + 1); 
	  multipole[index - 1] += cpz2*(ell + m - 1)*(ell + m)*scale;
	  multipole[index + 1] -= cpz2*(ell - m + 1)*(ell - m + 2)/scale;
	}

	cpz2 = cpz*cdn1/2.0;
	index = ell + (m + 1)*(g_pterms + 1); 
	multipole[index - 1] -= cpz2*scale;
	multipole[index + 1] += cpz2/scale;

	cpz2 = cpz*dn[3*i + 2];
	index = ell + m*(g_pterms + 1);
	multipole[index - 1] -= cpz2*(ell + m)*scale;
	multipole[index + 1] -= cpz2*(ell - m + 1)/scale;
      }

      m = ell-1;
      index = ell + m*(g_pterms + 1); 
      cp = charges[i]*bi[ell]*g_yukytop[index]*p[index]*g_beta/(2*ell + 1);
      cpz = cp*ephi[m - 1];
      cpz2 = cpz*cdn2;
      if (m == 1) {
	cp = creal(cpz2);
	multipole[ell - 1] += cp*ell*(ell + 1)*scale;
	multipole[ell + 1] -= cp*6.0/scale;
      } else {
	cpz2 = cpz2/2.0;
	index = ell + (m - 1)*(g_pterms + 1); 
	multipole[index - 1] += cpz2*(ell + m - 1)*(ell + m)*scale;
	multipole[index + 1] -= cpz2*6.0/scale;
      }

      cpz2 = cpz*cdn1/2.0;
      index = ell + 1 + (m + 1)*(g_pterms + 1);
      multipole[index] += cpz2/scale;

      cpz2 = cpz*dn[3*i + 2];
      index = ell + m*(g_pterms + 1); 
      multipole[index - 1] -= cpz2*(ell + m)*scale;
      multipole[index + 1] -= cpz2*2.0/scale;

      m = ell;
      index = ell + m*(g_pterms + 1);
      cp = charges[i]*bi[ell]*g_yukytop[index]*p[index]*g_beta/(2*ell + 1);
      cpz = cp*ephi[m - 1];

      cpz2 = cpz/2.0*cdn2;
      index = ell + (m - 1)*(g_pterms + 1);
      multipole[index - 1] += cpz2*(ell + m - 1)*(ell + m)*scale;
      multipole[index + 1] -= cpz2*2.0/scale;

      index = ell + 1 + (m + 1)*(g_pterms + 1);
      multipole[index] += cpz/(2.0*scale)*cdn1;

      index = ell + 1 + m*(g_pterms + 1);
      multipole[index] -= cpz*dn[3*i + 2]/scale;
    }

    ell = g_pterms;
    for (m = 1; m <= ell - 2; m++) {
      index = ell + m*(g_pterms + 1);
      cp = charges[i]*bi[ell]*g_yukytop[index]*p[index]*g_beta/(2*ell + 1)*scale;
      cpz = cp*ephi[m - 1];
      cpz2 = cpz*cdn2;
      if (m==1) {
	cp = creal(cpz2);
	multipole[ell - 1] += cp*ell*(ell + 1);
      } else {
	index = ell - 1 + (m - 1)*(g_pterms + 1);
	multipole[index] += cpz2*(ell + m - 1)*(ell + m)/2.0;
      }

      index = ell - 1 + (m + 1)*(g_pterms + 1); 
      multipole[index] -= cpz/2.0*cdn1; 

      index = ell - 1 + m*(g_pterms + 1);
      multipole[index] -= cpz*(ell + m)*dn[3*i + 2];
    }
 
    m = ell - 1;
    index = ell + m*(g_pterms + 1);
    cp = charges[i]*bi[ell]*g_yukytop[index]*p[index]*g_beta/(2*ell + 1)*scale;
    cpz = cp*ephi[m - 1];
    index = ell - 1 + (m - 1)*(g_pterms + 1);
    multipole[index] += cpz*(ell + m - 1)*(ell + m)/2.0*cdn2;
    index = ell - 1 + m*(g_pterms + 1); 
    multipole[index] -= cpz*(ell + m)*dn[3*i + 2];

    m = ell;
    index = ell + m*(g_pterms + 1);
    cp = charges[i]*bi[ell]*g_yukytop[index]*p[index]*g_beta/(2*ell + 1)*scale;
    cpz = cp*ephi[m - 1];
    index = ell - 1 + (m - 1)*(g_pterms + 1);
    multipole[index] += cpz*(ell + m - 1)*(ell + m)/2.0*cdn2;
  }
  free(p);
  free(bi);
  free(ephi);
}

void YukMultipoleToMultipole(const int pbox)
{
  dcomplex *pmultipole = &yuk_multipole[g_pgsz*pbox]; 
  int lev = g_sboxes[pbox].level; 
  double *dc = &g_yukdcu[g_yukdcpgsz*lev]; 
  static dcomplex var[5] = {0, -1 + _Complex_I, 1 + _Complex_I, 
			    1 - _Complex_I, -1 - _Complex_I};
  const double arg = sqrt(2.0)/2.0;
  dcomplex *mpolen = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *marray = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 2, sizeof(dcomplex));
  if (mpolen == 0 || marray == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int i = 0; i < 8; i++) {
    int cbox = g_sboxes[pbox].child[i]; 
    if (cbox > 0) {
      int ifl = g_yifl[i]; ;
      double *rd = (i < 4 ? g_yukrdsq3 : g_yukrdmsq3); 
      dcomplex *cmultipole = &yuk_multipole[g_pgsz*cbox]; 
      ephi[0] = 1.0;
      ephi[1] = arg*var[ifl];
      for (int ell = 1; ell <= g_pterms; ell++) 
	ephi[ell + 1] = ephi[ell]*ephi[1];

      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset;
	  mpolen[index] = conj(ephi[m])*cmultipole[index];
	}
      }

      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1); 
	int offset1 = (m + g_pterms)*g_pgsz; 
	int offset2 = (-m + g_pterms)*g_pgsz; 
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = offset + ell; 
	  marray[index] = mpolen[ell]*rd[ell + offset1]; 
	  for (int mp = 1; mp <= ell; mp++) {
	    int index1 = mp*(g_pterms + 1) + ell; 
	    marray[index] += mpolen[index1]*rd[index1 + offset1] + 
	      conj(mpolen[index1])*rd[index1 + offset2]; 
	  }
	}
      }

      for (int k = 0; k<= g_pterms; k++) {
	int offset = k*(g_pterms + 1); 
	for (int j = k; j <= g_pterms; j++) {
	  int index = offset + j; 
	  int offset1 = j*(g_pterms + 1) + k; 
	  mpolen[index] = 0; 
	  for (int ell = k; ell <= g_pterms; ell++) 
	    mpolen[index] += marray[ell + offset]*dc[offset1 + ell*g_pgsz]; 
	}
      }

      for (int m = 0; m <= g_pterms; m = m + 2) {
	int offset = m*(g_pterms + 1);
	int offset1 = (m + g_pterms)*g_pgsz; 
	int offset2 = (-m + g_pterms)*g_pgsz; 
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset; 
	  marray[index] = mpolen[ell]*rd[ell + offset1]; 
	  for (int mp = 1; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1); 
	    marray[index] -= mpolen[index1]*rd[offset1 + index1] + 
	      conj(mpolen[index1])*rd[offset2 + index1]; 
	  }
	  
	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1); 
	    marray[index] += mpolen[index1]*rd[offset1 + index1] + 
	      conj(mpolen[index1])*rd[offset2 + index1];
	  }
	}
      }

      for (int m = 1; m <= g_pterms; m = m + 2) {
	int offset = m*(g_pterms + 1);
	int offset1 = (m + g_pterms)*g_pgsz; 
	int offset2 = (-m + g_pterms)*g_pgsz; 
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset; 
	  marray[index] = -mpolen[ell]*rd[ell + offset1]; 
	  for (int mp = 1; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1); 
	    marray[index] += mpolen[index1]*rd[index1 + offset1] + 
	      conj(mpolen[index1])*rd[index1 + offset2]; 
	  }

	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1);
	    marray[index] -= mpolen[index1]*rd[index1 + offset1] + 
	      conj(mpolen[index1])*rd[index1 + offset2]; 
	  }
	}
      }

      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset; 
	  mpolen[index] = ephi[m]*marray[index]; 
	}
      }

      for (int m = 0; m < g_pgsz; m++) 
	pmultipole[m] += mpolen[m]; 
    }
  }
  free(mpolen);
  free(marray);
  free(ephi);
}

void YukMultipoleToExponential(const int ibox)
{
  dcomplex *mw = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex)); 
  dcomplex *mexpf1 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex)); 
  dcomplex *mexpf2 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex)); 
  if (mw == 0 || mexpf1 == 0 || mexpf2 == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  int level = g_sboxes[ibox].level; 
  YukMultipoleToExponentialPhase1(&yuk_multipole[g_pgsz*ibox], level, mexpf1, mexpf2); 
  YukMultipoleToExponentialPhase2(mexpf1, level, &yuk_expu[g_yukmnexptotp*ibox]); 
  YukMultipoleToExponentialPhase2(mexpf2, level, &yuk_expd[g_yukmnexptotp*ibox]);

  rotz2y(&yuk_multipole[g_pgsz*ibox], g_yukrdminus, mw);
  YukMultipoleToExponentialPhase1(mw, level, mexpf1, mexpf2);
  YukMultipoleToExponentialPhase2(mexpf1, level, &yuk_expn[g_yukmnexptotp*ibox]);
  YukMultipoleToExponentialPhase2(mexpf2, level, &yuk_exps[g_yukmnexptotp*ibox]); 

  rotz2x(&yuk_multipole[g_pgsz*ibox], g_yukrdplus, mw);
  YukMultipoleToExponentialPhase1(mw, level, mexpf1, mexpf2);
  YukMultipoleToExponentialPhase2(mexpf1, level, &yuk_expe[g_yukmnexptotp*ibox]);
  YukMultipoleToExponentialPhase2(mexpf2, level, &yuk_expw[g_yukmnexptotp*ibox]); 

  free(mw);
  free(mexpf1);
  free(mexpf2);
}

void YukMultipoleToExponentialPhase1(const dcomplex * const multipole, const int level, 
				     dcomplex * const mexpu, dcomplex * const mexpd)
{
  double *rlsc = &g_yukrlsc[g_pgsz*g_nlambs*(level - 1)];
  int ntot = 0;
  for (int nell = 0; nell < g_nlambs; nell++) {    
    double sgn = -1;
    for (int mth = 0; mth <= g_yuknumfour[nell] - 1; mth++) {
      int ncurrent = ntot + mth;
      dcomplex ztmp1 = 0;
      dcomplex ztmp2 = 0;
      sgn = -sgn;
      int offset = mth*(g_pterms + 1);
      int offset1 = offset + nell*g_pgsz; 
      for (int nm = mth; nm <= g_pterms; nm = nm + 2) 
	ztmp1 += multipole[nm + offset]*rlsc[nm + offset1];
      for (int nm = mth + 1; nm <= g_pterms; nm = nm + 2)
	ztmp2 += multipole[nm + offset]*rlsc[nm + offset1];
      mexpu[ncurrent] = ztmp1 + ztmp2;
      mexpd[ncurrent] = sgn*(ztmp1 - ztmp2);
    }
    ntot += g_yuknumfour[nell];
  }
}

void YukMultipoleToExponentialPhase2(const dcomplex * const mexpf, const int level, 
				     dcomplex * const mexpphys)
{
  dcomplex *fexpe = &yuk_fexpe[15000*(level - 1)]; 
  dcomplex *fexpo = &yuk_fexpo[15000*(level - 1)];  
  int nftot = 0;
  int nptot = 0;
  int nexte = 0;
  int nexto = 0;
  for (int i = 0; i < g_nlambs; i++) {
    for (int ival = 0; ival < g_yuknumphys[i]/2; ival++) {
      mexpphys[nptot + ival] = mexpf[nftot];
      double sgn = -2;        
      for (int nm = 1; nm < g_yuknumfour[i]; nm = nm + 2) {
	sgn = -sgn;
	double rtmp = sgn*(fexpe[nexte]*mexpf[nftot + nm]);
	nexte += 1;
	mexpphys[nptot + ival] += rtmp*_Complex_I;
      }
      
      sgn = 2;
      for (int nm = 2; nm < g_yuknumfour[i]; nm = nm + 2) {
	sgn = -sgn;
	double rtmp = sgn*(fexpo[nexto]*mexpf[nftot + nm]);
	nexto += 1;
	mexpphys[nptot + ival] += rtmp;
      }      
    }
    
    nftot += g_yuknumfour[i];
    nptot += g_yuknumphys[i]/2;
  }
}

void YukExponentialToLocal(const int ibox)
{
  if (g_tboxes[ibox].nchild) {
    int uall[36], nuall, xuall[36], yuall[36], 
      u1234[16], nu1234,x1234[16],y1234[16],
      dall[36], ndall, xdall[36], ydall[36], 
      d5678[16], nd5678, x5678[16], y5678[16], 
      nall[24], nnall, xnall[24], ynall[24], 
      n1256[8], nn1256, x1256[8], y1256[8], 
      n12[4], nn12, x12[4], y12[4], n56[4], nn56, x56[4], y56[4], 
      sall[24], nsall, xsall[24], ysall[24], 
      s3478[8], ns3478, x3478[8], y3478[8], 
      s34[4], ns34, x34[4], y34[4], s78[4], ns78, x78[4], y78[4], 
      eall[16], neall, xeall[16], yeall[16],
      e1357[4], ne1357, x1357[4], y1357[4], 
      e13[2], ne13, x13[2], y13[2], e57[2], ne57, x57[2], y57[2], 
      e1[3], ne1, x1[3], y1[3], e3[3], ne3, x3[3], y3[3], 
      e5[3], ne5, x5[3], y5[3], e7[3], ne7, x7[3], y7[3], 
      wall[16], nwall, xwall[16], ywall[16], 
      w2468[4], nw2468, x2468[4], y2468[4], 
      w24[2], nw24, x24[2], y24[2], w68[2], nw68, x68[2], y68[2], 
      w2[3], nw2, x2[3], y2[3], w4[3], nw4, x4[3], y4[3], 
      w6[3], nw6, x6[3], y6[3], w8[3], nw8, x8[3], y8[3];

    BuildMergedList2(ibox, uall, &nuall, xuall, yuall, 
		     u1234, &nu1234, x1234, y1234, 
		     dall, &ndall, xdall, ydall, 
		     d5678, &nd5678, x5678, y5678,
		     nall, &nnall, xnall, ynall, 
		     n1256, &nn1256, x1256, y1256, 
		     n12, &nn12, x12, y12, n56, &nn56, x56, y56, 
		     sall, &nsall, xsall, ysall,
		     s3478, &ns3478, x3478, y3478, 
		     s34, &ns34, x34, y34, s78, &ns78, x78, y78, 
		     eall, &neall, xeall, yeall, 
		     e1357, &ne1357, x1357, y1357, 
		     e13, &ne13, x13, y13, e57, &ne57, x57, y57, 
		     e1, &ne1, x1, y1, e3, &ne3, x3, y3, 
		     e5, &ne5, x5, y5, e7, &ne7, x7, y7, 
		     wall, &nwall, xwall, ywall, 
		     w2468, &nw2468, x2468, y2468, 
		     w24, &nw24, x24, y24, w68, &nw68, x68, y68, 
		     w2, &nw2, x2, y2, w4, &nw4, x4, y4, 
		     w6, &nw6, x6, y6, w8, &nw8, x8, y8);

    int *child = g_tboxes[ibox].child; 
    int level = g_tboxes[ibox].level;
    dcomplex *xs = &yuk_xs[3*g_yuknexpmax*level]; 
    dcomplex *ys = &yuk_ys[3*g_yuknexpmax*level]; 
    double *zs = &yuk_zs[3*g_yuknexpmax*level]; 
    int nexptotp = g_yukvnexptotp[level]; 

    dcomplex *mw1 = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
    dcomplex *mw2 = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
    dcomplex *temp = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexpf1 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexpf2 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));

    // Step 1: Process z-direction exponential expansions, and write
    // results to uall, u1234, dall, and d5678 lists. 
    dcomplex *mexuall = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexu1234 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexdall = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexd5678 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    if (mw1 == 0 || mw2 == 0 || temp == 0 || mexpf1 == 0 || mexpf2 == 0 ||
	mexuall == 0 || mexu1234 == 0 || mexdall == 0 || mexd5678 ==  0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    MakeUList(g_yukmnexptotp, yuk_expd, uall, nuall, xuall, yuall, xs, ys, mexuall); 
    MakeUList(g_yukmnexptotp, yuk_expd, u1234, nu1234, x1234, y1234, xs, ys, mexu1234); 
    MakeDList(g_yukmnexptotp, yuk_expu, dall, ndall, xdall, ydall, xs, ys, mexdall);
    MakeDList(g_yukmnexptotp, yuk_expu, d5678, nd5678, x5678, y5678, xs, ys, mexd5678);

    if (child[0]) {
      dcomplex *local = &yuk_local[child[0]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 2];
	iexp1++;
      }

      if (nu1234) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexu1234[j]*zs[3*j + 1];
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1); 
    
      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 1];
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }
      
      if (iexp1 + iexp2) {
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[1]) {
      dcomplex *local = &yuk_local[child[1]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 2]*conj(xs[3*j]);
	iexp1++;
      }
      
      if (nu1234) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexu1234[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 1]*xs[3*j];
	YukExponentialToLocalPhase1(temp, level, mexpf2);
	iexp2++;
      }

      if (iexp1 + iexp2) {
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[2]) {
      dcomplex *local = &yuk_local[child[2]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 2]*conj(ys[3*j]);
	iexp1++;
      }

      if (nu1234) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexu1234[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
         
      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;      
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[3]) {
      dcomplex *local = &yuk_local[child[3]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 2]*conj(xs[3*j]*ys[3*j]);
	iexp1++;
      }

      if (nu1234) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexu1234[j]*zs[3*j + 1]*conj(xs[3*j]*ys[3*j]);
	iexp1++;
      }
      
      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);

      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 1]*xs[3*j]*ys[3*j];
	iexp2++;      
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[4]) {
      dcomplex *local = &yuk_local[child[4]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 1];
	iexp1++;      
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }

      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 2];
	iexp2++;
      }

      if (nd5678) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexd5678[j]*zs[3*j + 1];
	iexp2++;
      }

      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
    
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[5]) {
      dcomplex *local = &yuk_local[child[5]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp1++;      
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }

      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 2]*xs[3*j];
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexd5678[j]*zs[3*j + 1]*xs[3*j];
	iexp2++;
      }

      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[6]) {
      dcomplex *local = &yuk_local[child[6]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }

      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (ndall)  {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 2]*ys[3*j];
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexd5678[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
      }
      
      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[7]) {
      dcomplex *local = &yuk_local[child[7]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (nuall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexuall[j]*zs[3*j + 1]*conj(xs[3*j]*ys[3*j]);
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }

      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (ndall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexdall[j]*zs[3*j + 2]*xs[3*j]*ys[3*j];
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexd5678[j]*zs[3*j + 1]*xs[3*j]*ys[3*j];
	iexp2++;
      }

      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }   

    free(mexuall);
    free(mexu1234);
    free(mexdall);
    free(mexd5678);

    // Step 2: Process y-direction exponential expansions, and write
    // results to nall, n1256, n12, n56, sall, s3478, s34, and s78
    // lists. 
    dcomplex *mexnall = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexn1256 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexn12 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexn56 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexsall = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexs3478 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexs34 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexs78 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    if (mexnall == 0 || mexn1256 == 0 || mexn12 == 0 || mexn56 == 0 ||
	mexsall == 0 || mexs3478 == 0 || mexs34 == 0 || mexs78 == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    MakeUList(g_yukmnexptotp, yuk_exps, nall, nnall, xnall, ynall, xs, ys, mexnall);
    MakeUList(g_yukmnexptotp, yuk_exps, n1256, nn1256, x1256, y1256, xs, ys, mexn1256);
    MakeUList(g_yukmnexptotp, yuk_exps, n12, nn12, x12, y12, xs, ys, mexn12);
    MakeUList(g_yukmnexptotp, yuk_exps, n56, nn56, x56, y56, xs, ys, mexn56);
    MakeDList(g_yukmnexptotp, yuk_expn, sall, nsall, xsall, ysall, xs, ys, mexsall);
    MakeDList(g_yukmnexptotp, yuk_expn, s3478, ns3478, x3478, y3478, xs, ys, mexs3478);
    MakeDList(g_yukmnexptotp, yuk_expn, s34, ns34, x34, y34, xs, ys, mexs34);
    MakeDList(g_yukmnexptotp, yuk_expn, s78, ns78, x78, y78, xs, ys, mexs78);

    if (child[0]) {
      dcomplex *local = &yuk_local[child[0]*g_pgsz]; 
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 2];
	iexp1++;
      }

      if (nn1256) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn1256[j]*zs[3*j + 1];
	iexp1++;
      }

      if (nn12) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn12[j]*zs[3*j + 1];
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);

      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 1];
	iexp2++; 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2); 
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[1]) {
      dcomplex *local = &yuk_local[child[1]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 2]*conj(ys[3*j]);
	iexp1++;
      }

      if (nn1256) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn1256[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
      }

      if (nn12) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn12[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[2]) {
      dcomplex *local = &yuk_local[child[2]*g_pgsz];
      int iexp1 = 0; 
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 1];
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }

      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 2];
	iexp2++;
      }

      if (ns3478) {
	for (int j = 0; j < nexptotp; j++)
	  temp[j] += mexs3478[j]*zs[3*j + 1];
	iexp2++;
      }

      if (ns34) {
	for (int j = 0; j < nexptotp; j++)
	  temp[j] += mexs34[j]*zs[3*j+1];
	iexp2++;
      }

      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[3]) {
      dcomplex *local = &yuk_local[child[3]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }
      
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 2]*ys[3*j];
	iexp2++;
      }

      if (ns3478) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexs3478[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
      }

      if (ns34) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexs34[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
      }

      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }

    if (child[4]) {
      dcomplex *local = &yuk_local[child[4]*g_pgsz]; 
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      int iexp1 = 0; 
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 2]*conj(xs[3*j]);
	iexp1++;
      }

      if (nn1256) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn1256[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp1++;
      }

      if (nn56) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn56[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 1]*xs[3*j];
	YukExponentialToLocalPhase1(temp, level, mexpf2);
	iexp2++;
      }
      
      if (iexp1 + iexp2) {
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[5]) {
      dcomplex *local = &yuk_local[child[5]*g_pgsz];
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      int iexp1 = 0;
      int iexp2 = 0;
    
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 2]*conj(xs[3*j]*ys[3*j]);
	iexp1++;
      }
      
      if (nn1256) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexn1256[j]*zs[3*j + 1]*conj(xs[3*j]*ys[3*j]);
	iexp1++;
      }

      if (nn56) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = temp[j] + mexn56[j]*zs[3*j+1]*conj(xs[3*j]*ys[3*j]);
	iexp1 = iexp1+1;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 1]*xs[3*j]*ys[3*j];
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }

      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[6]) {
      dcomplex *local = &yuk_local[child[6]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }

      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 2]*xs[3*j];
	iexp2++;
      }
      
      if (ns3478) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexs3478[j]*zs[3*j + 1]*xs[3*j];
	iexp2++;
      }

      if (ns78) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexs78[j]*zs[3*j + 1]*xs[3*j];
	iexp2++;
      }

      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    } 
    
    if (child[7]) {
      dcomplex *local = &yuk_local[child[7]*g_pgsz];
      int iexp1 = 0; 
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexnall[j]*zs[3*j + 1]*conj(xs[3*j]*ys[3*j]);
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);	    
      }
      
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nsall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexsall[j]*zs[3*j + 2]*xs[3*j]*ys[3*j];
	iexp2++;
      }
      
      if (ns3478) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexs3478[j]*zs[3*j + 1]*xs[3*j]*ys[3*j];
	iexp2++;
      }

      if (ns78) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexs78[j]*zs[3*j + 1]*xs[3*j]*ys[3*j];
	iexp2++;
      }

      if (iexp2 > 0) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	roty2z(mw1, g_yukrdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }
    
    free(mexnall);
    free(mexn1256);
    free(mexn12);
    free(mexn56);
    free(mexsall);
    free(mexs3478);
    free(mexs34);
    free(mexs78);

    // Step 3: Process x-direction exponential expansions, and write
    // results to eall, e1357, e13, e57, e1, e3, e5, e7, wall, w2468,
    // w24, w68, w2, w4, w6, and w8 lists. 
    dcomplex *mexeall = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe1357 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe13 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe57 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe1 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe3 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe5 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexe7 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexwall = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw2468 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw24 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw68 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw2 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw4= (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw6 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    dcomplex *mexw8 = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
    if (mexeall == 0 || mexe1357 == 0 || mexe13 == 0 || mexe57 == 0 || 
	mexe1 == 0 || mexe3 == 0 || mexe5 == 0 || mexe7 == 0 || mexwall == 0 || 
	mexw2468 == 0 || mexw24 == 0 || mexw68 == 0 || mexw2 == 0 || mexw4 == 0 
	|| mexw6 == 0 || mexw8 == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    MakeUList(g_yukmnexptotp, yuk_expw, eall, neall, xeall, yeall, xs, ys, mexeall);
    MakeUList(g_yukmnexptotp, yuk_expw, e1357, ne1357, x1357, y1357, xs, ys, mexe1357);
    MakeUList(g_yukmnexptotp, yuk_expw, e13, ne13, x13, y13, xs, ys, mexe13);
    MakeUList(g_yukmnexptotp, yuk_expw, e57, ne57, x57, y57, xs, ys, mexe57);
    MakeUList(g_yukmnexptotp, yuk_expw, e1, ne1, x1, y1, xs, ys, mexe1); 
    MakeUList(g_yukmnexptotp, yuk_expw, e3, ne3, x3, y3, xs, ys, mexe3);
    MakeUList(g_yukmnexptotp, yuk_expw, e5, ne5, x5, y5, xs, ys, mexe5);
    MakeUList(g_yukmnexptotp, yuk_expw, e7, ne7, x7, y7, xs, ys, mexe7);
    MakeDList(g_yukmnexptotp, yuk_expe, wall, nwall, xwall, ywall, xs, ys, mexwall);
    MakeDList(g_yukmnexptotp, yuk_expe, w2468, nw2468, x2468, y2468, xs, ys, mexw2468);
    MakeDList(g_yukmnexptotp, yuk_expe, w24, nw24, x24, y24, xs, ys, mexw24);
    MakeDList(g_yukmnexptotp, yuk_expe, w68, nw68, x68, y68, xs, ys, mexw68);
    MakeDList(g_yukmnexptotp, yuk_expe, w2, nw2, x2, y2, xs, ys, mexw2);
    MakeDList(g_yukmnexptotp, yuk_expe, w4, nw4, x4, y4, xs, ys, mexw4);
    MakeDList(g_yukmnexptotp, yuk_expe, w6, nw6, x6, y6, xs, ys, mexw6);
    MakeDList(g_yukmnexptotp, yuk_expe, w8, nw8, x8, y8, xs, ys, mexw8); 

    if (child[0]) {
      dcomplex *local = &yuk_local[child[0]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 2];
	  iexp1++;
      }

      if (ne1357) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe1357[j]*zs[3*j + 1];
	iexp1++;
      }

      if (ne13) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe13[j]*zs[3*j + 1];
	iexp1++;
      }
      
      if (ne1) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe1[j]*zs[3*j + 1];
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 1];
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }
      
      if (iexp1 + iexp2) {          
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
    
    if (child[1]) { 
      dcomplex *local = &yuk_local[child[1]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 1];
	  iexp1++;
	  YukExponentialToLocalPhase1(temp, level, mexpf1);
      }
      
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 2];
	iexp2++;
      }
      
      if (nw2468) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw2468[j]*zs[3*j + 1];
	iexp2++;
      }
      
      if (nw24) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw24[j]*zs[3*j + 1];
	iexp2++;
      }
      
      if (nw2) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw2[j]*zs[3*j + 1];
	iexp2++;
      }
      
      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }
    
    if (child[2]) {
      dcomplex *local = &yuk_local[child[2]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;

      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 2]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (ne1357) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe1357[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (ne13) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe13[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (ne3) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe3[j]*zs[3*j + 1]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
    
    if (child[3]) {
      dcomplex *local = &yuk_local[child[3]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 1]*conj(ys[3*j]);
	  iexp1++;
	  YukExponentialToLocalPhase1(temp, level, mexpf1);
      }
      
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 2]*ys[3*j];
	  iexp2++;
      }
      
      if (nw2468) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw2468[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
      }
      
      if (nw24) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw24[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
      }
      
      if (nw4) {
	for (int j = 0; j < nexptotp; j++) 
	    temp[j] += mexw4[j]*zs[3*j + 1]*ys[3*j];
	iexp2++;
      }
      
      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }
    
    if (child[4]) {
      dcomplex *local = &yuk_local[child[4]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
            
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	    temp[j] = mexeall[j]*zs[3*j + 2]*xs[3*j];
	iexp1++;
      }
      
      if (ne1357) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe1357[j]*zs[3*j + 1]*xs[3*j];
	iexp1++;
      }
      
      if (ne57) {
	for (int j = 0; j < nexptotp; j++) 
	    temp[j] += mexe57[j]*zs[3*j + 1]*xs[3*j];
	iexp1++;
      }

      if (ne5) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe5[j]*zs[3*j + 1]*xs[3*j];
	iexp1++;
      }
      
      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2); 
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }
    
    if (child[5]) {
      dcomplex *local = &yuk_local[child[5]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 1]*xs[3*j];
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }
      
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 2]*conj(xs[3*j]);
	iexp2++;
      }
      
      if (nw2468) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw2468[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp2++;
      }
      
      if (nw68) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw68[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp2++;
      }
      
      if (nw6) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] +=  mexw6[j]*zs[3*j + 1]*conj(xs[3*j]);
	iexp2++;
      }
      
      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
      
    if (child[6]) {
      dcomplex *local = &yuk_local[child[6]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 2]*xs[3*j]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (ne1357) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe1357[j]*zs[3*j + 1]*xs[3*j]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (ne57) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe57[j]*zs[3*j + 1]*xs[3*j]*conj(ys[3*j]);
	iexp1++;
      }
      
      if (ne7) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexe7[j]*zs[3*j + 1]*xs[3*j]*conj(ys[3*j]);
	iexp1++;
      }

      if (iexp1) 
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	    temp[j] = mexwall[j]*zs[3*j + 1]*conj(xs[3*j])*ys[3*j];
	iexp2++;
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      }
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
    
    if (child[7]) {
      dcomplex *local = &yuk_local[child[7]*g_pgsz]; 
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (neall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexeall[j]*zs[3*j + 1]*xs[3*j]*conj(ys[3*j]);
	iexp1++;
	YukExponentialToLocalPhase1(temp, level, mexpf1);
      }
      
      for (int j = 0; j < nexptotp; j++) 
	temp[j] = 0;
      
      if (nwall) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] = mexwall[j]*zs[3*j + 2]*conj(xs[3*j])*ys[3*j];
	iexp2++;
      }
      
      if (nw2468) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw2468[j]*zs[3*j + 1]*conj(xs[3*j])*ys[3*j];
	iexp2++;
      }
      
      if (nw68) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw68[j]*zs[3*j + 1]*conj(xs[3*j])*ys[3*j];
	iexp2++;
      }
      
      if (nw8) {
	for (int j = 0; j < nexptotp; j++) 
	  temp[j] += mexw8[j]*zs[3*j + 1]*conj(xs[3*j])*ys[3*j];
	iexp2++;
      }
      
      if (iexp2) 
	YukExponentialToLocalPhase1(temp, level, mexpf2);
      
      if (iexp1 + iexp2) {      
	YukExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, level, mw1);
	rotz2x(mw1, g_yukrdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
    
    free(mexeall);
    free(mexe1357);
    free(mexe13);
    free(mexe57);
    free(mexe1);
    free(mexe3);
    free(mexe5);
    free(mexe7);
    free(mexwall);
    free(mexw2468);
    free(mexw24);
    free(mexw68);
    free(mexw2);
    free(mexw4);
    free(mexw6);
    free(mexw8);
    
    free(mexpf1);
    free(mexpf2);
    free(temp);
    free(mw1);
    free(mw2);
  }
}

void YukExponentialToLocalPhase1(const dcomplex * const mexpphys, const int level, 
				 dcomplex * const mexpf)
{
  int *numphys = &g_yuknumphys[g_nlambs*level]; 
  dcomplex *fexpback = &yuk_fexpback[15000*level]; 
  int nftot = 0;
  int nptot = 0;
  int next = 0;
  for (int i = 0; i < g_nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha/2;
    mexpf[nftot] = 0;
    for (int ival = 0; ival < nalpha2; ival++) 
      mexpf[nftot] += 2.0*creal(mexpphys[nptot + ival]);
    mexpf[nftot] /= nalpha;

    for (int mm = 2; mm < g_yuknumfour[i]; mm = mm + 2) {
      mexpf[nftot + mm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
	double rtmp = 2*creal(mexpphys[nptot + ival]);
	mexpf[nftot + mm] += fexpback[next]*rtmp;
	next++;
      }
      mexpf[nftot + mm] /= nalpha;
    }

    for (int mm = 1; mm < g_yuknumfour[i]; mm = mm + 2) {
      mexpf[nftot + mm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
	dcomplex ztmp = 2*cimag(mexpphys[nptot + ival])*_Complex_I;
	mexpf[nftot + mm] += fexpback[next]*ztmp;
	next++;
      }
      mexpf[nftot + mm] /= nalpha;
    }
    nftot += g_yuknumfour[i];
    nptot += numphys[i]/2;
  }
}

void YukExponentialToLocalPhase2(const int iexpu, const dcomplex * const mexpu, 
				 const int iexpd, const dcomplex * const mexpd, 
				 const int level, dcomplex * const local)
{
  dcomplex *mexpplus = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
  dcomplex *mexpminus = (dcomplex *)calloc(g_yuknexpmax, sizeof(dcomplex));
  dcomplex *zeye = (dcomplex *)calloc(g_pterms + 1, sizeof(dcomplex)); 
  if (mexpplus == 0 || mexpminus == 0 || zeye == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }
  double beta = g_yukbetascale[level]; 
  double *rlsc = &g_yukrlsc[g_pgsz*g_nlambs*level]; 
  int nexptot = g_yukvnexptot[level];

  zeye[0] = 1.0; 
  for (int i = 1; i <= g_pterms; i++) 
    zeye[i] = zeye[i - 1]*_Complex_I; 

  for (int i = 0; i < g_pgsz; i++) 
    local[i] = 0.0; 

  if (iexpu <= 0) {
    for (int nm = 0; nm < nexptot; nm++) {
      mexpplus[nm] = mexpd[nm];
      mexpminus[nm] = mexpd[nm];
    }
  } else if (iexpd <= 0) {
    for (int nm = 0; nm < nexptot; nm++) {
      mexpplus[nm] = mexpu[nm];
      mexpminus[nm] = -mexpu[nm];
    }
  } else {
    for (int nm = 0; nm < nexptot; nm++) {
      mexpplus[nm] = mexpd[nm] + mexpu[nm];
      mexpminus[nm] = mexpd[nm] - mexpu[nm];
    }
  }

  int ntot = 1;

  for (int nl = 0; nl < g_nlambs; nl++) {    
    int mmax = g_yuknumfour[nl]-1;
    int offset1 = nl*g_pgsz;
    for (int mth = 0; mth <= mmax; mth = mth + 2) {
      int ncurrent = ntot + mth - 1;
      int offset = mth*(g_pterms + 1);
      dcomplex temp = g_yukwhts[nl]*mexpplus[ncurrent];
      for (int nm = mth; nm <= g_pterms; nm = nm + 2) {
	local[nm + offset]  += rlsc[nm + offset + offset1]*temp;
      }
      temp = g_yukwhts[nl]*mexpminus[ncurrent];
      for (int nm = mth + 1; nm <= g_pterms; nm = nm + 2) {
	local[nm + offset] += rlsc[nm + offset + offset1]*temp;
      }
    }

    for (int mth = 1; mth <= mmax; mth = mth + 2) {
      int ncurrent = ntot + mth - 1;
      int offset = mth*(g_pterms + 1);
      dcomplex temp = g_yukwhts[nl]*mexpminus[ncurrent];
      for (int nm = mth; nm <= g_pterms; nm = nm + 2) {
	local[nm + offset] += rlsc[nm + offset + offset1]*temp;
      }

      temp = g_yukwhts[nl]*mexpplus[ncurrent];
      for (int nm = mth + 1; nm <= g_pterms; nm = nm + 2) {
	local[nm + offset] += rlsc[nm + offset + offset1]*temp;
      }
    }
    ntot += g_yuknumfour[nl];
  }

  double rscale = M_PI/beta/2.0;

  for (int mth = 0; mth <= g_pterms; mth++) {
    int offset = mth*(g_pterms + 1);
    dcomplex temp = zeye[mth]*rscale;
    for (int nm = mth; nm <= g_pterms; nm++) {
      local[nm + offset] *= temp*g_yukytop[nm + offset];
    }
  }

  free(mexpplus);
  free(mexpminus);
  free(zeye);
}

void YukLocalToLocal(const int pbox)
{
  dcomplex *local = &yuk_local[g_pgsz*pbox]; 
  int level = g_tboxes[pbox].level; 
  double *dc = &g_yukdcd[g_yukdcpgsz*level]; 
  static dcomplex var[5]={0, 1 - _Complex_I, -1 - _Complex_I, 
			  -1 + _Complex_I, 1 + _Complex_I};
  const double arg = sqrt(2.0)/2.0;
  dcomplex *localn = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *marray = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 2, sizeof(dcomplex));
  if (localn == 0 || marray == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int i = 0; i < 8; i++) {   
    int cbox = g_tboxes[pbox].child[i];
    if (cbox) {
      int ifl = g_yifl[i]; 
      double *rd = (i < 4 ? g_yukrdmsq3 : g_yukrdsq3); 
      ephi[0] = 1.0;
      ephi[1] = arg*var[ifl];
      for (int ell = 1; ell <= g_pterms; ell++) 
	ephi[ell + 1] = ephi[ell]*ephi[1];

      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell=m; ell <= g_pterms; ell++) {
	  localn[ell + offset] = conj(ephi[m])*local[ell + offset];
	}
      }

      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1); 
	int offset1 = (m + g_pterms)*g_pgsz; 
	int offset2 = (-m + g_pterms)*g_pgsz; 
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset; 
	  marray[index] = localn[ell]*rd[ell + offset1]; 
	  for (int mp = 1; mp <= ell; mp++) {
	    int index1 = ell + mp*(g_pterms + 1); 
	    marray[index] += localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }
	}
      }
      
      for (int k = 0; k <= g_pterms; k++) {
	int offset = k*(g_pterms + 1); 
	for (int j = k; j <= g_pterms; j++) {
	  int index = j + offset; 
	  int offset1 = k + j*(g_pterms + 1);
	  localn[index] = 0;
	  for (int ell = k; ell <= g_pterms; ell++) 
	    localn[index] += marray[ell + offset]*dc[offset1 + ell*g_pgsz]; 
	}
      }

      for (int m = 0; m <= g_pterms; m = m + 2) {
	int offset = m*(g_pterms + 1); 
	int offset1 = (m + g_pterms)*g_pgsz; 
	int offset2 = (-m + g_pterms)*g_pgsz; 
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset; 
	  marray[index] = localn[ell]*rd[ell + offset1];
	  for (int mp = 1; mp <= ell; mp = mp + 2) {
	    int index1 = mp*(g_pterms + 1) + ell; 
	    marray[index] -= localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }

	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = mp*(g_pterms + 1) + ell; 
	    marray[index] += localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }
	}
      }

      for (int m = 1; m <= g_pterms; m = m + 2) {
	int offset = m*(g_pterms + 1);
	int offset1 = (m + g_pterms)*g_pgsz; 
	int offset2 = (-m + g_pterms)*g_pgsz; 
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = offset + ell; 
	  marray[index] = -localn[ell]*rd[ell + offset1]; 
	  for (int mp = 1; mp <= ell; mp = mp + 2) {
	    int index1 = mp*(g_pterms + 1) + ell; 
	    marray[index] += localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }

	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = mp*(g_pterms + 1) + ell; 
	    marray[index] -= localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2]; 
	  }
	}
      }

      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell = m; ell <= g_pterms; ell++) {
	  localn[ell + offset] = ephi[m]*marray[ell + offset];
	}
      }

      dcomplex *clocal = &yuk_local[g_pgsz*cbox];
      for (int m =0; m < g_pgsz; m++) 
	clocal[m] += localn[m];
    }
  }
  free(localn);
  free(marray);
  free(ephi);
}

void YukLocalToTarget(const int ibox)
{
  dcomplex *local = &yuk_local[g_pgsz*ibox]; 
  int level = g_tboxes[ibox].level; 
  double scale = g_yuksfactor[level]; 
  double h = g_size/pow(2, level + 1); 
  int ix = g_tboxes[ibox].idx; 
  int iy = g_tboxes[ibox].idy; 
  int iz = g_tboxes[ibox].idz; 
  double x0y0z0[3]; 
  x0y0z0[0] = g_bbcorner[0] + (2*ix + 1)*h; 
  x0y0z0[1] = g_bbcorner[1] + (2*iy + 1)*h;
  x0y0z0[2] = g_bbcorner[2] + (2*iz + 1)*h;

  double *p = (double *)calloc((g_pterms + 2)*(g_pterms + 2), sizeof(double));
  dcomplex *ptemp = (dcomplex *)calloc((g_pterms + 2)*(g_pterms + 2), sizeof(dcomplex));
  double *bi = (double *)calloc(g_pterms + 2, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 2, sizeof(dcomplex));
  if (p == 0 || ptemp == 0 || bi == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  const double precision = 1.0e-14;
  for (int i = 0; i < g_tboxes[ibox].npts; i++) {
    int ptr = g_tboxes[ibox].addr + i; 
    double *point = &g_fmmtargets[3*ptr]; 
    double *pot = &g_fmmpotential[ptr]; 
    double *field = &g_fmmfield[3*ptr]; 
    dcomplex cpz, fieldtemp[3] = {0}; 
    double rx = point[0] - x0y0z0[0];
    double ry = point[1] - x0y0z0[1];
    double rz = point[2] - x0y0z0[2];
    double proj = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 0.0: rz/d);
    ephi[0] = (proj <= precision*d ? 1.0 : rx/proj + _Complex_I*ry/proj);
    d *= g_beta; 
    for (int ell = 0; ell <= g_pterms; ell++) 
      ephi[ell + 1] = ephi[ell]*ephi[0];

    int ncalc;
    in(scale, d, g_pterms + 1, bi, &ncalc); 
    lgndr(g_pterms + 1, ctheta, p);

    *pot += creal(local[0])*bi[0];
    ptemp[0] = bi[0];

    for (int ell = 1; ell <= g_pterms; ell++) {
      double cp = bi[ell]*p[ell];
      ptemp[ell] = cp;
      *pot += creal(local[ell])*cp;
    }

    for (int m = 1; m <= g_pterms; m++) {
      int offset1 = m*(g_pterms + 2);
      int offset2 = m*(g_pterms + 1);
      for (int ell = m; ell <= g_pterms; ell++) {
	ptemp[ell + offset1] = bi[ell]*p[ell + offset1]*ephi[m - 1];
	*pot += 2.0*creal(local[ell + offset2]*ptemp[ell + offset1]);
      }
    }

    ptemp[g_pterms + 1] = 0;
    for (int m = 1; m <= g_pterms + 1; m++) {
      int offset1 = m*(g_pterms + 2);
      ptemp[g_pterms + 1 + offset1] = bi[g_pterms + 1]*p[g_pterms + 1 + offset1]*ephi[m - 1];
    }

    double rpotz = local[0]*scale;
    cpz = rpotz*ptemp[1 + g_pterms + 2]; 
    fieldtemp[0] = conj(cpz);
    fieldtemp[1] = rpotz*ptemp[1]; 
    fieldtemp[2] = -cpz;
    
    rpotz = local[1]/3.0;
    cpz = (rpotz*scale)*ptemp[2 + g_pterms + 2]; 
    fieldtemp[0] += conj(cpz);
    fieldtemp[1] += rpotz*(ptemp[0]/scale + 2.0*ptemp[2]*scale);
    fieldtemp[2] -= cpz;

    for (int ell = 2; ell <= g_pterms; ell++) {
      rpotz = local[ell]/(2*ell + 1);
      cpz = rpotz*(ptemp[ell + g_pterms + 1]/scale - ptemp[ell + g_pterms + 3]*scale);
      fieldtemp[0] -= conj(cpz);
      fieldtemp[1] += rpotz*
	(creal(ptemp[ell - 1])*ell/scale + (ell + 1)*creal(ptemp[ell + 1])*scale);
      fieldtemp[2] += cpz;
    }

    field[0] -= g_beta*creal(fieldtemp[0] - fieldtemp[2])/2;
    field[1] -= g_beta*creal((fieldtemp[0] + fieldtemp[2])*_Complex_I)/2.0;
    field[2] += g_beta*creal(fieldtemp[1]);
    
    cpz = local[g_pterms + 2]/3.0;
    fieldtemp[0] = cpz*((creal(ptemp[0])/scale - creal(ptemp[2])*scale)*2.0);
    cpz *= scale;
    fieldtemp[1] = cpz*ptemp[g_pterms + 4]; 
    fieldtemp[2] = cpz*ptemp[2 + 2*(g_pterms + 2)]; 

    for (int m = 1; m <= g_pterms - 2; m++) {
      int offset1 = m*(g_pterms + 1); 
      int offset2 = (m - 1)*(g_pterms + 2); 
      int offset3 = offset2 + (g_pterms + 2); 
      int offset4 = offset3 + (g_pterms + 2); 
      for (int ell = m + 2; ell <= g_pterms; ell++) {
	cpz = local[ell + offset1]/(2*ell + 1); 
	fieldtemp[0] += cpz*(ptemp[ell - 1 + offset2]*(ell + m - 1)*(ell + m)/scale - 
			     ptemp[ell + 1 + offset2]*(ell - m + 1)*(ell - m + 2)*scale); 
	fieldtemp[1] += cpz*(ptemp[ell - 1 + offset3]*(ell + m)/scale + 
			     ptemp[ell + 1 + offset3]*(ell - m + 1)*scale);
	fieldtemp[2] += cpz*(ptemp[ell - 1 + offset4]/scale - 
			     ptemp[ell + 1+ offset4]*scale); 
      }
    }

    for (int ell = 2; ell <= g_pterms; ell++) {    
      int offset1 = (ell - 1)*(g_pterms + 1);
      int offset2 = (ell - 2)*(g_pterms + 2) + ell;
      int offset3 = offset2 + (g_pterms + 2);
      int offset4 = offset3 + (g_pterms + 2);
      cpz = local[ell + offset1]/(2*ell + 1);
      fieldtemp[0] += cpz*(ptemp[-1 + offset2]*(ell + ell - 2)*(ell + ell - 1)/scale - 
			   ptemp[1 + offset2]*6*scale);
      fieldtemp[1] += cpz*(ptemp[-1 + offset3]*(ell + ell - 1)/scale + 
			   ptemp[1 + offset3]*2*scale);
      fieldtemp[2]  -= cpz*ptemp[1 + offset4]*scale;
    }

    for (int ell = 2; ell <= g_pterms; ell++) {
      int offset1 = ell*(g_pterms + 1);
      int offset2 = (ell - 1)*(g_pterms + 2) + ell;
      int offset3 = offset2 + (g_pterms + 2);
      int offset4 = offset3 + (g_pterms+2);
      cpz = local[ell + offset1]/(2*ell + 1);
      fieldtemp[0] += cpz*(ptemp[-1 + offset2]*(ell + ell - 1)*(ell + ell)/scale - 
			   ptemp[1 + offset2]*2.0*scale);
      fieldtemp[1]  += cpz*ptemp[1 + offset3]*scale;
      fieldtemp[2]  -= cpz*ptemp[1 + offset4]*scale;
    }
    
    field[0] -= g_beta*creal(fieldtemp[0] - fieldtemp[2]);
    field[1] -= g_beta*creal((fieldtemp[0] + fieldtemp[2])*_Complex_I);
    field[2] += 2.0*g_beta*creal(fieldtemp[1]);
  }
  free(p);
  free(ptemp);
  free(bi);
  free(ephi);
}

void sYukGreenFunction(const double *T, const double *S, const double charge, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz)
{
  double rx = T[0] - S[0];
  double ry = T[1] - S[1];
  double rz = T[2] - S[2]; 
  double rr = rx*rx + ry*ry + rz*rz;
  double rdis = sqrt(rr)*g_beta;
  double expr = cexp(-rdis);
  if ( rdis) {
    double term1 = charge/rdis*expr*M_PI_2;
    *pot += term1;
    double term2 = -term1*(1+rdis)/rr;
    *fx += rx*term2;
    *fy += ry*term2;
    *fz += rz*term2;
  }
}

void dYukGreenFunction(const double *T, const double *S, const double charge, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz)
{
  double rx = T[0] - S[0];
  double ry = T[1] - S[1];
  double rz = T[2] - S[2];
  double rr = rx*rx + ry*ry + rz*rz;
  double rdis = sqrt(rr)*g_beta;
  double expr = cexp(-rdis);
  if (rdis) {
    double term1 = charge/rdis*expr*M_PI_2;
    double term2 = -term1*(1.0 + rdis)/rr;
    double term3 = rx*dn[0] + ry*dn[1] + rz*dn[2];
    *pot += term2*term3;
    double term4 = term3*(3.0 + 3.0*rdis + rdis*rdis)/rr/rr;
    double term14 = term1*term4; 
    *fx = *fx + term2*dn[0] + term14*rx;
    *fy = *fy + term2*dn[1] + term14*ry;
    *fz = *fz + term2*dn[2] + term14*rz;
  }
}

void yhfrmini(void)
{
  double *fac = (double *)calloc(2*g_pterms + 2, sizeof(double)); 
  if (fac == 0) {
    ERRMSG("memory allocation failure");
    exit (-1);
  }

  fac[0] = 1; 
  for (int ell = 1; ell <= 2*g_pterms + 1; ell++) 
    fac[ell] = fac[ell - 1]*ell; 

  for (int m = 0; m <= g_pterms; m++) {
    int offset = m*(g_pterms + 1);
    for (int ell = m; ell <= g_pterms; ell++) {
      g_yukytop[ell + offset] = fac[ell - m]/fac[ell + m]*(2*ell + 1);
    }
  }

  free(fac);
}

void yhrotgen(void)
{
  //double *carray = (double *)calloc(pow(4*g_pterms + 1, 2), sizeof(double));
  double *carray = (double *)calloc((4*g_pterms + 1)*(4*g_pterms + 1), sizeof(double));
  if (carray==0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  bnlcft(carray, 4*g_pterms); 

  double theta = M_PI_2; 
  yhfstrtn(g_pterms, theta, carray, g_yukrdplus); 
  theta = -M_PI_2; 
  yhfstrtn(g_pterms, theta, carray, g_yukrdminus); 
  theta = acos(sqrt(3)/3); 
  yhfstrtn(g_pterms, theta, carray, g_yukrdsq3); 
  theta = acos(-sqrt(3)/3); 
  yhfstrtn(g_pterms, theta, carray, g_yukrdmsq3);

  free(carray);
}

void yhfstrtn(const int p, const double theta, const double *sqc, double *d)
{
  const double precision = 1.0e-19;
  const double ww = sqrt(2)/2; 
  double ctheta = cos(theta);
  ctheta = (fabs(ctheta) <= precision ? 0.0 : ctheta);
  double stheta = sin(-theta);
  stheta = (fabs(stheta) <= precision ? 0.0 : stheta);
  double hsthta = ww*stheta;
  double cthtap = ww*(1 + ctheta);
  double cthtan = -ww*(1 - ctheta);

  d[p*g_pgsz] = 1; 

  for (int ij = 1; ij <= p; ij++) {
    for (int im = -ij; im <=-1; im++) {
      int index = ij + (im + p)*g_pgsz; 
      
      d[index] = -sqc[ij - im + 2*(1 + 4*p)]*d[ij - 1 + (im + 1 + p)*g_pgsz]; 

      if (im > 1 - ij)
	d[index] += sqc[ij + im + 2*(1 + 4*p)]*d[ij - 1 + (im - 1 + p)*g_pgsz]; 

      d[index] *= hsthta;

      if (im > -ij) 
	d[index] += ctheta*sqc[ij + im + 4*p + 1]*sqc[ij - im + 4*p + 1]*d[index - 1]; 

      d[index] /= ij;
    }

    d[ij + p*g_pgsz] = d[ij - 1 + p*g_pgsz]*ctheta; 

    if (ij > 1) 
      d[ij + p*g_pgsz] += hsthta*sqc[ij + 2*(1 + 4*p)]*
	(d[ij - 1 + (-1 + p)*g_pgsz] + d[ij - 1 + (1 + p)*g_pgsz])/ij; 

    for (int im = 1; im <= ij; im++) {
      int index = ij + (im + p)*g_pgsz; 

      d[index] -= sqc[ij + im + 2*(1 + 4*p)]*d[ij - 1 + (im - 1 + p)*g_pgsz]; 

      if (im < ij - 1) 
	d[index] += sqc[ij - im + 2*(1 + 4*p)]*d[ij -1 + (im + 1 + p)*g_pgsz];

      d[index] *= hsthta;

      if (im < ij) 
	d[index] += ctheta*sqc[ij + im + 4*p + 1]*sqc[ij - im + 4*p + 1]*d[index - 1];
      
      d[index] /= ij;
    }

    for (int imp = 1; imp <= ij; imp++) {
      for (int im = -ij; im <= -1; im++) {
	int index1 = ij + imp*(p + 1) + (im + p)*g_pgsz; 
	int index2 = ij - 1 + (imp - 1)*(p + 1) + (im + p)*g_pgsz; 

	d[index1] = cthtan*sqc[ij - im + 2*(4*p + 1)]*d[index2 + g_pgsz]; 

	if (im > 1 - ij) 
	  d[index1] -= d[index2 - g_pgsz]*cthtap*sqc[ij + im + 2*(4*p + 1)]; 

	if (im > -ij) 
	  d[index1] += d[index2]*stheta*sqc[ij + im + 4*p + 1]*sqc[ij - im + 4*p + 1]; 

	d[index1] *= ww/sqc[ij + imp + 2*(4*p + 1)]; 
      }

      int index3 = ij + imp*(p + 1) + p*g_pgsz; 
      int index4 = ij - 1 + (imp - 1)*(p + 1) + p*g_pgsz; 
      d[index3] = ij*stheta*d[index4]; 

      if (ij > 1) 
	d[index3] -= sqc[ij + 2*(4*p + 1)]*
	  (d[index4 - g_pgsz]*cthtap + d[index4 + g_pgsz]*cthtan); 

      d[index3] *= ww/sqc[ij + imp + 2*(4*p + 1)]; 

      for (int im = 1; im <= ij; im++) {
	int index5 = ij + imp*(p + 1) + (im + p)*g_pgsz; 
	int index6 = ij - 1 + (imp - 1)*(p + 1) + (im + p)*g_pgsz; 

	d[index5] = d[index6 - g_pgsz]*cthtap*sqc[ij + im + 2*(4*p + 1)]; 
	
	if (im < ij - 1) 
	  d[index5] -= d[index6 + g_pgsz]*cthtan*sqc[ij - im + 2*(4*p + 1)]; 

	if (im < ij) 
	  d[index5] += d[index6]*stheta*sqc[ij + im + 4*p + 1]*sqc[ij - im + 4*p + 1];

	d[index5] *= ww/sqc[ij + imp + 2*(4*p + 1)]; 
      }
    }
  }
  
  double *fac = (double *)calloc(2*p + 1, sizeof(double));
  if (fac == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  fac[0] = 1;
  for (int ij = 1; ij <= 2*p; ij++) 
    fac[ij] = fac[ij - 1]*ij; 

  for (int ij = 0; ij <= p; ij++) {
    for (int im = 0; im <= ij; im++) {
      for (int imp = -ij; imp <= ij; imp++) {
	int impabs = fabs(imp); 
	d[ij + im*(p + 1) + (imp + p)*g_pgsz] *= 
	  sqrt(fac[ij + im]/fac[ij + impabs]*fac[ij - impabs]/fac[ij - im]);
      }
    }
  }

  free(fac);
}

void yukvwts(void)
{
  if (g_nlambs < 2 || g_nlambs > 39) {
    ERRMSG("wrong input for yukvwts()");
    exit (-1);
  }
  
  if (g_nlambs == 9) {
    g_yukrlams[0] = 0.99273996739714473469540223504736787e-01;
    g_yukrlams[1] = 0.47725674637049431137114652301534079e+00;
    g_yukrlams[2] = 0.10553366138218296388373573790886439e+01;
    g_yukrlams[3] = 0.17675934335400844688024335482623428e+01;
    g_yukrlams[4] = 0.25734262935147067530294862081063911e+01;
    g_yukrlams[5] = 0.34482433920158257478760788217186928e+01;
    g_yukrlams[6] = 0.43768098355472631055818055756390095e+01;
    g_yukrlams[7] = 0.53489575720546005399569367000367492e+01;
    g_yukrlams[8] = 0.63576578531337464283978988532908261e+01;
    g_yukwhts[0] = 0.24776441819008371281185532097879332e+00;
    g_yukwhts[1] = 0.49188566500464336872511239562300034e+00;
    g_yukwhts[2] = 0.65378749137677805158830324216978624e+00;
    g_yukwhts[3] = 0.76433038408784093054038066838984378e+00;
    g_yukwhts[4] = 0.84376180565628111640563702167128213e+00;
    g_yukwhts[5] = 0.90445883985098263213586733400006779e+00;
    g_yukwhts[6] = 0.95378613136833456653818075210438110e+00;
    g_yukwhts[7] = 0.99670261613218547047665651916759089e+00;
    g_yukwhts[8] = 0.10429422730252668749528766056755558e+01;
  } else if (g_nlambs == 18) { 
    g_yukrlams[0] = 0.52788527661177607475107009804560221e-01;
    g_yukrlams[1] = 0.26949859838931256028615734976483509e+00;
    g_yukrlams[2] = 0.63220353174689392083962502510985360e+00;
    g_yukrlams[3] = 0.11130756427760852833586113774799742e+01;
    g_yukrlams[4] = 0.16893949614021379623807206371566281e+01;
    g_yukrlams[5] = 0.23437620046953044905535534780938178e+01;
    g_yukrlams[6] = 0.30626998290780611533534738555317745e+01;
    g_yukrlams[7] = 0.38356294126529686394633245072327554e+01;
    g_yukrlams[8] = 0.46542473432156272750148673367220908e+01;
    g_yukrlams[9] = 0.55120938659358147404532246582675725e+01;
    g_yukrlams[10] = 0.64042126837727888499784967279992998e+01;
    g_yukrlams[11] = 0.73268800190617540124549122992902994e+01;
    g_yukrlams[12] = 0.82774009925823861522076185792684555e+01;
    g_yukrlams[13] = 0.92539718060248947750778825138695538e+01;
    g_yukrlams[14] = 0.10255602723746401139237605093512684e+02;
    g_yukrlams[15] = 0.11282088297877740146191172243561596e+02;
    g_yukrlams[16] = 0.12334067909676926788620221486780792e+02;
    g_yukrlams[17] = 0.13414920240172401477707353478763252e+02;
    g_yukwhts[0] = 0.13438265914335215112096477696468355e+00;
    g_yukwhts[1] = 0.29457752727395436487256574764614925e+00;
    g_yukwhts[2] = 0.42607819361148618897416895379137713e+00;
    g_yukwhts[3] = 0.53189220776549905878027857397682965e+00;
    g_yukwhts[4] = 0.61787306245538586857435348065337166e+00;
    g_yukwhts[5] = 0.68863156078905074508611505734734237e+00;
    g_yukwhts[6] = 0.74749099381426187260757387775811367e+00;
    g_yukwhts[7] = 0.79699192718599998208617307682288811e+00;
    g_yukwhts[8] = 0.83917454386997591964103548889397644e+00;
    g_yukwhts[9] = 0.87570092283745315508980411323136650e+00;
    g_yukwhts[10] = 0.90792943590067498593754180546966381e+00;
    g_yukwhts[11] = 0.93698393742461816291466902839601971e+00;
    g_yukwhts[12] = 0.96382546688788062194674921556725167e+00;
    g_yukwhts[13] = 0.98932985769673820186653756536543369e+00;
    g_yukwhts[14] = 0.10143828459791703888726033255807124e+01;
    g_yukwhts[15] = 0.10400365437416452252250564924906939e+01;
    g_yukwhts[16] = 0.10681548926956736522697610780596733e+01;
    g_yukwhts[17] = 0.11090758097553685690428437737864442e+01;
  }
}

void ympshftcoef(const double r0, const int lev)
{
  double scale = g_yuksfactor2[lev]; 
  double *c = &g_yukdcu[g_yukdcpgsz*lev]; 
  double *fac = (double *)calloc(2*g_pterms + 2, sizeof(double)); 
  double *bj = (double *)calloc(2*g_pterms + 3, sizeof(double)); 
  if (fac == 0 || bj == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  fac[0] = 1;
  for (int i = 1; i <= 2*g_pterms + 1; i++) 
    fac[i] = fac[i - 1]*i; 

  double r0k = r0*g_beta; 
  int ncalc; 
  in(r0k, r0k, 2*g_pterms + 2, bj, &ncalc);

  for (int mnew = 0; mnew <= g_pterms; mnew++) {
    for (int ellnew = mnew; ellnew <= g_pterms; ellnew++ ) {
      int offset = ellnew*(g_pterms + 1);
      for (int nn = mnew; nn <= g_pterms; nn++) {
	c[mnew + offset + nn*g_pgsz] = 0; 
	int np1 = MIN(nn, ellnew);
	for (int np = mnew; np <= np1; np++) {
	  c[mnew + offset + nn*g_pgsz] += pow(scale, nn - ellnew)*
	    pow(2, -ellnew - np)*pow(-1, ellnew + nn)*(2*ellnew + 1)*
	    fac[ellnew - mnew]/fac[np + mnew]*fac[nn + mnew]*fac[2*np]/fac[np]/
	    fac[np - mnew]/fac[ellnew - np]/fac[nn - np]*bj[ellnew + nn - np]*
	    pow(r0k, ellnew + nn - 2*np);
	}
      }
    }
  }
  free(fac);
  free(bj);
}

void ylcshftcoef(const double r0, const int lev)
{
  double scale = g_yuksfactor[lev]; 
  double *c = &g_yukdcd[g_yukdcpgsz*lev]; 
  double *fac = (double *)calloc(2*g_pterms + 2, sizeof(double)); 
  double *bj = (double *)calloc(2*g_pterms + 3, sizeof(double)); 
  if (fac == 0 || bj == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  fac[0] = 1; 
  for (int i = 1; i <= 2*g_pterms + 1; i++) 
    fac[i] = fac[i - 1]*i; 

  double r0k = r0*g_beta; 
  int ncalc; 
  in(r0k, r0k, 2*g_pterms + 2, bj, &ncalc); 

  for (int mnew = 0; mnew <= g_pterms; mnew++) {
    for (int ellnew = mnew; ellnew <= g_pterms; ellnew++) {
      int offset = ellnew*(g_pterms + 1);
      int nn; 
      for ( nn = mnew; nn <= g_pterms; nn++ ) {
	c[mnew + offset + nn*g_pgsz] = 0; 
	int np1 = MIN(nn, ellnew);
	for (int np = mnew; np <= np1; np++) {
	  c[mnew + offset + nn*g_pgsz] += pow(scale, ellnew - nn)*
	    pow(2, -ellnew - np)*(2*ellnew + 1)*fac[ellnew - mnew]/
	    fac[np + mnew]*fac[nn + mnew]*fac[2*np]/fac[np]/fac[np - mnew]/
	    fac[ellnew - np]/fac[nn - np]*bj[ellnew + nn - np]*
	    pow(r0k, ellnew + nn - 2*np);
	}
      }
    }
  }
  free(fac);
  free(bj);
}

void ymkfexp(const int lev)
{
  int *numphys = &g_yuknumphys[lev*g_nlambs]; 
  dcomplex *fexpe = &yuk_fexpe[15000*lev]; 
  dcomplex *fexpo = &yuk_fexpo[15000*lev]; 
  dcomplex *fexpback = &yuk_fexpback[15000*lev]; 

  int nexte = 0; 
  int nexto = 0;

  for (int i = 0; i < g_nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha/2;
    double halpha = 2.0*M_PI/nalpha;
    for (int j = 1; j <= nalpha2; j++) {
      double alpha = (j - 1)*halpha;
      for (int nm = 2; nm <= g_yuknumfour[i]; nm = nm + 2) {
	fexpe[nexte] = cexp((nm - 1)*_Complex_I*alpha);
	nexte++;
      }

      for (int nm = 3; nm <= g_yuknumfour[i]; nm = nm + 2) {
	fexpo[nexto] = cexp((nm - 1)*_Complex_I*alpha);
	nexto++;
      }
    }
  }

  int next = 0;
  for (int i = 0; i < g_nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha/2;
    double halpha = 2*M_PI/nalpha;
    for (int nm = 3; nm <= g_yuknumfour[i]; nm = nm + 2) {
      for (int j = 1; j <= nalpha2; j++) {
	double alpha = (j - 1)*halpha;
	fexpback[next] = cexp(-(nm - 1)*_Complex_I*alpha);
	next++;
      }
    }

    for (int nm=2; nm <= g_yuknumfour[i]; nm = nm + 2) {
      for (int j = 1; j <= nalpha2; j++) {
	double alpha = (j - 1)*halpha;
	fexpback[next] = cexp(-(nm - 1)*_Complex_I*alpha);
	next++;
      }
    }
  }
}

void yrlscini(const int lev)
{
  double scale = g_yuksfactor2[lev]; 
  double beta = g_yukbetascale[lev]; 
  double *rlsc = &g_yukrlsc[lev*g_pgsz*g_nlambs]; 

  for (int nell=1; nell <= g_nlambs; nell++) {
    double u1 = g_yukrlams[nell - 1]/beta + 1.0;
    lgndrgt1(scale, g_pterms, u1, &rlsc[g_pgsz*(nell - 1)]);
  }
}


void lgndrgt1(const double scale, const int nmax, const double x, double *y)
{
  double u = sqrt(x*x - 1)*scale; 
  double v = scale*x; 
  double w = scale*scale; 
  y[0] = 1; 
  y[1] = y[0]*v;
  for (int n = 2; n <= nmax; n++) 
    y[n] = ((2*n - 1)*v*y[n - 1] - (n - 1)*w*y[n - 2])/n;

  for (int m = 1; m < nmax; m++) {
    int offset = m*(nmax + 1);
    y[m + offset] = y[m - 1 + (m - 1)*(nmax + 1)]*u*(2*m - 1);
    y[m + 1 + offset] = y[m + offset]*(2*m + 1)*v;

    for (int n = m + 2; n <= nmax; n++) 
      y[n + offset] = ((2*n - 1)*v*y[n - 1 + offset] - (n + m - 1)*w*
		       y[n - 2 + offset])/(n - m);
  }

  y[nmax + nmax*(nmax + 1)] = y[nmax - 1 + (nmax - 1)*(nmax + 1)]*u*(2*nmax - 1);
}

void ymkexps(const int lev)
{
  double beta = g_yukbetascale[lev]; 
  int *numphys = &g_yuknumphys[g_nlambs*lev]; 
  dcomplex *xs = &yuk_xs[g_yuknexpmax*3*lev]; 
  dcomplex *ys = &yuk_ys[g_yuknexpmax*3*lev];
  double *zs = &yuk_zs[g_yuknexpmax*3*lev]; 

  int ntot = 0; 
  for (int nell = 0; nell < g_nlambs; nell++) {
    double w1 = g_yukrlams[nell] + beta;
    double w2 = sqrt(g_yukrlams[nell]*(g_yukrlams[nell] + beta*2));
    double hu = 2*M_PI/numphys[nell];
    for (int mth = 0; mth < numphys[nell]/2; mth++) {
      double u = mth*hu;
      int ncurrent = 3*(ntot + mth);
      zs[ncurrent] = exp(-w1);
      zs[ncurrent + 1] = zs[ncurrent]*zs[ncurrent];
      zs[ncurrent + 2] = zs[ncurrent + 1]*zs[ncurrent];
      xs[ncurrent] = cexp(w2*cos(u)*_Complex_I);
      xs[ncurrent + 1] = xs[ncurrent]*xs[ncurrent];
      xs[ncurrent + 2] = xs[ncurrent + 1]*xs[ncurrent];
      ys[ncurrent] = cexp(w2*sin(u)*_Complex_I);
      ys[ncurrent + 1] = ys[ncurrent]*ys[ncurrent];
      ys[ncurrent + 2] = ys[ncurrent + 1]*ys[ncurrent];
    }
    ntot += numphys[nell]/2;
  }
}
