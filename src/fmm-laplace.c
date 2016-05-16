#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "adap_fmm.h"

void LapFMMInit(void)
{
  g_iflu[0] = 3; g_ifld[0] = 1;
  g_iflu[1] = 4; g_ifld[1] = 2;
  g_iflu[2] = 2; g_ifld[2] = 4;
  g_iflu[3] = 1; g_ifld[3] = 3;
  g_iflu[4] = 3; g_ifld[4] = 1;
  g_iflu[5] = 4; g_ifld[5] = 2;
  g_iflu[6] = 2; g_ifld[6] = 4;
  g_iflu[7] = 1; g_ifld[7] = 3;

  g_lapnumphys = (int *)calloc(g_nlambs, sizeof(int));
  g_lapnumfour = (int *)calloc(g_nlambs, sizeof(int));
  g_lapwhts = (double *)calloc(g_nlambs, sizeof(double));
  g_laprlams = (double *)calloc(g_nlambs, sizeof(double));
  g_laprdplus = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  g_laprdminus = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  g_laprdsq3 = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  g_laprdmsq3 = (double *)calloc(g_pgsz*(2*g_pterms + 1), sizeof(double));
  //g_lapdc = (double *)calloc(pow(2*g_pterms + 1, 2), sizeof(double));
  g_lapdc = (double *)calloc((2*g_pterms + 1)*(2*g_pterms+1), sizeof(double));
  g_lapytopc = (double *)calloc(g_pgsz, sizeof(double));
  g_lapytopcs = (double *)calloc(g_pgsz, sizeof(double));
  g_lapytopcsinv = (double *)calloc(g_pgsz, sizeof(double));
  g_laprlsc = (double *)calloc(g_pgsz*g_nlambs, sizeof(double));
  if (g_lapnumphys == 0 || g_lapnumfour == 0 || g_lapwhts == 0 ||
      g_laprlams == 0 || g_laprdplus == 0 || g_laprdminus == 0 ||
      g_laprdsq3 == 0 || g_laprdmsq3 == 0 || g_lapdc == 0 || 
      g_lapytopc == 0 || g_lapytopcs == 0 || g_lapytopcsinv == 0 || 
      g_laprlsc == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  frmini(); 
  rotgen(); 
  lapvwts(); 
  numthetahalf(g_lapnumfour); 
  numthetafour(g_lapnumphys); 
  rlscini(); 

  g_lapnexptot = 0; 
  g_lapnthmax = 0; 
  g_lapnexptotp = 0; 
  for (int i = 1; i <= g_nlambs; i++) {
    g_lapnexptot += g_lapnumfour[i - 1];
    if (g_lapnumfour[i - 1] > g_lapnthmax) 
      g_lapnthmax = g_lapnumfour[i - 1];
    g_lapnexptotp += g_lapnumphys[i - 1];
  }
  g_lapnexptotp /= 2.0;
  g_lapnexpmax = MAX(g_lapnexptot, g_lapnexptotp) + 1;
  
  lap_xs = (dcomplex *)calloc(g_lapnexpmax*3, sizeof(dcomplex));
  lap_ys = (dcomplex *)calloc(g_lapnexpmax*3, sizeof(dcomplex));
  lap_zs = (double *)calloc(g_lapnexpmax*3, sizeof(double));
  lap_fexpe = (dcomplex *)calloc(15000, sizeof(dcomplex));
  lap_fexpo = (dcomplex *)calloc(15000, sizeof(dcomplex));
  lap_fexpback = (dcomplex *)calloc(15000, sizeof(dcomplex));
  lap_multipole = (dcomplex *)calloc((1 + g_nsboxes)*g_pgsz, sizeof(dcomplex));
  lap_local = (dcomplex *)calloc((1 + g_ntboxes)*g_pgsz, sizeof(dcomplex));
  lap_expu = (dcomplex *)calloc((1 + g_nsboxes)*g_lapnexpmax, sizeof(dcomplex));
  lap_expd = (dcomplex *)calloc((1 + g_nsboxes)*g_lapnexpmax, sizeof(dcomplex));
  lap_expn = (dcomplex *)calloc((1 + g_nsboxes)*g_lapnexpmax, sizeof(dcomplex));
  lap_exps = (dcomplex *)calloc((1 + g_nsboxes)*g_lapnexpmax, sizeof(dcomplex));
  lap_expe = (dcomplex *)calloc((1 + g_nsboxes)*g_lapnexpmax, sizeof(dcomplex));
  lap_expw = (dcomplex *)calloc((1 + g_nsboxes)*g_lapnexpmax, sizeof(dcomplex));
  g_lapscale = (double *)calloc(1 + g_nslev, sizeof(double));

  if (lap_xs == 0 || lap_ys == 0 || lap_zs == 0 || lap_fexpe == 0 ||
      lap_fexpo == 0 || lap_fexpback == 0 || lap_multipole == 0 || 
      lap_local == 0 || lap_expu == 0 || lap_expd == 0 || lap_exps == 0 || 
      lap_expe == 0 || lap_expw == 0 || g_lapscale == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  mkfexp();
  mkexps();
  g_lapscale[0] = 1.0/g_size; 
  for (int i = 1; i <= g_nslev; i++) 
    g_lapscale[i] = 2*g_lapscale[i - 1];  
}

void LapFMMClean(void)
{
  free(lap_xs);
  free(lap_ys);
  free(lap_zs);
  free(lap_fexpe);
  free(lap_fexpo);
  free(lap_fexpback);

  free(lap_multipole);
  free(lap_local);
  free(lap_expu);
  free(lap_expd);
  free(lap_expn);
  free(lap_exps);
  free(lap_expe);
  free(lap_expw);

  free(g_lapscale);
  free(g_lapnumphys);
  free(g_lapnumfour);
  free(g_lapwhts);
  free(g_laprlams);
  free(g_laprdplus);
  free(g_laprdminus);
  free(g_laprdsq3);
  free(g_laprdmsq3);
  free(g_lapdc);
  free(g_lapytopc);
  free(g_lapytopcs);
  free(g_lapytopcsinv);
  free(g_laprlsc);
}

void sLapSourceToMultipole(const int ibox)
{
  int ptr = g_sboxes[ibox].addr;
  int level = g_sboxes[ibox].level; 
  dcomplex *multipole = &lap_multipole[g_pgsz*ibox]; 
  double *sources = &g_fmmsources[3*ptr]; 
  double *charges = &g_fmmcharges[ptr]; 
  int nsources = g_sboxes[ibox].npts; 
  double scale = g_lapscale[level]; 
  double x0y0z0[3]; 
  double h = g_size/pow(2, level + 1);
  int ix = g_sboxes[ibox].idx; 
  int iy = g_sboxes[ibox].idy; 
  int iz = g_sboxes[ibox].idz; 
  x0y0z0[0] = g_bbcorner[0] + (2*ix + 1)*h; 
  x0y0z0[1] = g_bbcorner[1] + (2*iy + 1)*h;
  x0y0z0[2] = g_bbcorner[2] + (2*iz + 1)*h;

  double *powers = (double *)calloc(g_pterms + 1, sizeof(double));
  double *p = (double *)calloc(g_pgsz, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 1, sizeof(dcomplex));
  if (powers == 0 || p == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  const double precision = 1.0e-14;
  for (int i = 0; i < nsources; i++) {
    double rx = sources[3*i] - x0y0z0[0]; 
    double ry = sources[3*i + 1] - x0y0z0[1];
    double rz = sources[3*i + 2] - x0y0z0[2]; 
    double proj = rx*rx + ry*ry; 
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 1.0 : rz/d);
    ephi[0] = (proj <= precision*d ? 1.0 : rx/proj + _Complex_I*ry/proj);

    d *= scale; 
    powers[0] = 1.0;
    for (int ell = 1; ell <= g_pterms; ell++) {
      powers[ell] = powers[ell - 1]*d;
      ephi[ell] = ephi[ell - 1]*ephi[0];
    }

    multipole[0] += charges[i]; 

    lgndr(g_pterms, ctheta, p);

    for (int ell = 1; ell <= g_pterms; ell++) {
      double cp = charges[i]*powers[ell]*p[ell]; 
      multipole[ell] += cp;
    }

    for (int m = 1; m <= g_pterms; m++) {
      int offset = m*(g_pterms + 1);
      for (int ell = m; ell <= g_pterms; ell++) {
	double cp = charges[i]*powers[ell]*g_lapytopc[ell + offset]*p[ell + offset];
	multipole[ell + offset] += cp*conj(ephi[m - 1]);
      }
    }
  }

  free(powers);
  free(ephi);
  free(p);
}

void dLapSourceToMultipole(const int ibox)
{
  int ptr = g_sboxes[ibox].addr;
  int level = g_sboxes[ibox].level; 
  dcomplex *multipole = &lap_multipole[g_pgsz*ibox]; 
  double *sources = &g_fmmsources[3*ptr]; 
  double *charges = &g_fmmcharges[ptr];
  double *dn = &g_fmminnernormal[3*ptr]; 
  int nsources = g_sboxes[ibox].npts; 
  double scale = g_lapscale[level]; 
  double x0y0z0[3];
  double h = g_size/pow(2,level + 1);
  int ix = g_sboxes[ibox].idx;
  int iy = g_sboxes[ibox].idy;
  int iz = g_sboxes[ibox].idz;
  x0y0z0[0] = g_bbcorner[0] + (2*ix + 1)*h; 
  x0y0z0[1] = g_bbcorner[1] + (2*iy + 1)*h; 
  x0y0z0[2] = g_bbcorner[2] + (2*iz + 1)*h;

  double *powers = (double *)calloc(g_pterms + 3, sizeof(double));
  double *p = (double *)calloc(g_pgsz, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 2, sizeof(dcomplex));
  if (powers == 0 || p == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  const double precision = 1.0e-14;

  for (int i = 0; i < nsources; i++) {
    double rx = sources[3*i] - x0y0z0[0];
    double ry = sources[3*i + 1] - x0y0z0[1];
    double rz = sources[3*i + 2] - x0y0z0[2];
    double proj = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 1.0 : rz/d);
    ephi[0] = (proj <= precision*d ? 1.0 : rx/proj - _Complex_I*ry/proj);

    d = d*scale;
    powers[0] = 1.0;

    for (int ell = 1; ell <= g_pterms + 1; ell++) {
      powers[ell] = powers[ell - 1]*d;
      ephi[ell] = ephi[ell - 1]*ephi[0];
    }

    double cp = charges[i]*scale;
    dcomplex dnz1 = dn[3*i] - dn[3*i + 1]*_Complex_I;
    dcomplex dnz2 = conj(dnz1);
    multipole[1 + g_pterms + 1] += dnz1*cp*g_lapytopcs[1 + g_pterms + 1];
    multipole[1] -= cp*dn[3*i + 2];

    lgndr(g_pterms, ctheta, p);

    for (int ell = 1; ell <= g_pterms - 1; ell++) {
      cp = charges[i]*powers[ell]*p[ell]*scale;
      dcomplex cpz = cp*g_lapytopcsinv[ell + 1 + g_pterms + 1]*g_lapytopcs[ell]/2.0;
      multipole[ell + 1 + g_pterms + 1] += cpz*dnz1;
      multipole[ell + 1] -= cp*(ell + 1)*dn[3*i + 2];
    }

    for (int ell = 1; ell <= g_pterms - 1; ell++) {
      for (int m = 1; m <= ell; m++) {
	int offset = m*(g_pterms + 1); 
	cp = charges[i]*powers[ell]*g_lapytopc[ell + offset]*p[ell + offset]*scale;
	dcomplex cpz = cp*ephi[m - 1];
	double sr2 = g_lapytopcsinv[offset + ell - g_pterms]*g_lapytopcs[ell + offset];
	if (m == 1) {
	  cp = sr2*creal(cpz*dnz2);
	  multipole[ell + 1] -= cp;
	} else {
	  multipole[offset + ell - g_pterms] -= sr2/2.0*cpz*dnz2;
	}
	int index = offset + ell + g_pterms + 2; 
	sr2 = g_lapytopcsinv[index]*g_lapytopcs[ell + offset]/2.0;
	multipole[index] += sr2*cpz*dnz1;
	multipole[ell + 1 + offset] -= 
	  g_lapytopcsinv[index]*g_lapytopcs[ell + offset]*cpz*dn[3*i + 2];
      }
    }
  }
  free(powers);
  free(p);
  free(ephi);
}

void LapMultipoleToMultipole(const int pbox)
{
  dcomplex *pmultipole = &lap_multipole[g_pgsz*pbox]; 
  int lev = g_sboxes[pbox].level; 
  double sc1 = g_lapscale[lev + 1]; 
  double sc2 = g_lapscale[lev]; 

  static const dcomplex var[5]={1,-1 + _Complex_I, 1 + _Complex_I, 
				1 - _Complex_I, -1 - _Complex_I};
  const double arg = sqrt(2)/2.0;

  double *powers = (double *)calloc(g_pterms + 3,sizeof(double));
  dcomplex *mpolen = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *marray = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(g_pterms + 3, sizeof(dcomplex));

  if (powers == 0 || mpolen == 0 || marray == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int i = 0; i < 8; i++) {
    int cbox = g_sboxes[pbox].child[i]; 

    if (cbox) {
      int ifl = g_iflu[i]; 
      double *rd = (i < 4 ? g_laprdsq3 : g_laprdmsq3);
      dcomplex *cmultipole = &lap_multipole[g_pgsz*cbox];

      ephi[0] = 1.0; 
      ephi[1] = arg*var[ifl]; 
      double dd = -sqrt(3)/2.0; 
      powers[0] = 1.0; 

      for (int ell = 1; ell <= g_pterms + 1; ell++) {
	powers[ell] = powers[ell - 1]*dd; 
	ephi[ell + 1] = ephi[ell]*ephi[1]; 
      }

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
	    int index1 = ell + mp*(g_pterms + 1);
	    marray[index] += mpolen[index1]*rd[index1 + offset1] + 
	      conj(mpolen[index1])*rd[index1 + offset2];
	  }
	}
      }

      for (int k = 0; k <= g_pterms; k++) {
	int offset = k*(g_pterms + 1);
	for (int j = k; j <= g_pterms; j++) {
	  int index = offset + j;
	  mpolen[index] = marray[index]; 
	  for (int ell = 1; ell <= j - k; ell++) {
	    int index2 = j - k + ell*(2*g_pterms + 1);
	    int index3 = j + k + ell*(2*g_pterms + 1);
	    mpolen[index] += marray[index - ell]*powers[ell]*
	      g_lapdc[index2]*g_lapdc[index3]; 
	  }
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
	    marray[index] -= mpolen[index1]*rd[index1 + offset1] + 
	      conj(mpolen[index1])*rd[index1 + offset2];
	  }

	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1);
	    marray[index] += mpolen[index1]*rd[index1 + offset1] + 
	      conj(mpolen[index1])*rd[index1 + offset2];
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
      
      powers[0] = 1.0;
      dd = sc2/sc1;
      for (int ell = 1; ell <= g_pterms + 1; ell++) {
	powers[ell] = powers[ell-1]*dd;
      }
      
      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset;
	  mpolen[index] = ephi[m]*marray[index]*powers[ell];
	}
      }
      
      for (int m = 0; m < g_pgsz; m++) 
	pmultipole[m] += mpolen[m];
    }
  }
  free(ephi);
  free(powers);
  free(mpolen);
  free(marray);
}

void LapMultipoleToExponential(const int ibox)
{
  dcomplex *mw = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *mexpf1 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
  dcomplex *mexpf2 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
  if (mw == 0 || mexpf1 == 0 || mexpf2 == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  LapMultipoleToExponentialPhase1(&lap_multipole[g_pgsz*ibox], mexpf1, mexpf2);
  LapMultipoleToExponentialPhase2(mexpf1, &lap_expu[g_lapnexptotp*ibox]);
  LapMultipoleToExponentialPhase2(mexpf2, &lap_expd[g_lapnexptotp*ibox]); 

  rotz2y(&lap_multipole[g_pgsz*ibox], g_laprdminus, mw);
  LapMultipoleToExponentialPhase1(mw, mexpf1, mexpf2);
  LapMultipoleToExponentialPhase2(mexpf1, &lap_expn[g_lapnexptotp*ibox]);
  LapMultipoleToExponentialPhase2(mexpf2, &lap_exps[g_lapnexptotp*ibox]); 

  rotz2x(&lap_multipole[g_pgsz*ibox], g_laprdplus, mw); 
  LapMultipoleToExponentialPhase1(mw, mexpf1, mexpf2);
  LapMultipoleToExponentialPhase2(mexpf1, &lap_expe[g_lapnexptotp*ibox]);
  LapMultipoleToExponentialPhase2(mexpf2, &lap_expw[g_lapnexptotp*ibox]);

  free(mw);
  free(mexpf1);
  free(mexpf2);
}

void LapMultipoleToExponentialPhase1(const dcomplex * const multipole, 
				     dcomplex * const mexpu, dcomplex * const mexpd)
{
  int ntot = 0;
  for (int nell = 0; nell < g_nlambs; nell++) {
    double sgn = -1.0;
    dcomplex zeyep = 1.0;
    for (int mth = 0; mth <= g_lapnumfour[nell] - 1; mth++) {
      int ncurrent = ntot + mth ;
      dcomplex ztmp1 = 0.0;
      dcomplex ztmp2 = 0.0;
      sgn = -sgn;
      int offset = mth*(g_pterms + 1);
      int offset1 = offset + nell*g_pgsz;
      for (int nm = mth; nm <= g_pterms; nm = nm + 2) {
	ztmp1 += g_laprlsc[nm + offset1]*multipole[nm + offset];
      }
      for (int nm = mth + 1; nm <= g_pterms; nm = nm + 2) {
	ztmp2 += g_laprlsc[nm + offset1]*multipole[nm + offset];
      }

      mexpu[ncurrent] = (ztmp1 + ztmp2)*zeyep;
      mexpd[ncurrent] = sgn*(ztmp1 - ztmp2)*zeyep;
      zeyep = zeyep*_Complex_I;
    }
    ntot = ntot + g_lapnumfour[nell];
  }
}


void LapMultipoleToExponentialPhase2(const dcomplex * const mexpf, dcomplex * const mexpphys)
{
  int nftot = 0;
  int nptot = 0;
  int nexte = 0;
  int nexto = 0;
  for (int i = 0; i < g_nlambs; i++) {
    for (int ival = 0; ival < g_lapnumphys[i]/2; ival++) {
      mexpphys[nptot + ival] = mexpf[nftot];     
      for (int nm = 1; nm < g_lapnumfour[i]; nm = nm + 2) {
	double rt1 = cimag(lap_fexpe[nexte])*creal(mexpf[nftot + nm]);
	double rt2 = creal(lap_fexpe[nexte])*cimag(mexpf[nftot + nm]);
	double rtmp = 2*(rt1 + rt2);
	nexte++;
	mexpphys[nptot + ival] += rtmp*_Complex_I; 
      }
      
      for (int nm=2; nm < g_lapnumfour[i]; nm = nm + 2) {
	double rt1 = creal(lap_fexpo[nexto])*creal(mexpf[nftot + nm]);
	double rt2 = cimag(lap_fexpo[nexto])*cimag(mexpf[nftot + nm]);
	double rtmp = 2*(rt1 - rt2);
	nexto++;
	mexpphys[nptot + ival] += rtmp;
      }      
    }    
    nftot += g_lapnumfour[i]; 
    nptot += g_lapnumphys[i]/2; 
  }
}

void LapExponentialToLocal(const int ibox)
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
    int ilev = g_tboxes[ibox].level; 
    double scale = g_lapscale[ilev + 1];
    dcomplex *mw1 = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
    dcomplex *mw2 = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
    dcomplex *temp = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexpf1 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexpf2 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));

    // Process z-direction exponential expansions, and write results
    // to uall, u1234, dall, and d5678 lists. 
    dcomplex *mexuall = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexu1234 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexdall = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexd5678 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));

    if (mw1 == 0 || mw2 == 0 || temp == 0 || mexpf1 == 0 || mexpf2 == 0 ||
	mexuall == 0 || mexu1234 == 0 || mexdall == 0 || mexd5678 == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    MakeUList(g_lapnexptotp, lap_expd, uall, nuall, xuall, yuall, lap_xs, lap_ys, mexuall); 
    MakeUList(g_lapnexptotp, lap_expd, u1234, nu1234, x1234, y1234, lap_xs, lap_ys, mexu1234);
    MakeDList(g_lapnexptotp, lap_expu, dall, ndall, xdall, ydall, lap_xs, lap_ys, mexdall);
    MakeDList(g_lapnexptotp, lap_expu, d5678, nd5678, x5678, y5678, lap_xs, lap_ys, mexd5678);

    if (child[0]) {
      dcomplex *local = &lap_local[child[0]*g_pgsz];
      int iexp1 = 0; 
      int iexp2 = 0; 
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0; 

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 2]*scale; 
	iexp1++;
      }

      if (nu1234) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexu1234[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
      }

      if (iexp1)
	LapExponentialToLocalPhase1(temp, mexpf1);  
    
      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2); 
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[1]) {
      dcomplex *local = &lap_local[child[1]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 2]*conj(lap_xs[3*j])*scale;   
	iexp1++;
      }

      if (nu1234) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexu1234[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;   
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);

      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	LapExponentialToLocalPhase1(temp, mexpf2);
	iexp2++;
      }
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[2]) {
      dcomplex *local = &lap_local[child[2]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 2]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (nu1234) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexu1234[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);
    
      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[3]) {
      dcomplex *local = &lap_local[child[3]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 2]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
      }
      
      if (nu1234) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexu1234[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);
      
      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 1]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j];
      }
    }

    if (child[4]) {
      dcomplex *local = &lap_local[child[4]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;      

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1); 
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 2]*scale;
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexd5678[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
      }
      
      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[5]) {
      dcomplex *local = &lap_local[child[5]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }
      
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 2]*lap_xs[3*j]*scale;
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexd5678[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[6]) {
      dcomplex *local = &lap_local[child[6]*g_pgsz]; 
      int iexp1 = 0;
      int iexp2 = 0;      

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (ndall)  {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 2]*lap_ys[3*j]*scale;
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexd5678[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
      }
      
      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw1[j]; 
      }
    }

    if (child[7]) {
      dcomplex *local = &lap_local[child[7]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (nuall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexuall[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }
      
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (ndall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexdall[j]*lap_zs[3*j + 2]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
      }
      
      if (nd5678) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexd5678[j]*lap_zs[3*j + 1]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
      }
      
      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
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
    dcomplex *mexnall = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexn1256 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexn12 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexn56 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexsall = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexs3478 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexs34 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexs78 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));

    if (mexnall == 0 || mexn1256 == 0 || mexn12 == 0 || mexn56 == 0 ||
	mexsall == 0 || mexs3478 == 0 || mexs34 == 0 || mexs78 == 0 ) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    MakeUList(g_lapnexptotp, lap_exps, nall, nnall, xnall, ynall, lap_xs, lap_ys, mexnall); 
    MakeUList(g_lapnexptotp, lap_exps, n1256, nn1256, x1256, y1256, lap_xs, lap_ys, mexn1256);
    MakeUList(g_lapnexptotp, lap_exps, n12, nn12, x12, y12, lap_xs, lap_ys, mexn12); 
    MakeUList(g_lapnexptotp, lap_exps, n56, nn56, x56, y56, lap_xs, lap_ys, mexn56);
    MakeDList(g_lapnexptotp, lap_expn, sall, nsall, xsall, ysall, lap_xs, lap_ys, mexsall); 
    MakeDList(g_lapnexptotp, lap_expn, s3478, ns3478, x3478, y3478, lap_xs, lap_ys, mexs3478); 
    MakeDList(g_lapnexptotp, lap_expn, s34, ns34, x34, y34, lap_xs, lap_ys, mexs34); 
    MakeDList(g_lapnexptotp, lap_expn, s78, ns78, x78, y78, lap_xs, lap_ys, mexs78); 

    if (child[0]) {
      dcomplex *local = &lap_local[child[0]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 2]*scale;
	iexp1++;
      }
      
      if (nn1256) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn1256[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
      }
      
      if (nn12) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn12[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
      }
      
      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);
      
      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 1]*scale;
	iexp2++;  
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
    
    if (child[1]) {
      dcomplex *local = &lap_local[child[1]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 2]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }
      
      if (nn1256) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn1256[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (nn12) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn12[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);

      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[2]) {
      dcomplex *local = &lap_local[child[2]*g_pgsz];
      int iexp1 = 0; 
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }
      
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 2]*scale;
	iexp2++;
      }
      
      if (ns3478) {
	for (int j = 0; j < g_lapnexptotp; j++)
	  temp[j] += mexs3478[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
      }

      if (ns34) {
	for (int j = 0; j < g_lapnexptotp; j++)
	  temp[j] += mexs34[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
      }
    
      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);
         
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }
    
    if (child[3]) {
      dcomplex *local = &lap_local[child[3]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }
      
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 2]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (ns3478) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexs3478[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
      }
      
      if (ns34) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexs34[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (iexp2)
	LapExponentialToLocalPhase1(temp, mexpf2);
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }

    if (child[4]) {
      dcomplex *local = &lap_local[child[4]*g_pgsz];
      int iexp1 = 0; 
      int iexp2 = 0;

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 2]*conj(lap_xs[3*j])*scale;
	iexp1++;
      }

      if (nn1256) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn1256[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp1++;
      }

      if (nn56) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn56[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp1++;
      }
      
      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);
      
      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	LapExponentialToLocalPhase1(temp, mexpf2);
	iexp2++;
      }
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }
    
    if (child[5]) {
      dcomplex *local = &lap_local[child[5]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      for (int  j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 2]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
      }
      
      if (nn1256) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn1256[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
      }
      
      if (nn56) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexn56[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);
      
      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 1]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[6]) {
      dcomplex *local = &lap_local[child[6]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 2]*lap_xs[3*j]*scale;
	iexp2++;
      }

      if (ns3478) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexs3478[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp2++;
      }

      if (ns78) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexs78[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    } 

    if (child[7]) {
      dcomplex *local = &lap_local[child[7]*g_pgsz];
      int iexp1 = 0; 
      int iexp2 = 0;

      if (nnall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexnall[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j]*lap_ys[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }
      
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nsall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexsall[j]*lap_zs[3*j + 2]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (ns3478) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexs3478[j]*lap_zs[3*j + 1]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
      }
      
      if (ns78) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexs78[j]*lap_zs[3*j + 1]*lap_xs[3*j]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	roty2z(mw1, g_laprdplus, mw2);
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
    // w24, 268, w2, 24, w6, and w8 lists. 
    dcomplex *mexeall = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe1357 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe13 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe57 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe1 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe3 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe5 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexe7 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexwall = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw2468 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw24 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw68 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw2 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw4= (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw6 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    dcomplex *mexw8 = (dcomplex *)calloc(g_lapnexpmax, sizeof(dcomplex));
    if (mexeall == 0 || mexe1357 == 0 || mexe13 == 0 || mexe57 == 0 ||
	mexe1 == 0 || mexe3 == 0 || mexe5 == 0 || mexe7 == 0 ||
	mexwall == 0 || mexw2468 == 0 || mexw24 == 0 || mexw68 == 0 ||
	mexw2 == 0 || mexw4 == 0 || mexw6 == 0 || mexw8 == 0 ) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    MakeUList(g_lapnexptotp, lap_expw, eall, neall, xeall, yeall, lap_xs, lap_ys, mexeall);
    MakeUList(g_lapnexptotp, lap_expw, e1357, ne1357, x1357, y1357, lap_xs, lap_ys, mexe1357);
    MakeUList(g_lapnexptotp, lap_expw, e13, ne13, x13, y13, lap_xs, lap_ys, mexe13);
    MakeUList(g_lapnexptotp, lap_expw, e57, ne57, x57, y57, lap_xs, lap_ys, mexe57);
    MakeUList(g_lapnexptotp, lap_expw, e1, ne1, x1, y1, lap_xs, lap_ys, mexe1); 
    MakeUList(g_lapnexptotp, lap_expw, e3, ne3, x3, y3, lap_xs, lap_ys, mexe3);
    MakeUList(g_lapnexptotp, lap_expw, e5, ne5, x5, y5, lap_xs, lap_ys, mexe5);
    MakeUList(g_lapnexptotp, lap_expw, e7, ne7, x7, y7, lap_xs, lap_ys, mexe7);
    MakeDList(g_lapnexptotp, lap_expe, wall, nwall, xwall, ywall, lap_xs, lap_ys, mexwall);
    MakeDList(g_lapnexptotp, lap_expe, w2468, nw2468, x2468, y2468, lap_xs, lap_ys, mexw2468);
    MakeDList(g_lapnexptotp, lap_expe, w24, nw24, x24, y24, lap_xs, lap_ys, mexw24);
    MakeDList(g_lapnexptotp, lap_expe, w68, nw68, x68, y68, lap_xs, lap_ys, mexw68);
    MakeDList(g_lapnexptotp, lap_expe, w2, nw2, x2, y2, lap_xs, lap_ys, mexw2);
    MakeDList(g_lapnexptotp, lap_expe, w4, nw4, x4, y4, lap_xs, lap_ys, mexw4);
    MakeDList(g_lapnexptotp, lap_expe, w6, nw6, x6, y6, lap_xs, lap_ys, mexw6);
    MakeDList(g_lapnexptotp, lap_expe, w8, nw8, x8, y8, lap_xs, lap_ys, mexw8); 

    if (child[0]) {
      dcomplex *local = &lap_local[child[0]*g_pgsz]; 
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 2]*scale;
	iexp1++;
      }

      if (ne1357) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe1357[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
      }

      if (ne13) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe13[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
      }

      if (ne1) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe1[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);

      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2); 
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }

    if (child[1]) {  
      dcomplex *local = &lap_local[child[1]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 1]*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nwall) {      
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 2]*scale;
	iexp2++;
      }
      
      if (nw2468) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw2468[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
      }

      if (nw24) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw24[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
      }

      if (nw2) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw2[j]*lap_zs[3*j + 1]*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2); 
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[2]) {
      dcomplex *local = &lap_local[child[2]*g_pgsz]; 
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 2]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }
    
      if (ne1357) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe1357[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (ne13) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe13[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (ne3) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe3[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }
      
      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);

      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }

    if (child[3]) {
      dcomplex *local = &lap_local[child[3]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 1]*conj(lap_ys[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 2]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (nw2468) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw2468[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (nw24) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw24[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (nw4) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw4[j]*lap_zs[3*j + 1]*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[4]) {
      dcomplex *local = &lap_local[child[4]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;
      
      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 2]*lap_xs[3*j]*scale;
	iexp1++;
      }

      if (ne1357) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe1357[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp1++;
      }

      if (ne57) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe57[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp1++;
      }

      if (ne5) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe5[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);

      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j];
      }
    }

    if (child[5]) {
      dcomplex *local = &lap_local[child[5]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;

      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 1]*lap_xs[3*j]*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 2]*conj(lap_xs[3*j])*scale;
	iexp2++;
      }

      if (nw2468) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw2468[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp2++;
      }

      if (nw68) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw68[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp2++;
      }

      if (nw6) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw6[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[6]) {
      dcomplex *local = &lap_local[child[6]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 2]*lap_xs[3*j]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (ne1357) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe1357[j]*lap_zs[3*j + 1]*lap_xs[3*j]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (ne57) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe57[j]*lap_zs[3*j + 1]* lap_xs[3*j]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (ne7) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexe7[j]*lap_zs[3*j + 1]*lap_xs[3*j]*conj(lap_ys[3*j])*scale;
	iexp1++;
      }

      if (iexp1) 
	LapExponentialToLocalPhase1(temp, mexpf1);

      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*lap_ys[3*j]*scale;
	iexp2++;
	LapExponentialToLocalPhase1(temp, mexpf2);
      }

      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2);
	for (int j = 0; j < g_pgsz; j++) 
	  local[j] += mw2[j]; 
      }
    }

    if (child[7]) {
      dcomplex *local = &lap_local[child[7]*g_pgsz];
      int iexp1 = 0;
      int iexp2 = 0;
      
      if (neall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexeall[j]*lap_zs[3*j + 1]*lap_xs[3*j]*conj(lap_ys[3*j])*scale;
	iexp1++;
	LapExponentialToLocalPhase1(temp, mexpf1);
      }

      for (int j = 0; j < g_lapnexptotp; j++) 
	temp[j] = 0;

      if (nwall) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] = mexwall[j]*lap_zs[3*j + 2]*conj(lap_xs[3*j])*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (nw2468) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw2468[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (nw68) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw68[j]*lap_zs[3*j + 1]*conj(lap_xs[3*j])*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (nw8) {
	for (int j = 0; j < g_lapnexptotp; j++) 
	  temp[j] += mexw8[j]*lap_zs[3*j+1]*conj(lap_xs[3*j])*lap_ys[3*j]*scale;
	iexp2++;
      }

      if (iexp2) 
	LapExponentialToLocalPhase1(temp, mexpf2);
      
      if (iexp1 + iexp2) {
	LapExponentialToLocalPhase2(iexp2, mexpf2, iexp1, mexpf1, mw1);
	rotz2x(mw1, g_laprdminus, mw2);
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

void LapExponentialToLocalPhase1(const dcomplex *const mexpphys, dcomplex *const mexpf)
{
  int nftot = 0;
  int nptot = 0;
  int next = 0;
  for (int i = 0; i < g_nlambs; i++) {
    int nalpha = g_lapnumphys[i]; 
    int nalpha2 = nalpha/2;
    mexpf[nftot] = 0;
    for (int ival = 0; ival < nalpha2; ival++) {
      mexpf[nftot] += 2.0*creal(mexpphys[nptot + ival]);
    }
    mexpf[nftot] /= nalpha;

    for (int nm = 2; nm < g_lapnumfour[i]; nm = nm + 2) {
      mexpf[nftot + nm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
	double rtmp = 2*creal(mexpphys[nptot + ival]);
	mexpf[nftot + nm] += lap_fexpback[next]*rtmp;
	next += 1;
      }
      mexpf[nftot + nm] /= nalpha;
    }

    for (int nm = 1; nm < g_lapnumfour[i]; nm = nm + 2) {
      mexpf[nftot + nm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
	dcomplex ztmp = 2*cimag(mexpphys[nptot + ival])*_Complex_I;
	mexpf[nftot + nm] += lap_fexpback[next]*ztmp;
	next += 1;
      }
      mexpf[nftot + nm] /= nalpha;
    }
    nftot += g_lapnumfour[i]; 
    nptot += g_lapnumphys[i]/2; 
  }
}

void LapExponentialToLocalPhase2(const int iexpu, const dcomplex * const mexpu, 
				 const int iexpd, const dcomplex * const mexpd, 
				 dcomplex * const local)
{
  double *rlampow = (double *)calloc(g_pterms + 1, sizeof(double));
  dcomplex *zeye = (dcomplex *)calloc(g_pterms + 1, sizeof(dcomplex));
  dcomplex *mexpplus = (dcomplex *)calloc(g_lapnexptot, sizeof(dcomplex));
  dcomplex *mexpminus = (dcomplex *)calloc(g_lapnexptot, sizeof(dcomplex));
  if (rlampow == 0 || zeye == 0 || mexpplus == 0 || mexpminus == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  zeye[0] = 1.0;
  for (int i = 1; i <= g_pterms; i++) {
    zeye[i] = zeye[i - 1]*_Complex_I;
  }

  for (int i = 0; i < g_pgsz; i++) 
    local[i] = 0.0;

  for (int i = 0; i < g_lapnexptot; i++) {
    if (iexpu <= 0) {
      mexpplus[i] = mexpd[i];
      mexpminus[i] = mexpd[i];
    } else if (iexpd <= 0) {
      mexpplus[i] = mexpu[i];
      mexpminus[i] = -mexpu[i];
    } else {
      mexpplus[i] = mexpd[i] + mexpu[i];
      mexpminus[i] = mexpd[i] - mexpu[i];
    }
  }

  int ntot = 0;
  for (int nell = 0; nell < g_nlambs; nell++) {
    rlampow[0] = g_lapwhts[nell];
    double rmul = g_laprlams[nell];
    for (int j = 1; j <= g_pterms; j++) 
      rlampow[j] = rlampow[j - 1]*rmul;    

    int mmax = g_lapnumfour[nell]-1;
    for (int mth = 0; mth <= mmax; mth = mth + 2) {
      int offset = mth*(g_pterms + 1);
      for (int nm = mth; nm <= g_pterms; nm = nm + 2) {
	int index = offset + nm;
	int ncurrent = ntot + mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpplus[ncurrent];
      }

      for (int nm = mth + 1; nm <= g_pterms; nm = nm + 2) {
	int index = offset + nm;
	int ncurrent = ntot + mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpminus[ncurrent];
      }
    }

    for (int mth = 1; mth <= mmax; mth = mth + 2) {
      int offset = mth*(g_pterms + 1);
      for (int nm = mth + 1; nm <= g_pterms; nm = nm + 2) {
	int index = nm + offset;
	int ncurrent = ntot+mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpplus[ncurrent];
      }

      for (int nm = mth; nm <= g_pterms; nm = nm + 2) {
	int index = nm + offset;
	int ncurrent = ntot + mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpminus[ncurrent];
      }
    }
    ntot += g_lapnumfour[nell];
  }

  for (int mth = 0; mth <= g_pterms; mth++) {
    int offset = mth*(g_pterms + 1);
    for (int nm = mth; nm <= g_pterms; nm++) {
      int index = nm + offset;
      local[index] *= zeye[mth]*g_lapytopcs[index];
    }
  }
  free(rlampow);
  free(zeye);
  free(mexpplus);
  free(mexpminus);
}

void LapLocalToLocal(const int pbox)
{
  dcomplex *local = &lap_local[g_pgsz*pbox]; 
  int lev = g_tboxes[pbox].level; 
  double sc1 = g_lapscale[lev];
  double sc2 = g_lapscale[lev + 1];
  static dcomplex var[5]={1, 1 - _Complex_I, -1 - _Complex_I, 
			  -1 + _Complex_I, 1 + _Complex_I};
  const double arg=sqrt(2)/2.0;
  dcomplex *localn = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *marray = (dcomplex *)calloc(g_pgsz, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(1 + g_pterms, sizeof(dcomplex));
  double *powers = (double *)calloc(1 + g_pterms, sizeof(double));
  if (localn == 0 || marray == 0 || ephi == 0 || powers == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int i = 0; i < 8; i++) {
    int cbox = g_tboxes[pbox].child[i]; 
    if (cbox) {     
      int ifl = g_ifld[i]; 
      double *rd = (i < 4 ? g_laprdsq3 : g_laprdmsq3);
      ephi[0] = 1.0;
      ephi[1] = arg*var[ifl];
      double dd = -sqrt(3)/4.0;
      powers[0] = 1.0;

      for (int ell = 1; ell <= g_pterms; ell++) 
	powers[ell] = powers[ell - 1]*dd;
      
      for (int ell = 2; ell <= g_pterms; ell++) 
	ephi[ell] = ephi[ell-1]*ephi[1];
      
      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset;
	  localn[index] = conj(ephi[m])*local[index];
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
	for (int j = k; j<= g_pterms; j++) {
	  int index = j + offset;
	  localn[index] = marray[index];
	  for (int ell = 1; ell <= g_pterms - j; ell++) {
	    int index1 = ell + index;
	    int index2 = ell + j + k + ell*(2*g_pterms + 1);
	    int index3 = ell + j - k + ell*(2*g_pterms + 1);
	    localn[index] += marray[index1]*powers[ell]*g_lapdc[index2]*g_lapdc[index3];
	  }
	}
      }
      
      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	int offset1 = (m + g_pterms)*g_pgsz;
	int offset2 = (-m + g_pterms)*g_pgsz;
	for (int ell = m; ell <= g_pterms; ell++) {
	  int index = ell + offset;
	  marray[index] = localn[ell]*rd[ell + offset1];
	  for (int mp = 1; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1);
	    marray[index] -= localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }

	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1);
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
	  int index = ell + offset;
	  marray[index] = -localn[ell]*rd[ell + offset1];
	  for (int mp = 1; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1);
	    marray[index] += localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }

	  for (int mp = 2; mp <= ell; mp = mp + 2) {
	    int index1 = ell + mp*(g_pterms + 1);
	    marray[index] -= localn[index1]*rd[index1 + offset1] + 
	      conj(localn[index1])*rd[index1 + offset2];
	  }
	}
      }

      powers[0] = 1.0;
      dd = sc1/sc2;
      for (int ell = 1; ell <= g_pterms; ell++) 
	powers[ell] = powers[ell - 1]*dd;
           
      for (int m = 0; m <= g_pterms; m++) {
	int offset = m*(g_pterms + 1);
	for (int ell = m; ell <= g_pterms; ell++ ) {
	  int index = offset + ell;
	  localn[index] = ephi[m]*marray[index]*powers[ell];
	}
      }

      dcomplex *clocal = &lap_local[g_pgsz*cbox];
      for (int m = 0; m < g_pgsz; m++) 
	clocal[m] += localn[m];
    }
  }
  free(ephi);
  free(powers);
  free(localn);
  free(marray);
}

void LapLocalToTarget(const int ibox)
{
  dcomplex *local = &lap_local[g_pgsz*ibox]; 
  int level = g_tboxes[ibox].level; 
  double scale = g_lapscale[level]; 
  double h = g_size/pow(2, level + 1);
  double x0y0z0[3];
  int ix = g_tboxes[ibox].idx; 
  int iy = g_tboxes[ibox].idy; 
  int iz = g_tboxes[ibox].idz; 
  x0y0z0[0] = g_bbcorner[0] + (2*ix + 1)*h; 
  x0y0z0[1] = g_bbcorner[1] + (2*iy + 1)*h; 
  x0y0z0[2] = g_bbcorner[2] + (2*iz + 1)*h; 
  double *p = (double *)calloc(g_pgsz, sizeof(double));
  double *powers = (double *)calloc(1 + g_pterms, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(1 + g_pterms, sizeof(dcomplex));
  if (p == 0 || powers == 0 || ephi == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }
  const double precision = 1.0e-14;

  for (int i = 0; i < g_tboxes[ibox].npts; i++) {
    double field0, field1, field2, rloc, cp, rpotz = 0.0;
    dcomplex cpz, zs1 = 0.0, zs2 = 0.0, zs3 = 0.0;
    int ptr = g_tboxes[ibox].addr + i; 
    double *point = &g_fmmtargets[3*ptr]; 
    double *pot = &g_fmmpotential[ptr]; 
    double *field = &g_fmmfield[3*ptr]; 
    double rx = point[0] - x0y0z0[0];
    double ry = point[1] - x0y0z0[1];
    double rz = point[2] - x0y0z0[2];
    double proj  = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 0.0 : rz/d); 
    ephi[0] = (proj <= precision*d ? 1.0 : rx/proj + _Complex_I*ry/proj);
    d *= scale;
    double dd = d;

    powers[0] = 1.0;
    for (int ell = 1; ell <= g_pterms; ell++) {
      powers[ell] = dd; 
      dd *= d; 
      ephi[ell] = ephi[ell - 1]*ephi[0]; 
    }

    lgndr(g_pterms, ctheta, p); 
    *pot += creal(local[0]);
    
    field2 = 0.0;
    for (int ell = 1; ell <= g_pterms; ell++) {
      rloc = creal(local[ell]);
      cp = rloc*powers[ell]*p[ell];
      *pot += cp;
      cp = powers[ell - 1]*p[ell - 1]*g_lapytopcs[ell - 1];
      cpz = local[ell + g_pterms + 1]*cp*g_lapytopcsinv[ell + g_pterms + 1];
      zs2 = zs2 + cpz;
      cp = rloc*cp*g_lapytopcsinv[ell];
      field2 += cp;
    }

    for (int ell = 1; ell <= g_pterms; ell++) {
      for (int m = 1; m <= ell; m++) {
	int index = ell + m*(g_pterms + 1); 
	cpz = local[index]*ephi[m - 1];
	rpotz += creal(cpz)*powers[ell]*g_lapytopc[index]*p[index];
      }
    
      for (int m = 1; m <= ell - 1; m++) {
	int index1 = ell + m*(g_pterms + 1);
	int index2 = index1 - 1;
	zs3 += local[index1]*ephi[m - 1]*powers[ell - 1]*
	  g_lapytopc[index2]*p[index2]*g_lapytopcs[index2]*g_lapytopcsinv[index1];
      }

      for (int m = 2; m <= ell; m++) {
	int index1 = ell + m*(g_pterms+1); 
	int index2 = ell - 1 + (m - 1)*(g_pterms + 1);
	zs2 += local[index1]*ephi[m - 2]*g_lapytopcs[index2]*
	  g_lapytopcsinv[index1]*powers[ell - 1]*g_lapytopc[index2]*p[index2];
      }

      for (int m = 0; m <= ell - 2; m++) {
	int index1 = ell + m*(g_pterms + 1); 
	int index2 = ell - 1 + (m + 1)*(g_pterms + 1);
	zs1 += local[index1]*ephi[m]*g_lapytopcs[index2]*
	  g_lapytopcsinv[index1]*powers[ell - 1]*g_lapytopc[index2]*p[index2];
      }
    }

    *pot += 2.0*rpotz;
    field0 = creal(zs2 - zs1);
    field1 = -cimag(zs2 + zs1);
    field2 += 2.0*creal(zs3);

    field[0] += field0*scale;
    field[1] += field1*scale;
    field[2] -= field2*scale;
  } 
  free(powers);
  free(ephi);
  free(p);
}

void sLapGreenFunction(const double *T, const double *S, const double charge, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz)
{
  double rx = T[0] - S[0];
  double ry = T[1] - S[1];
  double rz = T[2] - S[2]; 
  double rr = rx*rx + ry*ry + rz*rz;
  double rdis = sqrt(rr);

  if (rr) {
    *pot += charge/rdis;
    double rmul = charge/(rdis*rr);
    *fx += rmul*rx;
    *fy += rmul*ry;
    *fz += rmul*rz;
  }

}

void dLapGreenFunction(const double *T, const double *S, const double charge, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz)
{
  double rx = T[0] - S[0];
  double ry = T[1] - S[1];
  double rz = T[2] - S[2];
  double rr = rx*rx + ry*ry + rz*rz;
  double rdis = sqrt(rr);

  if (rr) {
    double term1 = -charge/(rdis*rr);
    double term3 = rx*dn[0] + ry*dn[1] + rz*dn[2];
    double term2 = 3.0/rr*term1*term3;
    *pot += term1*term3;
    *fx = *fx - dn[0]*term1 + rx*term2;
    *fy = *fy - dn[1]*term1 + ry*term2;
    *fz = *fz - dn[2]*term1 + rz*term2; 
  }
}

void frmini(void)
{
  double *fac = (double *)calloc(1 + 2*g_pterms, sizeof(double));
  if (fac == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  double d = 1.0;
  fac[0] = d;
  for (int ell = 1; ell <= 2*g_pterms; ell++) {
    d *= sqrt(ell);
    fac[ell] = d;
  }

  g_lapytopcs[0] = 1.0;
  g_lapytopcsinv[0] = 1.0;
  for (int m = 0; m <= g_pterms; m++) {
    int offset = m*(g_pterms + 1);
    for (int ell = m; ell <= g_pterms; ell++) {
      g_lapytopc[ell + offset] = fac[ell - m]/fac[ell + m];
      g_lapytopcsinv[ell + offset] = fac[ell - m]*fac[ell + m];
      g_lapytopcs[ell + offset] = 1.0/g_lapytopcsinv[ell + offset];
    }
  }

  free(fac);
}

void rotgen(void)
{
  bnlcft(g_lapdc, 2*g_pterms);
  double theta = M_PI/2;
  fstrtn(g_pterms, g_laprdplus, g_lapdc, theta);
  theta = -theta;
  fstrtn(g_pterms, g_laprdminus, g_lapdc, theta);
  theta = acos(sqrt(3)/3);
  fstrtn(g_pterms, g_laprdsq3, g_lapdc, theta);
  theta = acos(-sqrt(3)/3);
  fstrtn(g_pterms, g_laprdmsq3, g_lapdc, theta);
}

void fstrtn(const int p, double *d, const double *sqc, const double theta)
{
  const double precision = 1.0e-19;
  const double ww = sqrt(2)/2;
  double ctheta = cos(theta);
  ctheta = (fabs(ctheta) <= precision ? 0.0 : ctheta);
  double stheta = sin(-theta);
  stheta = (fabs(stheta) <= precision ? 0.0 : stheta);
  double hsthta = ww*stheta;
  double cthtap = ww*(1.0+ctheta);
  double cthtan = -ww*(1.0-ctheta);

  d[p*g_pgsz] = 1.0;

  for (int ij = 1; ij <= p; ij++) {
    for (int im = -ij; im <= -1; im++) {
      int index = ij + (im + p)*g_pgsz; 
      d[index] = -sqc[ij - im + 2*(1 + 2*p)]*d[ij-1 + (im + 1 + p)*g_pgsz];

      if (im > 1 - ij) 
	d[index] += sqc[ij + im + 2*(1 + 2*p)]*d[ij - 1 + (im - 1 + p)*g_pgsz];      

      d[index] *= hsthta;

      if (im > -ij) 
	d[index] += d[ij - 1 + (im + p)*g_pgsz]*ctheta*
	  sqc[ij + im + 2*p + 1]*sqc[ij - im + 2*p + 1];
      
      d[index] /= ij;
    }

    d[ij + p*g_pgsz] = d[ij - 1 + p*g_pgsz]*ctheta;

    if (ij > 1) 
      d[ij + p*g_pgsz] += hsthta*sqc[ij + 2*(1 + 2*p)]*
	(d[ij - 1 + (-1 + p)*g_pgsz] + d[ij - 1 + (1 + p)*g_pgsz])/ij;
    
    for (int im = 1; im <= ij; im++) {
      int index = ij + (im + p)*g_pgsz; 
      d[index] = -sqc[ij + im + 2*(1 + 2*p)]*d[ij - 1 + (im - 1 + p)*g_pgsz];

      if (im < ij-1) 
	d[index] += sqc[ij - im + 2*(1 + 2*p)]*d[ij - 1 + (im + 1 + p)*g_pgsz];
      
      d[index] *= hsthta;

      if (im < ij) 
	d[index] += d[ij- 1 + (im + p)*g_pgsz]*ctheta*
	  sqc[ij + im + 2*p + 1]*sqc[ij - im + 2*p + 1];
      
      d[index] /= ij;
    }

    for (int imp = 1; imp <= ij; imp++) {
      for (int im = -ij; im <= -1; im++) {
	int index1 = ij + imp*(p + 1) + (im + p)*g_pgsz; 
	int index2 = ij - 1 + (imp - 1)*(p + 1) + (im + p)*g_pgsz; 
	
	d[index1] = d[index2 + g_pgsz]*cthtan*sqc[ij - im + 2*(2*p + 1)]; 

	if (im > 1 - ij) 
	  d[index1] -= d[index2 - g_pgsz]*cthtap*sqc[ij + im + 2*(2*p + 1)]; 

	if (im > -ij) 
	  d[index1] += d[index2]*stheta*sqc[ij + im + 2*p + 1]*sqc[ij - im + 2*p + 1];

	d[index1] *= ww/sqc[ij + imp + 2*(2*p + 1)];
      }      

      int index3 = ij + imp*(p + 1) + p*g_pgsz; 
      int index4 = ij - 1 + (imp - 1)*(p + 1) + p*g_pgsz; 

      d[index3] = ij*stheta*d[index4];
      if (ij > 1) 
	d[index3] -= sqc[ij + 2*(2*p + 1)]*
	  (d[index4 - g_pgsz]*cthtap + d[index4 + g_pgsz]*cthtan);

      d[index3] *= ww/sqc[ij + imp + 2*(2*p + 1)]; 

      for (int im = 1; im <= ij; im++) {
	int index5 = ij + imp*(p + 1) + (im + p)*g_pgsz; 
	int index6 = ij - 1 + (imp - 1)*(p + 1) + (im + p)*g_pgsz; 

	d[index5] = d[index6 - g_pgsz]*cthtap*sqc[ij + im + 2*(2*p + 1)]; 

	if (im < ij - 1) 
	  d[index5] -= d[index6 + g_pgsz]*cthtan*sqc[ij - im + 2*(2*p + 1)]; 
	
	if (im < ij) 
	  d[index5] += d[index6]*stheta*sqc[ij + im + 2*p + 1]*sqc[ij - im + 2*p + 1];

	d[index5] *= ww/sqc[ij + imp + 2*(2*p + 1)];
      }
    }
  }
}

void lapvwts(void)
{
  if (g_nlambs < 2 || g_nlambs > 39) {
    ERRMSG("wrong input for lapvwts()");
    exit (-1);
  }
  
  if (g_nlambs == 9) {
    g_laprlams[0] = 0.99273996739714473469540223504736787e-01;
    g_laprlams[1] = 0.47725674637049431137114652301534079e+00;
    g_laprlams[2] = 0.10553366138218296388373573790886439e+01;
    g_laprlams[3] = 0.17675934335400844688024335482623428e+01;
    g_laprlams[4] = 0.25734262935147067530294862081063911e+01;
    g_laprlams[5] = 0.34482433920158257478760788217186928e+01;
    g_laprlams[6] = 0.43768098355472631055818055756390095e+01;
    g_laprlams[7] = 0.53489575720546005399569367000367492e+01;
    g_laprlams[8] = 0.63576578531337464283978988532908261e+01;
    g_lapwhts[0] = 0.24776441819008371281185532097879332e+00;
    g_lapwhts[1] = 0.49188566500464336872511239562300034e+00;
    g_lapwhts[2] = 0.65378749137677805158830324216978624e+00;
    g_lapwhts[3] = 0.76433038408784093054038066838984378e+00;
    g_lapwhts[4] = 0.84376180565628111640563702167128213e+00;
    g_lapwhts[5] = 0.90445883985098263213586733400006779e+00;
    g_lapwhts[6] = 0.95378613136833456653818075210438110e+00;
    g_lapwhts[7] = 0.99670261613218547047665651916759089e+00;
    g_lapwhts[8] = 0.10429422730252668749528766056755558e+01;
  } else if (g_nlambs == 18) { 
    g_laprlams[0] = 0.52788527661177607475107009804560221e-01;
    g_laprlams[1] = 0.26949859838931256028615734976483509e+00;
    g_laprlams[2] = 0.63220353174689392083962502510985360e+00;
    g_laprlams[3] = 0.11130756427760852833586113774799742e+01;
    g_laprlams[4] = 0.16893949614021379623807206371566281e+01;
    g_laprlams[5] = 0.23437620046953044905535534780938178e+01;
    g_laprlams[6] = 0.30626998290780611533534738555317745e+01;
    g_laprlams[7] = 0.38356294126529686394633245072327554e+01;
    g_laprlams[8] = 0.46542473432156272750148673367220908e+01;
    g_laprlams[9] = 0.55120938659358147404532246582675725e+01;
    g_laprlams[10] = 0.64042126837727888499784967279992998e+01;
    g_laprlams[11] = 0.73268800190617540124549122992902994e+01;
    g_laprlams[12] = 0.82774009925823861522076185792684555e+01;
    g_laprlams[13] = 0.92539718060248947750778825138695538e+01;
    g_laprlams[14] = 0.10255602723746401139237605093512684e+02;
    g_laprlams[15] = 0.11282088297877740146191172243561596e+02;
    g_laprlams[16] = 0.12334067909676926788620221486780792e+02;
    g_laprlams[17] = 0.13414920240172401477707353478763252e+02;
    g_lapwhts[0] = 0.13438265914335215112096477696468355e+00;
    g_lapwhts[1] = 0.29457752727395436487256574764614925e+00;
    g_lapwhts[2] = 0.42607819361148618897416895379137713e+00;
    g_lapwhts[3] = 0.53189220776549905878027857397682965e+00;
    g_lapwhts[4] = 0.61787306245538586857435348065337166e+00;
    g_lapwhts[5] = 0.68863156078905074508611505734734237e+00;
    g_lapwhts[6] = 0.74749099381426187260757387775811367e+00;
    g_lapwhts[7] = 0.79699192718599998208617307682288811e+00;
    g_lapwhts[8] = 0.83917454386997591964103548889397644e+00;
    g_lapwhts[9] = 0.87570092283745315508980411323136650e+00;
    g_lapwhts[10] = 0.90792943590067498593754180546966381e+00;
    g_lapwhts[11] = 0.93698393742461816291466902839601971e+00;
    g_lapwhts[12] = 0.96382546688788062194674921556725167e+00;
    g_lapwhts[13] = 0.98932985769673820186653756536543369e+00;
    g_lapwhts[14] = 0.10143828459791703888726033255807124e+01;
    g_lapwhts[15] = 0.10400365437416452252250564924906939e+01;
    g_lapwhts[16] = 0.10681548926956736522697610780596733e+01;
    g_lapwhts[17] = 0.11090758097553685690428437737864442e+01;
  }
}

void rlscini(void)
{
  double *fac = (double *)calloc(2*g_pterms + 1, sizeof(double));
  double *rlampow = (double *)calloc(g_pterms + 1, sizeof(double));
  if (fac == 0 || rlampow == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  fac[0] = 1;
  for (int i = 1; i <= 2*g_pterms; i++)
    fac[i] = fac[i-1]*sqrt(i);
 
  for (int nell = 0; nell < g_nlambs; nell++) {
    double rmul = g_laprlams[nell];
    rlampow[0] = 1;
    for (int j = 1;  j <= g_pterms; j++)
      rlampow[j] = rlampow[j - 1]*rmul;      
    for (int j = 0; j <= g_pterms; j++) {
      for (int k = 0; k <= j; k++) {
	g_laprlsc[j + k*(g_pterms + 1) + nell*g_pgsz] = 
	  rlampow[j]/fac[j - k]/fac[j + k];
      }
    }    
  }
    
  free(fac);
  free(rlampow);
}

void mkfexp(void)
{
  int nexte = 0; 
  int nexto = 0; 
  for (int i = 0; i < g_nlambs; i++) {
    int nalpha = g_lapnumphys[i]; 
    int nalpha2 = nalpha/2; 
    double halpha = 2.0*M_PI/nalpha; 
    for (int j = 1; j <= nalpha2; j++) {
      double alpha = (j - 1)*halpha; 
      for (int nm = 2; nm <= g_lapnumfour[i]; nm = nm + 2) {
	lap_fexpe[nexte] = cexp((nm - 1)*_Complex_I*alpha);
	nexte++;
      }

      for (int nm = 3; nm <= g_lapnumfour[i]; nm = nm + 2) {
	lap_fexpo[nexto] = cexp((nm - 1)*_Complex_I*alpha);
	nexto++;
      }
    }
  }

  int next = 0; 
  for (int i = 0; i < g_nlambs; i++) {
    int nalpha = g_lapnumphys[i]; 
    int nalpha2 = nalpha/2; 
    double halpha = 2.0*M_PI/nalpha; 
    for (int nm = 3; nm <= g_lapnumfour[i]; nm = nm + 2) {
      for (int j = 1; j <= nalpha2; j++) {
	double alpha = (j - 1)*halpha; 
	lap_fexpback[next] = cexp(-(nm - 1)*_Complex_I*alpha);
	next++;
      }
    }

    for (int nm = 2; nm <= g_lapnumfour[i]; nm = nm + 2) {
      for (int j = 1; j <= nalpha2; j++) {
	double alpha = (j - 1)*halpha; 
	lap_fexpback[next] = cexp(-(nm - 1)*_Complex_I*alpha);
	next++;
      }
    }
  }
}

void mkexps(void)
{
  int ntot = 0; 
  for (int nell = 0; nell < g_nlambs; nell++) {
    double hu = 2.0*M_PI/g_lapnumphys[nell]; 
    for (int  mth = 0; mth < g_lapnumphys[nell]/2; mth++) {
      double u = mth*hu; 
      int ncurrent = 3*(ntot + mth);
      lap_zs[ncurrent] = exp(-g_laprlams[nell]);
      lap_zs[ncurrent + 1] = lap_zs[ncurrent]*lap_zs[ncurrent]; 
      lap_zs[ncurrent + 2] = lap_zs[ncurrent]*lap_zs[ncurrent + 1];
      lap_xs[ncurrent] = cexp(_Complex_I*cos(u)*g_laprlams[nell]);
      lap_xs[ncurrent + 1] = lap_xs[ncurrent]*lap_xs[ncurrent];
      lap_xs[ncurrent + 2] = lap_xs[ncurrent + 1]*lap_xs[ncurrent]; 
      lap_ys[ncurrent] = cexp(_Complex_I*sin(u)*g_laprlams[nell]);
      lap_ys[ncurrent + 1] = lap_ys[ncurrent]*lap_ys[ncurrent]; 
      lap_ys[ncurrent + 2] = lap_ys[ncurrent + 1]*lap_ys[ncurrent]; 
    }
    ntot += g_lapnumphys[nell]/2; 
  }
}

