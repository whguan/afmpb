#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "afmpb.h"
#include "adap_fmm.h"

void SetGNDGN(const double * const p, const double * const vn, const int nsources, 
	      const double * const sources, const double *charges, const int * const node, 
	      const int nelem, double * const gn, double * const dgn)
{
  static const double pi4 = M_PI*4;
  const double xi[] = {0.101286507323456, 0.797426958353087, 0.101286507323456,
		       0.470142064105115, 0.059715871789770, 0.470142064105115, 1.0/3.0};
  const double eta[] = {0.101286507323456, 0.101286507323456, 0.797426958353087,
			0.470142064105115, 0.470142064105115, 0.059715871789770, 1.0/3.0};

  double *targets = (double *)calloc(21*nelem, sizeof(double));
  double *targetsnd = (double *)calloc(21*nelem, sizeof(double));
  if (targets==NULL || targetsnd==NULL) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  int ntargets = 0;   
  for (int k = 0; k < nelem; k++) {
    int u0 = 3*k;
    int i1 = node[u0] - 1;
    int i2 = node[u0 + 1] - 1;
    int i3 = node[u0 + 2] - 1;
    for (int j = 0; j < 7; j++) {
      int u1 = 3*ntargets++;
      double zta = 1.0 - xi[j] - eta[j]; 	  
      targets[u1] = p[3*i1]*zta + p[3*i2]*xi[j] + p[3*i3]*eta[j];
      targets[u1 + 1] = p[3*i1 + 1]*zta + p[3*i2 + 1]*xi[j] + p[3*i3 + 1]*eta[j];
      targets[u1 + 2] = p[3*i1 + 2]*zta + p[3*i2 + 2]*xi[j] + p[3*i3 + 2]*eta[j];
      targetsnd[u1] = vn[u0];
      targetsnd[u1 + 1] = vn[u0 + 1];
      targetsnd[u1 + 2] = vn[u0 + 2];
    }
  }
  
  g_nsources = nsources; 
  g_ntargets = ntargets; 
  g_fmmsources = (double *)calloc(3*g_nsources, sizeof(double));
  g_fmmcharges = (double *)calloc(g_nsources, sizeof(double));
  g_fmmtargets = (double *)calloc(3*g_ntargets, sizeof(double));
  g_fmmpotential = (double *)calloc(g_ntargets, sizeof(double));
  g_fmmfield = (double *)calloc(3*g_ntargets, sizeof(double));
  g_fmmouternormal = (double *)calloc(3*g_ntargets, sizeof(double)); 
  g_mapsrc = (int *)calloc(g_nsources, sizeof(int));
  g_maptar = (int *)calloc(g_ntargets, sizeof(int));
  if (g_fmmsources == 0 || g_fmmcharges == 0 || g_fmmtargets == 0 ||
      g_fmmpotential == 0 || g_fmmfield == 0 || g_fmmouternormal == 0 ||
      g_mapsrc == 0 || g_maptar == 0){
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  BuildGraph(sources, nsources, targets, ntargets, 50);
  LapFMMInit(); 

  SourceToMultipole = sLapSourceToMultipole; 
  MultipoleToMultipole = LapMultipoleToMultipole; 
  MultipoleToExponential = LapMultipoleToExponential; 
  ExponentialToLocal = LapExponentialToLocal;
  LocalToLocal = LapLocalToLocal;
  LocalToTarget = LapLocalToTarget; 
  GreenFunction = sLapGreenFunction;   

  for (int i = 0; i < nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmsources[3*i] = sources[3*j]; 
    g_fmmsources[3*i + 1] = sources[3*j + 1];
    g_fmmsources[3*i + 2] = sources[3*j + 2]; 
    g_fmmcharges[i] = charges[j];
  }

  for (int i = 0; i < ntargets; i++) {
    int j = g_maptar[i];
    g_fmmtargets[3*i] = targets[3*j];
    g_fmmtargets[3*i + 1] = targets[3*j + 1];
    g_fmmtargets[3*i + 2] = targets[3*j + 2];
    g_fmmouternormal[3*i] = targetsnd[3*j]; 
    g_fmmouternormal[3*i + 1] = targetsnd[3*j + 1]; 
    g_fmmouternormal[3*i + 2] = targetsnd[3*j + 2]; 
  }

  AdapFMMCompute(1);

  for (int i = 0; i < ntargets; i++) {
    int j = g_maptar[i]; 
    gn[j] = g_fmmpotential[i]/pi4; 
    dgn[j] = -(g_fmmfield[3*i]*g_fmmouternormal[3*i] + 
	       g_fmmfield[3*i + 1]*g_fmmouternormal[3*i + 1] + 
	       g_fmmfield[3*i + 2]*g_fmmouternormal[3*i + 2])/pi4;
  }

  DestroyGraph();
  LapFMMClean();
  FMMClean();

  free(targets);
  free(targetsnd);
  return;
}

void CloseCoeff(const double * const rp0, const int jth, const double * const vnx0, 
		const double * const px, const int * const ippt, const double * const xnrm, 
		const double * const xwt, double * const ah, double * const bh, 
		double * const ch, double * const dh)
{
  static const double pi4 = M_PI*4;
  *ah = 0.0;
  *bh = 0.0;
  *ch = 0.0;
  *dh = 0.0;
  
  int i0 = ippt[jth]; 
  int i1 = ippt[jth + 1] - 1;
  
  for (int i = i0; i <= i1; i++) {
    double dx[3]; 
    dx[0] = px[3*i] - rp0[0];
    dx[1] = px[3*i + 1] - rp0[1]; 
    dx[2] = px[3*i + 2] - rp0[2]; 
    double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + sigma; 
    double r = sqrt(r2);
    double r3 = r2*r; 
    double r4 = r2*r2; 
    double r5 = r3*r2;
    double gr0 = 1.0/pi4; 
    double gr1 = gr0/r;
    double gr3 = gr0/r3;
    double gr5 = gr0/r5;
    double ur0 = cexp(-kap*r)/pi4;
    double ur1 = ur0/r;
    double ur2 = ur0/r2;
    double ur3 = ur0/r3;
    double ur4 = ur0/r4;
    double ur5 = ur0/r5;
    double ur3ur2 = ur3 + kap*ur2;
    double pur4ur3 = kap*(3*ur4 + kap*ur3) + 3*ur5;
    double gd[3], ud[3], gdd[3][3], udd[3][3];
    for (int j = 0; j < 3; j++) {
      gd[j] = -dx[j]*gr3;
      ud[j] = -dx[j]*(ur3 + kap*ur2);
      for (int t = 0; t < 3; t++) {
	int p0 = (t == j ? 1.0 : 0);
	double d0 = dx[j]*dx[t]; 
	gdd[j][t] = 3.0*d0*gr5 - p0*gr3;
	udd[j][t] = d0*pur4ur3 - p0*ur3ur2;
      }
    }

    double aht = gr1 - ur1; 
    double bht = 0.0;
    double cht = 0.0;
    double dht = 0.0;
    
    for (int j = 0; j < 3; j++) {
      int k0 = j + 3*i;
      bht += (gd[j]/dei - ud[j])*xnrm[k0];
      cht -= (gd[j] - ud[j]/dei)*vnx0[j];
      for (int k = 0; k < 3; k++) {
	dht -= vnx0[k]*(gdd[j][k] - udd[j][k])*xnrm[k0];
      }
    }

    *ah += aht*xwt[i];
    *bh += bht*xwt[i]; 
    *ch += cht*xwt[i]; 
    *dh += dht*xwt[i]; 
  }
  *dh /= dei;
}

void NSingularCoeff(const double * const rp0, const double * const rp1, const double arep,
		    const double arepn, const double * const arep3, const double * const vnp,
		    double * const ah, double * const bh, double * const ch, double * const dh)
{
  static const double pi4 = M_PI*4;
  double dx[3]; 
  dx[0] = rp1[0] - rp0[0];
  dx[1] = rp1[1] - rp0[1];
  dx[2] = rp1[2] - rp0[2];
  //dx[0:2] = rp1[0:2] - rp0[0:2];
  double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + sigma;
  double r = sqrt(r2);
  double r3 = r2*r;
  double r4 = r2*r2;
  double r5 = r3*r2;
  double gr0 = 1.0/pi4;
  double gr1 = gr0/r;
  double gr3 = gr0/r3;
  double gr5 = gr0/r5;
  double ur0 = cexp(-kap*r)/pi4;
  double ur1 = ur0/r;
  double ur2 = ur0/r2;
  double ur3 = ur0/r3;
  double ur4 = ur0/r4;
  double ur5 = ur0/r5;
  double ur3ur2 = ur3 + kap*ur2;
  double pur4ur3 = kap*(3*ur4 + kap*ur3) + 3*ur5;
  double gd[3], ud[3], gdd[3][3], udd[3][3]; 
  for (int j = 0; j < 3; j++) {
    gd[j] = -dx[j]*gr3;
    ud[j] = -dx[j]*(ur3 + kap*ur2);
    for (int t = 0; t < 3; t++) {
      int p0 = (t==j ? 1 : 0);
      double d0 = dx[j]*dx[t];
      gdd[j][t] = 3.0*d0*gr5 - p0*gr3;
      udd[j][t] = d0*pur4ur3 - p0*ur3ur2;
    }
  }

  *ah = gr1 - ur1;
  *bh = 0.0;
  *ch = 0.0;
  *dh = 0.0;

  for (int j = 0; j < 3; j++) {
    *bh += (gd[j]/dei - ud[j])*arep3[j]; 
    *ch -= (gd[j] - ud[j]/dei)*vnp[j];
    for (int k = 0; k < 3; k++) {
      *dh -= vnp[k]*(gdd[j][k] - udd[j][k])*arep3[j];
    }
  }

  *ah *= arep;
  *bh *= arepn;
  *ch *= arep;
  *dh *= arepn/dei;
}

double ComputeSolvationEnergy(const int nelem, const double * const gn, 
			      const double * const dgn, const double * const ff, 
			      const int npts, const double * const arel, 
			      const int * const node)
{
  const double xi[] = {0.101286507323456, 0.797426958353087, 0.101286507323456, 
		       0.470142064105115, 0.059715871789770, 0.470142064105115, 1.0/3.0};
  const double eta[] = {0.101286507323456, 0.101286507323456, 0.797426958353087,
			0.470142064105115, 0.470142064105115, 0.059715871789770, 1.0/3.0};
  const double w[] = {0.125939180544827, 0.125939180544827, 0.125939180544827, 
		      0.132394152788506, 0.132394152788506, 0.132394152788506, 0.225};
  double energy = 0; 
  int num = 0; 

  for (int k = 0; k < nelem; k++) {
    int u0 = 3*k;
    int i1 = node[u0] - 1; 
    int i2 = node[u0 + 1] - 1;
    int i3 = node[u0 + 2] - 1;
    double v = 0; 
    for (int j = 0; j < 7; j++) {
      double zta = 1.0 - xi[j] - eta[j]; 
      double f = ff[i1]*zta + ff[i2]*xi[j] + ff[i3]*eta[j];
      double h = ff[i1 + npts]*zta + ff[i2 + npts]*xi[j] + ff[i3 + npts]*eta[j];
      v += (gn[num]*h*dei - dgn[num]*f)*w[j]; 
      num++; 
    }
    energy += v*arel[k]; 
  }
  energy /= 2.0; 
  return energy; 
}


