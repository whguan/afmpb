#ifndef ADAP_FMM_H
#define ADAP_FMM_H

#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#endif

#define CAP 64
#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ERRMSG(str) printf("Error in %s, line %d: %s\n", __FILE__, __LINE__, str); 

typedef double complex dcomplex;
typedef void (*op)(const int boxid);
typedef void (*kernel)(const double *t, const double *s, const double q, 
		       const double *dn, double *p, double *fx, 
		       double *fy, double *fz);

typedef struct fmmnode {
  int level, boxid, parent, child[8], nchild, idx, idy, idz, npts, addr, 
    *list1, *list3, *list4, *list5; 
} fmmnode; 

typedef struct linklist {
  int ndata; 
  int data[CAP]; 
  struct linklist *next;
} linklist; 

op SourceToMultipole, MultipoleToMultipole, MultipoleToExponential, 
  ExponentialToLocal, LocalToLocal, LocalToTarget; 

kernel GreenFunction; 

int g_nslev, g_nsboxes, g_ntlev, g_ntboxes, *g_mapsrc, *g_maptar, 
  *g_scontent, *g_tcontent, g_nsources, g_ntargets, g_pterms, g_nlambs, g_pgsz;

double g_bbcenter[3], g_bbcorner[3], g_size, *g_fmmsources, *g_fmmcharges, 
  *g_fmminnernormal, *g_fmmtargets, *g_fmmpotential, *g_fmmfield, *g_fmmouternormal;  

fmmnode *g_sboxes, *g_tboxes; 

void BuildGraph(const double * const sources, const int nsources, 
		const double * const targets, const int ntargets, 
		const int s);

void DestroyGraph(void);

void PartitionBox(const fmmnode * const ibox, const double * const points, 
		  const double h, int * const counts, int * const addrs, 
		  int * const map);

void BuildList(const int ibox);

void BuildFinerList(fmmnode * const ibox, linklist *list1, int nlist1);

void PushStack(linklist * const head, int * const nelems, const int ibox);

void PopAll(linklist *head, int * const list);

int PopStack(linklist *head, int * const nelems);

int IfAdjacent(const fmmnode * const box1, const fmmnode * const box2);

void BuildMergedList2(const int tboxid, int * const uall, int * const nuall, 
		      int * const xuall, int * const yuall, 
		      int * const u1234, int * const nu1234, int * const x1234, int * const y1234,
		      int * const dall, int * const ndall, int * const xdall, int * const ydall, 
		      int * const d5678, int * const nd5678, int * const x5678, int * const y5678, 
		      int * const nall, int * const nnall, int * const xnall, int * const ynall, 
		      int * const n1256, int * const nn1256, int * const x1256, int * const y1256, 
		      int * const n12, int * const nn12, int * const x12, int * const y12, 
		      int * const n56, int * const nn56, int * const x56, int * const y56, 
		      int * const sall, int * const nsall, int * const xsall, int * const ysall, 
		      int * const s3478, int * const ns3478, int * const x3478, int * const y3478,
		      int * const s34, int * const ns34, int * const x34, int * const y34, 
		      int * const s78, int * const ns78, int * const x78, int * const y78, 
		      int * const eall, int * const neall, int * const xeall, int * const yeall, 
		      int * const e1357, int * const ne1357, int * const x1357, int * const y1357, 
		      int * const e13, int * const ne13, int * const x13, int * const y13,
		      int * const e57, int * const ne57, int * const x57, int * const y57, 
		      int * const e1, int * const ne1, int * const x1, int * const y1, 
		      int * const e3, int * const ne3, int * const x3, int * const y3, 
		      int * const e5, int * const ne5, int * const x5, int * const y5, 
		      int * const e7, int * const ne7, int * const x7, int * const y7, 
		      int * const wall, int * const nwall, int * const xwall, int * const ywall, 
		      int * const w2468, int * const nw2468, int * const x2468, int * const y2468, 
		      int * const w24, int * const nw24, int * const x24, int * const y24, 
		      int * const w68, int * const nw68, int * const x68, int * const y68, 
		      int * const w2, int * const nw2, int * const x2, int * const y2, 
		      int * const w4, int * const nw4, int * const x4, int * const y4, 
		      int * const w6, int * const nw6, int * const x6, int * const y6, 
		      int * const w8, int * const nw8, int * const x8, int * const y8);

void UpdateList(int * const list, int * const nlist, int * const xoff, 
		int * const yoff, int entry, int ix, int iy);

void MakeUList(const int nexpo, const dcomplex * const expo, 
	       const int * const list, const int nlist, const int * const xoff, 
	       const int * const yoff, const dcomplex * const xs, 
	       const dcomplex * const ys, dcomplex * const mexpo);

void MakeDList(const int nexpo, const dcomplex * const expo, 
	       const int * const list, const int nlist, const int * const xoff, 
	       const int * const yoff, const dcomplex * const xs, 
	       const dcomplex * const ys, dcomplex * const mexpo);

void BuildDirectList13(int **directlist);

void CleanDirectList13(int **directlist);

void AdapFMMCompute(const int mode);

void AggregateSweep(const int ibox);

void DisAggregateSweep(const int ibox);

void ProcessList13(const int ibox);

void ProcessList4(const int ibox);

void DirectEvaluation(const int tbox, const int sbox);

void FMMClean(void);

void numthetahalf(int *numfour);

void numthetafour(int *numphys);

void lgndr(const int n, const double x, double *y);

void rotz2y(const dcomplex *multipole, const double *rd, dcomplex *mrotate);

void roty2z(const dcomplex *multipole, const double *rd, dcomplex *mrotate);

void rotz2x(const dcomplex *multipole, const double *rd, dcomplex *mrotate);

void in(const double scal, const double x, const int nb, double *b, int *ncalc);

void ribesl_(const double *x, double *alpha, int *nb, int *ize, double *b, int *ncalc); 

void bnlcft(double *c, const int p);

int *g_lapnumphys, *g_lapnumfour, g_lapnexptot, g_lapnthmax, 
  g_lapnexptotp, g_lapnexpmax, g_iflu[8], g_ifld[8]; 

double *g_lapwhts, *g_laprlams, *g_laprdplus, *g_laprdminus, *g_laprdsq3, 
  *g_laprdmsq3, *g_lapdc, *g_lapytopc, *g_lapytopcs, *g_lapytopcsinv, 
  *g_laprlsc, *g_lapscale, *lap_zs; 

dcomplex *lap_xs, *lap_ys, *lap_fexpe, *lap_fexpo, *lap_fexpback, *lap_multipole, 
  *lap_local, *lap_expu, *lap_expd, *lap_expn, *lap_exps, *lap_expe, *lap_expw; 

void LapFMMInit(void);

void sLapSourceToMultipole(const int ibox);

void dLapSourceToMultipole(const int ibox);

void LapMultipoleToMultipole(const int pbox);

void LapMultipoleToExponential(const int ibox);

void LapMultipoleToExponentialPhase1(const dcomplex * const multipole, 
				     dcomplex * const mexpu, dcomplex * const mexpd);

void LapMultipoleToExponentialPhase2(const dcomplex * const mexpf, dcomplex * const mexpphys);

void LapExponentialToLocal(const int ibox);

void LapExponentialToLocalPhase1(const dcomplex *const mexpphys, dcomplex *const mexpf);

void LapExponentialToLocalPhase2(const int iexpu, const dcomplex * const mexpu, 
				 const int iexpd, const dcomplex * const mexpd, 
				 dcomplex * const local);

void LapLocalToLocal(const int pbox);

void LapLocalToTarget(const int ibox);

void sLapGreenFunction(const double *T, const double *S, const double Q, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz);

void dLapGreenFunction(const double *T, const double *S, const double Q, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz);

void LapFMMClean(void);

void frmini(void);

void rotgen(void);

void bnlcft(double *c, const int p);

void fstrtn(const int p, double *d, const double *sqc, const double theta);

void lapvwts(void); 

void rlscini(void);

void mkfexp(void);

void mkexps(void);

int *g_yuknumfour, *g_yuknumphys, g_yifl[8], g_yuknexptot, g_yuknthmax, 
  g_yukmnexptotp, g_yuknexpmax, *g_yukvnexptot, *g_yukvnexptotp, 
  *g_yukvnthmax, g_yukdcpgsz; 

double g_beta, *g_yukytop, *g_yukwhts, *g_yukrlams, *g_yukrdplus, 
  *g_yukrdminus, *g_yukrdsq3, *g_yukrdmsq3, *g_yukdcu, *g_yukdcd, 
  *g_yuksfactor, *g_yuksfactor2, *g_yukbetascale, *g_yukrlsc, 
  g_yukscale, *yuk_zs;

dcomplex *yuk_xs, *yuk_ys, *yuk_fexpe, *yuk_fexpo, *yuk_fexpback, *yuk_multipole, 
  *yuk_local, *yuk_expu, *yuk_expd, *yuk_expn, *yuk_exps, *yuk_expe, *yuk_expw; 

void YukFMMInit(const double beta);

void sYukSourceToMultipole(const int ibox);

void dYukSourceToMultipole(const int ibox);

void YukMultipoleToMultipole(const int pbox);

void YukMultipoleToExponential(const int ibox);

void YukMultipoleToExponentialPhase1(const dcomplex * const multipole, const int level, 
				     dcomplex * const mexpu, dcomplex * const mexpd);

void YukMultipoleToExponentialPhase2(const dcomplex * const mexpf, const int level, 
				     dcomplex * const mexpphys);

void YukExponentialToLocal(const int ibox);

void YukExponentialToLocalPhase1(const dcomplex * const mexpphys, const int level, 
				 dcomplex * const mexpf);

void YukExponentialToLocalPhase2(const int iexpu, const dcomplex * const mexpu, 
				 const int iexpd, const dcomplex * const mexpd, 
				 const int level, dcomplex * const local);

void YukLocalToLocal(const int pbox);

void YukLocalToTarget(const int ibox);

void sYukGreenFunction(const double *T, const double *S, const double Q, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz);

void dYukGreenFunction(const double *T, const double *S, const double Q, 
		       const double *dn, double *pot, double *fx, 
		       double *fy, double *fz);

void YukFMMClean(void);

void yhfrmini(void);

void yhrotgen(void);

void yukvwts(void);

void yhfstrtn(const int p, const double theta, const double *sqc, double *d);

void ympshftcoef(const double r0, const int lev); 

void ylcshftcoef(const double r0, const int lev);

void ribesl_(const double *x, double *alpha, int *nb, int *ize, double *b, int *ncalc); 

void ymkfexp(const int lev); 

void yrlscini(const int lev); 

void lgndrgt1(const double scale, const int nmax, const double x, double *y);

void ymkexps(const int lev); 
#endif
