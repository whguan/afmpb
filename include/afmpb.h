#ifndef AFMPB_H
#define AFMPB_H

#include "adap_fmm.h"
#define FMMACCURACY 3

typedef void(*IterativeSolver)(int *n, double *rhs, double *sol, int *ipar, 
			       double *fpar, double *w); 

IterativeSolver KrylovSpaceSolver; 

static const double unitfactor = 4171.8;

int nmol, ipotw, iforce, iintere, iselfe, ipotdx;

double di, de, concentration, temperature, dei, kap, sigma, cut1, cut2; 

int ParseCommandLine(int argc, char **argv, char * const input, char * const output, 
		     char * const surfp, char * const vmesh, char * const volp);

int SetParameters(const char * const input, int *meshfmt, char * const fpqr, 
		  char * const meshf);

int ProcessPQRFile(const char * const fname); 

void ReadPQRFile(const int ncharges, double * const xq, 
		 double * const q, double * const r); 

void ReadMeshFile(const char * const fname, double **p_, int **node_, double **vnp_, 
		  const int meshfmt, int *npts, int *nelems);

void RemoveIsolatedNodes(const int meshfmt, double * const p, int * const node, 
			 double * const vnp, const int nelems, int *npts);

void ProcessElementGeometry(const int meshfmt, const int npts, const int nelem, 
			    int *node, double *p, double *vn, double *arel, 
			    double *arep, double *arepn, double *arep3, 
			    double *vnp, int *ippt, double *px, double *xnrm, double *xwt);

void SetGNDGN(const double * const p, const double * const vn, const int nsources, 
	      const double * const sources, const double *charges, const int * const node, 
	      const int nelem, double * const gn, double * const dgn);

void SolvePoissonBoltzmann(const int nnodes, const double * const node, 
			   const double * const nodevn, const double * const arep, 
			   const double * const arepn, const double * const arep3, 
			   const int ncharges, const double * const charges, 
			   const double * const zcharges, const double * const px, 
			   const int * const ippt, const double * const xnrm, 
			   const double * const xwt, double * const fhiso, double * const rhs);

void SetRightHandSide(const int nsources, const double * const sources, 
		      const double * const charges, const int ntargets, 
		      const double * const targets, const double * const normal, 
		      double * const gn, double * const dgn);

void SetSelfMatrix(int **directlist, const double * const targets, 
		   const double * const normal, const double * const sources, 
		   const double * const arep, const double * const arepn, 
		   const double * const arep3, const double * const px, 
		   const int * const ippt, const double * const xnrm, 
		   const double * const xwt, double **A, double **B, 
		   double **C, double **D);

void CloseCoeff(const double * const rp0, const int jth, const double * const vnx0, 
		const double * const px, const int * const ippt, const double * const xnrm, 
		const double * const xwt, double * const ah, double * const bh, 
		double * const ch, double * const dh);

void NSingularCoeff(const double * const rp0, const double * const rp1, const double arep,
		    const double arepn, const double * const arep3, const double * const vnp, 
		    double * const ah, double * const bh, double * const ch, double * const dh);

void MatrixVectorMultiply(double * const U, const double * const arep, 
			  int **directlist, double **A, double **B, double **C, 
			  double **D, double * const V);

void gmres_(int *n, double *rhs, double *sol, int *ipar, double *fpar, double *w);

void bcgstab_(int *n, double *rhs, double *sol, int *ipar, double *fpar, double *w); 

void WritePotential(const char * const fname, const int npts, 
		    const int nelems, const int * const node, 
		    const double * const p, const double * const vnp, 
		    const double * const f);

double ComputeSolvationEnergy(const int nelem, const double * const gn, 
			      const double * const dgn, const double * const ff, 
			      const int npts, const double * const arel, 
			      const int * const node);
#endif

