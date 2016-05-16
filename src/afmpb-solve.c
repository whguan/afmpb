#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
//#include <cilk/reducer_opadd.h>
#include "afmpb.h"
#include <sys/time.h>

void SolvePoissonBoltzmann(const int nnodes, const double * const node, 
			   const double * const nodevn, const double * const arep, 
			   const double * const arepn, const double * const arep3, 
			   const int ncharges, const double * const charges, 
			   const double * const zcharges, const double * const px, 
			   const int * const ippt, const double * const xnrm, 
			   const double * const xwt, double * const fhiso, double * const rhs)
{
  // Step 1: Generate the right-hand side of the linear system. 
  struct timeval tic, toc;
  double elapsed;
   gettimeofday(&tic,0);
  SetRightHandSide(ncharges, zcharges, charges, nnodes, node, nodevn, &rhs[0], &rhs[nnodes]); 
  gettimeofday(&toc, 0);
  double right = (double)((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
//      printf("... SetRightHand  time:%f seconds \n", right);

  double pi4de = M_PI*4*de;

  for (int i = 0; i < nnodes; i++) {
    rhs[i] /= pi4de;
    rhs[i + nnodes] /= -pi4de;
  }

  // Step 2: Generate the adaptive FMM graph for matrix-vector
  // multiplication in the GMRES procedure. Different from the FMM
  // called in step 1, variable "node" is considered as the sources
  // and the targets. This step also initializes the FMM for Laplace
  // and Yukawa kernel functions.
  g_nsources = nnodes;
  g_ntargets = nnodes;
  g_mapsrc = (int *)calloc(g_nsources, sizeof(int));
  g_maptar = (int *)calloc(g_ntargets, sizeof(int));
  g_fmmsources = (double *)calloc(3*g_nsources, sizeof(double));
  g_fmmcharges = (double *)calloc(g_nsources, sizeof(double));
  g_fmminnernormal = (double *)calloc(3*g_nsources, sizeof(double));
  g_fmmtargets = (double *)calloc(3*g_ntargets, sizeof(double));
  g_fmmpotential = (double *)calloc(g_ntargets, sizeof(double));
  g_fmmfield = (double *)calloc(3*g_ntargets, sizeof(double));
  g_fmmouternormal = (double *)calloc(3*g_ntargets, sizeof(double));
  if (g_mapsrc == 0 || g_maptar == 0 || g_fmmsources == 0 || g_fmmcharges == 0 || 
      g_fmminnernormal == 0 || g_fmmtargets == 0 || g_fmmpotential == 0 || 
      g_fmmfield == 0 || g_fmmouternormal == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  BuildGraph(node, nnodes, node, nnodes, 50);
  LapFMMInit();
  YukFMMInit(kap);

  // Step 3: Generate self-matrix for near field interactions. 
  gettimeofday(&tic,0);
  int **directlist = (int **)calloc(1 + g_ntboxes, sizeof(int *));
  double **directA = (double **)calloc(g_ntargets, sizeof(double *));
  double **directB = (double **)calloc(g_ntargets, sizeof(double *));
  double **directC = (double **)calloc(g_ntargets, sizeof(double *));
  double **directD = (double **)calloc(g_ntargets, sizeof(double *));
  if (directlist == 0 || directA == 0 || directB == 0 || 
      directC == 0 || directD == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  BuildDirectList13(directlist);

  SetSelfMatrix(directlist, node, nodevn, node, arep, arepn, arep3, px, 
		ippt, xnrm, xwt, directA, directB, directC, directD);
  gettimeofday(&toc, 0);
  double time = (double)((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("... SetSelfMatrix  time:%f seconds \n", time);

  // When FMM is used in an iterative solver setting, sources and
  // targets are unchanged while charges are updated in each
  // step. Relocate the points of the same box into contiguous memory
  // space.   
  cilk_for (int i = 0; i < g_nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmsources[3*i] = node[3*j]; 
    g_fmmsources[3*i + 1] = node[3*j + 1];
    g_fmmsources[3*i + 2] = node[3*j + 2]; 
    g_fmminnernormal[3*i] = arep3[3*j]; 
    g_fmminnernormal[3*i + 1] = arep3[3*j + 1]; 
    g_fmminnernormal[3*i + 2] = arep3[3*j + 2];
  }
  
  cilk_for (int i = 0; i < g_ntargets; i++) {
    int j = g_maptar[i]; 
    g_fmmtargets[3*i] = node[3*j]; 
    g_fmmtargets[3*i + 1] = node[3*j + 1]; 
    g_fmmtargets[3*i + 2] = node[3*j + 2]; 
    g_fmmouternormal[3*i] = nodevn[3*j]; 
    g_fmmouternormal[3*i + 1] = nodevn[3*j + 1]; 
    g_fmmouternormal[3*i + 2] = nodevn[3*j + 2];
  }

  // Step 4: Based on the problem size, solve the linear system by
  // calling either GMRES or BCGStab from Sparskit. For a problem of
  // size n, the work space required is 8*n for BCGStab and
  // (n+3)*(m+2) + m*(m+1)/2 for GMRES, where m is the size of Krylov
  // space before it restarts. The work space is indexed using a
  // 4-byte integer, which can hold value up to 2^31 - 1. 

  // Using the current Sparskit, the program can handle problem size
  // up to 2^28 using BCGStab. The program prefers to use GMRES,
  // attempting to set the size of Krylov subspace to be 100. If this
  // request cannot be met, the program will choose BCGStab. 
  int ipar[16] = {0};
  double fpar[16] = {0};
  double *rwork; 
  long long unsigned int dynamic_iter_mem; 
  ipar[0] = 0; // To initialize the iterative solver. 
  ipar[1] = 0; // Right preconditioning. 
  ipar[2] = 2; // Stopping criteria. 
  fpar[0] = 1e-6; // Relative tolerance.
  fpar[1] = 1e-6; // Absolute tolerance.
  fpar[10] = 0.0; // Reset flops. 

  const int maxProblemSizeAllowed = 268435456; 
  if (nnodes*2 > maxProblemSizeAllowed) {
    ERRMSG("problem size too large to pass to Sparskit solvers");
    exit(-1);
  } else {
    if (nnodes*2*100 > maxProblemSizeAllowed) {
      // The program selects BCGStab solver. 
      KrylovSpaceSolver = bcgstab_; 
      rwork = (double *)calloc(16*nnodes, sizeof(double)); 
      dynamic_iter_mem = sizeof(double)*16*nnodes; 
    } else {
      // The program selects GMRES solver and the Krylov subspace can
      // be as large as 100. 
      KrylovSpaceSolver = gmres_; 
      ipar[4] = 100; 
      ipar[3] = (nnodes*2 + 3)*(ipar[4] + 2) + (ipar[4] + 1)*ipar[4]/2; 
      rwork = (double *)calloc(ipar[3], sizeof(double)); 
      dynamic_iter_mem = sizeof(double)*ipar[3]; 
    }
  } 

  double *uvec = (double *)calloc(nnodes*2, sizeof(double));
  if (rwork == 0 || uvec == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  int nnodes2 = nnodes*2, numiter = 0; 
  for (int i = 0; i < nnodes2; i++)
    fhiso[i] = rhs[i]; 

  long long unsigned dynamic_fmm_mem = sizeof(int)*(g_nsources + g_ntargets) + 
    sizeof(double)*(7*g_nsources + 7*g_ntargets) + 
    sizeof(dcomplex)*(g_pgsz*2*(1 + g_nsboxes) + g_pgsz*2*(1 + g_ntboxes) + 
		      g_lapnexpmax*(1 + g_nsboxes)*6 + 
		      g_yuknexpmax*(1 + g_nsboxes)*6); 

  // Print out program status. If ipar[3] = 0, then the program has
  // chosen BCGStab. Otherwise, the program has chosen GMRES.   
  printf("\n----------------------------------------------------------------------\n"
	 "AFMPB solver status:\n"
	 "%-50s%20s\n"
	 "%-50s%18.5eMB\n"
	 "%-33s%-17s%18.5eMB\n"
	 "%-50s%20.5e\n",
	 "... Iterative solver selected:", (ipar[3] ? "GMRES" : "BCGStab"), 
	 "... Dynamic memory allocated for FMM:", 1.0*dynamic_fmm_mem/1024/1024, 
	 "... Dynamic memory allocated for ",  (ipar[3] ? "GMRES:" : "BCGStab:"), 
	 1.0*dynamic_iter_mem/1024/1024, "... Iterative solver absolute tolerance:", 
	 fpar[1]);
  // Call the basic iterative solver with reverse communication. If
  // ipar[0] equals to 1, then another matrix-vector multiply is
  // performed. If ipar[0] > 1, other type of matrix-vector multiply
  // is requested by the GMRES subroutine, which is not supported
  // here. If ipar[0] becomes nonpositive, then the iterative solver
  // exits with a status specified in ipar[0] 
 // struct timeval tic, toc; 
 // double elapsed; 

  while (1) {
    KrylovSpaceSolver(&nnodes2, rhs, fhiso, ipar, fpar, rwork);
    
    if (ipar[0] == 1) {
      if (numiter) {
	printf("... Iteration %2d residual norm: %38.5e\n", numiter, fpar[4]);
      }
      numiter++;

      double *ptr = &rwork[ipar[7] - 1]; 
      for (int i = 0; i < nnodes2; i++) 
	uvec[i] = ptr[i]; 
      double *V = &rwork[ipar[8] - 1]; 
      
      gettimeofday(&tic, 0);
      MatrixVectorMultiply(uvec, arep, directlist, directA, directB, directC, directD, V);
      gettimeofday(&toc, 0);
      elapsed += (double)((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);

    } else if (ipar[0] >= 2) {
      ERRMSG("nonsupported matrix-vector multiply requested by the GMRES");
      exit(-1);
    } else {
      printf("... Iteration %2d residual norm: %38.5e\n"
	     "%-50s%20d\n"
	     "%-50s%20.5e\n"
	     "----------------------------------------------------------------------\n",
	     numiter, fpar[4], 	     
	     "... Iterative solver exits with status:", ipar[0], 
	     "... Total matrix-vector multiply time:", elapsed);
      break;
    }
  }
  
  // Step 5: Cleanup all the memory space.   
  free(rwork);
  free(uvec);
  CleanDirectList13(directlist);
  for (int i = 0; i < g_ntargets; i++) {
    free(directA[i]);
    free(directB[i]);
    free(directC[i]);
    free(directD[i]);
  }
  free(directA);
  free(directB);
  free(directC);
  free(directD);

  DestroyGraph();
  LapFMMClean();
  YukFMMClean();
  FMMClean();
}

void SetRightHandSide(const int nsources, const double * const sources, 
		      const double * const charges, const int ntargets, 
		      const double * const targets, const double * const normal, 
		      double * const gn, double * const dgn)
{
  // Function SetRightHandSide() computes the right-hand side of the
  // linear system via a single-layer Laplace kernel FMM. 
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

  cilk_for (int i = 0; i < nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmsources[3*i] = sources[3*j]; 
    g_fmmsources[3*i + 1] = sources[3*j + 1];
    g_fmmsources[3*i + 2] = sources[3*j + 2]; 
    g_fmmcharges[i] = charges[j];
  }

  cilk_for (int i = 0; i < ntargets; i++) {
    int j = g_maptar[i];
    g_fmmtargets[3*i] = targets[3*j];
    g_fmmtargets[3*i + 1] = targets[3*j + 1];
    g_fmmtargets[3*i + 2] = targets[3*j + 2];
    g_fmmouternormal[3*i] = normal[3*j]; 
    g_fmmouternormal[3*i + 1] = normal[3*j + 1]; 
    g_fmmouternormal[3*i + 2] = normal[3*j + 2]; 
  }

  AdapFMMCompute(1); 

  cilk_for (int i = 0; i < ntargets; i++) {
    int j = g_maptar[i]; 
    gn[j] = g_fmmpotential[i]; 
    dgn[j] = g_fmmfield[3*i]*g_fmmouternormal[3*i] + 
      g_fmmfield[3*i + 1]*g_fmmouternormal[3*i + 1] + 
      g_fmmfield[3*i + 2]*g_fmmouternormal[3*i + 2];
  }

  DestroyGraph();
  LapFMMClean();
  FMMClean();
}

void SetSelfMatrix(int **directlist, const double * const targets, 
		   const double * const normal, const double * const sources, 
		   const double * const arep, const double * const arepn, 
		   const double * const arep3, const double * const px, 
		   const int * const ippt, const double * const xnrm, 
		   const double * const xwt, double **A, double **B, 
		   double **C, double **D)		   
{
  double gfactor1 = 0.5;
  double gfactor2 = 1.0 - gfactor1;
  double gf1 = gfactor1/dei + gfactor2; 
  double gf2 = gfactor1 + gfactor2/dei; 
 // printf("cut1=%f,cut2=%f,g_nthboxes=%d\n",cut1,cut2, g_ntboxes);

  cilk_for (int ibox = 1; ibox <= g_ntboxes; ibox++) {
    if (directlist[ibox] != 0 ) {
      int pid1 = g_tboxes[ibox].addr;
      int pid2 = pid1 + g_tboxes[ibox].npts - 1;
      int nlist = directlist[ibox][0]; 
      int nelems = directlist[ibox][1]; 

      // Allocate memory to hold the A, B, C, and D lists for each
      // point contained in the box. 
      for (int i = pid1; i <= pid2; i++) {
	A[i] = (double *)calloc(nlist, sizeof(double));
	B[i] = (double *)calloc(nlist, sizeof(double));
	C[i] = (double *)calloc(nlist, sizeof(double));
	D[i] = (double *)calloc(nlist, sizeof(double));
	if (A[i] == 0 || B[i] == 0 || C[i] == 0 || D[i] == 0) {
	  ERRMSG("memory allocation failure");
	  exit(-1);
	}
      }

      // Fill out the contents of the A, B, C, and D lists.
      cilk_for (int i = pid1; i <= pid2; i++) {
	int j = g_maptar[i]; 
	double rp0[3], vnx0[3]; 
	rp0[0] = targets[3*j];
	rp0[1] = targets[3*j + 1];
	rp0[2] = targets[3*j + 2];
	vnx0[0] = normal[3*j];
	vnx0[1] = normal[3*j + 1];
	vnx0[2] = normal[3*j + 2];
	int it = 0; 
	for (int j = 1; j <= nelems; j++) {
	  int pid3 = directlist[ibox][2*j];
	  int pid4 = directlist[ibox][2*j + 1]; 
	  for (int k = pid3; k <= pid4; k++) {
	    int ell = g_mapsrc[k]; 
	    double rp1[3], vnx[3]; 
	    rp1[0] = sources[3*ell];
	    rp1[1] = sources[3*ell + 1];
	    rp1[2] = sources[3*ell + 2];
	    vnx[0] = arep3[3*ell];
	    vnx[1] = arep3[3*ell + 1];
	    vnx[2] = arep3[3*ell + 2];
	    double dx = rp0[0] - rp1[0];
	    double dy = rp0[1] - rp1[1];
	    double dz = rp0[2] - rp1[2];
	    double dst = sqrt(dx*dx + dy*dy + dz*dz);
	    double ah, bh, ch, dh;
	    if (dst == 0) {
	      A[i][it] = 0.0;
	      B[i][it] = gf1; 
	      C[i][it] = gf2; 
	      D[i][it] = 0.0;
	    } else if (dst < cut2) {
		    CloseCoeff(rp0, ell, vnx0, px, ippt, xnrm, xwt, &ah, &bh, &ch, &dh);
	      A[i][it] = -ah;
	      B[i][it] = bh;
	      C[i][it] = -ch;
	      D[i][it] = dh;
	    } else if (dst < cut1) {
	      //CloseCoeff(rp0, ell, vnx0, px, ippt, xnrm, xwt, &ah, &bh, &ch, &dh);
	      A[i][it] = 0.0;
	      B[i][it] = 0.0;
	      C[i][it] = 0.0;
	      D[i][it] = 0.0;
	    } else {
	      NSingularCoeff(rp0, rp1, arep[ell], arepn[ell], vnx, vnx0, 
			     &ah, &bh, &ch, &dh);
	      A[i][it] = -ah;
	      B[i][it] = bh;
	      C[i][it] = -ch;
	      D[i][it] = dh;
	    }
	    it++;
	  }
	}
      }
    }
  }
}

void MatrixVectorMultiply(double * const U, const double * const arep, 
			  int **directlist, double **A, double **B, double **C, 
			  double **D, double * const V)
{
  // Function MatrixVectorMultiply() computes the matrix-vector
  // multiply required by the GMRES iterative solver. U is the input
  // vector and V is the output result. The length of U is
  // 2*g_nsources, and the length of V is 2*g_ntargets. 

  // The matrix-vector multiply takes two steps to complete. In the
  // first step, near-field contribution is processed with direct
  // interaction using the coefficients stored in A, B, C, and D lists
  // for each target point. In the second step, far-field contribution
  // is processed via four FMM calls. 

  // Variable arep, arep3, vn are the assigned area, inner normal
  // derivative, and outer normal derivate of each node. 

  const double scalar1 = g_beta*2.0/M_PI; 
  const double scalar2 = 1.0/dei*g_beta*2.0/M_PI; 

  for (int i = 0; i < 2*g_ntargets; i++)
    V[i] = 0.0; 

  // Step 1: Perform near-field computation. 
  cilk_for (int ibox = 1; ibox <= g_ntboxes; ibox++) {
    if (directlist[ibox] != 0) {
      int pid1 = g_tboxes[ibox].addr; 
      int pid2 = pid1 + g_tboxes[ibox].npts - 1; 
      int nelems = directlist[ibox][1]; 
      cilk_for (int i = pid1; i <= pid2; i++) {
	int t = g_maptar[i]; 
	double vt = 0, vtn = 0; 
	int it = 0; 
	for (int j = 1; j <= nelems; j++) {
	  int pid3 = directlist[ibox][2*j]; 
	  int pid4 = directlist[ibox][2*j + 1]; 
	    for (int k = pid3; k <= pid4; k++) {
	    int ell = g_mapsrc[k]; 
	    vt += A[i][it]*U[ell + g_nsources] + B[i][it]*U[ell]; 
	    vtn += C[i][it]*U[ell + g_nsources] + D[i][it]*U[ell];
	    it++;
	  }
	}
	V[t] += vt;
	V[t + g_ntargets] += vtn;
      }
    }
  }
  // Step 2: Perform far-field computation. It consists of four
  // successive FMM calls, in the order of single-layer Laplace
  // potential, double-layer Laplace potential, single-layer Yukawa
  // potential, and double-layer Yukawa potential. 
  double * const ff = &U[0]; 
  double * const fh = &U[g_nsources]; 
  const double fac1 = 0.79577472e-1; 
  
  for (int i = 0; i < g_nsources; i++) {
    double temp = fac1*arep[i]; 
    // Variable ff and fh are the charges in the context of FMM for
    // double and single layer potentials, respectively. 
    ff[i] *= temp; 
    fh[i] *= temp; 
  }

  // Perform single-layer Laplace potential computation. 
  cilk_for (int i = 0; i < g_nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmcharges[i] = fh[j]; 
  }

  SourceToMultipole = sLapSourceToMultipole; 
  MultipoleToMultipole = LapMultipoleToMultipole; 
  MultipoleToExponential = LapMultipoleToExponential; 
  ExponentialToLocal = LapExponentialToLocal;
  LocalToLocal = LapLocalToLocal;
  LocalToTarget = LapLocalToTarget; 
  GreenFunction = sLapGreenFunction;   

  AdapFMMCompute(0); 

  cilk_for (int i = 0; i < g_ntargets; i++) {
    int j = g_maptar[i]; 
    V[j] -= g_fmmpotential[i]; 
    V[j + g_ntargets] += g_fmmfield[3*i]*g_fmmouternormal[3*i] + 
      g_fmmfield[3*i + 1]*g_fmmouternormal[3*i + 1] + 
      g_fmmfield[3*i + 2]*g_fmmouternormal[3*i + 2]; 
  }

  memset(g_fmmpotential, 0, g_ntargets*sizeof(double));
  memset(g_fmmfield, 0, 3*g_ntargets*sizeof(double));
  memset(lap_multipole, 0, (1 + g_nsboxes)*g_pgsz*sizeof(dcomplex)); 
  memset(lap_local, 0, (1 + g_ntboxes)*g_pgsz*sizeof(dcomplex));
  memset(lap_expu, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expd, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expn, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_exps, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expe, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expw, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  
  // Perform double-layer Laplace potential computation. 
   cilk_for (int i = 0; i < g_nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmcharges[i] = ff[j];
  }

  SourceToMultipole = dLapSourceToMultipole; 
  GreenFunction = dLapGreenFunction; 
  
  AdapFMMCompute(0); 

  for (int i = 0; i < g_ntargets; i++) {
    int j = g_maptar[i]; 
    V[j] -= g_fmmpotential[i]/dei; 
    V[j + g_ntargets] -= (g_fmmfield[3*i]*g_fmmouternormal[3*i] + 
			  g_fmmfield[3*i + 1]*g_fmmouternormal[3*i + 1] + 
			  g_fmmfield[3*i + 2]*g_fmmouternormal[3*i + 2])/dei; 
  }

  memset(g_fmmpotential, 0, g_ntargets*sizeof(double));
  memset(g_fmmfield, 0, 3*g_ntargets*sizeof(double));
  memset(lap_multipole, 0, (1 + g_nsboxes)*g_pgsz*sizeof(dcomplex)); 
  memset(lap_local, 0, (1 + g_ntboxes)*g_pgsz*sizeof(dcomplex));
  memset(lap_expu, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expd, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expn, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_exps, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expe, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  memset(lap_expw, 0, (1 + g_nsboxes)*g_lapnexpmax*sizeof(dcomplex));
  
  // Perform single-layer Yukawa potential computation.   
  cilk_for (int i = 0; i < g_nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmcharges[i] = fh[j];
  }

  SourceToMultipole = sYukSourceToMultipole; 
  MultipoleToMultipole = YukMultipoleToMultipole; 
  MultipoleToExponential = YukMultipoleToExponential; 
  ExponentialToLocal = YukExponentialToLocal;
  LocalToLocal = YukLocalToLocal;
  LocalToTarget = YukLocalToTarget; 
  GreenFunction = sYukGreenFunction;   

  AdapFMMCompute(0);

  cilk_for (int i = 0; i < g_ntargets; i++) {
    int j = g_maptar[i]; 
    V[j] += g_fmmpotential[i]*scalar1;
    V[j + g_ntargets] += (g_fmmfield[3*i]*g_fmmouternormal[3*i] + 
			  g_fmmfield[3*i + 1]*g_fmmouternormal[3*i + 1] + 
			  g_fmmfield[3*i + 2]*g_fmmouternormal[3*i + 2])*scalar2;
  }

  memset(g_fmmpotential, 0, g_ntargets*sizeof(double));
  memset(g_fmmfield, 0, 3*g_ntargets*sizeof(double));
  memset(yuk_multipole, 0, (1 + g_nsboxes)*g_pgsz*sizeof(dcomplex));
  memset(yuk_local, 0, (1 + g_ntboxes)*g_pgsz*sizeof(dcomplex));
  memset(yuk_expu, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expd, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expn, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_exps, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expe, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expw, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  
  // Perform double-layer Yukawa potential computation. 
  cilk_for (int i = 0; i < g_nsources; i++) {
    int j = g_mapsrc[i]; 
    g_fmmcharges[i] = ff[j];
  }

  SourceToMultipole = dYukSourceToMultipole;
  GreenFunction = dYukGreenFunction; 

  AdapFMMCompute(0);
  
  cilk_for (int i = 0; i < g_ntargets; i++) {
    int j = g_maptar[i]; 
    V[j] += g_fmmpotential[i]*scalar1; 
    V[j + g_ntargets] -= (g_fmmfield[3*i]*g_fmmouternormal[3*i] + 
			  g_fmmfield[3*i + 1]*g_fmmouternormal[3*i + 1] + 
			  g_fmmfield[3*i + 2]*g_fmmouternormal[3*i + 2])*scalar2;
  }

  memset(g_fmmpotential, 0, g_ntargets*sizeof(double));
  memset(g_fmmfield, 0, 3*g_ntargets*sizeof(double));
  memset(yuk_multipole, 0, (1 + g_nsboxes)*g_pgsz*sizeof(dcomplex));
  memset(yuk_local, 0, (1 + g_ntboxes)*g_pgsz*sizeof(dcomplex));
  memset(yuk_expu, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expd, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expn, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_exps, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expe, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
  memset(yuk_expw, 0, (1 + g_nsboxes)*g_yuknexpmax*sizeof(dcomplex));
}
