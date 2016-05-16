#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "afmpb.h"
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

double wctime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec+1E-6*tv.tv_usec);
}
//enum BOOL {FALSE = 0,TRUE = !FALSE};

int main(int argc, char **argv)
{
  int ncharges, npts, nelems, *node, meshfmt, *ippt; 

  char input[80], output[80], surfp[80], vmesh[80], volp[80], fpqr[80], fmesh[80];

  double *xq, *q, *r, *p, *vnp, *ff, *rh, *gn, *dgn, 
    *vn, *arel, *arep, *arepn, *arep3, *px, *xnrm, *xwt;

  printf("\n----------------------------------------------------------------------\n"
	 "*          Adaptive Fast Multipole Poisson Boltzmann Solver          *\n"
	 "----------------------------------------------------------------------\n\n");
  
  if (ParseCommandLine(argc, argv, input, output, surfp, vmesh, volp)) {
    ERRMSG("ParseCommandLine() failure");
    return -1;
  }

  if (SetParameters(input, &meshfmt, fpqr, fmesh)) {
    ERRMSG("SetParameters() failure");
    return -1;
  }

  ncharges = ProcessPQRFile(fpqr); 
  if (ncharges <= 0) {
    ERRMSG("number of charges in the PQR file is nonpositive");
    return -1;
  }

  xq = (double *)calloc(3*ncharges, sizeof(double)); 
  q = (double *)calloc(ncharges, sizeof(double));
  r = (double *)calloc(ncharges, sizeof(double)); 
  if (xq == NULL || q == NULL || r == NULL) {
    ERRMSG("memory allocation failure");
    return -1;
  }

  ReadPQRFile(ncharges, xq, q, r); 

  ReadMeshFile(fmesh, &p, &node, &vnp, meshfmt, &npts, &nelems); 

  RemoveIsolatedNodes(meshfmt, p, node, vnp, nelems, &npts);

  printf("----------------------------------------------------------------------\n"
	 "Problem parameters:\n"
	 "%-50s%20d\n"
	 "%-50s%20d\n"
	 "%-50s%20d\n"
	 "%-50s%20d\n",
	 "... Number of molecules:", nmol, 
	 "... Number of point charges inside the molecule:", ncharges, 
	 "... Number of elements on molecular boundary:", nelems,
	 "... Number of nodes on molecular boundary:", npts); 

  // Allocate memories. 
  vn    = (double *)calloc(3*nelems, sizeof(double));
  arel  = (double *)calloc(nelems, sizeof(double));
  arep  = (double *)calloc(npts, sizeof(double));
  arepn = (double *)calloc(npts, sizeof(double));
  arep3 = (double *)calloc(3*npts, sizeof(double));
  px    = (double *)calloc(18*nelems, sizeof(double));
  xnrm  = (double *)calloc(18*nelems, sizeof(double));
  xwt   = (double *)calloc(6*nelems, sizeof(double));
  ippt  = (int    *)calloc((npts+1), sizeof(int   ));
  gn    = (double *)calloc(7*nelems, sizeof(double));
  dgn   = (double *)calloc(7*nelems, sizeof(double));
  ff    = (double *)calloc(2*npts, sizeof(double));
  rh    = (double *)calloc(2*npts, sizeof(double));
  if (vn == NULL || arel == NULL || arep == NULL || arepn == NULL ||
      arep3 == NULL || px == NULL || xnrm == NULL || xwt == NULL ||
      ippt == NULL || gn == NULL || dgn == NULL || ff == NULL || rh == NULL) {
    ERRMSG("memory allocation failure");
    return -1;
  }

  ProcessElementGeometry(meshfmt, npts, nelems, node, p, vn, arel, 
			 arep, arepn, arep3, vnp, ippt, px, xnrm, xwt); 

  // Construct potential and its normal derivative for future energy
  // computation. 
  SetGNDGN(p, vn, ncharges, xq, q, node, nelems, gn, dgn); 
  printf("%-50s%20.5e\n"
	 "%-50s%20.5e\n"
	 "%-50s%20d\n"
	 "----------------------------------------------------------------------\n\n",
	 "... Kap:", kap, "... Screening length:", 1.0/kap, 
	 "... Number of cilk workers:", __cilkrts_get_nworkers());
  // Timing
  double start = wctime();

  // Solve linearized Poisson Boltzmann equation with GMRES and FMM. 
  SolvePoissonBoltzmann(npts, p, vnp, arep, arepn, arep3, ncharges, q, 
			xq, px, ippt, xnrm, xwt, ff, rh); 
  double end = wctime();
  double t_elapsed = end - start;
  printf("The elapsed time is %f seconds.\n", t_elapsed);

  for (int i = 0; i < 2*npts; i++)
    ff[i] *= unitfactor; 

  WritePotential(surfp, npts, nelems, node, p, vnp, ff); 

  // Compute solvation energy. 
  printf("\n----------------------------------------------------------------------\n"
	 "%-50s%20.5e\n"
	 "----------------------------------------------------------------------\n", 
	 "Solvation energy:",  
	 ComputeSolvationEnergy(nelems, gn, dgn, ff, npts, arel, node));
   
 // Cleanup the memory space
  free(rh); 
  free(ff);
  free(dgn); 
  free(gn);
  free(ippt); 
  free(xwt);
  free(xnrm);
  free(px);
  free(vn);
  free(arep);
  free(arel); 
  free(arepn);
  free(node);
  free(vnp);  
  free(p);     
  free(arep3);
  free(xq);   
  free(q); 
  free(r);

  return 0;
}
