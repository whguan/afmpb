#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "adap_fmm.h"

void AdapFMMCompute(const int mode)
{
  // Function AdapFMMCompute() traverses the adaptive FMM graph and
  // completes related computation. If mode equals to 1, lists 1 and 3
  // for leaf target boxes are processed. Otherwise, they are
  // skipped. 
  AggregateSweep(1); 
  DisAggregateSweep(1); 
  if (mode) {
    cilk_for (int ibox = 1; ibox <= g_ntboxes; ibox++) {
      if (g_tboxes[ibox].nchild == 0)
	ProcessList13(ibox);
    }
  }
}

void AggregateSweep(const int ibox)
{
  // Function AggregateSweep() generates the multipole and six directional
  // exponential expansions for each source box of the subtree rooted
  // at ibox. 
  if (g_sboxes[ibox].nchild == 0) {
    SourceToMultipole(ibox);
  } else {
    cilk_for (int i = 0; i < 8; i++) {
      int child = g_sboxes[ibox].child[i]; 
      if (child)
	AggregateSweep(child);
    }
    MultipoleToMultipole(ibox);
  }

  // This version of the FMM considers free-space boundary
  // condition. This means, the exponential expansions of box 1 is not
  // used in any computation. Hence, they are not generated. 
  if (ibox > 1)
    MultipoleToExponential(ibox); 
} 


void DisAggregateSweep(const int ibox)
{
  // Function DisAggregateSweep() performs a top-down
  // sweep. Afterwards, each target box of the subtree rooted at ibox
  // has received the far-field potential result. 
  ProcessList4(ibox); 
  if (g_tboxes[ibox].nchild == 0) {
    LocalToTarget(ibox);
  } else {
    ExponentialToLocal(ibox); 
    LocalToLocal(ibox); 
    cilk_for (int i = 0; i < 8; i++) {
      int child = g_tboxes[ibox].child[i];
      if (child)
	 DisAggregateSweep(child);
	
    }
  }
}

void ProcessList13(const int ibox)
{
  if (g_tboxes[ibox].list3) {
    int nlist = g_tboxes[ibox].list3[0]; 
    for (int it = 1; it <= nlist; it++)
      DirectEvaluation(ibox, g_tboxes[ibox].list3[it]);
  }

  if (g_tboxes[ibox].list1) {
    int nlist = g_tboxes[ibox].list1[0];
    for (int it = 1; it <= nlist; it++)
      DirectEvaluation(ibox, g_tboxes[ibox].list1[it]);
  }
}

void ProcessList4(const int ibox)
{
  if (g_tboxes[ibox].list4) {
    int nlist = g_tboxes[ibox].list4[0];
    for (int it = 1; it <= nlist; it++)
      DirectEvaluation(ibox, g_tboxes[ibox].list4[it]);
  }
}

void DirectEvaluation(const int tbox, const int sbox)
{
  // Function DirectInteraction computes the pairwise interaction
  // between a source box and a target box. 

  int start1 = g_tboxes[tbox].addr;
  int num1 = g_tboxes[tbox].npts; 
  int end1 = start1 + num1 - 1;
  int start2 = g_sboxes[sbox].addr;
  int num2 = g_sboxes[sbox].npts; 
  int end2 = start2 + num2 - 1;

  for (int it1 = start1; it1 <= end1; it1++) {
    double pot = 0, fx = 0, fy = 0, fz = 0; 
    int ptr1 = 3*it1; 
    for (int it2 = start2; it2 <= end2; it2++) {
      int ptr2 = 3*it2; 
      GreenFunction(&g_fmmtargets[ptr1], &g_fmmsources[ptr2], g_fmmcharges[it2], 
		    &g_fmminnernormal[ptr2], &pot, &fx, &fy, &fz); 
    }
    g_fmmpotential[it1] += pot; 
    g_fmmfield[ptr1] += fx; 
    g_fmmfield[ptr1 + 1] += fy; 
    g_fmmfield[ptr1 + 2] += fz;
  }
}
