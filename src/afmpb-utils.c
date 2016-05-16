#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "afmpb.h"

int ParseCommandLine(int argc, char **argv, char * const input, char * const output, 
		     char * const surfp, char * const vmesh, char * const volp)
{
  // Setup default parameters
  strcpy(input, "input.txt");
  strcpy(output, "output.txt");
  strcpy(surfp, "surfp.dat");
  strcpy(vmesh, "vmesh.dat");
  strcpy(volp, "volp.dat");

  // Parse available command line options if argc > 1. 
  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
	switch (argv[i][1]) {
	case 'i':
	  strcpy(input, argv[++i]);
	  break;       	  
	case 'o':
	  strcpy(output, argv[++i]);
	  break;		  
	case 's':
	  strcpy(surfp, argv[++i]);
	  break;		  
	case 'm':
	  strcpy(vmesh, argv[++i]);
	  break;		  
	case 'p':
	  strcpy(volp, argv[++i]);
	  break;		  
	default:
	  break;
	}
      }
    }   
  }
  return 0; 
}

int SetParameters(const char * const input, int *meshfmt, char * const fpqr, 
		  char * const meshf)
{
  FILE *pf; 
  pf = fopen(input, "r"); 
  if (pf == NULL) {
    ERRMSG("cannot open input file");
    return -1;
  }

  // The input file consists of 10 lines. A line starting with the `#'
  // sign stands for a comment line. 
  char temp[100]; 

  // Skip comment line and get number of molecules. 
  fgets(temp, 100, pf);  
  fgets(temp, 100, pf); 
  sscanf(temp, "%d", &nmol); 

  // Skip comment line and get values for di, dei, ion concentration,
  // and temperature. 
  fgets(temp, 100, pf);
  fgets(temp, 100, pf);
  sscanf(temp, "%lf%lf%lf%lf", &di, &de, &concentration, &temperature);

  // Skip comment line and get mesh format. 
  fgets(temp, 100, pf);
  fgets(temp, 100, pf);
  sscanf(temp, "%d", meshfmt);

  // Skip comment line and get output keys. 
  fgets(temp, 100, pf);
  fgets(temp, 100, pf);
  sscanf(temp, "%d%d%d%d%d", &ipotw, &iforce, &iintere, &iselfe, &ipotdx);

  // Skip comment line and get pqr and mesh file names. 
  fgets(temp, 100, pf);
  fgets(temp, 100, pf);
  sscanf(temp, "%s%s", fpqr, meshf);

  fclose(pf);

  concentration = MAX(concentration, 1e-10); 
  dei = de/di; 
  kap = sqrt(2.528639884*concentration/de/temperature); 

  if (*meshfmt == 2) {
    cut1 = 0.4;
    cut2 = 0.4;
  } else {
    cut1 = 0.0;
    cut2 = 0.3;
  } 

  sigma = 0.001; 

  if (FMMACCURACY == 3) {
    g_pterms = 9;
    g_nlambs = 9;
    g_pgsz = pow(g_pterms + 1, 2);
  } else if (FMMACCURACY == 6) {
    g_pterms = 18;
    g_nlambs = 18;
    g_pgsz = pow(g_pterms + 1, 2);
  } else {
    ERRMSG("wrong FMM accuracy");
    return -1;
  }

  return 0;
}

int ProcessPQRFile(const char * const fname)
{
  // Function ProcessPQRFile parses the input PQR file and extracts
  // all the lines starting with either ATOM or HETATM into a
  // processed file called processed.pqr. The number of point charges
  // inside the molecule is returned by counting the number of lines
  // in the processed PQR file. 
  FILE *pf; 
  char command[200]; 
  int ncharge; 

  sprintf(command, "grep \'ATOM\\|HETATM\' %s > %s; wc -l < %s", 
	  fname, "processed.pqr", "processed.pqr"); 

  pf = popen(command, "r"); 

  if (pf == NULL) {
    ERRMSG("cannot process pqr file");
    return -1;
  } else {
    while (1) {
      char *line; 
      char buf[1000];
      line = fgets(buf, sizeof buf, pf); 
      if (line == NULL) {
	break;
      } else {
	// We only expect 1 line output from the command. 
	sscanf(line, "%d", &ncharge); 
      }
    }
    pclose(pf); 
  }

  return ncharge; 
}

void ReadPQRFile(const int ncharges, double * const xq, 
		 double * const q, double * const r)
{
  FILE *pf; 
  pf = fopen("processed.pqr", "r"); 
  if (pf == NULL) {
    ERRMSG("cannot open processed PQR file");
    exit(-1);
  }

  char temp[100], atom[10], atom1[10], res[10]; 
  for (int k = 0; k < ncharges; k++) {
    int j = k*3; 
    int num, nres; 
    fgets(temp, 100, pf); 
    sscanf(temp, "%s%d%s%s%d%lf%lf%lf%lf%lf", atom, &num, atom1, res, &nres, 
	   xq + j, xq + j + 1, xq + j + 2, q + k, r + k);
  }

  fclose(pf); 
}

void ReadMeshFile(const char * const fname, double **p_, int **node_, double **vnp_, 
		  const int meshfmt, int *npts, int *nelems)
{
  FILE *pf; 
  pf = fopen(fname, "r"); 
  if (pf == NULL) {
    ERRMSG("cannot open mesh file");
    exit(-1);
  }

  // Read the first line of the mesh file, and get the number of
  // points and the number of elements. If the line starts with OFF,
  // then it is a comment line. Skip it and then get the desired
  // information from the next line. 
  char temp[100]; 
  while (1) {
    fgets(temp, 100, pf); 
    char *comment; 
    comment = strstr(temp, "OFF"); 
    if (comment == NULL) {
      // The word "OFF" does not appear in the buffer and we retrieve
      // values for npts and nelems. 
      sscanf(temp, "%d%d", npts, nelems); 
      break;
    }
  }

  double *p   = (double *)calloc(3*(*npts), sizeof(double)); 
  double *vnp = (double *)calloc(3*(*npts), sizeof(double)); 
  int *node   = (int *)calloc(3*(*nelems), sizeof(int)); 
  if (p == NULL || vnp == NULL || node == NULL) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  if (meshfmt == 2) {
    for (int i = 0; i < *npts; i++) {
      int j = 3*i;
      fgets(temp, 100, pf);
      sscanf(temp, "%lf%lf%lf%lf%lf%lf", p + j, p + j + 1, p + j + 2, 
	     vnp + j, vnp + j + 1, vnp + j + 2);	
    }
  }
  
  if(meshfmt == 5 || meshfmt == 3) {
    for (int i = 0; i < *npts; i++) {
      int j = 3*i;
      fgets(temp, 100, pf);
      sscanf(temp, "%lf%lf%lf", p + j, p + j + 1, p + j + 2);
    }
  }

  if(meshfmt == 5) {
    for (int i = 0; i < *nelems; i++) {
      int j = 3*i;
      int k; 
      fgets(temp, 100, pf);
      sscanf(temp,"%d%d%d%d", &k, node + j, node + j + 1, node + j + 2);
    }
  } else {
    for (int i = 0; i < *nelems; i++) {
      int j = 3*i;
      fgets(temp, 100, pf);
      sscanf(temp,"%d%d%d", node + j, node + j + 1, node + j + 2);
    }
  }

  if(meshfmt == 5) {
    for (int i = 0; i < *nelems; i++) {
      int j = 3*i;
      node[j]++; 
      node[j + 1]++; 
      node[j + 2]++;
    }
  }
  
  fclose(pf);
  *p_ = p;
  *node_ = node; 
  *vnp_ = vnp; 
}

void RemoveIsolatedNodes(const int meshfmt, double * const p, int * const node, 
			 double * const vnp, const int nelems, int *npts)
{
  int *isolated = (int *)calloc(*npts + 1, sizeof(int)); 
  int *nnbr = (int *)calloc(*npts + 1, sizeof(int)); 
  if (isolated == 0 || nnbr == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int k = 0; k < nelems; k++) {
    int i1 = node[3*k]; 
    int i2 = node[3*k + 1]; 
    int i3 = node[3*k + 2]; 
    nnbr[i1]++;
    nnbr[i2]++;
    nnbr[i3]++;
  }
  
  int count = 0; 

  for (int k = 1; k <= *npts; k++) {
    if (nnbr[k] == 0) {
      isolated[count++] = k;
    }
  }

  isolated[count] = *npts + 1;
  if (count) {
    for (int k = 0; k < 3*nelems; k++) {
      int m0 = count - 1;
      while (m0 >=0 && node[k] < isolated[m0])
	m0--;
      node[k] -= m0 + 1;
    }

    for (int m0 = 0; m0 < count; m0++) {
      for (int j = isolated[m0] + 1; j < isolated[m0 + 1]; j++) {
	int u0 = 3*(j - m0 - 2);
	int u1 = 3*(j - 1);
	p[u0] = p[u1]; 
	p[u0 + 1] = p[u1 + 1]; 
	p[u0 + 2] = p[u1 + 2]; 
	if (meshfmt == 2) {
	  vnp[u0] = vnp[u1]; 
	  vnp[u0 + 1] = vnp[u1 + 1]; 
	  vnp[u0 + 2] = vnp[u1 + 2]; 
	}
      }
    }
  }

  free(isolated);
  free(nnbr); 
  *npts -= count; 	

}

void ProcessElementGeometry(const int meshfmt, const int npts, const int nelem, 
			    int *node, double *p, double *vn, double *arel, 
			    double *arep, double *arepn, double *arep3, 
			    double *vnp, int *ippt, double *px, double *xnrm, double *xwt)
{
  double rp1[3], rp2[3], rp3[3], rpm[3]; 
  double area = 0.0;
  double vol = 0.0;
  int *nne = (int *)calloc(npts + 1, sizeof(int));
  int *idfcl = (int *)calloc(npts + 1, sizeof(int)); 
  if (nne == NULL || idfcl == NULL) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  if (meshfmt != 2) {
    for (int i = 0; i < npts; i++) {
      int j = 3*i; 
      vnp[j] = 0.0;
      vnp[j + 1] = 0.0; 
      vnp[j + 2] = 0.0;
    }
  }

  for (int i = 0; i <= npts; i++) 
    nne[i] = 0; 

  for (int k = 0; k < nelem; k++) {
    int i1 = node[3*k] - 1; 
    int i2 = node[3*k + 1] - 1;
    int i3 = node[3*k + 2] - 1;
    double x1 = p[3*i1]; 
    double y1 = p[3*i1 + 1]; 
    double z1 = p[3*i1 + 2]; 
    double x2 = p[3*i2]; 
    double y2 = p[3*i2 + 1]; 
    double z2 = p[3*i2 + 2]; 
    double x3 = p[3*i3]; 
    double y3 = p[3*i3 + 1]; 
    double z3 = p[3*i3 + 2]; 
    double x21 = x2 - x1;
    double y21 = y2 - y1; 
    double z21 = z2 - z1;
    double x31 = x3 - x1; 
    double y31 = y3 - y1; 
    double z31 = z3 - z1;
    double x32 = x3 - x2; 
    double y32 = y3 - y2; 
    double z32 = z3 - z2; 
    int u0 = 3*k;
    vn[u0] = y21*z31 - y31*z21; 
    vn[u0 + 1] = z21*x31 - z31*x21; 
    vn[u0 + 2] = x21*y31 - x31*y21;       
    vol += ((x1 + x2 + x3)*vn[u0] + (y1 + y2 + y3)*vn[u0 + 1] +
	    (z1 + z2 + z3)*vn[u0 + 2])/18.0; 
    double hs = sqrt(vn[u0]*vn[u0] + vn[u0 + 1]*vn[u0 + 1] + vn[u0 + 2]*vn[u0 + 2]);
    arel[k] = hs/2.0; 

    if (hs > 0.00000001) {
      vn[u0] /= hs; 
      vn[u0 + 1] /= hs; 
      vn[u0 + 2] /= hs; 

      if (meshfmt != 2) {
	double b = sqrt(x21*x21 + y21*y21 + z21*z21);
	double a = sqrt(x31*x31 + y31*y31 + z31*z31);
	double c = sqrt(x32*x32 + y32*y32 + z32*z32);	  
	vnp[3*i1] += vn[u0]/a/b;
	vnp[3*i1 + 1] += vn[u0 + 1]/a/b;
	vnp[3*i1 + 2] += vn[u0 + 2]/a/b;
    
	vnp[3*i2] += vn[u0]/c/b;
	vnp[3*i2 + 1] += vn[u0 + 1]/c/b;
	vnp[3*i2 + 2] += vn[u0 + 2]/c/b;

	vnp[3*i3] += vn[u0]/a/c;
	vnp[3*i3 + 1] += vn[u0 + 1]/a/c;
	vnp[3*i3 + 2] += vn[u0 + 2]/a/c;
      }
    }      
    nne[i1 + 1]++;
    nne[i2 + 1]++;
    nne[i3 + 1]++;    
    area += arel[k]; 
  }

  if (vol < 0) {
    for (int k = 0; k < nelem; k++) {
      int u0 = 3*k;
      vn[u0] = -vn[u0]; 
      vn[u0 + 1] = -vn[u0 + 1]; 
      vn[u0 + 2] = -vn[u0 + 2];
    }
    
    vol = -vol;
  }

  printf("%-50s%20.3f\n"
	 "%-50s%20.3f\n", 
	 "... Total area of molecular boundary:", area,
	 "... Total volume of the molecule:", vol);

  idfcl[0] = 1; 
  for (int i = 1; i <= npts; i++) {
    idfcl[i] = idfcl[i - 1] + nne[i]; 
    nne[i] = 0; 
  }

  int *ne0 = (int *)calloc(idfcl[npts], sizeof(int));
  int *ne1 = (int *)calloc(idfcl[npts], sizeof(int));
  int *ne2 = (int *)calloc(idfcl[npts], sizeof(int)); 
  if (ne0 == NULL || ne1 == NULL || ne2 == NULL) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int k = 0; k < nelem; k++) {
    int i1 = node[3*k]; 
    int i2 = node[3*k + 1]; 
    int i3 = node[3*k + 2]; 

    nne[i1]++;
    nne[i2]++;
    nne[i3]++;

    int n1 = nne[i1] + idfcl[i1 - 1] - 1;
    int n2 = nne[i2] + idfcl[i2 - 1] - 1;
    int n3 = nne[i3] + idfcl[i3 - 1] - 1;

    ne0[n1] = k;
    ne0[n2] = k;
    ne0[n3] = k;
    
    ne1[n1] = 2;
    ne1[n2] = 3; 
    ne1[n3] = 1;
    
    ne2[n1] = 3;
    ne2[n2] = 1;
    ne2[n3] = 2;
  }

  int igenus = -(npts - nelem/2 - 2)/2;
  printf("%-50s%20d\n", "... The genus of surface (holes):", igenus); 

  if(meshfmt != 2) {
    for(int k = 0; k < npts; k++) {
      double vnpl = sqrt(vnp[3*k]*vnp[3*k] + vnp[3*k + 1]*vnp[3*k + 1] 
			 + vnp[3*k + 2]*vnp[3*k + 2]);
      int u0 = 3*k;
      vnp[u0] /= vnpl;  
      vnp[u0 + 1] /= vnpl;  
      vnp[u0 + 2] /= vnpl; 
    }
  }
 
  for(int k = 0; k < npts; k++) {
    int u0 = 3*k;      
    arep3[u0] = 0.0;
    arep3[u0 + 1] = 0.0; 
    arep3[u0 + 2] = 0.0; 
    arep[k] = 0.0;      
    for(int j = 1; j <= nne[k + 1]; j++) {  
      int k0 = ne0[idfcl[k] + j - 1];
      arep3[u0] += arel[k0]*vn[3*k0];
      arep3[u0 + 1] += arel[k0]*vn[3*k0 + 1];
      arep3[u0 + 2] += arel[k0]*vn[3*k0 + 2];      
      arep [k] += arel[k0];
    }
    
    arepn[k] = sqrt(arep3[u0]*arep3[u0] + arep3[u0 + 1]*arep3[u0 + 1] 
		    + arep3[u0 + 2]*arep3[u0 + 2]);
   
    arep3[u0] /= arepn[k];
    arep3[u0 + 1] /= arepn[k];
    arep3[u0 + 2] /= arepn[k];
    
    arep[k] /= 3.0;
    arepn[k] /= 3.0;
  }

  int num = 0;
  ippt[0] = 0;
  for (int i = 1; i <= npts; i++) { 
    int u0 = 3*i - 3;        
    rp1[0] = p[u0]; 
    rp1[1] = p[u0 + 1];
    rp1[2] = p[u0 + 2];    
    for (int j = idfcl[i - 1]; j < idfcl[i]; j++) {
      int np2 = ne1[j]; 
      int np3 = ne2[j];
      int np1 = ne0[j];      
      np2 = node[np2 - 1 + 3*np1];
      np3 = node[np3 - 1 + 3*np1];      
      np2 = 3*np2 - 3;
      np3 = 3*np3 - 3;
      rp2[0] = (rp1[0] + p[np2])/2; 
      rp2[1] = (rp1[1] + p[np2 + 1])/2;
      rp2[2] = (rp1[2] + p[np2 + 2])/2; 

      rp3[0] = (rp1[0] + p[np3])/2; 
      rp3[1] = (rp1[1] + p[np3 + 1])/2;
      rp3[2] = (rp1[2] + p[np3 + 2])/2;
      
      rpm[0] = (rp1[0] + p[np2] + p[np3])/3;
      rpm[1] = (rp1[1] + p[np2 + 1] + p[np3 + 1])/3;
      rpm[2] = (rp1[2] + p[np2 + 2] + p[np3 + 2])/3;
      
      int k = 3*num;      
      px[k] = (rp1[0] + rp2[0] + rp3[0])/3.0;
      px[k + 1] = (rp1[1] + rp2[1] + rp3[1])/3.0;
      px[k + 2] = (rp1[2] + rp2[2] + rp3[2])/3.0;
      
      xnrm[k] = vn[3*np1];
      xnrm[k + 1] = vn[3*np1 + 1];
      xnrm[k + 2] = vn[3*np1 + 2];
      
      xwt[num] = arel[np1]/4;
	  
      num++;
      k = 3*num;

      px[k] = (rpm[0] + rp2[0] + rp3[0])/3.0;
      px[k + 1] = (rpm[1] + rp2[1] + rp3[1])/3.0;
      px[k + 2] = (rpm[2] + rp2[2] + rp3[2])/3.0;
      
      xnrm[k] = vn[3*np1];
      xnrm[k + 1] = vn[3*np1 + 1];
      xnrm[k + 2] = vn[3*np1 + 2];
      
      xwt[num] = arel[np1]/12;
      num++;
    }
      
    ippt[i] = num;
  }

  free(ne0);
  free(ne1);
  free(ne2);
  free(idfcl);
  free(nne);

  return;
}

void WritePotential(const char * const fname, const int npts, 
		    const int nelems, const int * const node, 
		    const double * const p, const double * const vnp, 
		    const double * const f)
{
  FILE *pf; 
  pf = fopen(fname, "w");
  if (pf == NULL) {
    ERRMSG("cannot open surfp file");
    exit(-1);
  }

  for (int i = 0; i < npts; i++) {
    int j = 3*i; 
    fprintf(pf, "%10.4lf %10.4lf %10.4lf %8.4lf %8.4lf %8.4lf %13.8lf %13.8lf\n", 
	    p[j], p[j + 1], p[j + 2], vnp[j], vnp[j + 1], vnp[j + 2], f[i], f[i + npts]);
  }

  for (int i = 0; i < nelems; i++) {
    int j = 3*i; 
    fprintf(pf, "%8d%8d%8d\n", node[j], node[j + 1], node[j + 2]);
  }

  fclose(pf);
}

