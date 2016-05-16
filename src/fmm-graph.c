#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "adap_fmm.h"

void BuildGraph(const double * const sources, const int nsources, 
		const double * const targets, const int ntargets, 
		const int s)
{
  // BuildGraph() function generates the adaptive graph for the fast
  // multipole method. Source and target ensembles are allowed to be
  // different. 

  // The partition criteria are: (1) A (target or source) box
  // containing more than s points, and (2) The box is adajcent to a
  // box of the other type that satisfies (1). 

  // Each target box is associated with up to five lists of source
  // boxes. 

  // The construction takes two phases. In the first phase, source and
  // target trees are generated level by level. At the same time, the
  // list 5 for each target box (consisting of same level adajcent
  // source boxes) is formed. In the second phase, a top-down sweep is
  // performed to form the lists 1, 3, and 4 for each target box if
  // they exist. 
  g_nslev = 0; 
  g_nsboxes = 0;
  g_ntlev = 0;
  g_ntboxes = 0; 
  const int xoff[] = {0, 1, 0, 1, 0, 1, 0, 1};
  const int yoff[] = {0, 0, 1, 1, 0, 0, 1, 1};
  const int zoff[] = {0, 0, 0, 0, 1, 1, 1, 1};

  // Phase 1, step 1: Determine the smallest bounding box enclosing
  // all the source and target points. This box is the root node of
  // the source and target trees. 
  double xmin = DBL_MAX;
  double xmax = -DBL_MAX;
  double ymin = xmin;
  double ymax = xmax;
  double zmin = xmin;
  double zmax = xmax;

  for (int it = 0; it < nsources; it++) {
    double sx = sources[3*it];
    double sy = sources[3*it + 1];
    double sz = sources[3*it + 2];
    xmin = MIN(xmin, sx); 
    xmax = MAX(xmax, sx); 
    ymin = MIN(ymin, sy);
    ymax = MAX(ymax, sy); 
    zmin = MIN(zmin, sz); 
    zmax = MAX(zmax, sz); 
  }
  
  for (int it = 0; it < ntargets; it++) {
    double tx = targets[3*it];
    double ty = targets[3*it + 1];
    double tz = targets[3*it + 2];
    xmin = MIN(xmin, tx); 
    xmax = MAX(xmax, tx); 
    ymin = MIN(ymin, ty); 
    ymax = MAX(ymax, ty); 
    zmin = MIN(zmin, tz); 
    zmax = MAX(zmax, tz); 
  }
   
  g_bbcenter[0] = (xmax + xmin)*0.5;
  g_bbcenter[1] = (ymax + ymin)*0.5;
  g_bbcenter[2] = (zmax + zmin)*0.5;
  
  double sidex = xmax - xmin;
  double sidey = ymax - ymin;
  double sidez = zmax - zmin;  
  g_size = MAX(MAX(sidex, sidey), sidez); 

  g_bbcorner[0] = g_bbcenter[0] - g_size/2;
  g_bbcorner[1] = g_bbcenter[1] - g_size/2;
  g_bbcorner[2] = g_bbcenter[2] - g_size/2;

  // Phase 1, step 2: Allocate temporary work space which can hold up
  // to 200 levels of partition attempts. 
  for (int it = 0; it < nsources; it++)
    g_mapsrc[it] = it;

  for (int it = 0; it < ntargets; it++)
    g_maptar[it] = it;

  const int maxlev = 200;
  fmmnode **snodes = (fmmnode **)calloc(maxlev, sizeof(fmmnode *));
  fmmnode **tnodes = (fmmnode **)calloc(maxlev, sizeof(fmmnode *));
  int *nsnodes = (int *)calloc(maxlev, sizeof(int));
  int *ntnodes = (int *)calloc(maxlev, sizeof(int));
  if (snodes == 0 || tnodes == 0 || nsnodes == 0 || ntnodes == 0) {
    ERRMSG("memory allocation failure"); 
    exit(-1); 
  }

  snodes[0] = (fmmnode *)calloc(1, sizeof(fmmnode));
  tnodes[0] = (fmmnode *)calloc(1, sizeof(fmmnode));
  if (snodes[0]==0 || tnodes[0]==0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  snodes[0][0] = (fmmnode) {0, ++g_nsboxes, 0, {0, 0, 0, 0, 0, 0, 0, 0},
			    0, 0, 0, 0, nsources, 0, 0, 0, 0, 0};
  tnodes[0][0] = (fmmnode) {0, ++g_ntboxes, 0, {0, 0, 0, 0, 0, 0, 0, 0},
			    0, 0, 0, 0, ntargets, 0, 0, 0, 0, 0};

  tnodes[0][0].list5 = (int *)calloc(2, sizeof(int));
  if (tnodes[0][0].list5 == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }
  tnodes[0][0].list5[0] = 1;
  tnodes[0][0].list5[1] = 1;

  nsnodes[0] = 1;
  ntnodes[0] = 1;

  // Phase 1, step 3: Generate source and target trees. 
  for (int L = 0; L < maxlev; L++) {
    if (L == maxlev) {
      ERRMSG("Too many partition levels attempted");
      exit(-1);
    }

    int ns = nsnodes[L];
    int nt = ntnodes[L];
    int nspart = 0;
    int ntpart = 0; 
    int sbox0 = snodes[L][0].boxid;
    double h = g_size/pow(2, L + 1);

    // Step 3.1: Label boxes that meet partition criteria. 
    for (int it = 0; it < nt; it++) {
      if (tnodes[L][it].npts > s) {
	int *list5 = &tnodes[L][it].list5[1];
	int nlist5 = tnodes[L][it].list5[0];
	for (int it1 = 0; it1 < nlist5; it1++) {
	  int index = list5[it1] - sbox0;
	  if (snodes[L][index].npts > s) {
	    tnodes[L][it].nchild = 1;
	    snodes[L][index].nchild = 1;
	  }
	}
      }
    }

    for (int it = 0; it < nt; it++)
      ntpart += tnodes[L][it].nchild;

    for (int it = 0; it < ns; it++)
      nspart += snodes[L][it].nchild; 

    // Step 3.2: Partition level L source boxes if needed. 
    if (nspart) {
      int *work = (int *)calloc(ns*32, sizeof(int));
      if (work == 0) {
	ERRMSG("memory allocation failure");
	exit(-1);
      }
      
      int *counts = work;
      int *addrs = counts + ns*8;
      int *ords = addrs + ns*8; 
      int *ids = ords + ns*8;
      
      for (int it = 0; it < ns; it++) {
	if (snodes[L][it].nchild) 
	  PartitionBox(&snodes[L][it], sources, h, &counts[8*it], &addrs[8*it], 
		       g_mapsrc); 
      }

      int boxcreated = 0;
      for (int it = 0; it < ns*8; it++) {
	if (counts[it]) {
	  ords[it] = boxcreated++;
	  ids[it] = ++g_nsboxes;
	}
      }

      snodes[L + 1] = (fmmnode *)calloc(boxcreated, sizeof(fmmnode));
      if (snodes[L + 1] == 0) {
	ERRMSG("memory allocation failure");
	exit(-1);
      }
      nsnodes[L + 1] = boxcreated; 

      for (int it = 0; it < ns; it++) {
	if (snodes[L][it].nchild) {
	  snodes[L][it].nchild = 0;
	  int *count = &counts[8*it]; 
	  int *addr = &addrs[8*it];
	  int *ord = &ords[8*it];
	  int *id = &ids[8*it]; 
	  for (int it1 = 0; it1 < 8; it1++) {
	    if (count[it1]) {
	      int index = ord[it1];
	      snodes[L + 1][index].level = L + 1;
	      snodes[L + 1][index].boxid = id[it1];
	      snodes[L + 1][index].parent = snodes[L][it].boxid;
	      snodes[L][it].child[it1] = id[it1];
	      snodes[L][it].nchild++;
	      snodes[L + 1][index].idx = 2*snodes[L][it].idx + xoff[it1];
	      snodes[L + 1][index].idy = 2*snodes[L][it].idy + yoff[it1];
	      snodes[L + 1][index].idz = 2*snodes[L][it].idz + zoff[it1];
	      snodes[L + 1][index].npts = count[it1];
	      snodes[L + 1][index].addr = snodes[L][it].addr + addr[it1];
	    }
	  }
	}
      }

      free(work);
    } else {
      g_nslev = L;
    }

    // Step 3.3: Partition level L target boxes if needed. 
    if (ntpart) {
      int *work = (int *)calloc(nt*32, sizeof(int));
      if (work == 0) {
	ERRMSG("memory allocation failure"); 
	exit(-1);
      }

      int *counts = work;
      int *addrs = counts + nt*8;
      int *ords = addrs + nt*8;
      int *ids = ords + nt*8;

      for (int it = 0; it < nt; it++) {
	if (tnodes[L][it].nchild) 
	  PartitionBox(&tnodes[L][it], targets, h, &counts[8*it], &addrs[8*it],
		       g_maptar);
      }

      int boxcreated = 0;
      for (int it = 0; it < 8*nt; it++) {
	if (counts[it]) {
	  ords[it] = boxcreated++;
	  ids[it] = ++g_ntboxes;
	}
      }

      tnodes[L + 1] = (fmmnode *)calloc(boxcreated, sizeof(fmmnode));
      if (tnodes[L + 1] == 0) {
	ERRMSG("memory allocation failure");
	exit(-1);
      }
      ntnodes[L + 1] = boxcreated;

      for (int it = 0; it < nt; it++) {
	if (tnodes[L][it].nchild) {
	  tnodes[L][it].nchild = 0;
	  int *count = &counts[8*it];
	  int *addr = &addrs[8*it];
	  int *ord = &ords[8*it];
	  int *id = &ids[8*it];
	  for (int it1 = 0; it1 < 8; it1++) {
	    if (count[it1]) {
	      int index = ord[it1];
	      tnodes[L + 1][index].level = L + 1;
	      tnodes[L + 1][index].boxid = id[it1]; 
	      tnodes[L + 1][index].parent = tnodes[L][it].boxid;
	      tnodes[L][it].child[it1] = id[it1];  
	      tnodes[L][it].nchild++;
	      tnodes[L + 1][index].idx = 2*tnodes[L][it].idx + xoff[it1];
	      tnodes[L + 1][index].idy = 2*tnodes[L][it].idy + yoff[it1];
	      tnodes[L + 1][index].idz = 2*tnodes[L][it].idz + zoff[it1];
	      tnodes[L + 1][index].npts = count[it1];
	      tnodes[L + 1][index].addr = tnodes[L][it].addr + addr[it1];
	    }
	  }
	}
      }

      for (int it = 0; it < boxcreated; it++) {
	int nlist5 = 0;
	int temp5[27] = {0};
	int tidx = tnodes[L + 1][it].idx;
	int tidy = tnodes[L + 1][it].idy;
	int tidz = tnodes[L + 1][it].idz; 
	int ord = tnodes[L + 1][it].parent - tnodes[L][0].boxid; 
	int *plist5 = &tnodes[L][ord].list5[1];
	int nplist5 = tnodes[L][ord].list5[0];
	for (int it1 = 0; it1 < nplist5; it1++) {
	  int ord1 = plist5[it1] - sbox0; 
	  for (int it2 = 0; it2 < 8; it2++) {
	    int child = snodes[L][ord1].child[it2];
	    if (child) {
	      int ord2 = child - snodes[L + 1][0].boxid;
	      int sidx = snodes[L + 1][ord2].idx;
	      int sidy = snodes[L + 1][ord2].idy;
	      int sidz = snodes[L + 1][ord2].idz;
	      int diffx = fabs(tidx - sidx) <= 1;
	      int diffy = fabs(tidy - sidy) <= 1;
	      int diffz = fabs(tidz - sidz) <= 1;
	      if (diffx*diffy*diffz)
		temp5[nlist5++] = child;
	    }
	  }
	}

	tnodes[L + 1][it].list5 = (int *)calloc(1 + nlist5, sizeof(int));
	if (tnodes[L + 1][it].list5 == 0) {
	  ERRMSG("memory allocation failure");
	  exit(-1);
	}
	
	tnodes[L + 1][it].list5[0] = nlist5;
	for (int it1 = 0; it1 < nlist5; it1++)
	  tnodes[L + 1][it].list5[it1 + 1] = temp5[it1];
      }

      free(work);
    } else {
      g_ntlev = L;
      break;
    }
  }

  // Phase 1, step 4: write levelwise stored source and target trees
  // into contiguous one-dimensional arrays. 
  g_sboxes = (fmmnode *)calloc(g_nsboxes + 1, sizeof(fmmnode));
  g_tboxes = (fmmnode *)calloc(g_ntboxes + 1, sizeof(fmmnode));
  g_scontent = (int *)calloc((g_nslev + 1)*2, sizeof(int));
  g_tcontent = (int *)calloc((g_ntlev + 1)*2, sizeof(int));

  if (g_sboxes == 0 || g_tboxes == 0 || g_scontent == 0 || g_tcontent == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int L = 0; L <= g_nslev; L++) {
    int ns = nsnodes[L];
    for (int it1 = 0; it1 < ns; it1++) {
      fmmnode *ibox = &snodes[L][it1];
      int boxid = ibox->boxid;
      g_sboxes[boxid].level = ibox->level;
      g_sboxes[boxid].boxid = boxid;
      g_sboxes[boxid].parent = ibox->parent;
      for (int it2 = 0; it2 < 8; it2++) 
	g_sboxes[boxid].child[it2] = ibox->child[it2];
      g_sboxes[boxid].nchild = ibox->nchild;
      g_sboxes[boxid].idx = ibox->idx;
      g_sboxes[boxid].idy = ibox->idy;
      g_sboxes[boxid].idz = ibox->idz;
      g_sboxes[boxid].npts = ibox->npts;
      g_sboxes[boxid].addr = ibox->addr;
    }
    g_scontent[2*L] = snodes[L][0].boxid;
    g_scontent[2*L + 1] = ns;
  }

  for (int L = 0; L <= g_ntlev; L++) {
    int nt = ntnodes[L];
    for (int it1 = 0; it1 < nt; it1++) {
      fmmnode *ibox = &tnodes[L][it1];
      int boxid = ibox->boxid;
      g_tboxes[boxid].level = ibox->level;
      g_tboxes[boxid].boxid = boxid;
      g_tboxes[boxid].parent = ibox->parent;
      for (int it2 = 0; it2 < 8; it2++) 
	g_tboxes[boxid].child[it2] = ibox->child[it2];
      g_tboxes[boxid].nchild = ibox->nchild;
      g_tboxes[boxid].idx = ibox->idx;
      g_tboxes[boxid].idy = ibox->idy;
      g_tboxes[boxid].idz = ibox->idz;
      g_tboxes[boxid].npts = ibox->npts;
      g_tboxes[boxid].addr = ibox->addr;
      int nlist5 = ibox->list5[0];
      g_tboxes[boxid].list5 = (int *)calloc(1 + nlist5, sizeof(int)); 
      if (g_tboxes[boxid].list5 == 0) {
	ERRMSG("memory allocation failure");
	exit(-1);
      }
      for (int it2 = 0; it2 <= nlist5; it2++)
	g_tboxes[boxid].list5[it2] = ibox->list5[it2];
      free(ibox->list5);
    }
    g_tcontent[2*L] = tnodes[L][0].boxid;
    g_tcontent[2*L + 1] = nt;
  }

  for (int it = 0; it < maxlev; it++) {  
    free(snodes[it]);
    free(tnodes[it]);
  }
  free(snodes);
  free(tnodes);
  free(nsnodes);
  free(ntnodes);

  // Phase 2: Perform a top-down sweep to form the lists 1, 3, and 4
  // for each target box if they exist. 
  for (int it = 0; it < 8; it++) {
    int child = g_tboxes[1].child[it]; 
    if (child)
      BuildList(child);
  }
}


void PartitionBox(const fmmnode * const ibox, const double * const points, 
		  const double h, int * const counts, int * const addrs, 
		  int * const map)
{
  int npoints = ibox->npts; 
  int start = ibox->addr; 
  int *imap = &map[start]; 
  double center0 = g_bbcorner[0] + (2*ibox->idx + 1)*h; 
  double center1 = g_bbcorner[1] + (2*ibox->idy + 1)*h;
  double center2 = g_bbcorner[2] + (2*ibox->idz + 1)*h; 

  int *record = (int *)calloc(npoints, sizeof(int)); 
  if (record == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int i = 0; i < npoints; i++) {
    int j = 3*imap[i]; 
    int bin = 4*(points[j + 2] > center2) + 2*(points[j + 1] > center1) 
      + (points[j] > center0);  
    counts[bin]++; 
    record[i] = bin; 
  }

  for (int i = 1; i < 8; i++) 
    addrs[i] = addrs[i - 1] + counts[i - 1]; 

  int *swap = (int *)calloc(npoints, sizeof(int)); 
  if (swap == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  int assigned[8] = {0}; 
  for (int i = 0; i < npoints; i++) {
    int bin = record[i]; 
    int offset = addrs[bin] + assigned[bin]; 
    swap[offset] = imap[i]; 
    assigned[bin]++; 
  }

  for (int i = 0; i < npoints; i++) 
    imap[i] = swap[i]; 

  free(record);
  free(swap);
}

void BuildList(const int ibox)
{
  // BuildList() function generates the lists for all the target
  // boxes of the subtree rooted at ibox. For each target box
  // processed, the construction takes two steps. In the first step,
  // it forms the lists that connect ibox with coarser level source
  // boxes. In the second step, it forms the lists that connect ibox
  // with the same or finer level source boxes. Afterwards,
  // BuildList() is called recursively on each of ibox's child boxes. 

  // BuildList() function slightly modifies the definition of list
  // 1. In the original paper from Greengard and Rokhlin, list 1 is
  // defined only for leaf target box. Here, list 1 is extended to
  // nonleaf target box as well, consisting of the adjacent leaf
  // source boxes that are at a coarser level than ibox. 

  linklist *list1 = (linklist *)calloc(1, sizeof(linklist));
  linklist *list4 = (linklist *)calloc(1, sizeof(linklist));
  if (list1 == 0 || list4 == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  int nlist1 = 0;
  int nlist4 = 0;
  int parent = g_tboxes[ibox].parent;
  int *plist1 = (g_tboxes[parent].list1 ? &g_tboxes[parent].list1[1] : 0);
  int nplist1 = (g_tboxes[parent].list1 ? g_tboxes[parent].list1[0] : 0);

  // Go over the list 1 of ibox's parent. For each source box in the
  // list, it belongs to ibox's list 1 if it is adjacent to
  // ibox. Otherwise, it belongs to ibox's list 4. 
  for (int it = 0; it < nplist1; it++) {
    int sbox = plist1[it]; 
    (IfAdjacent(&g_sboxes[sbox], &g_tboxes[ibox]) ? 
     PushStack(list1, &nlist1, sbox) : PushStack(list4, &nlist4, sbox));
  }

  // The list 4 of ibox is now completely constructed. Write its
  // content from linked list form into array storage. 
  if (nlist4) {
    g_tboxes[ibox].list4 = (int *)calloc(nlist4 + 1, sizeof(int));
    if (g_tboxes[ibox].list4 == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    g_tboxes[ibox].list4[0] = nlist4;
    PopAll(list4, &g_tboxes[ibox].list4[1]);
  } else {
    // Release the temporary working space. 
    free(list4);
  }

  if (g_tboxes[ibox].nchild) {
    // If ibox is a nonleaf box, go over its own list 4 and put each
    // leaf source box that is adjacent to ibox to its list
    // 1. Afterwards, invoke BuildList() on each of its child boxes. 
    int *list5 = &g_tboxes[ibox].list5[1];
    int nlist5 = g_tboxes[ibox].list5[0]; 

    for (int it = 0; it < nlist5; it++) {
      int sbox = list5[it];
      if (g_sboxes[sbox].nchild == 0) 
	PushStack(list1, &nlist1, sbox);
    }

    if (nlist1) {
      g_tboxes[ibox].list1 = (int *)calloc(nlist1 + 1, sizeof(int));
      if (g_tboxes[ibox].list1==0) {
	ERRMSG("memory allocation failure");
	exit(-1);
      }

      g_tboxes[ibox].list1[0] = nlist1;
      PopAll(list1, &g_tboxes[ibox].list1[1]);
    } else {
      // Release the working space
      free(list1);
    }

    for (int it = 0; it < 8; it++ ) {
      int child = g_tboxes[ibox].child[it]; 
      if (child) 
	BuildList(child);
    }
  } else {
    // If ibox is a leaf box, go over its own list 5, and complete the
    // configuration of the lists 1 and 3. 
    BuildFinerList(&g_tboxes[ibox], list1, nlist1);
  }
}

void PushStack(linklist * const head, int * const nelems, const int ibox)
{
  // Each linklist object has a capacity specified by CAP. 
  if (head->ndata < CAP) {
    head->data[head->ndata++] = ibox;
  } else {
    linklist *push = (linklist *)calloc(1, sizeof(linklist));
    if (push == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }
    push->data[push->ndata++] = ibox;
    push->next = head->next;
    head->next = push;
  }
  (*nelems)++;
}

void PopAll(linklist *head, int * const list)
{
  // The PopAll() function writes the content of a link list into
  // array storage. 
  int offset = 0; 
  while (head) {
    for (int it = 0; it < head->ndata; it++ )
      list[offset + it] = head->data[it];
    offset += head->ndata;
    linklist *pop = head; 
    head = head->next;
    free(pop);
  }
}

int PopStack(linklist *head, int * const nelems)
{
  int sbox = head->data[--head->ndata];
  (*nelems)--;
  if (*nelems && head->ndata < 0) {
    linklist *pop = head->next;
    head = head->next;
    free(pop);
  }
  return sbox;
}

void BuildFinerList(fmmnode * const ibox, linklist *list1, int nlist1)
{
  int *list5 = &ibox->list5[1];
  int nlist5 = ibox->list5[0]; 
  int nlist3 = 0; 
  int nstack = 0;

  linklist *list3 = (linklist *)calloc(1, sizeof(linklist));
  linklist *stack = (linklist *)calloc(1, sizeof(linklist));
  if (list3 == 0 || stack == 0) {
    ERRMSG("memory allocation failure");
    exit(-1);
  }

  for (int it = 0; it < nlist5; it++) {
    int member = list5[it]; 
    PushStack(stack, &nstack, member);
    while (nstack) {
      // Pop the first element from the stack for examination. 
      int sbox = PopStack(stack, &nstack);
      const fmmnode *snode = &g_sboxes[sbox];
      if (IfAdjacent(ibox, snode)) {
	// If sbox is a leaf box, then sbox belongs to the list 1 of
	// ibox. If sbox is a nonleaf box, need to further examine
	// each of its child. 
	if (snode->nchild) {
	  for (int it1 = 0; it1 < 8; it1++ ) {
	    int child = snode->child[it1]; 
	    if (child) 
	      PushStack(stack, &nstack, child);
	  }
	} else {
	  PushStack(list1, &nlist1, sbox);
	}
      } else {
	PushStack(list3, &nlist3, sbox);
      }
    }      
  }
  free(stack);

  if (nlist1) {
    ibox->list1 = (int *)calloc(nlist1 + 1, sizeof(int));
    if (ibox->list1 == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    ibox->list1[0] = nlist1;
    PopAll(list1, &ibox->list1[1]);
  } else {
    free(list1);
  }

  if (nlist3) {
    ibox->list3 = (int *)calloc(nlist3 + 1, sizeof(int));
    if (ibox->list3 == 0) {
      ERRMSG("memory allocation failure");
      exit(-1);
    }

    ibox->list3[0] = nlist3;
    PopAll(list3, &ibox->list3[1]);
  } else {
    free(list3);
  }
}

int IfAdjacent(const fmmnode * const box1, const fmmnode * const box2)
{
  // Function IfAdjacent() checks whether box1 and box2 are adjacent
  // or not. The function assumes that box1 is at a coarser or the
  // same level as box2. 
  int dim = pow(2, box2->level - box1->level);
  int idxL = dim*box1->idx;
  int idyL = dim*box1->idy;
  int idzL = dim*box1->idz;
  int idxH = idxL + dim - 1;
  int idyH = idyL + dim - 1;
  int idzH = idzL + dim - 1;
  int adjx = (box2->idx >= idxL - 1) && (box2->idx <= idxH + 1);
  int adjy = (box2->idy >= idyL - 1) && (box2->idy <= idyH + 1);
  int adjz = (box2->idz >= idzL - 1) && (box2->idz <= idzH + 1);
  return adjx*adjy*adjz; 
}

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
		      int * const w8, int * const nw8, int * const x8, int * const y8)
{
  int *list5 = &g_tboxes[tboxid].list5[1];
  int nlist5 = g_tboxes[tboxid].list5[0];

  *nuall = 0;  
  *nu1234 = 0;
  *ndall = 0; 
  *nd5678 = 0;
  *nnall = 0; 
  *nn1256 = 0; 
  *nn12 = 0; 
  *nn56 = 0;
  *nsall = 0; 
  *ns3478 = 0; 
  *ns34 = 0; 
  *ns78 = 0;
  *neall = 0; 
  *ne1357 = 0; 
  *ne13 = 0; 
  *ne57 = 0; 
  *ne1 = 0;
  *ne3 = 0; 
  *ne5 = 0; 
  *ne7 = 0;
  *nwall = 0; 
  *nw2468 = 0; 
  *nw24 = 0; 
  *nw68 = 0; 
  *nw2 = 0; 
  *nw4 = 0; 
  *nw6 = 0; 
  *nw8 = 0;

  int tidx = g_tboxes[tboxid].idx;
  int tidy = g_tboxes[tboxid].idy;
  int tidz = g_tboxes[tboxid].idz;

  for (int it = 0; it < nlist5; it++ ) {
    int member = list5[it];
    int sidx = g_sboxes[member].idx;
    int sidy = g_sboxes[member].idy;
    int sidz = g_sboxes[member].idz;    

    // There are up to 27 members in the list5, [x][y][z], where each
    // index is from -1 to 1. 
    int offset = 9*(sidz - tidz + 1) + 3*(sidy - tidy + 1) + sidx - tidx + 1;

    switch (offset) {
    case 0:
      // [-1][-1][-1], update dall, sall, wall, d5678, s34, w2 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], -2, -2);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], -1,-2);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], -2, -1);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], -1,-1);
      if (g_sboxes[member].child[4]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[4], -1, -2);
      if (g_sboxes[member].child[5]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[5], -1, -1);
      if (g_sboxes[member].child[6]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[6], 1, -1);
      if (g_sboxes[member].child[7]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[7], -1, -1);
	UpdateList(s34, ns34, x34, y34, g_sboxes[member].child[7], -1, -1);
	UpdateList(w2, nw2, x2, y2, g_sboxes[member].child[7], 1, -1);
      }
      break;
    case 1://[0][-1][-1], update dall, sall, d5678, s34 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], 0, -2);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], 1, -2);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], 0, -1);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], 1,-1);
      if (g_sboxes[member].child[4]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[4], -1, 0);
      if (g_sboxes[member].child[5]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[5], -1, 1);     
      if (g_sboxes[member].child[6]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[6], 0, -1);
	UpdateList(s34, ns34, x34, y34, g_sboxes[member].child[6], -1,0);
      }
      if (g_sboxes[member].child[7]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[7], 1, -1);
	UpdateList(s34, ns34, x34, y34, g_sboxes[member].child[7], -1, 1);
      }
      break;
    case 2: //[1][-1][-1], update dall, sall, eall, d5678, s34, e1 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], 2, -2);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], 3, -2);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], 2, -1);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], 3, -1);
      if (g_sboxes[member].child[4]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[4], -1, 2);
      if (g_sboxes[member].child[5]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[5], -1, 3);
      if (g_sboxes[member].child[6]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[6], 2, -1);
	UpdateList(s34, ns34, x34, y34, g_sboxes[member].child[6], -1, 2);
	UpdateList(e1, ne1, x1, y1, g_sboxes[member].child[6], 1, -1);
      }
      if (g_sboxes[member].child[7]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[7], 1, -1);
      break;
    case 3: //[-1][0][-1], update dall, wall, d5678, w24 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], -2, 0);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], -1, 0);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], -2, 1);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], -1, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[4], 1, 0);
      if (g_sboxes[member].child[5]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[5], -1, 0);
	UpdateList(w24, nw24, x24, y24, g_sboxes[member].child[5], 1, 0);
      }
      if (g_sboxes[member].child[6]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[6], 1, 1);
      if (g_sboxes[member].child[7]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[7], -1, 1);
	UpdateList(w24, nw24, x24, y24, g_sboxes[member].child[7], 1, 1);
      }
      break;
    case 4://[0][0][-1], update dall and d5678 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], 0, 0);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], 1, 0);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], 0, 1);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], 1, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[4], 0, 0);
      if (g_sboxes[member].child[5]) 
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[5], 1, 0);
      if (g_sboxes[member].child[6]) 
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[6], 0, 1);
      if (g_sboxes[member].child[7]) 
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[7], 1, 1);
      break;
    case 5: // [1][0][-1], update dall, d5678, eall, e13 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], 2, 0);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], 3, 0);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], 2, 1);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], 3, 1);
      if (g_sboxes[member].child[4]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[4], 2, 0);
	UpdateList(e13, ne13, x13, y13, g_sboxes[member].child[4], 1, 0);
      }
      if (g_sboxes[member].child[5]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[5], 1, 0);
      if (g_sboxes[member].child[6]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[6], 2, 1);
	UpdateList(e13, ne13, x13, y13, g_sboxes[member].child[6], 1, 1);
      }
      if (g_sboxes[member].child[7]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[7], 1, 1);
      break;
    case 6://[-1][1][-1], update dall, nall, wall, d5678, n12, w4 list3
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], -2, 2);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], -1, 2);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], -2, 3);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], -1, 3);
      if (g_sboxes[member].child[4]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[4], 1, 2);
      if (g_sboxes[member].child[5]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[5], -1, 2);
	UpdateList(n12, nn12, x12, y12, g_sboxes[member].child[5], -1, -1);
	UpdateList(w4, nw4, x4, y4, g_sboxes[member].child[5], 1, 2);
      }
      if (g_sboxes[member].child[6]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[6], -1, -2);
      if (g_sboxes[member].child[7]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[7], -1, -1);
      break;
    case 7: //[0][1][-1], update dallm d5678, nall, n12 lists 
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], 0, 2);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], 1, 2);
      if (g_sboxes[member].child[2]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], 0, 3);
      if (g_sboxes[member].child[3]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], 1, 3);
      if (g_sboxes[member].child[4]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[4], 0, 2);
	UpdateList(n12, nn12, x12, y12, g_sboxes[member].child[4], -1, 0);
      }
      if (g_sboxes[member].child[5]) {
	UpdateList(d5678, nd5678,x5678, y5678, g_sboxes[member].child[5], 1, 2);
	UpdateList(n12, nn12, x12, y12, g_sboxes[member].child[5], -1, 1);
      }
      if (g_sboxes[member].child[6]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[6], -1, 0);
      if (g_sboxes[member].child[7])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[7], -1, 1);
      break;
    case 8: //[1][1][-1], update dall, d5678, nall, eall, n12, e3 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[0], 2,2);
      if (g_sboxes[member].child[1]) 
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[1], 3, 2);
      if (g_sboxes[member].child[2])
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[2], 2, 3);
      if (g_sboxes[member].child[3])
	UpdateList(dall, ndall, xdall, ydall, g_sboxes[member].child[3], 3,3);
      if (g_sboxes[member].child[4]) {
	UpdateList(d5678, nd5678, x5678, y5678, g_sboxes[member].child[4], 2, 2);
	UpdateList(n12, nn12, x12, y12, g_sboxes[member].child[4], -1, 2);
	UpdateList(e3, ne3, x3, y3, g_sboxes[member].child[4], 1, 2);
      }
      if (g_sboxes[member].child[5]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[5], 1, 2);
      if (g_sboxes[member].child[6]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[6], -1, 2);
      if (g_sboxes[member].child[7]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[7], -1, 3);
      break;
    case 9: // [-1][-1][0], update sall, wall, s3478, w2, w6 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[0], 0, -2);
      if (g_sboxes[member].child[1])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[1], 0, -1);
      if (g_sboxes[member].child[2])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[2], 0, -1);
      if (g_sboxes[member].child[3]) {
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[3], 0, -1);
	UpdateList(w2, nw2, x2, y2, g_sboxes[member].child[3], 0, -1);
	UpdateList(w6, nw6, x6, y6, g_sboxes[member].child[3], 0, -1);
      }
      if (g_sboxes[member].child[4]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[4], 1, -2);
      if (g_sboxes[member].child[5])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[5], 1, -1);
      if (g_sboxes[member].child[6]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[6], -1, -1);
      if (g_sboxes[member].child[7]) {
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[7], 1, -1);
	UpdateList(w2, nw2, x2, y2, g_sboxes[member].child[7], -1, -1);
	UpdateList(w6, nw6, x6, y6, g_sboxes[member].child[7], -1, -1);
      }
      break;
    case 10://[0][-1][0], update sall, s3478 lists
      if (g_sboxes[member].child[0]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[0], 0, 0);
      if (g_sboxes[member].child[1])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[1], 0, 1);
      if (g_sboxes[member].child[2])
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[2], 0, 0);
      if (g_sboxes[member].child[3])
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[3], 0, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[4], 1, 0);
      if (g_sboxes[member].child[5])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[5], 1, 1);
      if (g_sboxes[member].child[6])
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[6], 1, 0);
      if (g_sboxes[member].child[7])
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[7], 1, 1);
      break;
    case 11: //[1][-1][0], update eall, sall, s3478, e1, e5 lists
      if (g_sboxes[member].child[0])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[0], 0, 2);
      if (g_sboxes[member].child[1])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[1], 0, 3);
      if (g_sboxes[member].child[2]) {
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[2], 0, 2);
	UpdateList(e1, ne1, x1, y1, g_sboxes[member].child[2], 0, -1);
	UpdateList(e5, ne5, x5, y5, g_sboxes[member].child[2], 0, -1);
      }
      if (g_sboxes[member].child[3]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[3], 0, -1);
      if (g_sboxes[member].child[4])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[4], 1, 2);
      if (g_sboxes[member].child[5]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[5], 1, 3);
      if (g_sboxes[member].child[6]) {
	UpdateList(s3478, ns3478, x3478, y3478, g_sboxes[member].child[6], 1, 2);
	UpdateList(e1, ne1, x1, y1, g_sboxes[member].child[6], -1, -1);
	UpdateList(e5, ne5, x5, y5, g_sboxes[member].child[6], -1, -1);
      }
      if (g_sboxes[member].child[7]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[7], -1, -1);
      break;
    case 12: // [-1][0][0], update wall, w2468 lists
      if (g_sboxes[member].child[0])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[0], 0, 0);
      if (g_sboxes[member].child[1])
	UpdateList(w2468, nw2468, x2468, y2468, g_sboxes[member].child[1], 0, 0);
      if (g_sboxes[member].child[2])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[2], 0, 1);
      if (g_sboxes[member].child[3])
	UpdateList(w2468, nw2468, x2468, y2468, g_sboxes[member].child[3], 0, 1);
      if (g_sboxes[member].child[4])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[4], -1, 0);
      if (g_sboxes[member].child[5])
	UpdateList(w2468, nw2468, x2468, y2468, g_sboxes[member].child[5], -1, 0);
      if (g_sboxes[member].child[6])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[6], -1, 1);
      if (g_sboxes[member].child[7])
	UpdateList(w2468, nw2468, x2468, y2468, g_sboxes[member].child[7], -1, 1);
      break;
    case 13: //[0][0][0], nothing here
      break;
    case 14: //[1][0][0], update eall, e1357 lists
      if (g_sboxes[member].child[0])
	UpdateList(e1357, ne1357, x1357, y1357, g_sboxes[member].child[0], 0, 0);
      if (g_sboxes[member].child[1])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[1], 0, 0);
      if (g_sboxes[member].child[2])
	UpdateList(e1357, ne1357, x1357, y1357, g_sboxes[member].child[2], 0, 1);
      if (g_sboxes[member].child[3])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[3], 0, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(e1357, ne1357, x1357, y1357, g_sboxes[member].child[4], -1, 0);
      if (g_sboxes[member].child[5])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[5], -1, 0);
      if (g_sboxes[member].child[6])
	UpdateList(e1357, ne1357, x1357, y1357, g_sboxes[member].child[6], -1, 1);
      if (g_sboxes[member].child[7])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[7], -1, 1);
      break;
    case 15://[-1][1][0], update wall, nall, n1256, w4, w8 lists
      if (g_sboxes[member].child[0])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[0], 0, 2);
      if (g_sboxes[member].child[1]) {
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[1], 0, -1);
	UpdateList(w4, nw4, x4, y4, g_sboxes[member].child[1], 0, 2);
	UpdateList(w8, nw8, x8, y8, g_sboxes[member].child[1], 0, 2);
      }
      if (g_sboxes[member].child[2])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[2], 0, -2);
      if (g_sboxes[member].child[3])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[3], 0, -1);
      if (g_sboxes[member].child[4])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[4], -1, 2);
      if (g_sboxes[member].child[5]) {
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[5], 1, -1);
	UpdateList(w4, nw4, x4, y4, g_sboxes[member].child[5], -1, 2);
	UpdateList(w8, nw8, x8, y8, g_sboxes[member].child[5], -1, 2);
      }
      if (g_sboxes[member].child[6]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[6], 1, -2);
      if (g_sboxes[member].child[7])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[7], 1, -1);
      break;
    case 16: //[0][1][0], update nall, n1256 lists
      if (g_sboxes[member].child[0])
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[0], 0, 0);
      if (g_sboxes[member].child[1])
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[1], 0, 1);
      if (g_sboxes[member].child[2])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[2], 0, 0);
      if (g_sboxes[member].child[3])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[3], 0, 1);
      if (g_sboxes[member].child[4])
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[4], 1, 0);
      if (g_sboxes[member].child[5])
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[5], 1, 1);
      if (g_sboxes[member].child[6])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[6], 1, 0);
      if (g_sboxes[member].child[7])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[7], 1, 1);
      break;
    case 17: //[1][1][0], update nall, n1256, eall, e3, e7 lists
      if (g_sboxes[member].child[0]) {
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[0], 0, 2);
	UpdateList(e3, ne3, x3, y3, g_sboxes[member].child[0], 0, 2);
	UpdateList(e7, ne7, x7, y7, g_sboxes[member].child[0], 0, 2);
      }
      if (g_sboxes[member].child[1]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[1], 0, 2);
      if (g_sboxes[member].child[2])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[2], 0, 2);
      if (g_sboxes[member].child[3])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[3], 0, 3);
      if (g_sboxes[member].child[4]) {
	UpdateList(n1256, nn1256, x1256, y1256, g_sboxes[member].child[4], 1, 2);
	UpdateList(e3, ne3, x3, y3, g_sboxes[member].child[4], -1, 2);
	UpdateList(e7, ne7, x7, y7, g_sboxes[member].child[4], -1, 2);
      }
      if (g_sboxes[member].child[5])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[5], -1, 2);
      if (g_sboxes[member].child[6])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[6], 1, 2);
      if (g_sboxes[member].child[7])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[7], 1, 3);
      break;
    case 18: //[-1][-1][1], update sall, wall, u1234, s78, w6, uall lists
      if (g_sboxes[member].child[0]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[0], 2, -2);
      if (g_sboxes[member].child[1]) 
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[1], 2, -1);
      if (g_sboxes[member].child[2]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[2], -2, -1);
      if (g_sboxes[member].child[3]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[3], -1,-1);
	UpdateList(s78, ns78, x78, y78, g_sboxes[member].child[3], 2, -1);
	UpdateList(w6, nw6, x6, y6, g_sboxes[member].child[3], -2, -1);
      }
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], -2, -2);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], -1, -2);
      if (g_sboxes[member].child[6])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], -2, -1);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], -1, -1);
      break;
    case 19: //[0][-1][1], update sall, u1234, s78, uall lists
      if (g_sboxes[member].child[0])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[0], 2, 0);
      if (g_sboxes[member].child[1])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[1], 2, 1);
      if (g_sboxes[member].child[2]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[2], 0, -1);
	UpdateList(s78, ns78, x78, y78, g_sboxes[member].child[2], 2, 0);
      }
      if (g_sboxes[member].child[3]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[3], 1, -1);
	UpdateList(s78, ns78, x78, y78, g_sboxes[member].child[3], 2, 1);
      }
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], 0, -2);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], 1, -2);
      if (g_sboxes[member].child[6])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], 0, -1);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], 1, -1);
      break;
    case 20: // [1][-1][1], update sall, eall, u1234, s78, e5, uall lists
      if (g_sboxes[member].child[0])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[0], 2, 2);
      if (g_sboxes[member].child[1])
	UpdateList(sall, nsall, xsall, ysall, g_sboxes[member].child[1], 2, 3);
      if (g_sboxes[member].child[2]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[2], 2, -1);
	UpdateList(s78, ns78, x78, y78, g_sboxes[member].child[2], 2, 2);
	UpdateList(e5, ne5, x5, y5, g_sboxes[member].child[2], -2, -1);
      }
      if (g_sboxes[member].child[3])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[3], -2, -1);
      if (g_sboxes[member].child[4])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], 2, -2);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], 3, -2);
      if (g_sboxes[member].child[6])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], 2, -1);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], 3, -1);
      break;
    case 21: // [-1][0][1], update wall, u1234, w68, uall lists
      if (g_sboxes[member].child[0])
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[0], -2, 0);
      if (g_sboxes[member].child[1]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[1], -1, 0);
	UpdateList(w68, nw68, x68, y68, g_sboxes[member].child[1], -2, 0);
      }
      if (g_sboxes[member].child[2]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[2], -2, 1);
      if (g_sboxes[member].child[3]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[3], -1, 1);
	UpdateList(w68, nw68, x68, y68, g_sboxes[member].child[3], -2, 1);
      }
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], -2, 0);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], -1, 0);
      if (g_sboxes[member].child[6])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], -2, 1);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], -1, 1);
      break;
    case 22: //[0][0][1], update u1234, uall lists
      if (g_sboxes[member].child[0]) 
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[0], 0, 0);
      if (g_sboxes[member].child[1])
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[1], 1, 0);
      if (g_sboxes[member].child[2]) 
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[2], 0, 1);
      if (g_sboxes[member].child[3])
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[3], 1, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], 0, 0);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], 1, 0);
      if (g_sboxes[member].child[6]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], 0, 1);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], 1, 1);
      break;
    case 23: // [1][0][1], update u1234, e57, eall, uall lists
      if (g_sboxes[member].child[0]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[0], 2, 0);
	UpdateList(e57, ne57, x57, y57, g_sboxes[member].child[0], -2, 0);
      }
      if (g_sboxes[member].child[1])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[1], -2, 0);
      if (g_sboxes[member].child[2]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[2], 2, 1);
	UpdateList(e57, ne57, x57, y57, g_sboxes[member].child[2], -2, 1);
      }
      if (g_sboxes[member].child[3])
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[3], -2, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], 2, 0);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], 3, 0);
      if (g_sboxes[member].child[6]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], 2, 1);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], 3, 1);	
      break;
    case 24: // [-1][1][1], update nall, wall, u1234, n56, w8, uall lists
      if (g_sboxes[member].child[0]) 
	UpdateList(wall, nwall, xwall, ywall, g_sboxes[member].child[0], -2, 2);
      if (g_sboxes[member].child[1]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[1], -1, 2);
	UpdateList(n56, nn56, x56, y56, g_sboxes[member].child[1], 2, -1);
	UpdateList(w8, nw8, x8, y8, g_sboxes[member].child[1], -2, 2);
      }
      if (g_sboxes[member].child[2]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[2], 2, -2);
      if (g_sboxes[member].child[3])
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[3], 2, -1);
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], -2, 2);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], -1, 2);
      if (g_sboxes[member].child[6]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], -2, 3);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], -1, 3);	
      break;
    case 25: // [0][1][1], update u1234, nall, n56, nall lists
      if (g_sboxes[member].child[0]) {
	UpdateList(u1234, nu1234, x1234, y1234,  g_sboxes[member].child[0], 0, 2);
	UpdateList(n56, nn56, x56, y56,  g_sboxes[member].child[0], 2, 0);
      }
      if (g_sboxes[member].child[1]) {
	UpdateList(u1234, nu1234, x1234, y1234,  g_sboxes[member].child[1],1, 2);
	UpdateList(n56, nn56, x56, y56,  g_sboxes[member].child[1], 2, 1);
      }
      if (g_sboxes[member].child[2]) 
	UpdateList(nall, nnall, xnall, ynall,  g_sboxes[member].child[2], 2, 0);
      if (g_sboxes[member].child[3])
	UpdateList(nall, nnall, xnall, ynall,  g_sboxes[member].child[3], 2, 1);
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], 0, 2);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], 1, 2);
      if (g_sboxes[member].child[6]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], 0, 3);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], 1, 3);	
      break;
    case 26: // [1][1][1], update u1234, n56, e7, eall, nall, uall lists 
      if (g_sboxes[member].child[0]) {
	UpdateList(u1234, nu1234, x1234, y1234, g_sboxes[member].child[0], 2, 2);
	UpdateList(n56, nn56, x56, y56, g_sboxes[member].child[0], 2, 2);
	UpdateList(e7, ne7, x7, y7, g_sboxes[member].child[0], -2, 2);
      }
      if (g_sboxes[member].child[1]) 
	UpdateList(eall, neall, xeall, yeall, g_sboxes[member].child[1], -2, 2);
      if (g_sboxes[member].child[2]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[2], 2, 2);
      if (g_sboxes[member].child[3]) 
	UpdateList(nall, nnall, xnall, ynall, g_sboxes[member].child[3], 2, 3);
      if (g_sboxes[member].child[4]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[4], 2, 2);
      if (g_sboxes[member].child[5])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[5], 3, 2);
      if (g_sboxes[member].child[6]) 
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[6], 2, 3);
      if (g_sboxes[member].child[7])
	UpdateList(uall, nuall, xuall, yuall, g_sboxes[member].child[7], 3, 3);	
      break;
    default:
      break;
    }
  }
}

void UpdateList(int * const list, int * const nlist, int * const xoff, 
		int * const yoff, int entry, int ix, int iy)
{
  list[*nlist] = entry; 
  xoff[*nlist] = ix; 
  yoff[*nlist] = iy; 
  *nlist += 1;
}

void BuildDirectList13(int **directlist)
{
  // BuildDirectList13() function is specially developed for the AFMPB
  // solver. For each target point in a leaf target box, the function
  // forms the set of source points enclosed in the lists 1 and 3 of
  // the target box. 

  for (int ibox = 1; ibox <= g_ntboxes; ibox++) {
    int nbrpts = 0, nnbr1 = 0, nnbr3 = 0; 
    if (g_tboxes[ibox].list1 != 0 && g_tboxes[ibox].nchild == 0) {
      nnbr1 = g_tboxes[ibox].list1[0]; 
      for (int i = 1; i <= nnbr1; i++) {
	int sbox = g_tboxes[ibox].list1[i]; 
	nbrpts += g_sboxes[sbox].npts; 
      }
    }

    if (g_tboxes[ibox].list3 != 0) {
      nnbr3 = g_tboxes[ibox].list3[0]; 
      for (int i = 1; i <= nnbr3; i++) {
	int sbox = g_tboxes[ibox].list3[i]; 
	nbrpts += g_sboxes[sbox].npts; 
      }
    }

    if (nbrpts) {
      directlist[ibox] = (int *)calloc(2*(1 + nnbr1 + nnbr3), sizeof(int)); 
      if (directlist[ibox] == 0) {
	ERRMSG("memory allocation failure");
	exit(-1);
      }

      directlist[ibox][0] = nbrpts; 
      directlist[ibox][1] = nnbr1 + nnbr3; 

      int j = 1; 
      for (int i = 1; i <= nnbr1; i++) {
	int sbox = g_tboxes[ibox].list1[i]; 
	int first = g_sboxes[sbox].addr;
	int last = first + g_sboxes[sbox].npts - 1; 
	directlist[ibox][2*j] = first; 
	directlist[ibox][2*j + 1] = last; 
	j++;
      }

      for (int i = 1; i <= nnbr3; i++) {
	int sbox = g_tboxes[ibox].list3[i]; 
	int first = g_sboxes[sbox].addr;
	int last = first + g_sboxes[sbox].npts - 1; 
	directlist[ibox][2*j] = first; 
	directlist[ibox][2*j + 1] = last; 
	j++;
      }	
    }
  }
}

void CleanDirectList13(int **directlist)
{
  for (int ibox = 0; ibox <= g_ntboxes; ibox++) {
    if (directlist[ibox] != 0)
      free(directlist[ibox]);
  }
  free(directlist);
}

void DestroyGraph(void)
{
  for (int i = 0; i <= g_ntboxes; i++) {
    free(g_tboxes[i].list1);
    free(g_tboxes[i].list3);
    free(g_tboxes[i].list4);
    free(g_tboxes[i].list5);
  }
  free(g_sboxes);
  free(g_scontent);
  free(g_tboxes);
  free(g_tcontent);
  free(g_mapsrc);
  free(g_maptar);
}
