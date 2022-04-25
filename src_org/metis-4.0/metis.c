/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * balance.c
 *
 * This file contains code that is used to forcefully balance either
 * bisections or k-sections
 *
 * Started 7/29/97
 * George
 *
 * $Id$
 *
 */

#include <metis.h>

/*************************************************************************
* This function is the entry point of the bisection balancing algorithms.
**************************************************************************/
void Balance2Way(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  int mindiff;

  /* Return right away if the balance is OK */
  mindiff = abs(tpwgts[0]-graph->pwgts[0]);
  if (mindiff < 3*(graph->pwgts[0]+graph->pwgts[1])/graph->nvtxs)
    return;
  if (graph->pwgts[0] > tpwgts[0] && graph->pwgts[0] < (int)(ubfactor*tpwgts[0]))
    return;
  if (graph->pwgts[1] > tpwgts[1] && graph->pwgts[1] < (int)(ubfactor*tpwgts[1]))
    return;

  if (graph->nbnd > 0)
    Bnd2WayBalance(ctrl, graph, tpwgts);
  else
    General2WayBalance(ctrl, graph, tpwgts);

}



/*************************************************************************
* This function balances two partitions by moving boundary nodes
* from the domain that is overweight to the one that is underweight.
**************************************************************************/
void Bnd2WayBalance(CtrlType *ctrl, GraphType *graph, int *tpwgts)
{
  int i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, tmp;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind, *pwgts;
  idxtype *moved, *perm;
  PQueueType parts;
  int higain, oldgain, mincut, mindiff;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  pwgts = graph->pwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);

  /* Determine from which domain you will be moving data */
  mindiff = abs(tpwgts[0]-pwgts[0]);
  from = (pwgts[0] < tpwgts[0] ? 1 : 0);
  to = (from+1)%2;

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d] T[%6d %6d], Nv-Nb[%6d %6d]. ICut: %6d [B]\n",
             pwgts[0], pwgts[1], tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  tmp = graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)];
  PQueueInit(ctrl, &parts, nvtxs, tmp);

  idxset(nvtxs, -1, moved);

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  /* Insert the boundary nodes of the proper partition whose size is OK in the priority queue */
  nbnd = graph->nbnd;
  RandomPermute(nbnd, perm, 1);
  for (ii=0; ii<nbnd; ii++) {
    i = perm[ii];
    ASSERT(ed[bndind[i]] > 0 || id[bndind[i]] == 0);
    ASSERT(bndptr[bndind[i]] != -1);
    if (where[bndind[i]] == from && vwgt[bndind[i]] <= mindiff)
      PQueueInsert(&parts, bndind[i], ed[bndind[i]]-id[bndind[i]]);
  }

  mincut = graph->mincut;
  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if ((higain = PQueueGetMax(&parts)) == -1)
      break;
    ASSERT(bndptr[higain] != -1);

    if (pwgts[to]+vwgt[higain] > tpwgts[to])
      break;

    mincut -= (ed[higain]-id[higain]);
    INC_DEC(pwgts[to], pwgts[from], vwgt[higain]);

    where[higain] = to;
    moved[higain] = nswaps;

    IFSET(ctrl->dbglvl, DBG_MOVEINFO,
      printf("Moved %6d from %d. [%3d %3d] %5d [%4d %4d]\n", higain, from, ed[higain]-id[higain], vwgt[higain], mincut, pwgts[0], pwgts[1]));

    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update its boundary information and queue position */
      if (bndptr[k] != -1) { /* If k was a boundary vertex */
        if (ed[k] == 0) { /* Not a boundary vertex any more */
          BNDDelete(nbnd, bndind, bndptr, k);
          if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff)  /* Remove it if in the queues */
            PQueueDelete(&parts, k, oldgain);
        }
        else { /* If it has not been moved, update its position in the queue */
          if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff)
            PQueueUpdate(&parts, k, oldgain, ed[k]-id[k]);
        }
      }
      else {
        if (ed[k] > 0) {  /* It will now become a boundary vertex */
          BNDInsert(nbnd, bndind, bndptr, k);
          if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff)
            PQueueInsert(&parts, k, ed[k]-id[k]);
        }
      }
    }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("\tMinimum cut: %6d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, pwgts[0], pwgts[1], nbnd));

  graph->mincut = mincut;
  graph->nbnd = nbnd;

  PQueueFree(ctrl, &parts);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function balances two partitions by moving the highest gain
* (including negative gain) vertices to the other domain.
* It is used only when tha unbalance is due to non contigous
* subdomains. That is, the are no boundary vertices.
* It moves vertices from the domain that is overweight to the one that
* is underweight.
**************************************************************************/
void General2WayBalance(CtrlType *ctrl, GraphType *graph, int *tpwgts)
{
  int i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, tmp;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind, *pwgts;
  idxtype *moved, *perm;
  PQueueType parts;
  int higain, oldgain, mincut, mindiff;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  pwgts = graph->pwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);

  /* Determine from which domain you will be moving data */
  mindiff = abs(tpwgts[0]-pwgts[0]);
  from = (pwgts[0] < tpwgts[0] ? 1 : 0);
  to = (from+1)%2;

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d] T[%6d %6d], Nv-Nb[%6d %6d]. ICut: %6d [B]\n",
             pwgts[0], pwgts[1], tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  tmp = graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)];
  PQueueInit(ctrl, &parts, nvtxs, tmp);

  idxset(nvtxs, -1, moved);

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  /* Insert the nodes of the proper partition whose size is OK in the priority queue */
  RandomPermute(nvtxs, perm, 1);
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    if (where[i] == from && vwgt[i] <= mindiff)
      PQueueInsert(&parts, i, ed[i]-id[i]);
  }

  mincut = graph->mincut;
  nbnd = graph->nbnd;
  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if ((higain = PQueueGetMax(&parts)) == -1)
      break;

    if (pwgts[to]+vwgt[higain] > tpwgts[to])
      break;

    mincut -= (ed[higain]-id[higain]);
    INC_DEC(pwgts[to], pwgts[from], vwgt[higain]);

    where[higain] = to;
    moved[higain] = nswaps;

    IFSET(ctrl->dbglvl, DBG_MOVEINFO,
      printf("Moved %6d from %d. [%3d %3d] %5d [%4d %4d]\n", higain, from, ed[higain]-id[higain], vwgt[higain], mincut, pwgts[0], pwgts[1]));

    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff)
        PQueueUpdate(&parts, k, oldgain, ed[k]-id[k]);

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("\tMinimum cut: %6d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, pwgts[0], pwgts[1], nbnd));

  graph->mincut = mincut;
  graph->nbnd = nbnd;

  PQueueFree(ctrl, &parts);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * bucketsort.c
 *
 * This file contains code that implement a variety of counting sorting
 * algorithms
 *
 * Started 7/25/97
 * George
 *
 * $Id$
 *
 */





/*************************************************************************
* This function uses simple counting sort to return a permutation array
* corresponding to the sorted order. The keys are assumed to start from
* 0 and they are positive.  This sorting is used during matching.
**************************************************************************/
void BucketSortKeysInc(int n, int max, idxtype *keys, idxtype *tperm, idxtype *perm)
{
  int i, ii;
  idxtype *counts;

  counts = idxsmalloc(max+2, 0, "BucketSortKeysInc: counts");

  for (i=0; i<n; i++)
    counts[keys[i]]++;
  MAKECSR(i, max+1, counts);

  for (ii=0; ii<n; ii++) {
    i = tperm[ii];
    perm[counts[keys[i]]++] = i;
  }

  free(counts);
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ccgraph.c
 *
 * This file contains the functions that create the coarse graph
 *
 * Started 8/11/97
 * George
 *
 * $Id$
 *
 */





/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, jj, k, kk, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask, dovsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);

  mask = HTLENGTH;
  if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) {
    CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
    return;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore;
  memcpy(auxadj, adjncy, iend*sizeof(idxtype));
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1));

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs)
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      scopy(ncon, nvwgt+v*ncon, cnvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m] += adjwgt[j];
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj] += adjwgt[j];
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = adjwgt[j];
        }
      }
    }

    if (v != u) {
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        saxpy(ncon, 1.0, nvwgt+u*ncon, 1, cnvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m] += adjwgt[j];
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj] += adjwgt[j];
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = adjwgt[j];
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs)
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) { /* This 2nd check is needed for non-adjacent matchings */
        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraphNoMask(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, k, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, dovsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;


  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  htable = idxset(cnvtxs, -1, idxwspacemalloc(ctrl, cnvtxs));

  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore;
  memcpy(auxadj, adjncy, iend*sizeof(idxtype));
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs)
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      scopy(ncon, nvwgt+v*ncon, cnvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      if ((m = htable[k]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[k] = nedges++;
      }
      else {
        cadjwgt[m] += adjwgt[j];
      }
    }

    if (v != u) {
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        saxpy(ncon, 1.0, nvwgt+u*ncon, 1, cnvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        if ((m = htable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }

      /* Remove the contracted adjacency weight */
      if ((j = htable[cnvtxs]) != -1) {
        ASSERT(cadjncy[j] == cnvtxs);
        cadjwgtsum[cnvtxs] -= cadjwgt[j];
        cadjncy[j] = cadjncy[--nedges];
        cadjwgt[j] = cadjwgt[nedges];
        htable[cnvtxs] = -1;
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d\n", cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt)));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]] = -1;  /* Zero out the htable */

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, cnvtxs);
}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph_NVW(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, jj, k, kk, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask;
  idxtype *xadj, *adjncy, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, 0);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore;
  memcpy(auxadj, adjncy, iend*sizeof(idxtype));
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  mask = HTLENGTH;
  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1));

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs)
      continue;

    u = match[v];
    cvwgt[cnvtxs] = 1;
    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = 1;
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m]++;
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj]++;
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = 1;
        }
      }
    }

    if (v != u) {
      cvwgt[cnvtxs]++;
      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = 1;
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m]++;
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj]++;
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = 1;
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs)
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) { /* This 2nd check is needed for non-adjacent matchings */
        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, 0);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}


/*************************************************************************
* Setup the various arrays for the coarse graph
**************************************************************************/
GraphType *SetUpCoarseGraph(GraphType *graph, int cnvtxs, int dovsize)
{
  GraphType *cgraph;

  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->ncon = graph->ncon;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /* Allocate memory for the coarser graph */
  if (graph->ncon == 1) {
    if (dovsize) {
      cgraph->gdata = idxmalloc(5*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->vwgt 		= cgraph->gdata + cnvtxs+1;
      cgraph->vsize 		= cgraph->gdata + 2*cnvtxs+1;
      cgraph->adjwgtsum 	= cgraph->gdata + 3*cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 4*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 5*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 5*cnvtxs+1 + graph->nedges;
    }
    else {
      cgraph->gdata = idxmalloc(4*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->vwgt 		= cgraph->gdata + cnvtxs+1;
      cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 3*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 4*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 4*cnvtxs+1 + graph->nedges;
    }
  }
  else {
    if (dovsize) {
      cgraph->gdata = idxmalloc(4*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->vsize 		= cgraph->gdata + cnvtxs+1;
      cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 3*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 4*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 4*cnvtxs+1 + graph->nedges;
    }
    else {
      cgraph->gdata = idxmalloc(3*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->adjwgtsum 	= cgraph->gdata + cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 2*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 3*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 3*cnvtxs+1 + graph->nedges;
    }

    cgraph->nvwgt 	= fmalloc(graph->ncon*cnvtxs, "SetUpCoarseGraph: nvwgt");
  }

  return cgraph;
}


/*************************************************************************
* This function re-adjusts the amount of memory that was allocated if
* it will lead to significant savings
**************************************************************************/
void ReAdjustMemory(GraphType *graph, GraphType *cgraph, int dovsize)
{

  if (cgraph->nedges > 100000 && graph->nedges < 0.7*graph->nedges) {
    idxcopy(cgraph->nedges, cgraph->adjwgt, cgraph->adjncy+cgraph->nedges);

    if (graph->ncon == 1) {
      if (dovsize) {
        cgraph->gdata = realloc(cgraph->gdata, (5*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 		= cgraph->gdata;
        cgraph->vwgt 		= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->vsize 		= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->adjwgtsum	= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->cmap 		= cgraph->gdata + 4*cgraph->nvtxs+1;
        cgraph->adjncy 		= cgraph->gdata + 5*cgraph->nvtxs+1;
        cgraph->adjwgt 		= cgraph->gdata + 5*cgraph->nvtxs+1 + cgraph->nedges;
      }
      else {
        cgraph->gdata = realloc(cgraph->gdata, (4*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 	= cgraph->gdata;
        cgraph->vwgt 	= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->adjwgtsum	= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->cmap 	= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->adjncy 	= cgraph->gdata + 4*cgraph->nvtxs+1;
        cgraph->adjwgt 	= cgraph->gdata + 4*cgraph->nvtxs+1 + cgraph->nedges;
      }
    }
    else {
      if (dovsize) {
        cgraph->gdata = realloc(cgraph->gdata, (4*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 		= cgraph->gdata;
        cgraph->vsize		= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->adjwgtsum	= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->cmap 		= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->adjncy 		= cgraph->gdata + 4*cgraph->nvtxs+1;
        cgraph->adjwgt 		= cgraph->gdata + 4*cgraph->nvtxs+1 + cgraph->nedges;
      }
      else {
        cgraph->gdata = realloc(cgraph->gdata, (3*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 		= cgraph->gdata;
        cgraph->adjwgtsum	= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->cmap 		= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->adjncy 		= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->adjwgt 		= cgraph->gdata + 3*cgraph->nvtxs+1 + cgraph->nedges;
      }
    }
  }

}
/*
 * coarsen.c
 *
 * This file contains the driving routines for the coarsening process
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *Coarsen2Way(CtrlType *ctrl, GraphType *graph)
{
  int clevel;
  GraphType *cgraph;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->CoarsenTmr));

  cgraph = graph;

  /* The following is ahack to allow the multiple bisections to go through with correct
     coarsening */
  if (ctrl->CType > 20) {
    clevel = 1;
    ctrl->CType -= 20;
  }
  else
    clevel = 0;

  do {
    IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d [%d] [%d %d]\n",
          cgraph->nvtxs, cgraph->nedges, ctrl->CoarsenTo, ctrl->maxvwgt,
          (cgraph->vwgt ? idxsum(cgraph->nvtxs, cgraph->vwgt) : cgraph->nvtxs)));

    if (cgraph->adjwgt) {
      switch (ctrl->CType) {
        case MATCH_RM:
          Match_RM(ctrl, cgraph);
          break;
        case MATCH_HEM:
          if (clevel < 1)
            Match_RM(ctrl, cgraph);
          else
            Match_HEM(ctrl, cgraph);
          break;
        case MATCH_SHEM:
          if (clevel < 1)
            Match_RM(ctrl, cgraph);
          else
            Match_SHEM(ctrl, cgraph);
          break;
        case MATCH_SHEMKWAY:
          Match_SHEM(ctrl, cgraph);
          break;
        default:
          errexit("Unknown CType: %d\n", ctrl->CType);
      }
    }
    else {
      Match_RM_NVW(ctrl, cgraph);
    }

    cgraph = cgraph->coarser;
    clevel++;

  } while (cgraph->nvtxs > ctrl->CoarsenTo && cgraph->nvtxs < COARSEN_FRACTION2*cgraph->finer->nvtxs && cgraph->nedges > cgraph->nvtxs/2);

  IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d [%d] [%d %d]\n",
        cgraph->nvtxs, cgraph->nedges, ctrl->CoarsenTo, ctrl->maxvwgt,
        (cgraph->vwgt ? idxsum(cgraph->nvtxs, cgraph->vwgt) : cgraph->nvtxs)));

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->CoarsenTmr));

  return cgraph;
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * compress.c
 *
 * This file contains code for compressing nodes with identical adjacency
 * structure and for prunning dense columns
 *
 * Started 9/17/97
 * George
 *
 * $Id$
 */



/*************************************************************************
* This function compresses a graph by merging identical vertices
* The compression should lead to at least 10% reduction.
**************************************************************************/
void CompressGraph(CtrlType *ctrl, GraphType *graph, int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *cptr, idxtype *cind)
{
  int i, ii, iii, j, jj, k, l, cnvtxs, cnedges;
  idxtype *cxadj, *cadjncy, *cvwgt, *mark, *map;
  KeyValueType *keys;

  mark = idxsmalloc(nvtxs, -1, "CompressGraph: mark");
  map = idxsmalloc(nvtxs, -1, "CompressGraph: map");
  keys = (KeyValueType *)GKmalloc(nvtxs*sizeof(KeyValueType), "CompressGraph: keys");

  /* Compute a key for each adjacency list */
  for (i=0; i<nvtxs; i++) {
    k = 0;
    for (j=xadj[i]; j<xadj[i+1]; j++)
      k += adjncy[j];
    keys[i].key = k+i; /* Add the diagonal entry as well */
    keys[i].val = i;
  }

  ikeysort(nvtxs, keys);

  l = cptr[0] = 0;
  for (cnvtxs=i=0; i<nvtxs; i++) {
    ii = keys[i].val;
    if (map[ii] == -1) {
      mark[ii] = i;  /* Add the diagonal entry */
      for (j=xadj[ii]; j<xadj[ii+1]; j++)
        mark[adjncy[j]] = i;

      cind[l++] = ii;
      map[ii] = cnvtxs;

      for (j=i+1; j<nvtxs; j++) {
        iii = keys[j].val;

        if (keys[i].key != keys[j].key || xadj[ii+1]-xadj[ii] != xadj[iii+1]-xadj[iii])
          break; /* Break if keys or degrees are different */

        if (map[iii] == -1) { /* Do a comparison if iii has not been mapped */
          for (jj=xadj[iii]; jj<xadj[iii+1]; jj++) {
            if (mark[adjncy[jj]] != i)
              break;
          }

          if (jj == xadj[iii+1]) { /* Identical adjacency structure */
            map[iii] = cnvtxs;
            cind[l++] = iii;
          }
        }
      }

      cptr[++cnvtxs] = l;
    }
  }

  /* printf("Original: %6d, Compressed: %6d\n", nvtxs, cnvtxs); */


  InitGraph(graph);

  if (cnvtxs >= COMPRESSION_FRACTION*nvtxs) {
    graph->nvtxs = nvtxs;
    graph->nedges = xadj[nvtxs];
    graph->ncon = 1;
    graph->xadj = xadj;
    graph->adjncy = adjncy;

    graph->gdata = idxmalloc(3*nvtxs+graph->nedges, "CompressGraph: gdata");
    graph->vwgt    	= graph->gdata;
    graph->adjwgtsum    = graph->gdata+nvtxs;
    graph->cmap		= graph->gdata+2*nvtxs;
    graph->adjwgt	= graph->gdata+3*nvtxs;

    idxset(nvtxs, 1, graph->vwgt);
    idxset(graph->nedges, 1, graph->adjwgt);
    for (i=0; i<nvtxs; i++)
      graph->adjwgtsum[i] = xadj[i+1]-xadj[i];

    graph->label = idxmalloc(nvtxs, "CompressGraph: label");
    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }
  else { /* Ok, form the compressed graph  */
    cnedges = 0;
    for (i=0; i<cnvtxs; i++) {
      ii = cind[cptr[i]];
      cnedges += xadj[ii+1]-xadj[ii];
    }

    /* Allocate memory for the compressed graph*/
    graph->gdata = idxmalloc(4*cnvtxs+1 + 2*cnedges, "CompressGraph: gdata");
    cxadj = graph->xadj		= graph->gdata;
    cvwgt = graph->vwgt         = graph->gdata + cnvtxs+1;
    graph->adjwgtsum        	= graph->gdata + 2*cnvtxs+1;
    graph->cmap                 = graph->gdata + 3*cnvtxs+1;
    cadjncy = graph->adjncy     = graph->gdata + 4*cnvtxs+1;
    graph->adjwgt            	= graph->gdata + 4*cnvtxs+1 + cnedges;

    /* Now go and compress the graph */
    idxset(nvtxs, -1, mark);
    l = cxadj[0] = 0;
    for (i=0; i<cnvtxs; i++) {
      cvwgt[i] = cptr[i+1]-cptr[i];
      mark[i] = i;  /* Remove any dioganal entries in the compressed graph */
      for (j=cptr[i]; j<cptr[i+1]; j++) {
        ii = cind[j];
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
          k = map[adjncy[jj]];
          if (mark[k] != i)
            cadjncy[l++] = k;
          mark[k] = i;
        }
      }
      cxadj[i+1] = l;
    }

    graph->nvtxs = cnvtxs;
    graph->nedges = l;
    graph->ncon = 1;

    idxset(graph->nedges, 1, graph->adjwgt);
    for (i=0; i<cnvtxs; i++)
      graph->adjwgtsum[i] = cxadj[i+1]-cxadj[i];

    graph->label = idxmalloc(cnvtxs, "CompressGraph: label");
    for (i=0; i<cnvtxs; i++)
      graph->label[i] = i;

  }

  GKfree((void **) &keys, (void **) &map, (void **) &mark, LTERM);
}



/*************************************************************************
* This function prunes all the vertices in a graph with degree greater
* than factor*average
**************************************************************************/
void PruneGraph(CtrlType *ctrl, GraphType *graph, int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *iperm, float factor)
{
  int i, j, k, l, nlarge, pnvtxs, pnedges;
  idxtype *pxadj, *padjncy;
  idxtype *perm;

  perm = idxmalloc(nvtxs, "PruneGraph: perm");

  factor = factor*xadj[nvtxs]/nvtxs;

  pnvtxs = pnedges = nlarge = 0;
  for (i=0; i<nvtxs; i++) {
    if (xadj[i+1]-xadj[i] < factor) {
      perm[i] = pnvtxs;
      iperm[pnvtxs++] = i;
      pnedges += xadj[i+1]-xadj[i];
    }
    else {
      perm[i] = nvtxs - ++nlarge;
      iperm[nvtxs-nlarge] = i;
    }
  }

  /* printf("Pruned %d vertices\n", nlarge); */

  InitGraph(graph);

  if (nlarge == 0) { /* No prunning */
    graph->nvtxs = nvtxs;
    graph->nedges = xadj[nvtxs];
    graph->ncon = 1;
    graph->xadj = xadj;
    graph->adjncy = adjncy;

    graph->gdata = idxmalloc(3*nvtxs+graph->nedges, "CompressGraph: gdata");
    graph->vwgt    	= graph->gdata;
    graph->adjwgtsum    = graph->gdata+nvtxs;
    graph->cmap		= graph->gdata+2*nvtxs;
    graph->adjwgt	= graph->gdata+3*nvtxs;

    idxset(nvtxs, 1, graph->vwgt);
    idxset(graph->nedges, 1, graph->adjwgt);
    for (i=0; i<nvtxs; i++)
      graph->adjwgtsum[i] = xadj[i+1]-xadj[i];

    graph->label = idxmalloc(nvtxs, "CompressGraph: label");
    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }
  else { /* Prune the graph */
    /* Allocate memory for the compressed graph*/
    graph->gdata = idxmalloc(4*pnvtxs+1 + 2*pnedges, "PruneGraph: gdata");
    pxadj = graph->xadj		= graph->gdata;
    graph->vwgt         	= graph->gdata + pnvtxs+1;
    graph->adjwgtsum        	= graph->gdata + 2*pnvtxs+1;
    graph->cmap                 = graph->gdata + 3*pnvtxs+1;
    padjncy = graph->adjncy     = graph->gdata + 4*pnvtxs+1;
    graph->adjwgt            	= graph->gdata + 4*pnvtxs+1 + pnedges;

    pxadj[0] = pnedges = l = 0;
    for (i=0; i<nvtxs; i++) {
      if (xadj[i+1]-xadj[i] < factor) {
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = perm[adjncy[j]];
          if (k < pnvtxs)
            padjncy[pnedges++] = k;
        }
        pxadj[++l] = pnedges;
      }
    }

    graph->nvtxs = pnvtxs;
    graph->nedges = pnedges;
    graph->ncon = 1;

    idxset(pnvtxs, 1, graph->vwgt);
    idxset(pnedges, 1, graph->adjwgt);
    for (i=0; i<pnvtxs; i++)
      graph->adjwgtsum[i] = pxadj[i+1]-pxadj[i];

    graph->label = idxmalloc(pnvtxs, "CompressGraph: label");
    for (i=0; i<pnvtxs; i++)
      graph->label[i] = i;
  }

  free(perm);

}









/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * debug.c
 *
 * This file contains code that performs self debuging
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function computes the cut given the graph and a where vector
**************************************************************************/
int ComputeCut(GraphType *graph, idxtype *where)
{
  int i, j, cut;

  if (graph->adjwgt == NULL) {
    for (cut=0, i=0; i<graph->nvtxs; i++) {
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
        if (where[i] != where[graph->adjncy[j]])
          cut++;
    }
  }
  else {
    for (cut=0, i=0; i<graph->nvtxs; i++) {
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
        if (where[i] != where[graph->adjncy[j]])
          cut += graph->adjwgt[j];
    }
  }

  return cut/2;
}


/*************************************************************************
* This function checks whether or not the boundary information is correct
**************************************************************************/
int CheckBnd(GraphType *graph)
{
  int i, j, nvtxs, nbnd;
  idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  for (nbnd=0, i=0; i<nvtxs; i++) {
    if (xadj[i+1]-xadj[i] == 0)
      nbnd++;   /* Islands are considered to be boundary vertices */

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        nbnd++;
        ASSERT(bndptr[i] != -1);
        ASSERT(bndind[bndptr[i]] == i);
        break;
      }
    }
  }

  ASSERTP(nbnd == graph->nbnd, ("%d %d\n", nbnd, graph->nbnd));

  return 1;
}



/*************************************************************************
* This function checks whether or not the boundary information is correct
**************************************************************************/
int CheckBnd2(GraphType *graph)
{
  int i, j, nvtxs, nbnd, id, ed;
  idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  for (nbnd=0, i=0; i<nvtxs; i++) {
    id = ed = 0;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]])
        ed += graph->adjwgt[j];
      else
        id += graph->adjwgt[j];
    }
    if (ed - id >= 0 && xadj[i] < xadj[i+1]) {
      nbnd++;
      ASSERTP(bndptr[i] != -1, ("%d %d %d\n", i, id, ed));
      ASSERT(bndind[bndptr[i]] == i);
    }
  }

  ASSERTP(nbnd == graph->nbnd, ("%d %d\n", nbnd, graph->nbnd));

  return 1;
}

/*************************************************************************
* This function checks whether or not the boundary information is correct
**************************************************************************/
int CheckNodeBnd(GraphType *graph, int onbnd)
{
  int i, nvtxs, nbnd;
  idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  for (nbnd=0, i=0; i<nvtxs; i++) {
    if (where[i] == 2)
      nbnd++;
  }

  ASSERTP(nbnd == onbnd, ("%d %d\n", nbnd, onbnd));

  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2) {
      ASSERTP(bndptr[i] == -1, ("%d %d\n", i, bndptr[i]));
    }
    else {
      ASSERTP(bndptr[i] != -1, ("%d %d\n", i, bndptr[i]));
    }
  }

  return 1;
}



/*************************************************************************
* This function checks whether or not the rinfo of a vertex is consistent
**************************************************************************/
int CheckRInfo(RInfoType *rinfo)
{
  int i, j;

  for (i=0; i<rinfo->ndegrees; i++) {
    for (j=i+1; j<rinfo->ndegrees; j++)
      ASSERTP(rinfo->edegrees[i].pid != rinfo->edegrees[j].pid, ("%d %d %d %d\n", i, j, rinfo->edegrees[i].pid, rinfo->edegrees[j].pid));
  }

  return 1;
}



/*************************************************************************
* This function checks the correctness of the NodeFM data structures
**************************************************************************/
int CheckNodePartitionParams(GraphType *graph)
{
  int i, j, nvtxs, me, other;
  idxtype *xadj, *adjncy, *adjwgt, *vwgt, *where;
  idxtype edegrees[2], pwgts[3];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;

  /*------------------------------------------------------------
  / Compute now the separator external degrees
  /------------------------------------------------------------*/
  pwgts[0] = pwgts[1] = pwgts[2] = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];

    if (me == 2) { /* If it is on the separator do some computations */
      edegrees[0] = edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (other != 2)
          edegrees[other] += vwgt[adjncy[j]];
      }
      if (edegrees[0] != graph->nrinfo[i].edegrees[0] || edegrees[1] != graph->nrinfo[i].edegrees[1]) {
        printf("Something wrong with edegrees: %d %d %d %d %d\n", i, edegrees[0], edegrees[1], graph->nrinfo[i].edegrees[0], graph->nrinfo[i].edegrees[1]);
        return 0;
      }
    }
  }

  if (pwgts[0] != graph->pwgts[0] || pwgts[1] != graph->pwgts[1] || pwgts[2] != graph->pwgts[2])
    printf("Something wrong with part-weights: %d %d %d %d %d %d\n", pwgts[0], pwgts[1], pwgts[2], graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]);

  return 1;
}


/*************************************************************************
* This function checks if the separator is indeed a separator
**************************************************************************/
int IsSeparable(GraphType *graph)
{
  int i, j, nvtxs, other;
  idxtype *xadj, *adjncy, *where;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  for (i=0; i<nvtxs; i++) {
    if (where[i] == 2)
      continue;
    other = (where[i]+1)%2;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ASSERTP(where[adjncy[j]] != other, ("%d %d %d %d %d %d\n", i, where[i], adjncy[j], where[adjncy[j]], xadj[i+1]-xadj[i], xadj[adjncy[j]+1]-xadj[adjncy[j]]));
    }
  }

  return 1;
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * estmem.c
 *
 * This file contains code for estimating the amount of memory required by
 * the various routines in METIS
 *
 * Started 11/4/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function computes how much memory will be required by the various
* routines in METIS
**************************************************************************/
void METIS_EstimateMemory(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *optype, int *nbytes)
{
  int nedges, nlevels;
  float vfraction, efraction, vmult, emult;
  int coresize, gdata, rdata;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  nedges = xadj[*nvtxs];

  InitRandom(-1);
  EstimateCFraction(*nvtxs, xadj, adjncy, &vfraction, &efraction);

  /* Estimate the amount of memory for coresize */
  if (*optype == 2)
    coresize = nedges;
  else
    coresize = 0;
  coresize += nedges + 11*(*nvtxs) + 4*1024 + 2*(NEG_GAINSPAN+PLUS_GAINSPAN+1)*(sizeof(ListNodeType *)/sizeof(idxtype));
  coresize += 2*(*nvtxs);  /* add some more fore other vectors */

  gdata = nedges;   /* Assume that the user does not pass weights */

  nlevels = (int)(log(100.0/(*nvtxs))/log(vfraction) + .5);
  vmult = 0.5 + (1.0 - pow(vfraction, nlevels))/(1.0 - vfraction);
  emult = 1.0 + (1.0 - pow(efraction, nlevels+1))/(1.0 - efraction);

  gdata += vmult*4*(*nvtxs) + emult*2*nedges;
  if ((vmult-1.0)*4*(*nvtxs) + (emult-1.0)*2*nedges < 5*(*nvtxs))
    rdata = 0;
  else
    rdata = 5*(*nvtxs);

  *nbytes = sizeof(idxtype)*(coresize+gdata+rdata+(*nvtxs));

  if (*numflag == 1)
    Change2FNumbering2(*nvtxs, xadj, adjncy);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void EstimateCFraction(int nvtxs, idxtype *xadj, idxtype *adjncy, float *vfraction, float *efraction)
{
  int i, ii, j, cnvtxs, cnedges, maxidx;
  idxtype *match, *cmap, *perm;

  cmap = idxmalloc(nvtxs, "cmap");
  match = idxsmalloc(nvtxs, UNMATCHED, "match");
  perm = idxmalloc(nvtxs, "perm");
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED) {
          maxidx = adjncy[j];
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  cnedges = ComputeCoarseGraphSize(nvtxs, xadj, adjncy, cnvtxs, cmap, match, perm);

  *vfraction = (1.0*cnvtxs)/(1.0*nvtxs);
  *efraction = (1.0*cnedges)/(1.0*xadj[nvtxs]);

  GKfree((void **) &cmap, (void **) &match, (void **) &perm, LTERM);
}




/*************************************************************************
* This function computes the size of the coarse graph
**************************************************************************/
int ComputeCoarseGraphSize(int nvtxs, idxtype *xadj, idxtype *adjncy, int cnvtxs, idxtype *cmap, idxtype *match, idxtype *perm)
{
  int i, j, k, istart, iend, cnedges, v, u;
  idxtype *htable;

  htable = idxsmalloc(cnvtxs, -1, "htable");

  cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs)
      continue;

    htable[cnvtxs] = cnvtxs;

    u = match[v];

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = cmap[adjncy[j]];
      if (htable[k] != cnvtxs) {
        htable[k] = cnvtxs;
        cnedges++;
      }
    }

    if (v != u) {
      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = cmap[adjncy[j]];
        if (htable[k] != cnvtxs) {
          htable[k] = cnvtxs;
          cnedges++;
        }
      }
    }
    cnvtxs++;
  }

  GKfree((void **) &htable, LTERM);

  return cnedges;
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * fm.c
 *
 * This file contains code that implements the edge-based FM refinement
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void FM_2WayEdgeRefine(CtrlType *ctrl, GraphType *graph, int *tpwgts, int npasses)
{
  int i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, pass, limit, tmp;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind, *pwgts;
  idxtype *moved, *swaps, *perm;
  PQueueType parts[2];
  int higain, oldgain, mincut, mindiff, origdiff, initcut, newcut, mincutorder, avgvwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  pwgts = graph->pwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);

  limit = amin(amax(0.01*nvtxs, 15), 100);
  avgvwgt = amin((pwgts[0]+pwgts[1])/20, 2*(pwgts[0]+pwgts[1])/nvtxs);

  tmp = graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)];
  PQueueInit(ctrl, &parts[0], nvtxs, tmp);
  PQueueInit(ctrl, &parts[1], nvtxs, tmp);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d] T[%6d %6d], Nv-Nb[%6d %6d]. ICut: %6d\n",
             pwgts[0], pwgts[1], tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  origdiff = abs(tpwgts[0]-pwgts[0]);
  idxset(nvtxs, -1, moved);
  for (pass=0; pass<npasses; pass++) { /* Do a number of passes */
    PQueueReset(&parts[0]);
    PQueueReset(&parts[1]);

    mincutorder = -1;
    newcut = mincut = initcut = graph->mincut;
    mindiff = abs(tpwgts[0]-pwgts[0]);

    ASSERT(ComputeCut(graph, where) == graph->mincut);
    ASSERT(CheckBnd(graph));

    /* Insert boundary nodes in the priority queues */
    nbnd = graph->nbnd;
    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = perm[ii];
      ASSERT(ed[bndind[i]] > 0 || id[bndind[i]] == 0);
      ASSERT(bndptr[bndind[i]] != -1);
      PQueueInsert(&parts[where[bndind[i]]], bndind[i], ed[bndind[i]]-id[bndind[i]]);
    }

    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      from = (tpwgts[0]-pwgts[0] < tpwgts[1]-pwgts[1] ? 0 : 1);
      to = (from+1)%2;

      if ((higain = PQueueGetMax(&parts[from])) == -1)
        break;
      ASSERT(bndptr[higain] != -1);

      newcut -= (ed[higain]-id[higain]);
      INC_DEC(pwgts[to], pwgts[from], vwgt[higain]);

      if ((newcut < mincut && abs(tpwgts[0]-pwgts[0]) <= origdiff+avgvwgt) ||
          (newcut == mincut && abs(tpwgts[0]-pwgts[0]) < mindiff)) {
        mincut = newcut;
        mindiff = abs(tpwgts[0]-pwgts[0]);
        mincutorder = nswaps;
      }
      else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
        newcut += (ed[higain]-id[higain]);
        INC_DEC(pwgts[from], pwgts[to], vwgt[higain]);
        break;
      }

      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
        printf("Moved %6d from %d. [%3d %3d] %5d [%4d %4d]\n", higain, from, ed[higain]-id[higain], vwgt[higain], newcut, pwgts[0], pwgts[1]));

      /**************************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***************************************************************/
      SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);

      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        oldgain = ed[k]-id[k];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        /* Update its boundary information and queue position */
        if (bndptr[k] != -1) { /* If k was a boundary vertex */
          if (ed[k] == 0) { /* Not a boundary vertex any more */
            BNDDelete(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)  /* Remove it if in the queues */
              PQueueDelete(&parts[where[k]], k, oldgain);
          }
          else { /* If it has not been moved, update its position in the queue */
            if (moved[k] == -1)
              PQueueUpdate(&parts[where[k]], k, oldgain, ed[k]-id[k]);
          }
        }
        else {
          if (ed[k] > 0) {  /* It will now become a boundary vertex */
            BNDInsert(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)
              PQueueInsert(&parts[where[k]], k, ed[k]-id[k]);
          }
        }
      }

    }


    /****************************************************************
    * Roll back computations
    *****************************************************************/
    for (i=0; i<nswaps; i++)
      moved[swaps[i]] = -1;  /* reset moved array */
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      to = where[higain] = (where[higain]+1)%2;
      SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);
      else if (ed[higain] > 0 && bndptr[higain] == -1)
        BNDInsert(nbnd, bndind,  bndptr, higain);

      INC_DEC(pwgts[to], pwgts[(to+1)%2], vwgt[higain]);
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        if (bndptr[k] != -1 && ed[k] == 0)
          BNDDelete(nbnd, bndind, bndptr, k);
        if (bndptr[k] == -1 && ed[k] > 0)
          BNDInsert(nbnd, bndind, bndptr, k);
      }
    }

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum cut: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, mincutorder, pwgts[0], pwgts[1], nbnd));

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut == initcut)
      break;
  }

  PQueueFree(ctrl, &parts[0]);
  PQueueFree(ctrl, &parts[1]);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * fortran.c
 *
 * This file contains code for the fortran to C interface
 *
 * Started 8/19/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function changes the numbering to start from 0 instead of 1
**************************************************************************/
void Change2CNumbering(int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, nedges;

  for (i=0; i<=nvtxs; i++)
    xadj[i]--;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]--;
}

/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumbering(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vector)
{
  int i, nedges;

  for (i=0; i<nvtxs; i++)
    vector[i]++;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;
}

/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumbering2(int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, nedges;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;
}



/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumberingOrder(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *v1, idxtype *v2)
{
  int i, nedges;

  for (i=0; i<nvtxs; i++) {
    v1[i]++;
    v2[i]++;
  }

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}



/*************************************************************************
* This function changes the numbering to start from 0 instead of 1
**************************************************************************/
void ChangeMesh2CNumbering(int n, idxtype *mesh)
{
  int i;

  for (i=0; i<n; i++)
    mesh[i]--;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering(int n, idxtype *mesh, int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, nedges;

  for (i=0; i<n; i++)
    mesh[i]++;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering2(int n, idxtype *mesh, int ne, int nn, idxtype *epart, idxtype *npart)
{
  int i;

  for (i=0; i<n; i++)
    mesh[i]++;

  for (i=0; i<ne; i++)
    epart[i]++;

  for (i=0; i<nn; i++)
    npart[i]++;

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * frename.c
 *
 * This file contains some renaming routines to deal with different Fortran compilers
 *
 * Started 9/15/97
 * George
 *
 * $Id$
 *
 */




void METIS_PARTGRAPHRECURSIVE(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}


void METIS_WPARTGRAPHRECURSIVE(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}



void METIS_PARTGRAPHKWAY(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}



void METIS_WPARTGRAPHKWAY(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}



void METIS_EDGEND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend_(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend__(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}



void METIS_NODEND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend_(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend__(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}



void METIS_NODEWND(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}



void METIS_PARTMESHNODAL(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal_(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal__(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}


void METIS_PARTMESHDUAL(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshdual(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshdual_(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshdual__(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}


void METIS_MESHTONODAL(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal_(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal__(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}


void METIS_MESHTODUAL(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual_(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual__(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}


void METIS_ESTIMATEMEMORY(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *optype, int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *optype, int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory_(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *optype, int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory__(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *optype, int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}



void METIS_MCPARTGRAPHRECURSIVE(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive_(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive__(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}


void METIS_MCPARTGRAPHKWAY(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *rubvec, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *rubvec, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway_(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *rubvec, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway__(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *rubvec, int *options, int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}


void METIS_PARTGRAPHVKWAY(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, int *options, int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, int *options, int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkway_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, int *options, int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkway__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, int *options, int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}

void METIS_WPARTGRAPHVKWAY(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkway_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkway__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag, int *numflag, int *nparts, float *tpwgts, int *options, int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}



/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * graph.c
 *
 * This file contains functions that deal with setting up the graphs
 * for METIS.
 *
 * Started 7/25/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void SetUpGraph(GraphType *graph, int OpType, int nvtxs, int ncon,
       idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int wgtflag)
{
  int i, j, sum, gsize;
  float *nvwgt;
  idxtype tvwgt[MAXNCON];

  if (OpType == OP_KMETIS && ncon == 1 && (wgtflag&2) == 0 && (wgtflag&1) == 0) {
    SetUpGraphKway(graph, nvtxs, xadj, adjncy);
    return;
  }

  InitGraph(graph);

  graph->nvtxs = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon = ncon;
  graph->xadj = xadj;
  graph->adjncy = adjncy;

  if (ncon == 1) { /* We are in the non mC mode */
    gsize = 0;
    if ((wgtflag&2) == 0)
      gsize += nvtxs;
    if ((wgtflag&1) == 0)
      gsize += graph->nedges;

    gsize += 2*nvtxs;

    graph->gdata = idxmalloc(gsize, "SetUpGraph: gdata");

    /* Create the vertex/edge weight vectors if they are not supplied */
    gsize = 0;
    if ((wgtflag&2) == 0) {
      vwgt = graph->vwgt = idxset(nvtxs, 1, graph->gdata);
      gsize += nvtxs;
    }
    else
      graph->vwgt = vwgt;

    if ((wgtflag&1) == 0) {
      adjwgt = graph->adjwgt = idxset(graph->nedges, 1, graph->gdata+gsize);
      gsize += graph->nedges;
    }
    else
      graph->adjwgt = adjwgt;


    /* Compute the initial values of the adjwgtsum */
    graph->adjwgtsum = graph->gdata + gsize;
    gsize += nvtxs;

    for (i=0; i<nvtxs; i++) {
      sum = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++)
        sum += adjwgt[j];
      graph->adjwgtsum[i] = sum;
    }

    graph->cmap = graph->gdata + gsize;
    gsize += nvtxs;

  }
  else {  /* Set up the graph in MOC mode */
    gsize = 0;
    if ((wgtflag&1) == 0)
      gsize += graph->nedges;

    gsize += 2*nvtxs;

    graph->gdata = idxmalloc(gsize, "SetUpGraph: gdata");
    gsize = 0;

    for (i=0; i<ncon; i++)
      tvwgt[i] = idxsum_strd(nvtxs, vwgt+i, ncon);

    nvwgt = graph->nvwgt = fmalloc(ncon*nvtxs, "SetUpGraph: nvwgt");

    for (i=0; i<nvtxs; i++) {
      for (j=0; j<ncon; j++)
        nvwgt[i*ncon+j] = (1.0*vwgt[i*ncon+j])/(1.0*tvwgt[j]);
    }


    /* Create the edge weight vectors if they are not supplied */
    if ((wgtflag&1) == 0) {
      adjwgt = graph->adjwgt = idxset(graph->nedges, 1, graph->gdata+gsize);
      gsize += graph->nedges;
    }
    else
      graph->adjwgt = adjwgt;

    /* Compute the initial values of the adjwgtsum */
    graph->adjwgtsum = graph->gdata + gsize;
    gsize += nvtxs;

    for (i=0; i<nvtxs; i++) {
      sum = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++)
        sum += adjwgt[j];
      graph->adjwgtsum[i] = sum;
    }

    graph->cmap = graph->gdata + gsize;
    gsize += nvtxs;

  }

  if (OpType != OP_KMETIS && OpType != OP_KVMETIS) {
    graph->label = idxmalloc(nvtxs, "SetUpGraph: label");

    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }

}


/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void SetUpGraphKway(GraphType *graph, int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i;

  InitGraph(graph);

  graph->nvtxs = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon = 1;
  graph->xadj = xadj;
  graph->vwgt = NULL;
  graph->adjncy = adjncy;
  graph->adjwgt = NULL;

  graph->gdata = idxmalloc(2*nvtxs, "SetUpGraph: gdata");
  graph->adjwgtsum = graph->gdata;
  graph->cmap = graph->gdata + nvtxs;

  /* Compute the initial values of the adjwgtsum */
  for (i=0; i<nvtxs; i++)
    graph->adjwgtsum[i] = xadj[i+1]-xadj[i];

}



/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void SetUpGraph2(GraphType *graph, int nvtxs, int ncon, idxtype *xadj,
       idxtype *adjncy, float *nvwgt, idxtype *adjwgt)
{
  int i, j, sum;

  InitGraph(graph);

  graph->nvtxs = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon = ncon;
  graph->xadj = xadj;
  graph->adjncy = adjncy;
  graph->adjwgt = adjwgt;

  graph->nvwgt = fmalloc(nvtxs*ncon, "SetUpGraph2: graph->nvwgt");
  scopy(nvtxs*ncon, nvwgt, graph->nvwgt);

  graph->gdata = idxmalloc(2*nvtxs, "SetUpGraph: gdata");

  /* Compute the initial values of the adjwgtsum */
  graph->adjwgtsum = graph->gdata;
  for (i=0; i<nvtxs; i++) {
    sum = 0;
    for (j=xadj[i]; j<xadj[i+1]; j++)
      sum += adjwgt[j];
    graph->adjwgtsum[i] = sum;
  }

  graph->cmap = graph->gdata+nvtxs;

  graph->label = idxmalloc(nvtxs, "SetUpGraph: label");
  for (i=0; i<nvtxs; i++)
    graph->label[i] = i;

}


/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void VolSetUpGraph(GraphType *graph, int OpType, int nvtxs, int ncon, idxtype *xadj,
                   idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int wgtflag)
{
  int i, j, sum, gsize;
  idxtype *adjwgt;
  float *nvwgt;
  idxtype tvwgt[MAXNCON];

  InitGraph(graph);

  graph->nvtxs = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon = ncon;
  graph->xadj = xadj;
  graph->adjncy = adjncy;

  if (ncon == 1) { /* We are in the non mC mode */
    gsize = graph->nedges;  /* This is for the edge weights */
    if ((wgtflag&2) == 0)
      gsize += nvtxs; /* vwgts */
    if ((wgtflag&1) == 0)
      gsize += nvtxs; /* vsize */

    gsize += 2*nvtxs;

    graph->gdata = idxmalloc(gsize, "SetUpGraph: gdata");

    /* Create the vertex/edge weight vectors if they are not supplied */
    gsize = 0;
    if ((wgtflag&2) == 0) {
      vwgt = graph->vwgt = idxset(nvtxs, 1, graph->gdata);
      gsize += nvtxs;
    }
    else
      graph->vwgt = vwgt;

    if ((wgtflag&1) == 0) {
      vsize = graph->vsize = idxset(nvtxs, 1, graph->gdata);
      gsize += nvtxs;
    }
    else
      graph->vsize = vsize;

    /* Allocate memory for edge weights and initialize them to the sum of the vsize */
    adjwgt = graph->adjwgt = graph->gdata+gsize;
    gsize += graph->nedges;

    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        adjwgt[j] = 1+vsize[i]+vsize[adjncy[j]];
    }


    /* Compute the initial values of the adjwgtsum */
    graph->adjwgtsum = graph->gdata + gsize;
    gsize += nvtxs;

    for (i=0; i<nvtxs; i++) {
      sum = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++)
        sum += adjwgt[j];
      graph->adjwgtsum[i] = sum;
    }

    graph->cmap = graph->gdata + gsize;
    gsize += nvtxs;

  }
  else {  /* Set up the graph in MOC mode */
    gsize = graph->nedges;
    if ((wgtflag&1) == 0)
      gsize += nvtxs;

    gsize += 2*nvtxs;

    graph->gdata = idxmalloc(gsize, "SetUpGraph: gdata");
    gsize = 0;

    /* Create the normalized vertex weights along each constrain */
    if ((wgtflag&2) == 0)
      vwgt = idxsmalloc(nvtxs, 1, "SetUpGraph: vwgt");

    for (i=0; i<ncon; i++)
      tvwgt[i] = idxsum_strd(nvtxs, vwgt+i, ncon);

    nvwgt = graph->nvwgt = fmalloc(ncon*nvtxs, "SetUpGraph: nvwgt");

    for (i=0; i<nvtxs; i++) {
      for (j=0; j<ncon; j++)
        nvwgt[i*ncon+j] = (1.0*vwgt[i*ncon+j])/(1.0*tvwgt[j]);
    }
    if ((wgtflag&2) == 0)
      free(vwgt);


    /* Create the vsize vector if it is not supplied */
    if ((wgtflag&1) == 0) {
      vsize = graph->vsize = idxset(nvtxs, 1, graph->gdata);
      gsize += nvtxs;
    }
    else
      graph->vsize = vsize;

    /* Allocate memory for edge weights and initialize them to the sum of the vsize */
    adjwgt = graph->adjwgt = graph->gdata+gsize;
    gsize += graph->nedges;

    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        adjwgt[j] = 1+vsize[i]+vsize[adjncy[j]];
    }

    /* Compute the initial values of the adjwgtsum */
    graph->adjwgtsum = graph->gdata + gsize;
    gsize += nvtxs;

    for (i=0; i<nvtxs; i++) {
      sum = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++)
        sum += adjwgt[j];
      graph->adjwgtsum[i] = sum;
    }

    graph->cmap = graph->gdata + gsize;
    gsize += nvtxs;

  }

  if (OpType != OP_KVMETIS) {
    graph->label = idxmalloc(nvtxs, "SetUpGraph: label");

    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }

}


/*************************************************************************
* This function randomly permutes the adjacency lists of a graph
**************************************************************************/
void RandomizeGraph(GraphType *graph)
{
  int i, j, k, l, tmp, nvtxs;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  for (i=0; i<nvtxs; i++) {
    l = xadj[i+1]-xadj[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = xadj[i] + RandomInRange(l);
      SWAP(adjncy[j], adjncy[k], tmp);
      SWAP(adjwgt[j], adjwgt[k], tmp);
    }
  }
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
int IsConnectedSubdomain(CtrlType *ctrl, GraphType *graph, int pid, int report)
{
  int i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;
  idxtype *cptr;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");
  cptr = idxmalloc(nvtxs, "IsConnected: cptr");

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] == pid)
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] == pid)
      break;
  }

  touched[i] = 1;
  queue[0] = i;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (where[i] == pid && !touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] == pid && !touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  if (ncmps > 1 && report) {
    printf("The graph has %d connected components in partition %d:\t", ncmps, pid);
    for (i=0; i<ncmps; i++) {
      wgt = 0;
      for (j=cptr[i]; j<cptr[i+1]; j++)
        wgt += graph->vwgt[queue[j]];
      printf("[%5d %5d] ", cptr[i+1]-cptr[i], wgt);
      /*
      if (cptr[i+1]-cptr[i] == 1)
        printf("[%d %d] ", queue[cptr[i]], xadj[queue[cptr[i]]+1]-xadj[queue[cptr[i]]]);
      */
    }
    printf("\n");
  }

  GKfree((void **) &touched, (void **) &queue, (void **) &cptr, LTERM);

  return (ncmps == 1 ? 1 : 0);
}


/*************************************************************************
* This function checks whether a graph is contigous or not
**************************************************************************/
int IsConnected(CtrlType *ctrl, GraphType *graph, int report)
{
  int i, j, k, nvtxs, first, last;
  idxtype *xadj, *adjncy, *touched, *queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");

  touched[0] = 1;
  queue[0] = 0;
  first = 0; last = 1;

  while (first < last) {
    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }

  if (first != nvtxs && report)
    printf("The graph is not connected. It has %d disconnected vertices!\n", nvtxs-first);

  return (first == nvtxs ? 1 : 0);
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
int IsConnected2(GraphType *graph, int report)
{
  int i, j, k, nvtxs, first, last, nleft, ncmps;
  idxtype *xadj, *adjncy, *where, *touched, *queue;
  idxtype *cptr;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");
  cptr = idxmalloc(nvtxs, "IsConnected: cptr");

  nleft = nvtxs;
  touched[0] = 1;
  queue[0] = 0;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  if (ncmps > 1 && report) {
    printf("%d connected components:\t", ncmps);
    for (i=0; i<ncmps; i++) {
      if (cptr[i+1]-cptr[i] > 200)
        printf("[%5d] ", cptr[i+1]-cptr[i]);
    }
    printf("\n");
  }

  GKfree((void **) &touched, (void **) &queue, (void **) &cptr, LTERM);

  return (ncmps == 1 ? 1 : 0);
}


/*************************************************************************
* This function returns the number of connected components in cptr,cind
* The separator of the graph is used to split it and then find its components.
**************************************************************************/
int FindComponents(CtrlType *ctrl, GraphType *graph, idxtype *cptr, idxtype *cind)
{
  int i, j, k, nvtxs, first, last, nleft, ncmps;
  idxtype *xadj, *adjncy, *where, *touched, *queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: queue");

  for (i=0; i<graph->nbnd; i++)
    touched[graph->bndind[i]] = 1;

  queue = cind;

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2)
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2)
      break;
  }

  touched[i] = 1;
  queue[0] = i;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  free(touched);

  return ncmps;
}



/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that performs the initial partition of the
 * coarsest graph
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function computes the initial bisection of the coarsest graph
**************************************************************************/
void Init2WayPartition(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  int dbglvl;

  dbglvl = ctrl->dbglvl;
  IFSET(ctrl->dbglvl, DBG_REFINE, ctrl->dbglvl -= DBG_REFINE);
  IFSET(ctrl->dbglvl, DBG_MOVEINFO, ctrl->dbglvl -= DBG_MOVEINFO);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  switch (ctrl->IType) {
    case IPART_GGPKL:
      GrowBisection(ctrl, graph, tpwgts, ubfactor);
      break;
    case 3:
      RandomBisection(ctrl, graph, tpwgts, ubfactor);
      break;
    default:
      errexit("Unknown initial partition type: %d\n", ctrl->IType);
  }

  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial Cut: %d\n", graph->mincut));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  ctrl->dbglvl = dbglvl;

/*
  IsConnectedSubdomain(ctrl, graph, 0);
  IsConnectedSubdomain(ctrl, graph, 1);
*/
}

/*************************************************************************
* This function computes the initial bisection of the coarsest graph
**************************************************************************/
void InitSeparator(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  int dbglvl;

  dbglvl = ctrl->dbglvl;
  IFSET(ctrl->dbglvl, DBG_REFINE, ctrl->dbglvl -= DBG_REFINE);
  IFSET(ctrl->dbglvl, DBG_MOVEINFO, ctrl->dbglvl -= DBG_MOVEINFO);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  GrowBisectionNode(ctrl, graph, ubfactor);
  Compute2WayNodePartitionParams(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial Sep: %d\n", graph->mincut));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

  ctrl->dbglvl = dbglvl;

}



/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void GrowBisection(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  int i, j, k, nvtxs, drain, nleft, first, last, pwgts[2], minpwgt[2], maxpwgt[2], bestcut,  nbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where;
  idxtype *queue, *touched, *bestwhere;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  queue = idxmalloc(nvtxs, "BisectGraph: queue");
  touched = idxmalloc(nvtxs, "BisectGraph: touched");

  ASSERTP(tpwgts[0]+tpwgts[1] == idxsum(nvtxs, vwgt), ("%d %d\n", tpwgts[0]+tpwgts[1], idxsum(nvtxs, vwgt)));

  maxpwgt[0] = ubfactor*tpwgts[0];
  maxpwgt[1] = ubfactor*tpwgts[1];
  minpwgt[0] = (1.0/ubfactor)*tpwgts[0];
  minpwgt[1] = (1.0/ubfactor)*tpwgts[1];

  nbfs = (nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = idxsum(nvtxs, graph->adjwgtsum)+1;  /* The +1 is for the 0 edges case */
  for (; nbfs>0; nbfs--) {
    idxset(nvtxs, 0, touched);

    pwgts[1] = tpwgts[0]+tpwgts[1];
    pwgts[0] = 0;

    idxset(nvtxs, 1, where);

    queue[0] = RandomInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    for (;;) {
      if (first == last) { /* Empty. Disconnected graph! */
        if (nleft == 0 || drain)
          break;

        k = RandomInRange(nleft);
        for (i=0; i<nvtxs; i++) {
          if (touched[i] == 0) {
            if (k == 0)
              break;
            else
              k--;
          }
        }

        queue[0] = i;
        touched[i] = 1;
        first = 0; last = 1;;
        nleft--;
      }

      i = queue[first++];
      if (pwgts[0] > 0 && pwgts[1]-vwgt[i] < minpwgt[1]) {
        drain = 1;
        continue;
      }

      where[i] = 0;
      INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
      if (pwgts[1] <= maxpwgt[1])
        break;

      drain = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (touched[k] == 0) {
          queue[last++] = k;
          touched[k] = 1;
          nleft--;
        }
      }
    }

    /* Check to see if we hit any bad limiting cases */
    if (pwgts[1] == 0) {
      i = RandomInRange(nvtxs);
      where[i] = 1;
      INC_DEC(pwgts[1], pwgts[0], vwgt[i]);
    }

    /*************************************************************
    * Do some partition refinement
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    /*printf("IPART: %3d [%5d %5d] [%5d %5d] %5d\n", graph->nvtxs, pwgts[0], pwgts[1], graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    Balance2Way(ctrl, graph, tpwgts, ubfactor);
    /*printf("BPART: [%5d %5d] %5d\n", graph->pwgts[0], graph->pwgts[1], graph->mincut);*/

    FM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);
    /*printf("RPART: [%5d %5d] %5d\n", graph->pwgts[0], graph->pwgts[1], graph->mincut);*/

    if (bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  GKfree((void **) &bestwhere, (void **) &queue, (void **) &touched, LTERM);
}




/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void GrowBisectionNode(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  int i, j, k, nvtxs, drain, nleft, first, last, pwgts[2], tpwgts[2], minpwgt[2], maxpwgt[2], bestcut, nbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where, *bndind;
  idxtype *queue, *touched, *bestwhere;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  queue = idxmalloc(nvtxs, "BisectGraph: queue");
  touched = idxmalloc(nvtxs, "BisectGraph: touched");

  tpwgts[0] = idxsum(nvtxs, vwgt);
  tpwgts[1] = tpwgts[0]/2;
  tpwgts[0] -= tpwgts[1];

  maxpwgt[0] = ubfactor*tpwgts[0];
  maxpwgt[1] = ubfactor*tpwgts[1];
  minpwgt[0] = (1.0/ubfactor)*tpwgts[0];
  minpwgt[1] = (1.0/ubfactor)*tpwgts[1];

  /* Allocate memory for graph->rdata. Allocate sufficient memory for both edge and node */
  graph->rdata = idxmalloc(5*nvtxs+3, "GrowBisectionNode: graph->rdata");
  graph->pwgts    = graph->rdata;
  graph->where    = graph->rdata + 3;
  graph->bndptr   = graph->rdata + nvtxs + 3;
  graph->bndind   = graph->rdata + 2*nvtxs + 3;
  graph->nrinfo   = (NRInfoType *)(graph->rdata + 3*nvtxs + 3);
  graph->id       = graph->rdata + 3*nvtxs + 3;
  graph->ed       = graph->rdata + 4*nvtxs + 3;

  where = graph->where;
  bndind = graph->bndind;

  nbfs = (nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = tpwgts[0]+tpwgts[1];
  for (nbfs++; nbfs>0; nbfs--) {
    idxset(nvtxs, 0, touched);

    pwgts[1] = tpwgts[0]+tpwgts[1];
    pwgts[0] = 0;

    idxset(nvtxs, 1, where);

    queue[0] = RandomInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    if (nbfs >= 1) {
      for (;;) {
        if (first == last) { /* Empty. Disconnected graph! */
          if (nleft == 0 || drain)
            break;

          k = RandomInRange(nleft);
          for (i=0; i<nvtxs; i++) {
            if (touched[i] == 0) {
              if (k == 0)
                break;
              else
                k--;
            }
          }

          queue[0] = i;
          touched[i] = 1;
          first = 0; last = 1;;
          nleft--;
        }

        i = queue[first++];
        if (pwgts[1]-vwgt[i] < minpwgt[1]) {
          drain = 1;
          continue;
        }

        where[i] = 0;
        INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
        if (pwgts[1] <= maxpwgt[1])
          break;

        drain = 0;
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (touched[k] == 0) {
            queue[last++] = k;
            touched[k] = 1;
            nleft--;
          }
        }
      }
    }

    /*************************************************************
    * Do some partition refinement
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    Balance2Way(ctrl, graph, tpwgts, ubfactor);
    FM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);

    /* Construct and refine the vertex separator */
    for (i=0; i<graph->nbnd; i++)
      where[bndind[i]] = 2;

    Compute2WayNodePartitionParams(ctrl, graph);
    FM_2WayNodeRefine(ctrl, graph, ubfactor, 6);

    /* printf("ISep: [%d %d %d] %d\n", graph->pwgts[0], graph->pwgts[1], graph->pwgts[2], bestcut); */

    if (bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  Compute2WayNodePartitionParams(ctrl, graph);

  GKfree((void **) &bestwhere, (void **) &queue, (void **) &touched, LTERM);
}


/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void RandomBisection(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  int i, ii, nvtxs, pwgts[2], minpwgt[2], maxpwgt[2], bestcut, nbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where;
  idxtype *perm, *bestwhere;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  perm = idxmalloc(nvtxs, "BisectGraph: queue");

  ASSERTP(tpwgts[0]+tpwgts[1] == idxsum(nvtxs, vwgt), ("%d %d\n", tpwgts[0]+tpwgts[1], idxsum(nvtxs, vwgt)));

  maxpwgt[0] = ubfactor*tpwgts[0];
  maxpwgt[1] = ubfactor*tpwgts[1];
  minpwgt[0] = (1.0/ubfactor)*tpwgts[0];
  minpwgt[1] = (1.0/ubfactor)*tpwgts[1];

  nbfs = (nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = idxsum(nvtxs, graph->adjwgtsum)+1;  /* The +1 is for the 0 edges case */
  for (; nbfs>0; nbfs--) {
    RandomPermute(nvtxs, perm, 1);

    idxset(nvtxs, 1, where);
    pwgts[1] = tpwgts[0]+tpwgts[1];
    pwgts[0] = 0;


    if (nbfs != 1) {
      for (ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        if (pwgts[0]+vwgt[i] < maxpwgt[0]) {
          where[i] = 0;
          pwgts[0] += vwgt[i];
          pwgts[1] -= vwgt[i];
          if (pwgts[0] > minpwgt[0])
            break;
        }
      }
    }

    /*************************************************************
    * Do some partition refinement
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    /* printf("IPART: %3d [%5d %5d] [%5d %5d] %5d\n", graph->nvtxs, pwgts[0], pwgts[1], graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    Balance2Way(ctrl, graph, tpwgts, ubfactor);
    /* printf("BPART: [%5d %5d] %5d\n", graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    FM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);
    /* printf("RPART: [%5d %5d] %5d\n", graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    if (bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  GKfree((void **) &bestwhere, (void **) &perm, LTERM);
}




/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This file contains the top level routines for the multilevel k-way partitioning
 * algorithm KMETIS.
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void METIS_PartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                         idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                         int *options, int *edgecut, idxtype *part)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++)
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts,
                       tpwgts, options, edgecut, part);

  free(tpwgts);
}


/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_WPartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                          idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                          float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = KMETIS_CTYPE;
    ctrl.IType = KMETIS_ITYPE;
    ctrl.RType = KMETIS_RTYPE;
    ctrl.dbglvl = KMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = amax((*nvtxs)/(40*log2(*nparts)), 20*(*nparts));
  ctrl.maxvwgt = 1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt) : (*nvtxs))/ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MlevelKWayPartitioning(&ctrl, &graph, *nparts, part, tpwgts, 1.03);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MlevelKWayPartitioning(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part, float *tpwgts, float ubfactor)
{
  GraphType *cgraph;
  int wgtflag=3, numflag=0, options[10], edgecut;

  cgraph = Coarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
  AllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  options[0] = 1;
  options[OPTION_CTYPE] = MATCH_SHEMKWAY;
  options[OPTION_ITYPE] = IPART_GGPKL;
  options[OPTION_RTYPE] = RTYPE_FM;
  options[OPTION_DBGLVL] = 0;

  METIS_WPartGraphRecursive(&cgraph->nvtxs, cgraph->xadj, cgraph->adjncy, cgraph->vwgt,
                            cgraph->adjwgt, &wgtflag, &numflag, &nparts, tpwgts, options,
                            &edgecut, cgraph->where);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial %d-way partitioning cut: %d\n", nparts, edgecut));

  IFSET(ctrl->dbglvl, DBG_KWAYPINFO, ComputePartitionInfo(cgraph, nparts, cgraph->where));

  RefineKWay(ctrl, graph, cgraph, nparts, tpwgts, ubfactor);

  idxcopy(graph->nvtxs, graph->where, part);

  GKfree((void **) &graph->gdata, (void **) &graph->rdata, LTERM);

  return graph->mincut;

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kvmetis.c
 *
 * This file contains the top level routines for the multilevel k-way partitioning
 * algorithm KMETIS.
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void METIS_PartGraphVKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                         idxtype *vsize, int *wgtflag, int *numflag, int *nparts,
                         int *options, int *volume, idxtype *part)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++)
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts,
                       tpwgts, options, volume, part);

  free(tpwgts);
}


/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_WPartGraphVKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                          idxtype *vsize, int *wgtflag, int *numflag, int *nparts,
                          float *tpwgts, int *options, int *volume, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  VolSetUpGraph(&graph, OP_KVMETIS, *nvtxs, 1, xadj, adjncy, vwgt, vsize, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = KVMETIS_CTYPE;
    ctrl.IType = KVMETIS_ITYPE;
    ctrl.RType = KVMETIS_RTYPE;
    ctrl.dbglvl = KVMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KVMETIS;
  ctrl.CoarsenTo = amax((*nvtxs)/(40*log2(*nparts)), 20*(*nparts));
  ctrl.maxvwgt = 1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt) : (*nvtxs))/ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *volume = MlevelVolKWayPartitioning(&ctrl, &graph, *nparts, part, tpwgts, 1.03);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MlevelVolKWayPartitioning(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part,
                              float *tpwgts, float ubfactor)
{
  GraphType *cgraph;
  int wgtflag=3, numflag=0, options[10], edgecut;

  cgraph = Coarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
  AllocateVolKWayPartitionMemory(ctrl, cgraph, nparts);

  options[0] = 1;
  options[OPTION_CTYPE] = MATCH_SHEMKWAY;
  options[OPTION_ITYPE] = IPART_GGPKL;
  options[OPTION_RTYPE] = RTYPE_FM;
  options[OPTION_DBGLVL] = 0;

  METIS_WPartGraphRecursive(&cgraph->nvtxs, cgraph->xadj, cgraph->adjncy, cgraph->vwgt,
                            cgraph->adjwgt, &wgtflag, &numflag, &nparts, tpwgts, options,
                            &edgecut, cgraph->where);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial %d-way partitioning cut: %d\n", nparts, edgecut));

  IFSET(ctrl->dbglvl, DBG_KWAYPINFO, ComputePartitionInfo(cgraph, nparts, cgraph->where));

  RefineVolKWay(ctrl, graph, cgraph, nparts, tpwgts, ubfactor);

  idxcopy(graph->nvtxs, graph->where, part);

  GKfree((void **) &graph->gdata, (void **) &graph->rdata, LTERM);

  return graph->minvol;

}

/*
 * kwayfm.c
 *
 * This file contains code that implements the multilevel k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Random_KWayEdgeRefine(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor, int npasses, int ffactor)
{
  int i, ii, iii, j, k, pass, nvtxs, nmoves, nbnd, tvwgt, myndegrees;
  int from, me, to, oldcut, vwgt, gain;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *itpwgts;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d]-[%6d %6d], Balance: %5.3f, Nv-Nb[%6d %6d]. Cut: %6d\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (nmoves=iii=0; iii<graph->nbnd; iii++) {
      ii = perm[iii];
      if (ii >= nbnd)
        continue;
      i = bndind[ii];

      myrinfo = graph->rinfo+i;

      if (myrinfo->ed >= myrinfo->id) { /* Total ED is too high */
        from = where[i];
        vwgt = graph->vwgt[i];

        if (myrinfo->id > 0 && pwgts[from]-vwgt < minwgt[from])
          continue;   /* This cannot be moved! */

        myedegrees = myrinfo->edegrees;
        myndegrees = myrinfo->ndegrees;

        j = myrinfo->id;
        for (k=0; k<myndegrees; k++) {
          to = myedegrees[k].pid;
          gain = myedegrees[k].ed-j; /* j = myrinfo->id. Allow good nodes to move */
          if (pwgts[to]+vwgt <= maxwgt[to]+ffactor*gain && gain >= 0)
            break;
        }
        if (k == myndegrees)
          continue;  /* break out if you did not find a candidate */

        for (j=k+1; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          if ((myedegrees[j].ed > myedegrees[k].ed && pwgts[to]+vwgt <= maxwgt[to]) ||
              (myedegrees[j].ed == myedegrees[k].ed &&
               itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid]))
            k = j;
        }

        to = myedegrees[k].pid;

        j = 0;
        if (myedegrees[k].ed-myrinfo->id > 0)
          j = 1;
        else if (myedegrees[k].ed-myrinfo->id == 0) {
          if ((iii&7) == 0 || pwgts[from] >= maxwgt[from] || itpwgts[from]*(pwgts[to]+vwgt) < itpwgts[to]*pwgts[from])
            j = 1;
        }
        if (j == 0)
          continue;

        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to'
        *======================================================================*/
        graph->mincut -= myedegrees[k].ed-myrinfo->id;

        IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

        /* Update where, weight, and ID/ED information of the vertex you moved */
        where[i] = to;
        INC_DEC(pwgts[to], pwgts[from], vwgt);
        myrinfo->ed += myrinfo->id-myedegrees[k].ed;
        SWAP(myrinfo->id, myedegrees[k].ed, j);
        if (myedegrees[k].ed == 0)
          myedegrees[k] = myedegrees[--myrinfo->ndegrees];
        else
          myedegrees[k].pid = from;

        if (myrinfo->ed-myrinfo->id < 0)
          BNDDelete(nbnd, bndind, bndptr, i);

        /* Update the degrees of adjacent vertices */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          ii = adjncy[j];
          me = where[ii];

          myrinfo = graph->rinfo+ii;
          if (myrinfo->edegrees == NULL) {
            myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
            ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
          }
          myedegrees = myrinfo->edegrees;

          ASSERT(CheckRInfo(myrinfo));

          if (me == from) {
            INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

            if (myrinfo->ed-myrinfo->id >= 0 && bndptr[ii] == -1)
              BNDInsert(nbnd, bndind, bndptr, ii);
          }
          else if (me == to) {
            INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

            if (myrinfo->ed-myrinfo->id < 0 && bndptr[ii] != -1)
              BNDDelete(nbnd, bndind, bndptr, ii);
          }

          /* Remove contribution from the .ed of 'from' */
          if (me != from) {
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (myedegrees[k].pid == from) {
                if (myedegrees[k].ed == adjwgt[j])
                  myedegrees[k] = myedegrees[--myrinfo->ndegrees];
                else
                  myedegrees[k].ed -= adjwgt[j];
                break;
              }
            }
          }

          /* Add contribution to the .ed of 'to' */
          if (me != to) {
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (myedegrees[k].pid == to) {
                myedegrees[k].ed += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->ndegrees) {
              myedegrees[myrinfo->ndegrees].pid = to;
              myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
            }
          }

          ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
          ASSERT(CheckRInfo(myrinfo));

        }
        nmoves++;
      }
    }

    graph->nbnd = nbnd;

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d, Vol: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut, ComputeVolume(graph, where)));

    if (graph->mincut == oldcut)
      break;
  }

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
}






/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Greedy_KWayEdgeRefine(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor, int npasses)
{
  int i, ii, iii, j, k, pass, nvtxs, nbnd, tvwgt, myndegrees, oldgain, gain;
  int from, me, to, oldcut, vwgt;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *moved, *itpwgts;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;
  PQueueType queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxwspacemalloc(ctrl, nvtxs);

  PQueueInit(ctrl, &queue, nvtxs, graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)]);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d]-[%6d %6d], Balance: %5.3f, Nv-Nb[%6d %6d]. Cut: %6d\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    PQueueReset(&queue);
    idxset(nvtxs, -1, moved);

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      PQueueInsert(&queue, i, graph->rinfo[i].ed - graph->rinfo[i].id);
      moved[i] = 2;
    }

    for (iii=0;;iii++) {
      if ((i = PQueueGetMax(&queue)) == -1)
        break;
      moved[i] = 1;

      myrinfo = graph->rinfo+i;
      from = where[i];
      vwgt = graph->vwgt[i];

      if (pwgts[from]-vwgt < minwgt[from])
        continue;   /* This cannot be moved! */

      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      j = myrinfo->id;
      for (k=0; k<myndegrees; k++) {
        to = myedegrees[k].pid;
        gain = myedegrees[k].ed-j; /* j = myrinfo->id. Allow good nodes to move */
        if (pwgts[to]+vwgt <= maxwgt[to]+gain && gain >= 0)
          break;
      }
      if (k == myndegrees)
        continue;  /* break out if you did not find a candidate */

      for (j=k+1; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        if ((myedegrees[j].ed > myedegrees[k].ed && pwgts[to]+vwgt <= maxwgt[to]) ||
            (myedegrees[j].ed == myedegrees[k].ed &&
             itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid]))
          k = j;
      }

      to = myedegrees[k].pid;

      j = 0;
      if (myedegrees[k].ed-myrinfo->id > 0)
        j = 1;
      else if (myedegrees[k].ed-myrinfo->id == 0) {
        if ((iii&7) == 0 || pwgts[from] >= maxwgt[from] || itpwgts[from]*(pwgts[to]+vwgt) < itpwgts[to]*pwgts[from])
          j = 1;
      }
      if (j == 0)
        continue;

      /*=====================================================================
      * If we got here, we can now move the vertex from 'from' to 'to'
      *======================================================================*/
      graph->mincut -= myedegrees[k].ed-myrinfo->id;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

      /* Update where, weight, and ID/ED information of the vertex you moved */
      where[i] = to;
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      myrinfo->ed += myrinfo->id-myedegrees[k].ed;
      SWAP(myrinfo->id, myedegrees[k].ed, j);
      if (myedegrees[k].ed == 0)
        myedegrees[k] = myedegrees[--myrinfo->ndegrees];
      else
        myedegrees[k].pid = from;

      if (myrinfo->ed < myrinfo->id)
        BNDDelete(nbnd, bndind, bndptr, i);

      /* Update the degrees of adjacent vertices */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        me = where[ii];

        myrinfo = graph->rinfo+ii;
        if (myrinfo->edegrees == NULL) {
          myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
          ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
        }
        myedegrees = myrinfo->edegrees;

        ASSERT(CheckRInfo(myrinfo));

        oldgain = (myrinfo->ed-myrinfo->id);

        if (me == from) {
          INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

          if (myrinfo->ed-myrinfo->id >= 0 && bndptr[ii] == -1)
            BNDInsert(nbnd, bndind, bndptr, ii);
        }
        else if (me == to) {
          INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

          if (myrinfo->ed-myrinfo->id < 0 && bndptr[ii] != -1)
            BNDDelete(nbnd, bndind, bndptr, ii);
        }

        /* Remove contribution from the .ed of 'from' */
        if (me != from) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == from) {
              if (myedegrees[k].ed == adjwgt[j])
                myedegrees[k] = myedegrees[--myrinfo->ndegrees];
              else
                myedegrees[k].ed -= adjwgt[j];
              break;
            }
          }
        }

        /* Add contribution to the .ed of 'to' */
        if (me != to) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == to) {
              myedegrees[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].pid = to;
            myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
          }
        }

        /* Update the queue */
        if (me == to || me == from) {
          gain = myrinfo->ed-myrinfo->id;
          if (moved[ii] == 2) {
            if (gain >= 0)
              PQueueUpdate(&queue, ii, oldgain, gain);
            else {
              PQueueDelete(&queue, ii, oldgain);
              moved[ii] = -1;
            }
          }
          else if (moved[ii] == -1 && gain >= 0) {
            PQueueInsert(&queue, ii, gain);
            moved[ii] = 2;
          }
        }

        ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
        ASSERT(CheckRInfo(myrinfo));

      }
    }

    graph->nbnd = nbnd;

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Cut: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, graph->mincut));

    if (graph->mincut == oldcut)
      break;
  }

  PQueueFree(ctrl, &queue);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Greedy_KWayEdgeBalance(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor, int npasses)
{
  int i, ii, j, k, pass, nvtxs, nbnd, tvwgt, myndegrees, oldgain, gain, nmoves;
  int from, me, to, oldcut, vwgt;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *moved, *itpwgts;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;
  PQueueType queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxwspacemalloc(ctrl, nvtxs);

  PQueueInit(ctrl, &queue, nvtxs, graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)]);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d]-[%6d %6d], Balance: %5.3f, Nv-Nb[%6d %6d]. Cut: %6d [B]\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    /* Check to see if things are out of balance, given the tolerance */
    for (i=0; i<nparts; i++) {
      if (pwgts[i] > maxwgt[i])
        break;
    }
    if (i == nparts) /* Things are balanced. Return right away */
      break;

    PQueueReset(&queue);
    idxset(nvtxs, -1, moved);

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      PQueueInsert(&queue, i, graph->rinfo[i].ed - graph->rinfo[i].id);
      moved[i] = 2;
    }

    nmoves = 0;
    for (;;) {
      if ((i = PQueueGetMax(&queue)) == -1)
        break;
      moved[i] = 1;

      myrinfo = graph->rinfo+i;
      from = where[i];
      vwgt = graph->vwgt[i];

      if (pwgts[from]-vwgt < minwgt[from])
        continue;   /* This cannot be moved! */

      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      for (k=0; k<myndegrees; k++) {
        to = myedegrees[k].pid;
        if (pwgts[to]+vwgt <= maxwgt[to] || itpwgts[from]*(pwgts[to]+vwgt) <= itpwgts[to]*pwgts[from])
          break;
      }
      if (k == myndegrees)
        continue;  /* break out if you did not find a candidate */

      for (j=k+1; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        if (itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid])
          k = j;
      }

      to = myedegrees[k].pid;

      if (pwgts[from] < maxwgt[from] && pwgts[to] > minwgt[to] && myedegrees[k].ed-myrinfo->id < 0)
        continue;

      /*=====================================================================
      * If we got here, we can now move the vertex from 'from' to 'to'
      *======================================================================*/
      graph->mincut -= myedegrees[k].ed-myrinfo->id;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

      /* Update where, weight, and ID/ED information of the vertex you moved */
      where[i] = to;
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      myrinfo->ed += myrinfo->id-myedegrees[k].ed;
      SWAP(myrinfo->id, myedegrees[k].ed, j);
      if (myedegrees[k].ed == 0)
        myedegrees[k] = myedegrees[--myrinfo->ndegrees];
      else
        myedegrees[k].pid = from;

      if (myrinfo->ed == 0)
        BNDDelete(nbnd, bndind, bndptr, i);

      /* Update the degrees of adjacent vertices */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        me = where[ii];

        myrinfo = graph->rinfo+ii;
        if (myrinfo->edegrees == NULL) {
          myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
          ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
        }
        myedegrees = myrinfo->edegrees;

        ASSERT(CheckRInfo(myrinfo));

        oldgain = (myrinfo->ed-myrinfo->id);

        if (me == from) {
          INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

          if (myrinfo->ed > 0 && bndptr[ii] == -1)
            BNDInsert(nbnd, bndind, bndptr, ii);
        }
        else if (me == to) {
          INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

          if (myrinfo->ed == 0 && bndptr[ii] != -1)
            BNDDelete(nbnd, bndind, bndptr, ii);
        }

        /* Remove contribution from the .ed of 'from' */
        if (me != from) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == from) {
              if (myedegrees[k].ed == adjwgt[j])
                myedegrees[k] = myedegrees[--myrinfo->ndegrees];
              else
                myedegrees[k].ed -= adjwgt[j];
              break;
            }
          }
        }

        /* Add contribution to the .ed of 'to' */
        if (me != to) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == to) {
              myedegrees[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].pid = to;
            myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
          }
        }

        /* Update the queue */
        if (me == to || me == from) {
          gain = myrinfo->ed-myrinfo->id;
          if (moved[ii] == 2) {
            if (myrinfo->ed > 0)
              PQueueUpdate(&queue, ii, oldgain, gain);
            else {
              PQueueDelete(&queue, ii, oldgain);
              moved[ii] = -1;
            }
          }
          else if (moved[ii] == -1 && myrinfo->ed > 0) {
            PQueueInsert(&queue, ii, gain);
            moved[ii] = 2;
          }
        }

        ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
        ASSERT(CheckRInfo(myrinfo));
      }
      nmoves++;
    }

    graph->nbnd = nbnd;

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut));
  }

  PQueueFree(ctrl, &queue);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kwayrefine.c
 *
 * This file contains the driving routines for multilevel k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void RefineKWay(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int nparts, float *tpwgts, float ubfactor)
{
  int i, nlevels, mustfree=0;
  GraphType *ptr;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  ComputeKWayPartitionParams(ctrl, graph, nparts);

  /* Take care any non-contiguity */
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->AuxTmr1));
  if (ctrl->RType == RTYPE_KWAYRANDOM_MCONN) {
    EliminateComponents(ctrl, graph, nparts, tpwgts, 1.25);
    EliminateSubDomainEdges(ctrl, graph, nparts, tpwgts);
    EliminateComponents(ctrl, graph, nparts, tpwgts, 1.25);
  }
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->AuxTmr1));

  /* Determine how many levels are there */
  for (ptr=graph, nlevels=0; ptr!=orggraph; ptr=ptr->finer, nlevels++);

  for (i=0; ;i++) {
    /* PrintSubDomainGraph(graph, nparts, graph->where); */
    if (ctrl->RType == RTYPE_KWAYRANDOM_MCONN && (i == nlevels/2 || i == nlevels/2+1))
      EliminateSubDomainEdges(ctrl, graph, nparts, tpwgts);

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));

    if (2*i >= nlevels && !IsBalanced(graph->pwgts, nparts, tpwgts, 1.04*ubfactor)) {
      ComputeKWayBalanceBoundary(ctrl, graph, nparts);
      if (ctrl->RType == RTYPE_KWAYRANDOM_MCONN)
        Greedy_KWayEdgeBalanceMConn(ctrl, graph, nparts, tpwgts, ubfactor, 1);
      else
        Greedy_KWayEdgeBalance(ctrl, graph, nparts, tpwgts, ubfactor, 1);
      ComputeKWayBoundary(ctrl, graph, nparts);
    }

    switch (ctrl->RType) {
      case RTYPE_KWAYRANDOM:
        Random_KWayEdgeRefine(ctrl, graph, nparts, tpwgts, ubfactor, 10, 1);
        break;
      case RTYPE_KWAYGREEDY:
        Greedy_KWayEdgeRefine(ctrl, graph, nparts, tpwgts, ubfactor, 10);
        break;
      case RTYPE_KWAYRANDOM_MCONN:
        Random_KWayEdgeRefineMConn(ctrl, graph, nparts, tpwgts, ubfactor, 10, 1);
        break;
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    GKfree((void **) &graph->gdata, LTERM);  /* Deallocate the graph related arrays */

    graph = graph->finer;

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    if (graph->vwgt == NULL) {
      graph->vwgt = idxsmalloc(graph->nvtxs, 1, "RefineKWay: graph->vwgt");
      graph->adjwgt = idxsmalloc(graph->nedges, 1, "RefineKWay: graph->adjwgt");
      mustfree = 1;
    }
    ProjectKWayPartition(ctrl, graph, nparts);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  if (!IsBalanced(graph->pwgts, nparts, tpwgts, ubfactor)) {
    ComputeKWayBalanceBoundary(ctrl, graph, nparts);
    if (ctrl->RType == RTYPE_KWAYRANDOM_MCONN) {
      Greedy_KWayEdgeBalanceMConn(ctrl, graph, nparts, tpwgts, ubfactor, 8);
      Random_KWayEdgeRefineMConn(ctrl, graph, nparts, tpwgts, ubfactor, 10, 0);
    }
    else {
      Greedy_KWayEdgeBalance(ctrl, graph, nparts, tpwgts, ubfactor, 8);
      Random_KWayEdgeRefine(ctrl, graph, nparts, tpwgts, ubfactor, 10, 0);
    }
  }

  /* Take care any trivial non-contiguity */
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->AuxTmr2));
  EliminateComponents(ctrl, graph, nparts, tpwgts, ubfactor);
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->AuxTmr2));

  if (mustfree)
    GKfree((void **) &graph->vwgt, (void **) &graph->adjwgt, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}


/*************************************************************************
* This function allocates memory for k-way edge refinement
**************************************************************************/
void AllocateKWayPartitionMemory(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int nvtxs, pad64;

  nvtxs = graph->nvtxs;

  pad64 = (3*nvtxs+nparts)%2;

  graph->rdata = idxmalloc(3*nvtxs+nparts+(sizeof(RInfoType)/sizeof(idxtype))*nvtxs+pad64, "AllocateKWayPartitionMemory: rdata");
  graph->pwgts          = graph->rdata;
  graph->where          = graph->rdata + nparts;
  graph->bndptr         = graph->rdata + nvtxs + nparts;
  graph->bndind         = graph->rdata + 2*nvtxs + nparts;
  graph->rinfo          = (RInfoType *)(graph->rdata + 3*nvtxs+nparts + pad64);

/*
  if (ctrl->wspace.edegrees != NULL)
    free(ctrl->wspace.edegrees);
  ctrl->wspace.edegrees = (EDegreeType *)GKmalloc(graph->nedges*sizeof(EDegreeType), "AllocateKWayPartitionMemory: edegrees");
*/
}


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void ComputeKWayPartitionParams(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, j, k, nvtxs, nbnd, mincut, me, other;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *pwgts, *where, *bndind, *bndptr;
  RInfoType *rinfo, *myrinfo;
  EDegreeType *myedegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = idxset(nparts, 0, graph->pwgts);
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  rinfo = graph->rinfo;


  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  ctrl->wspace.cdegree = 0;
  nbnd = mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];

    myrinfo = rinfo+i;
    myrinfo->id = myrinfo->ed = myrinfo->ndegrees = 0;
    myrinfo->edegrees = NULL;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me != where[adjncy[j]])
        myrinfo->ed += adjwgt[j];
    }
    myrinfo->id = graph->adjwgtsum[i] - myrinfo->ed;

    if (myrinfo->ed > 0)
      mincut += myrinfo->ed;

    if (myrinfo->ed-myrinfo->id >= 0)
      BNDInsert(nbnd, bndind, bndptr, i);

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      myedegrees = myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += xadj[i+1]-xadj[i];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == other) {
              myedegrees[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].pid = other;
            myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
          }
        }
      }

      ASSERT(myrinfo->ndegrees <= xadj[i+1]-xadj[i]);
    }
  }

  graph->mincut = mincut/2;
  graph->nbnd = nbnd;

}



/*************************************************************************
* This function projects a partition, and at the same time computes the
* parameters for refinement.
**************************************************************************/
void ProjectKWayPartition(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, j, k, nvtxs, nbnd, me, other, istart, iend, ndegrees;
  idxtype *xadj, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *cmap, *where, *bndptr, *bndind;
  idxtype *cwhere;
  GraphType *cgraph;
  RInfoType *crinfo, *rinfo, *myrinfo;
  EDegreeType *myedegrees;
  idxtype *htable;

  cgraph = graph->coarser;
  cwhere = cgraph->where;
  crinfo = cgraph->rinfo;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  AllocateKWayPartitionMemory(ctrl, graph, nparts);
  where = graph->where;
  rinfo = graph->rinfo;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);

  /* Go through and project partition and compute id/ed for the nodes */
  for (i=0; i<nvtxs; i++) {
    k = cmap[i];
    where[i] = cwhere[k];
    cmap[i] = crinfo[k].ed;  /* For optimization */
  }

  htable = idxset(nparts, -1, idxwspacemalloc(ctrl, nparts));

  ctrl->wspace.cdegree = 0;
  for (nbnd=0, i=0; i<nvtxs; i++) {
    me = where[i];

    myrinfo = rinfo+i;
    myrinfo->id = myrinfo->ed = myrinfo->ndegrees = 0;
    myrinfo->edegrees = NULL;

    myrinfo->id = adjwgtsum[i];

    if (cmap[i] > 0) { /* If it is an interface node. Note cmap[i] = crinfo[cmap[i]].ed */
      istart = xadj[i];
      iend = xadj[i+1];

      myedegrees = myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += iend-istart;

      ndegrees = 0;
      for (j=istart; j<iend; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          myrinfo->ed += adjwgt[j];
          if ((k = htable[other]) == -1) {
            htable[other] = ndegrees;
            myedegrees[ndegrees].pid = other;
            myedegrees[ndegrees++].ed = adjwgt[j];
          }
          else {
            myedegrees[k].ed += adjwgt[j];
          }
        }
      }
      myrinfo->id -= myrinfo->ed;

      /* Remove space for edegrees if it was interior */
      if (myrinfo->ed == 0) {
        myrinfo->edegrees = NULL;
        ctrl->wspace.cdegree -= iend-istart;
      }
      else {
        if (myrinfo->ed-myrinfo->id >= 0)
          BNDInsert(nbnd, bndind, bndptr, i);

        myrinfo->ndegrees = ndegrees;

        for (j=0; j<ndegrees; j++)
          htable[myedegrees[j].pid] = -1;
      }
    }
  }

  idxcopy(nparts, cgraph->pwgts, graph->pwgts);
  graph->mincut = cgraph->mincut;
  graph->nbnd = nbnd;

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

  idxwspacefree(ctrl, nparts);

  ASSERT(CheckBnd2(graph));

}



/*************************************************************************
* This function checks if the partition weights are within the balance
* contraints
**************************************************************************/
int IsBalanced(idxtype *pwgts, int nparts, float *tpwgts, float ubfactor)
{
  int i, tvwgt;

  tvwgt = idxsum(nparts, pwgts);
  for (i=0; i<nparts; i++) {
    if (pwgts[i] > tpwgts[i]*tvwgt*(ubfactor+0.005))
      return 0;
  }

  return 1;
}


/*************************************************************************
* This function computes the boundary definition for balancing
**************************************************************************/
void ComputeKWayBoundary(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, nvtxs, nbnd;
  idxtype *bndind, *bndptr;

  nvtxs = graph->nvtxs;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);


  /*------------------------------------------------------------
  / Compute the new boundary
  /------------------------------------------------------------*/
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    if (graph->rinfo[i].ed-graph->rinfo[i].id >= 0)
      BNDInsert(nbnd, bndind, bndptr, i);
  }

  graph->nbnd = nbnd;
}

/*************************************************************************
* This function computes the boundary definition for balancing
**************************************************************************/
void ComputeKWayBalanceBoundary(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, nvtxs, nbnd;
  idxtype *bndind, *bndptr;

  nvtxs = graph->nvtxs;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);


  /*------------------------------------------------------------
  / Compute the new boundary
  /------------------------------------------------------------*/
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    if (graph->rinfo[i].ed > 0)
      BNDInsert(nbnd, bndind, bndptr, i);
  }

  graph->nbnd = nbnd;
}

/*
 * kwayvolfm.c
 *
 * This file contains code that implements the multilevel k-way refinement
 *
 * Started 7/8/98
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Random_KWayVolRefine(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts,
                          float ubfactor, int npasses, int ffactor)
{
  int i, ii, iii, j, k, pass, nvtxs, nmoves, tvwgt, myndegrees, xgain;
  int from, to, oldcut, oldvol, vwgt;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *itpwgts, *updind, *marker, *phtable;
  VEDegreeType *myedegrees;
  VRInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  updind = idxmalloc(nvtxs, "Random_KWayVolRefine: updind");
  marker = idxsmalloc(nvtxs, 0, "Random_KWayVolRefine: marker");
  phtable = idxsmalloc(nparts, -1, "Random_KWayVolRefine: phtable");

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("VolPart: [%5d %5d]-[%5d %5d], Balance: %3.2f, Nv-Nb[%5d %5d]. Cut: %5d, Vol: %5d\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut, graph->minvol));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    oldcut = graph->mincut;
    oldvol = graph->minvol;

    RandomPermute(graph->nbnd, perm, 1);
    for (nmoves=iii=0; iii<graph->nbnd; iii++) {
      ii = perm[iii];
      if (ii >= graph->nbnd)
        continue;
      i = bndind[ii];
      myrinfo = graph->vrinfo+i;

      if (myrinfo->gv >= 0) { /* Total volume gain is too high */
        from = where[i];
        vwgt = graph->vwgt[i];

        if (myrinfo->id > 0 && pwgts[from]-vwgt < minwgt[from])
          continue;   /* This cannot be moved! */

        xgain = (myrinfo->id == 0 && myrinfo->ed > 0 ? graph->vsize[i] : 0);

        myedegrees = myrinfo->edegrees;
        myndegrees = myrinfo->ndegrees;

        for (k=0; k<myndegrees; k++) {
          to = myedegrees[k].pid;
          if (pwgts[to]+vwgt <= maxwgt[to]+ffactor*myedegrees[k].gv && xgain+myedegrees[k].gv >= 0)
            break;
        }
        if (k == myndegrees)
          continue;  /* break out if you did not find a candidate */

        for (j=k+1; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          if (pwgts[to]+vwgt > maxwgt[to])
            continue;
          if (myedegrees[j].gv > myedegrees[k].gv ||
              (myedegrees[j].gv == myedegrees[k].gv && myedegrees[j].ed > myedegrees[k].ed) ||
              (myedegrees[j].gv == myedegrees[k].gv && myedegrees[j].ed == myedegrees[k].ed &&
               itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid]))
            k = j;
        }

        to = myedegrees[k].pid;

        j = 0;
        if (xgain+myedegrees[k].gv > 0 || myedegrees[k].ed-myrinfo->id > 0)
          j = 1;
        else if (myedegrees[k].ed-myrinfo->id == 0) {
          if ((iii&5) == 0 || pwgts[from] >= maxwgt[from] || itpwgts[from]*(pwgts[to]+vwgt) < itpwgts[to]*pwgts[from])
            j = 1;
        }
        if (j == 0)
          continue;

        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to'
        *======================================================================*/
        INC_DEC(pwgts[to], pwgts[from], vwgt);
        graph->mincut -= myedegrees[k].ed-myrinfo->id;
        graph->minvol -= (xgain+myedegrees[k].gv);
        where[i] = to;

        IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d from %3d to %3d. Gain: [%4d %4d]. Cut: %6d, Vol: %6d\n",
              i, from, to, xgain+myedegrees[k].gv, myedegrees[k].ed-myrinfo->id, graph->mincut, graph->minvol));

        KWayVolUpdate(ctrl, graph, i, from, to, marker, phtable, updind);

        nmoves++;

        /* CheckVolKWayPartitionParams(ctrl, graph, nparts); */
      }
    }

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d, Vol: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut,
               graph->minvol));

    if (graph->minvol == oldvol && graph->mincut == oldcut)
      break;
  }

  GKfree((void **) &marker, (void **) &updind, (void **) &phtable, LTERM);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Random_KWayVolRefineMConn(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts,
            float ubfactor, int npasses, int ffactor)
{
  int i, ii, iii, j, k, l, pass, nvtxs, nmoves, tvwgt, myndegrees, xgain;
  int from, me, to, oldcut, oldvol, vwgt, nadd, maxndoms;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *itpwgts, *updind, *marker, *phtable;
  idxtype *pmat, *pmatptr, *ndoms;
  VEDegreeType *myedegrees;
  VRInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  updind = idxmalloc(nvtxs, "Random_KWayVolRefine: updind");
  marker = idxsmalloc(nvtxs, 0, "Random_KWayVolRefine: marker");
  phtable = idxsmalloc(nparts, -1, "Random_KWayVolRefine: phtable");

  pmat = ctrl->wspace.pmat;
  ndoms = idxwspacemalloc(ctrl, nparts);

  ComputeVolSubDomainGraph(graph, nparts, pmat, ndoms);

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("VolPart: [%5d %5d]-[%5d %5d], Balance: %3.2f, Nv-Nb[%5d %5d]. Cut: %5d, Vol: %5d\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut, graph->minvol));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    maxndoms = ndoms[idxamax(nparts, ndoms)];

    oldcut = graph->mincut;
    oldvol = graph->minvol;

    RandomPermute(graph->nbnd, perm, 1);
    for (nmoves=iii=0; iii<graph->nbnd; iii++) {
      ii = perm[iii];
      if (ii >= graph->nbnd)
        continue;
      i = bndind[ii];
      myrinfo = graph->vrinfo+i;

      if (myrinfo->gv >= 0) { /* Total volume gain is too high */
        from = where[i];
        vwgt = graph->vwgt[i];

        if (myrinfo->id > 0 && pwgts[from]-vwgt < minwgt[from])
          continue;   /* This cannot be moved! */

        xgain = (myrinfo->id == 0 && myrinfo->ed > 0 ? graph->vsize[i] : 0);

        myedegrees = myrinfo->edegrees;
        myndegrees = myrinfo->ndegrees;

        /* Determine the valid domains */
        for (j=0; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          phtable[to] = 1;
          pmatptr = pmat + to*nparts;
          for (nadd=0, k=0; k<myndegrees; k++) {
            if (k == j)
              continue;

            l = myedegrees[k].pid;
            if (pmatptr[l] == 0) {
              if (ndoms[l] > maxndoms-1) {
                phtable[to] = 0;
                nadd = maxndoms;
                break;
              }
              nadd++;
            }
          }
          if (ndoms[to]+nadd > maxndoms)
            phtable[to] = 0;
          if (nadd == 0)
            phtable[to] = 2;
        }

        for (k=0; k<myndegrees; k++) {
          to = myedegrees[k].pid;
          if (!phtable[to])
            continue;
          if (pwgts[to]+vwgt <= maxwgt[to]+ffactor*myedegrees[k].gv && xgain+myedegrees[k].gv >= 0)
            break;
        }
        if (k == myndegrees)
          continue;  /* break out if you did not find a candidate */

        for (j=k+1; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          if (!phtable[to] || pwgts[to]+vwgt > maxwgt[to])
            continue;
          if (myedegrees[j].gv > myedegrees[k].gv ||
              (myedegrees[j].gv == myedegrees[k].gv && myedegrees[j].ed > myedegrees[k].ed) ||
              (myedegrees[j].gv == myedegrees[k].gv && myedegrees[j].ed == myedegrees[k].ed &&
               itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid]))
            k = j;
        }

        to = myedegrees[k].pid;

        j = 0;
        if (xgain+myedegrees[k].gv > 0 || myedegrees[k].ed-myrinfo->id > 0)
          j = 1;
        else if (myedegrees[k].ed-myrinfo->id == 0) {
          if ((iii&5) == 0 || phtable[myedegrees[k].pid] == 2 || pwgts[from] >= maxwgt[from] || itpwgts[from]*(pwgts[to]+vwgt) < itpwgts[to]*pwgts[from])
            j = 1;
        }

        if (j == 0)
          continue;

        for (j=0; j<myndegrees; j++)
          phtable[myedegrees[j].pid] = -1;


        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to'
        *======================================================================*/
        INC_DEC(pwgts[to], pwgts[from], vwgt);
        graph->mincut -= myedegrees[k].ed-myrinfo->id;
        graph->minvol -= (xgain+myedegrees[k].gv);
        where[i] = to;

        IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d from %3d to %3d. Gain: [%4d %4d]. Cut: %6d, Vol: %6d\n",
              i, from, to, xgain+myedegrees[k].gv, myedegrees[k].ed-myrinfo->id, graph->mincut, graph->minvol));

        /* Update pmat to reflect the move of 'i' */
        pmat[from*nparts+to] += (myrinfo->id-myedegrees[k].ed);
        pmat[to*nparts+from] += (myrinfo->id-myedegrees[k].ed);
        if (pmat[from*nparts+to] == 0) {
          ndoms[from]--;
          if (ndoms[from]+1 == maxndoms)
            maxndoms = ndoms[idxamax(nparts, ndoms)];
        }
        if (pmat[to*nparts+from] == 0) {
          ndoms[to]--;
          if (ndoms[to]+1 == maxndoms)
            maxndoms = ndoms[idxamax(nparts, ndoms)];
        }

        for (j=xadj[i]; j<xadj[i+1]; j++) {
          ii = adjncy[j];
          me = where[ii];

          /* Update pmat to reflect the move of 'i' for domains other than 'from' and 'to' */
          if (me != from && me != to) {
            pmat[me*nparts+from] -= adjwgt[j];
            pmat[from*nparts+me] -= adjwgt[j];
            if (pmat[me*nparts+from] == 0) {
              ndoms[me]--;
              if (ndoms[me]+1 == maxndoms)
                maxndoms = ndoms[idxamax(nparts, ndoms)];
            }
            if (pmat[from*nparts+me] == 0) {
              ndoms[from]--;
              if (ndoms[from]+1 == maxndoms)
                maxndoms = ndoms[idxamax(nparts, ndoms)];
            }

            if (pmat[me*nparts+to] == 0) {
              ndoms[me]++;
              if (ndoms[me] > maxndoms) {
                printf("You just increased the maxndoms: %d %d\n", ndoms[me], maxndoms);
                maxndoms = ndoms[me];
              }
            }
            if (pmat[to*nparts+me] == 0) {
              ndoms[to]++;
              if (ndoms[to] > maxndoms) {
                printf("You just increased the maxndoms: %d %d\n", ndoms[to], maxndoms);
                maxndoms = ndoms[to];
              }
            }
            pmat[me*nparts+to] += adjwgt[j];
            pmat[to*nparts+me] += adjwgt[j];
          }
        }

        KWayVolUpdate(ctrl, graph, i, from, to, marker, phtable, updind);

        nmoves++;

        /* CheckVolKWayPartitionParams(ctrl, graph, nparts); */
      }
    }

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d, Vol: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut,
               graph->minvol));

    if (graph->minvol == oldvol && graph->mincut == oldcut)
      break;
  }

  GKfree((void **) &marker, (void **) &updind, (void **) &phtable, LTERM);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
}




/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Greedy_KWayVolBalance(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts,
                           float ubfactor, int npasses)
{
  int i, ii, j, k, pass, nvtxs, nmoves, tvwgt, myndegrees, xgain;
  int from, to, vwgt;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *moved, *bndptr, *bndind, *minwgt, *maxwgt, *itpwgts, *updind, *marker, *phtable;
  VEDegreeType *myedegrees;
  VRInfoType *myrinfo;
  PQueueType queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  updind = idxmalloc(nvtxs, "Random_KWayVolRefine: updind");
  marker = idxsmalloc(nvtxs, 0, "Random_KWayVolRefine: marker");
  phtable = idxsmalloc(nparts, -1, "Random_KWayVolRefine: phtable");

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxwspacemalloc(ctrl, nvtxs);

  PQueueInit(ctrl, &queue, nvtxs, graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)]);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("VolPart: [%5d %5d]-[%5d %5d], Balance: %3.2f, Nv-Nb[%5d %5d]. Cut: %5d, Vol: %5d [B]\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut, graph->minvol));


  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);
    /* Check to see if things are out of balance, given the tolerance */
    for (i=0; i<nparts; i++) {
      if (pwgts[i] > maxwgt[i])
        break;
    }
    if (i == nparts) /* Things are balanced. Return right away */
      break;

    PQueueReset(&queue);
    idxset(nvtxs, -1, moved);

    RandomPermute(graph->nbnd, perm, 1);
    for (ii=0; ii<graph->nbnd; ii++) {
      i = bndind[perm[ii]];
      PQueueInsert(&queue, i, graph->vrinfo[i].gv);
      moved[i] = 2;
    }

    for (nmoves=0;;) {
      if ((i = PQueueGetMax(&queue)) == -1)
        break;
      moved[i] = 1;

      myrinfo = graph->vrinfo+i;
      from = where[i];
      vwgt = graph->vwgt[i];

      if (pwgts[from]-vwgt < minwgt[from])
        continue;   /* This cannot be moved! */

      xgain = (myrinfo->id == 0 && myrinfo->ed > 0 ? graph->vsize[i] : 0);

      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      for (k=0; k<myndegrees; k++) {
        to = myedegrees[k].pid;
        if (pwgts[to]+vwgt <= maxwgt[to] ||
            itpwgts[from]*(pwgts[to]+vwgt) <= itpwgts[to]*pwgts[from])
          break;
      }
      if (k == myndegrees)
        continue;  /* break out if you did not find a candidate */

      for (j=k+1; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        if (itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid])
          k = j;
      }

      to = myedegrees[k].pid;

      if (pwgts[from] < maxwgt[from] && pwgts[to] > minwgt[to] &&
          (xgain+myedegrees[k].gv < 0 ||
           (xgain+myedegrees[k].gv == 0 &&  myedegrees[k].ed-myrinfo->id < 0))
         )
        continue;


      /*=====================================================================
      * If we got here, we can now move the vertex from 'from' to 'to'
      *======================================================================*/
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      graph->mincut -= myedegrees[k].ed-myrinfo->id;
      graph->minvol -= (xgain+myedegrees[k].gv);
      where[i] = to;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d from %3d to %3d. Gain: [%4d %4d]. Cut: %6d, Vol: %6d\n",
            i, from, to, xgain+myedegrees[k].gv, myedegrees[k].ed-myrinfo->id, graph->mincut, graph->minvol));

      KWayVolUpdate(ctrl, graph, i, from, to, marker, phtable, updind);

      nmoves++;

      /*CheckVolKWayPartitionParams(ctrl, graph, nparts); */
    }

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d, Vol: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut,
               graph->minvol));

  }

  GKfree((void **) &marker, (void **) &updind, (void **) &phtable, LTERM);

  PQueueFree(ctrl, &queue);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Greedy_KWayVolBalanceMConn(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts,
                                float ubfactor, int npasses)
{
  int i, ii, j, k, l, pass, nvtxs, nmoves, tvwgt, myndegrees, xgain;
  int from, me, to, vwgt, maxndoms, nadd;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *moved, *bndptr, *bndind, *minwgt, *maxwgt, *itpwgts, *updind, *marker, *phtable;
  idxtype *pmat, *pmatptr, *ndoms;
  VEDegreeType *myedegrees;
  VRInfoType *myrinfo;
  PQueueType queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  pwgts = graph->pwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  updind = idxmalloc(nvtxs, "Random_KWayVolRefine: updind");
  marker = idxsmalloc(nvtxs, 0, "Random_KWayVolRefine: marker");
  phtable = idxsmalloc(nparts, -1, "Random_KWayVolRefine: phtable");

  pmat = ctrl->wspace.pmat;
  ndoms = idxwspacemalloc(ctrl, nparts);

  ComputeVolSubDomainGraph(graph, nparts, pmat, ndoms);

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxwspacemalloc(ctrl, nvtxs);

  PQueueInit(ctrl, &queue, nvtxs, graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)]);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("VolPart: [%5d %5d]-[%5d %5d], Balance: %3.2f, Nv-Nb[%5d %5d]. Cut: %5d, Vol: %5d [B]\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut, graph->minvol));


  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);
    /* Check to see if things are out of balance, given the tolerance */
    for (i=0; i<nparts; i++) {
      if (pwgts[i] > maxwgt[i])
        break;
    }
    if (i == nparts) /* Things are balanced. Return right away */
      break;

    PQueueReset(&queue);
    idxset(nvtxs, -1, moved);

    RandomPermute(graph->nbnd, perm, 1);
    for (ii=0; ii<graph->nbnd; ii++) {
      i = bndind[perm[ii]];
      PQueueInsert(&queue, i, graph->vrinfo[i].gv);
      moved[i] = 2;
    }

    maxndoms = ndoms[idxamax(nparts, ndoms)];

    for (nmoves=0;;) {
      if ((i = PQueueGetMax(&queue)) == -1)
        break;
      moved[i] = 1;

      myrinfo = graph->vrinfo+i;
      from = where[i];
      vwgt = graph->vwgt[i];

      if (pwgts[from]-vwgt < minwgt[from])
        continue;   /* This cannot be moved! */

      xgain = (myrinfo->id == 0 && myrinfo->ed > 0 ? graph->vsize[i] : 0);

      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      /* Determine the valid domains */
      for (j=0; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        phtable[to] = 1;
        pmatptr = pmat + to*nparts;
        for (nadd=0, k=0; k<myndegrees; k++) {
          if (k == j)
            continue;

          l = myedegrees[k].pid;
          if (pmatptr[l] == 0) {
            if (ndoms[l] > maxndoms-1) {
              phtable[to] = 0;
              nadd = maxndoms;
              break;
            }
            nadd++;
          }
        }
        if (ndoms[to]+nadd > maxndoms)
          phtable[to] = 0;
      }

      for (k=0; k<myndegrees; k++) {
        to = myedegrees[k].pid;
        if (!phtable[to])
          continue;
        if (pwgts[to]+vwgt <= maxwgt[to] ||
            itpwgts[from]*(pwgts[to]+vwgt) <= itpwgts[to]*pwgts[from])
          break;
      }
      if (k == myndegrees)
        continue;  /* break out if you did not find a candidate */

      for (j=k+1; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        if (!phtable[to])
          continue;
        if (itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid])
          k = j;
      }

      to = myedegrees[k].pid;

      for (j=0; j<myndegrees; j++)
        phtable[myedegrees[j].pid] = -1;

      if (pwgts[from] < maxwgt[from] && pwgts[to] > minwgt[to] &&
          (xgain+myedegrees[k].gv < 0 ||
           (xgain+myedegrees[k].gv == 0 &&  myedegrees[k].ed-myrinfo->id < 0))
         )
        continue;


      /*=====================================================================
      * If we got here, we can now move the vertex from 'from' to 'to'
      *======================================================================*/
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      graph->mincut -= myedegrees[k].ed-myrinfo->id;
      graph->minvol -= (xgain+myedegrees[k].gv);
      where[i] = to;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d from %3d to %3d. Gain: [%4d %4d]. Cut: %6d, Vol: %6d\n",
            i, from, to, xgain+myedegrees[k].gv, myedegrees[k].ed-myrinfo->id, graph->mincut, graph->minvol));

      /* Update pmat to reflect the move of 'i' */
      pmat[from*nparts+to] += (myrinfo->id-myedegrees[k].ed);
      pmat[to*nparts+from] += (myrinfo->id-myedegrees[k].ed);
      if (pmat[from*nparts+to] == 0) {
        ndoms[from]--;
        if (ndoms[from]+1 == maxndoms)
          maxndoms = ndoms[idxamax(nparts, ndoms)];
      }
      if (pmat[to*nparts+from] == 0) {
        ndoms[to]--;
        if (ndoms[to]+1 == maxndoms)
          maxndoms = ndoms[idxamax(nparts, ndoms)];
      }

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        me = where[ii];

        /* Update pmat to reflect the move of 'i' for domains other than 'from' and 'to' */
        if (me != from && me != to) {
          pmat[me*nparts+from] -= adjwgt[j];
          pmat[from*nparts+me] -= adjwgt[j];
          if (pmat[me*nparts+from] == 0) {
            ndoms[me]--;
            if (ndoms[me]+1 == maxndoms)
              maxndoms = ndoms[idxamax(nparts, ndoms)];
          }
          if (pmat[from*nparts+me] == 0) {
            ndoms[from]--;
            if (ndoms[from]+1 == maxndoms)
              maxndoms = ndoms[idxamax(nparts, ndoms)];
          }

          if (pmat[me*nparts+to] == 0) {
            ndoms[me]++;
            if (ndoms[me] > maxndoms) {
              printf("You just increased the maxndoms: %d %d\n", ndoms[me], maxndoms);
              maxndoms = ndoms[me];
            }
          }
          if (pmat[to*nparts+me] == 0) {
            ndoms[to]++;
            if (ndoms[to] > maxndoms) {
              printf("You just increased the maxndoms: %d %d\n", ndoms[to], maxndoms);
              maxndoms = ndoms[to];
            }
          }
          pmat[me*nparts+to] += adjwgt[j];
          pmat[to*nparts+me] += adjwgt[j];
        }
      }

      KWayVolUpdate(ctrl, graph, i, from, to, marker, phtable, updind);

      nmoves++;

      /*CheckVolKWayPartitionParams(ctrl, graph, nparts); */
    }

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d, Vol: %6d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut,
               graph->minvol));

  }

  GKfree((void **) &marker, (void **) &updind, (void **) &phtable, LTERM);

  PQueueFree(ctrl, &queue);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}




/*************************************************************************
* This function updates the edge and volume gains as a result of moving
* v from 'from' to 'to'.
* The working arrays marker and phtable are assumed to be initialized to
* -1, and they left to -1 upon return
**************************************************************************/
void KWayVolUpdate(CtrlType *ctrl, GraphType *graph, int v, int from, int to,
                   idxtype *marker, idxtype *phtable, idxtype *updind)
{
  int ii, j, jj, k, kk, u, nupd, other, me, myidx;
  idxtype *xadj, *vsize, *adjncy, *adjwgt, *where;
  VEDegreeType *myedegrees, *oedegrees;
  VRInfoType *myrinfo, *orinfo;

  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  vsize = graph->vsize;
  where = graph->where;

  myrinfo = graph->vrinfo+v;
  myedegrees = myrinfo->edegrees;


  /*======================================================================
   * Remove the contributions on the gain made by 'v'.
   *=====================================================================*/
  for (k=0; k<myrinfo->ndegrees; k++)
    phtable[myedegrees[k].pid] = k;
  phtable[from] = k;

  myidx = phtable[to];  /* Keep track of the index in myedegrees of the 'to' domain */

  for (j=xadj[v]; j<xadj[v+1]; j++) {
    ii = adjncy[j];
    other = where[ii];
    orinfo = graph->vrinfo+ii;
    oedegrees = orinfo->edegrees;

    if (other == from) {
      for (k=0; k<orinfo->ndegrees; k++) {
        if (phtable[oedegrees[k].pid] == -1)
          oedegrees[k].gv += vsize[v];
      }
    }
    else {
      ASSERT(phtable[other] != -1);

      if (myedegrees[phtable[other]].ned > 1) {
        for (k=0; k<orinfo->ndegrees; k++) {
          if (phtable[oedegrees[k].pid] == -1)
            oedegrees[k].gv += vsize[v];
        }
      }
      else { /* There is only one connection */
        for (k=0; k<orinfo->ndegrees; k++) {
          if (phtable[oedegrees[k].pid] != -1)
            oedegrees[k].gv -= vsize[v];
        }
      }
    }
  }

  for (k=0; k<myrinfo->ndegrees; k++)
    phtable[myedegrees[k].pid] = -1;
  phtable[from] = -1;


  /*======================================================================
   * Update the id/ed of vertex 'v'
   *=====================================================================*/
  myrinfo->ed += myrinfo->id-myedegrees[myidx].ed;
  SWAP(myrinfo->id, myedegrees[myidx].ed, j);
  SWAP(myrinfo->nid, myedegrees[myidx].ned, j);
  if (myedegrees[myidx].ed == 0)
    myedegrees[myidx] = myedegrees[--myrinfo->ndegrees];
  else
    myedegrees[myidx].pid = from;

  /*======================================================================
   * Update the degrees of adjacent vertices and their volume gains
   *=====================================================================*/
  marker[v] = 1;
  updind[0] = v;
  nupd = 1;
  for (j=xadj[v]; j<xadj[v+1]; j++) {
    ii = adjncy[j];
    me = where[ii];

    if (!marker[ii]) {  /* The marking is done for boundary and max gv calculations */
      marker[ii] = 2;
      updind[nupd++] = ii;
    }

    myrinfo = graph->vrinfo+ii;
    if (myrinfo->edegrees == NULL) {
      myrinfo->edegrees = ctrl->wspace.vedegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
    }
    myedegrees = myrinfo->edegrees;

    if (me == from) {
      INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);
      myrinfo->nid--;
    }
    else if (me == to) {
      INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);
      myrinfo->nid++;
    }

    /* Remove the edgeweight from the 'pid == from' entry of the vertex */
    if (me != from) {
      for (k=0; k<myrinfo->ndegrees; k++) {
        if (myedegrees[k].pid == from) {
          if (myedegrees[k].ned == 1) {
            myedegrees[k] = myedegrees[--myrinfo->ndegrees];
            marker[ii] = 1;  /* You do a complete .gv calculation */

            /* All vertices adjacent to 'ii' need to be updated */
            for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
              u = adjncy[jj];
              other = where[u];
              orinfo = graph->vrinfo+u;
              oedegrees = orinfo->edegrees;

              for (kk=0; kk<orinfo->ndegrees; kk++) {
                if (oedegrees[kk].pid == from) {
                  oedegrees[kk].gv -= vsize[ii];
                  break;
                }
              }
            }
          }
          else {
            myedegrees[k].ed -= adjwgt[j];
            myedegrees[k].ned--;

            /* Update the gv due to single 'ii' connection to 'from' */
            if (myedegrees[k].ned == 1) {
              /* find the vertex 'u' that 'ii' was connected into 'from' */
              for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
                u = adjncy[jj];
                other = where[u];
                orinfo = graph->vrinfo+u;
                oedegrees = orinfo->edegrees;

                if (other == from) {
                  for (kk=0; kk<orinfo->ndegrees; kk++)
                    oedegrees[kk].gv += vsize[ii];
                  break;
                }
              }
            }
          }

          break;
        }
      }
    }

    /* Add the edgeweight to the 'pid == to' entry of the vertex */
    if (me != to) {
      for (k=0; k<myrinfo->ndegrees; k++) {
        if (myedegrees[k].pid == to) {
          myedegrees[k].ed += adjwgt[j];
          myedegrees[k].ned++;

          /* Update the gv due to non-single 'ii' connection to 'to' */
          if (myedegrees[k].ned == 2) {
            /* find the vertex 'u' that 'ii' was connected into 'to' */
            for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
              u = adjncy[jj];
              other = where[u];
              orinfo = graph->vrinfo+u;
              oedegrees = orinfo->edegrees;

              if (u != v && other == to) {
                for (kk=0; kk<orinfo->ndegrees; kk++)
                  oedegrees[kk].gv -= vsize[ii];
                break;
              }
            }
          }
          break;
        }
      }

      if (k == myrinfo->ndegrees) {
        myedegrees[myrinfo->ndegrees].pid = to;
        myedegrees[myrinfo->ndegrees].ed = adjwgt[j];
        myedegrees[myrinfo->ndegrees++].ned = 1;
        marker[ii] = 1;  /* You do a complete .gv calculation */

        /* All vertices adjacent to 'ii' need to be updated */
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
          u = adjncy[jj];
          other = where[u];
          orinfo = graph->vrinfo+u;
          oedegrees = orinfo->edegrees;

          for (kk=0; kk<orinfo->ndegrees; kk++) {
            if (oedegrees[kk].pid == to) {
              oedegrees[kk].gv += vsize[ii];
              if (!marker[u]) { /* Need to update boundary etc */
                marker[u] = 2;
                updind[nupd++] = u;
              }
              break;
            }
          }
        }
      }
    }

    ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
  }

  /*======================================================================
   * Add the contributions on the volume gain due to 'v'
   *=====================================================================*/
  myrinfo = graph->vrinfo+v;
  myedegrees = myrinfo->edegrees;
  for (k=0; k<myrinfo->ndegrees; k++)
    phtable[myedegrees[k].pid] = k;
  phtable[to] = k;

  for (j=xadj[v]; j<xadj[v+1]; j++) {
    ii = adjncy[j];
    other = where[ii];
    orinfo = graph->vrinfo+ii;
    oedegrees = orinfo->edegrees;

    if (other == to) {
      for (k=0; k<orinfo->ndegrees; k++) {
        if (phtable[oedegrees[k].pid] == -1)
          oedegrees[k].gv -= vsize[v];
      }
    }
    else {
      ASSERT(phtable[other] != -1);

      if (myedegrees[phtable[other]].ned > 1) {
        for (k=0; k<orinfo->ndegrees; k++) {
          if (phtable[oedegrees[k].pid] == -1)
            oedegrees[k].gv -= vsize[v];
        }
      }
      else { /* There is only one connection */
        for (k=0; k<orinfo->ndegrees; k++) {
          if (phtable[oedegrees[k].pid] != -1)
            oedegrees[k].gv += vsize[v];
        }
      }
    }
  }
  for (k=0; k<myrinfo->ndegrees; k++)
    phtable[myedegrees[k].pid] = -1;
  phtable[to] = -1;


  /*======================================================================
   * Recompute the volume information of the 'hard' nodes, and update the
   * max volume gain for all the update vertices
   *=====================================================================*/
  ComputeKWayVolume(graph, nupd, updind, marker, phtable);


  /*======================================================================
   * Maintain a consistent boundary
   *=====================================================================*/
  for (j=0; j<nupd; j++) {
    k = updind[j];
    marker[k] = 0;
    myrinfo = graph->vrinfo+k;

    if ((myrinfo->gv >= 0 || myrinfo->ed-myrinfo->id >= 0) && graph->bndptr[k] == -1)
      BNDInsert(graph->nbnd, graph->bndind, graph->bndptr, k);

    if (myrinfo->gv < 0 && myrinfo->ed-myrinfo->id < 0 && graph->bndptr[k] != -1)
      BNDDelete(graph->nbnd, graph->bndind, graph->bndptr, k);
  }

}




/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void ComputeKWayVolume(GraphType *graph, int nupd, idxtype *updind, idxtype *marker, idxtype *phtable)
{
  int ii, iii, i, j, k, kk, nvtxs, me, other;
  idxtype *xadj, *vsize, *adjncy, *adjwgt, *where;
  VRInfoType *rinfo, *myrinfo, *orinfo;
  VEDegreeType *myedegrees, *oedegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  rinfo = graph->vrinfo;


  /*------------------------------------------------------------
  / Compute now the iv/ev degrees
  /------------------------------------------------------------*/
  for (iii=0; iii<nupd; iii++) {
    i = updind[iii];
    me = where[i];

    myrinfo = rinfo+i;
    myedegrees = myrinfo->edegrees;

    if (marker[i] == 1) {  /* Only complete gain updates go through */
      for (k=0; k<myrinfo->ndegrees; k++)
        myedegrees[k].gv = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        other = where[ii];
        orinfo = rinfo+ii;
        oedegrees = orinfo->edegrees;

        for (kk=0; kk<orinfo->ndegrees; kk++)
          phtable[oedegrees[kk].pid] = kk;
        phtable[other] = 1;

        if (me == other) {
          /* Find which domains 'i' is connected and 'ii' is not and update their gain */
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (phtable[myedegrees[k].pid] == -1)
              myedegrees[k].gv -= vsize[ii];
          }
        }
        else {
          ASSERT(phtable[me] != -1);

          /* I'm the only connection of 'ii' in 'me' */
          if (oedegrees[phtable[me]].ned == 1) {
            /* Increase the gains for all the common domains between 'i' and 'ii' */
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (phtable[myedegrees[k].pid] != -1)
                myedegrees[k].gv += vsize[ii];
            }
          }
          else {
            /* Find which domains 'i' is connected and 'ii' is not and update their gain */
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (phtable[myedegrees[k].pid] == -1)
                myedegrees[k].gv -= vsize[ii];
            }
          }
        }

        for (kk=0; kk<orinfo->ndegrees; kk++)
          phtable[oedegrees[kk].pid] = -1;
        phtable[other] = -1;

      }
    }

    myrinfo->gv = -MAXIDX;
    for (k=0; k<myrinfo->ndegrees; k++) {
      if (myedegrees[k].gv > myrinfo->gv)
        myrinfo->gv = myedegrees[k].gv;
    }
    if (myrinfo->ed > 0 && myrinfo->id == 0)
      myrinfo->gv += vsize[i];

  }

}



/*************************************************************************
* This function computes the total volume
**************************************************************************/
int ComputeVolume(GraphType *graph, idxtype *where)
{
  int i, j, k, nvtxs, nparts, totalv;
  idxtype *xadj, *adjncy, *vsize, *marker;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vsize = (graph->vsize == NULL ? graph->vwgt : graph->vsize);

  nparts = where[idxamax(nvtxs, where)]+1;
  marker = idxsmalloc(nparts, -1, "ComputeVolume: marker");

  totalv = 0;

  for (i=0; i<nvtxs; i++) {
    marker[where[i]] = i;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = where[adjncy[j]];
      if (marker[k] != i) {
        marker[k] = i;
        totalv += vsize[i];
      }
    }
  }

  free(marker);

  return totalv;
}





/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void CheckVolKWayPartitionParams(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, ii, j, k, kk, nvtxs, me, other, pid;
  idxtype *xadj, *vsize, *adjncy, *adjwgt, *where;
  VRInfoType *rinfo, *myrinfo, *orinfo, tmprinfo;
  VEDegreeType *myedegrees, *oedegrees, *tmpdegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  rinfo = graph->vrinfo;

  tmpdegrees = (VEDegreeType *)GKmalloc(nparts*sizeof(VEDegreeType), "CheckVolKWayPartitionParams: tmpdegrees");

  /*------------------------------------------------------------
  / Compute now the iv/ev degrees
  /------------------------------------------------------------*/
  for (i=0; i<nvtxs; i++) {
    me = where[i];

    myrinfo = rinfo+i;
    myedegrees = myrinfo->edegrees;

    for (k=0; k<myrinfo->ndegrees; k++)
      tmpdegrees[k] = myedegrees[k];

    tmprinfo.ndegrees = myrinfo->ndegrees;
    tmprinfo.id = myrinfo->id;
    tmprinfo.ed = myrinfo->ed;

    myrinfo = &tmprinfo;
    myedegrees = tmpdegrees;


    for (k=0; k<myrinfo->ndegrees; k++)
      myedegrees[k].gv = 0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ii = adjncy[j];
      other = where[ii];
      orinfo = rinfo+ii;
      oedegrees = orinfo->edegrees;

      if (me == other) {
        /* Find which domains 'i' is connected and 'ii' is not and update their gain */
        for (k=0; k<myrinfo->ndegrees; k++) {
          pid = myedegrees[k].pid;
          for (kk=0; kk<orinfo->ndegrees; kk++) {
            if (oedegrees[kk].pid == pid)
              break;
          }
          if (kk == orinfo->ndegrees)
            myedegrees[k].gv -= vsize[ii];
        }
      }
      else {
        /* Find the orinfo[me].ed and see if I'm the only connection */
        for (k=0; k<orinfo->ndegrees; k++) {
          if (oedegrees[k].pid == me)
            break;
        }

        if (oedegrees[k].ned == 1) { /* I'm the only connection of 'ii' in 'me' */
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == other) {
              myedegrees[k].gv += vsize[ii];
              break;
            }
          }

          /* Increase the gains for all the common domains between 'i' and 'ii' */
          for (k=0; k<myrinfo->ndegrees; k++) {
            if ((pid = myedegrees[k].pid) == other)
              continue;
            for (kk=0; kk<orinfo->ndegrees; kk++) {
              if (oedegrees[kk].pid == pid) {
                myedegrees[k].gv += vsize[ii];
                break;
              }
            }
          }

        }
        else {
          /* Find which domains 'i' is connected and 'ii' is not and update their gain */
          for (k=0; k<myrinfo->ndegrees; k++) {
            if ((pid = myedegrees[k].pid) == other)
              continue;
            for (kk=0; kk<orinfo->ndegrees; kk++) {
              if (oedegrees[kk].pid == pid)
                break;
            }
            if (kk == orinfo->ndegrees)
              myedegrees[k].gv -= vsize[ii];
          }
        }
      }
    }

    myrinfo = rinfo+i;
    myedegrees = myrinfo->edegrees;

    for (k=0; k<myrinfo->ndegrees; k++) {
      pid = myedegrees[k].pid;
      for (kk=0; kk<tmprinfo.ndegrees; kk++) {
        if (tmpdegrees[kk].pid == pid) {
          if (tmpdegrees[kk].gv != myedegrees[k].gv)
            printf("[%d %d %d %d]\n", i, pid, myedegrees[k].gv, tmpdegrees[kk].gv);
          break;
        }
      }
    }

  }

  free(tmpdegrees);

}


/*************************************************************************
* This function computes the subdomain graph
**************************************************************************/
void ComputeVolSubDomainGraph(GraphType *graph, int nparts, idxtype *pmat, idxtype *ndoms)
{
  int i, j, k, me, nvtxs, ndegrees;
  idxtype *xadj, *adjncy, *adjwgt, *where;
  VRInfoType *rinfo;
  VEDegreeType *edegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  rinfo = graph->vrinfo;

  idxset(nparts*nparts, 0, pmat);

  for (i=0; i<nvtxs; i++) {
    if (rinfo[i].ed > 0) {
      me = where[i];
      ndegrees = rinfo[i].ndegrees;
      edegrees = rinfo[i].edegrees;

      k = me*nparts;
      for (j=0; j<ndegrees; j++)
        pmat[k+edegrees[j].pid] += edegrees[j].ed;
    }
  }

  for (i=0; i<nparts; i++) {
    ndoms[i] = 0;
    for (j=0; j<nparts; j++) {
      if (pmat[i*nparts+j] > 0)
        ndoms[i]++;
    }
  }
}



/*************************************************************************
* This function computes the subdomain graph
**************************************************************************/
void EliminateVolSubDomainEdges(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts)
{
  int i, ii, j, k, me, other, nvtxs, total, max, avg, totalout, nind, ncand, ncand2, target, target2, nadd;
  int min, move, cpwgt, tvwgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt, *pwgts, *where, *maxpwgt, *pmat, *ndoms, *mypmat, *otherpmat, *ind;
  KeyValueType *cand, *cand2;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = idxset(nparts, 0, graph->pwgts);

  maxpwgt = idxwspacemalloc(ctrl, nparts);
  ndoms = idxwspacemalloc(ctrl, nparts);
  otherpmat = idxwspacemalloc(ctrl, nparts);
  ind = idxwspacemalloc(ctrl, nvtxs);
  pmat = idxset(nparts*nparts, 0, ctrl->wspace.pmat);

  cand = (KeyValueType *)GKmalloc(nparts*sizeof(KeyValueType), "EliminateSubDomainEdges: cand");
  cand2 = (KeyValueType *)GKmalloc(nparts*sizeof(KeyValueType), "EliminateSubDomainEdges: cand");

  /* Compute the pmat matrix */
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] != me)
        pmat[me*nparts+where[k]] += adjwgt[j];
    }
  }

  /* Compute the maximum allowed weight for each domain */
  tvwgt = idxsum(nparts, pwgts);
  for (i=0; i<nparts; i++)
    maxpwgt[i] = 1.25*tpwgts[i]*tvwgt;

  /* Determine the domain connectivity */
  for (i=0; i<nparts; i++) {
    for (k=0, j=0; j<nparts; j++) {
      if (pmat[i*nparts+j] > 0)
        k++;
    }
    ndoms[i] = k;
  }

  /* Get into the loop eliminating subdomain connections */
  for (;;) {
    total = idxsum(nparts, ndoms);
    avg = total/nparts;
    max = ndoms[idxamax(nparts, ndoms)];

    /* printf("Adjacent Subdomain Stats: Total: %3d, Max: %3d, Avg: %3d\n", total, max, avg); */

    if (max < 1.5*avg)
      break;

    me = idxamax(nparts, ndoms);
    mypmat = pmat + me*nparts;
    totalout = idxsum(nparts, mypmat);

    /*printf("Me: %d, TotalOut: %d,\n", me, totalout);*/

    /* Sort the connections according to their cut */
    for (ncand2=0, i=0; i<nparts; i++) {
      if (mypmat[i] > 0) {
        cand2[ncand2].key = mypmat[i];
        cand2[ncand2++].val = i;
      }
    }
    ikeysort(ncand2, cand2);

    move = 0;
    for (min=0; min<ncand2; min++) {
      if (cand2[min].key > totalout/(2*ndoms[me]))
        break;

      other = cand2[min].val;

      /*printf("\tMinOut: %d to %d\n", mypmat[other], other);*/

      idxset(nparts, 0, otherpmat);

      /* Go and find the vertices in 'other' that are connected in 'me' */
      for (nind=0, i=0; i<nvtxs; i++) {
        if (where[i] == other) {
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            if (where[adjncy[j]] == me) {
              ind[nind++] = i;
              break;
            }
          }
        }
      }

      /* Go and construct the otherpmat to see where these nind vertices are connected to */
      for (cpwgt=0, ii=0; ii<nind; ii++) {
        i = ind[ii];
        cpwgt += vwgt[i];

        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (where[k] != other)
            otherpmat[where[k]] += adjwgt[j];
        }
      }

      for (ncand=0, i=0; i<nparts; i++) {
        if (otherpmat[i] > 0) {
          cand[ncand].key = -otherpmat[i];
          cand[ncand++].val = i;
        }
      }
      ikeysort(ncand, cand);

      /*
       * Go through and the select the first domain that is common with 'me', and
       * does not increase the ndoms[target] higher than my ndoms, subject to the
       * maxpwgt constraint. Traversal is done from the mostly connected to the least.
       */
      target = target2 = -1;
      for (i=0; i<ncand; i++) {
        k = cand[i].val;

        if (mypmat[k] > 0) {
          if (pwgts[k] + cpwgt > maxpwgt[k])  /* Check if balance will go off */
            continue;

          for (j=0; j<nparts; j++) {
            if (otherpmat[j] > 0 && ndoms[j] >= ndoms[me]-1 && pmat[nparts*j+k] == 0)
              break;
          }
          if (j == nparts) { /* No bad second level effects */
            for (nadd=0, j=0; j<nparts; j++) {
              if (otherpmat[j] > 0 && pmat[nparts*k+j] == 0)
                nadd++;
            }

            /*printf("\t\tto=%d, nadd=%d, %d\n", k, nadd, ndoms[k]);*/
            if (target2 == -1 && ndoms[k]+nadd < ndoms[me]) {
              target2 = k;
            }
            if (nadd == 0) {
              target = k;
              break;
            }
          }
        }
      }
      if (target == -1 && target2 != -1)
        target = target2;

      if (target == -1) {
        /* printf("\t\tCould not make the move\n");*/
        continue;
      }

      /*printf("\t\tMoving to %d\n", target);*/

      /* Update the partition weights */
      INC_DEC(pwgts[target], pwgts[other], cpwgt);

      /* Set all nind vertices to belong to 'target' */
      for (ii=0; ii<nind; ii++) {
        i = ind[ii];
        where[i] = target;

        /* First remove any contribution that this vertex may have made */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (where[k] != other) {
            if (pmat[nparts*other + where[k]] == 0)
              printf("Something wrong\n");
            pmat[nparts*other + where[k]] -= adjwgt[j];
            if (pmat[nparts*other + where[k]] == 0)
              ndoms[other]--;

            if (pmat[nparts*where[k] + other] == 0)
              printf("Something wrong\n");
            pmat[nparts*where[k] + other] -= adjwgt[j];
            if (pmat[nparts*where[k] + other] == 0)
              ndoms[where[k]]--;
          }
        }

        /* Next add the new contributions as a result of the move */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (where[k] != target) {
            if (pmat[nparts*target + where[k]] == 0)
              ndoms[target]++;
            pmat[nparts*target + where[k]] += adjwgt[j];

            if (pmat[nparts*where[k] + target] == 0)
              ndoms[where[k]]++;
            pmat[nparts*where[k] + target] += adjwgt[j];
          }
        }
      }

      move = 1;
      break;
    }

    if (move == 0)
      break;
  }

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);

  GKfree((void **) &cand, (void **) &cand2, LTERM);
}



/*************************************************************************
* This function finds all the connected components induced by the
* partitioning vector in wgraph->where and tries to push them around to
* remove some of them
**************************************************************************/
void EliminateVolComponents(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor)
{
  int i, ii, j, jj, k, me, nvtxs, tvwgt, first, last, nleft, ncmps, cwgt, ncand, other, target, deltawgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt, *where, *pwgts, *maxpwgt;
  idxtype *cpvec, *touched, *perm, *todo, *cind, *cptr, *npcmps;
  KeyValueType *cand;
  int recompute=0;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = idxset(nparts, 0, graph->pwgts);

  touched = idxset(nvtxs, 0, idxwspacemalloc(ctrl, nvtxs));
  cptr = idxwspacemalloc(ctrl, nvtxs);
  cind = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  todo = idxwspacemalloc(ctrl, nvtxs);
  maxpwgt = idxwspacemalloc(ctrl, nparts);
  cpvec = idxwspacemalloc(ctrl, nparts);
  npcmps = idxset(nparts, 0, idxwspacemalloc(ctrl, nparts));

  for (i=0; i<nvtxs; i++)
    perm[i] = todo[i] = i;

  /* Find the connected componends induced by the partition */
  ncmps = -1;
  first = last = 0;
  nleft = nvtxs;
  while (nleft > 0) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      ASSERT(touched[todo[0]] == 0);
      i = todo[0];
      cind[last++] = i;
      touched[i] = 1;
      me = where[i];
      npcmps[me]++;
    }

    i = cind[first++];
    k = perm[i];
    j = todo[k] = todo[--nleft];
    perm[j] = k;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] == me && !touched[k]) {
        cind[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  /* printf("I found %d components, for this %d-way partition\n", ncmps, nparts); */

  if (ncmps > nparts) { /* There are more components than processors */
    cand = (KeyValueType *)GKmalloc(nparts*sizeof(KeyValueType), "EliminateSubDomainEdges: cand");

    /* First determine the partition sizes and max allowed load imbalance */
    for (i=0; i<nvtxs; i++)
      pwgts[where[i]] += vwgt[i];
    tvwgt = idxsum(nparts, pwgts);
    for (i=0; i<nparts; i++)
      maxpwgt[i] = ubfactor*tpwgts[i]*tvwgt;

    deltawgt = tvwgt/(100*nparts);
    deltawgt = 5;

    for (i=0; i<ncmps; i++) {
      me = where[cind[cptr[i]]];  /* Get the domain of this component */
      if (npcmps[me] == 1)
        continue;  /* Skip it because it is contigous */

      /*printf("Trying to move %d from %d\n", i, me); */

      /* Determine the connectivity */
      idxset(nparts, 0, cpvec);
      for (cwgt=0, j=cptr[i]; j<cptr[i+1]; j++) {
        ii = cind[j];
        cwgt += vwgt[ii];
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
          other = where[adjncy[jj]];
          if (me != other)
            cpvec[other] += adjwgt[jj];
        }
      }

      /*printf("\tCmp weight: %d\n", cwgt);*/

      if (cwgt > .30*pwgts[me])
        continue;  /* Skip the component if it is over 30% of the weight */

      for (ncand=0, j=0; j<nparts; j++) {
        if (cpvec[j] > 0) {
          cand[ncand].key = -cpvec[j];
          cand[ncand++].val = j;
        }
      }
      if (ncand == 0)
        continue;

      ikeysort(ncand, cand);

      target = -1;
      for (j=0; j<ncand; j++) {
        k = cand[j].val;
        if (cwgt < deltawgt || pwgts[k] + cwgt < maxpwgt[k]) {
          target = k;
          break;
        }
      }

      /*printf("\tMoving it to %d [%d]\n", target, cpvec[target]);*/

      if (target != -1) {
        /* Assign all the vertices of 'me' to 'target' and update data structures */
        pwgts[me] -= cwgt;
        pwgts[target] += cwgt;
        npcmps[me]--;

        for (j=cptr[i]; j<cptr[i+1]; j++)
          where[cind[j]] = target;

        graph->mincut -= cpvec[target];
        recompute = 1;
      }
    }

    free(cand);
  }

  if (recompute) {
    int ttlv;
    idxtype *marker;

    marker = idxset(nparts, -1, cpvec);
    for (ttlv=0, i=0; i<nvtxs; i++) {
      marker[where[i]] = i;
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (marker[where[adjncy[j]]] != i) {
          ttlv += graph->vsize[i];
          marker[where[adjncy[j]]] = i;
        }
      }
    }
    graph->minvol = ttlv;
  }

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kwayvolrefine.c
 *
 * This file contains the driving routines for multilevel k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void RefineVolKWay(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int nparts,
                   float *tpwgts, float ubfactor)
{
  int i, nlevels;
  GraphType *ptr;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  /* Take care any non-contiguity */
  if (ctrl->RType == RTYPE_KWAYRANDOM_MCONN) {
    ComputeVolKWayPartitionParams(ctrl, graph, nparts);
    EliminateVolComponents(ctrl, graph, nparts, tpwgts, 1.25);
    EliminateVolSubDomainEdges(ctrl, graph, nparts, tpwgts);
    EliminateVolComponents(ctrl, graph, nparts, tpwgts, 1.25);
  }


  /* Determine how many levels are there */
  for (ptr=graph, nlevels=0; ptr!=orggraph; ptr=ptr->finer, nlevels++);

  /* Compute the parameters of the coarsest graph */
  ComputeVolKWayPartitionParams(ctrl, graph, nparts);

  for (i=0; ;i++) {
    /*PrintSubDomainGraph(graph, nparts, graph->where);*/
    MALLOC_CHECK(NULL);
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));

    if (2*i >= nlevels && !IsBalanced(graph->pwgts, nparts, tpwgts, 1.04*ubfactor)) {
      ComputeVolKWayBalanceBoundary(ctrl, graph, nparts);
      switch (ctrl->RType) {
        case RTYPE_KWAYRANDOM:
          Greedy_KWayVolBalance(ctrl, graph, nparts, tpwgts, ubfactor, 1);
          break;
        case RTYPE_KWAYRANDOM_MCONN:
          Greedy_KWayVolBalanceMConn(ctrl, graph, nparts, tpwgts, ubfactor, 1);
          break;
      }
      ComputeVolKWayBoundary(ctrl, graph, nparts);
    }

    switch (ctrl->RType) {
      case RTYPE_KWAYRANDOM:
        Random_KWayVolRefine(ctrl, graph, nparts, tpwgts, ubfactor, 10, 0);
        break;
      case RTYPE_KWAYRANDOM_MCONN:
        Random_KWayVolRefineMConn(ctrl, graph, nparts, tpwgts, ubfactor, 10, 0);
        break;
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    GKfree((void **) &graph->gdata, LTERM);  /* Deallocate the graph related arrays */

    graph = graph->finer;

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    ProjectVolKWayPartition(ctrl, graph, nparts);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  if (!IsBalanced(graph->pwgts, nparts, tpwgts, ubfactor)) {
    ComputeVolKWayBalanceBoundary(ctrl, graph, nparts);
    switch (ctrl->RType) {
      case RTYPE_KWAYRANDOM:
        Greedy_KWayVolBalance(ctrl, graph, nparts, tpwgts, ubfactor, 8);
        Random_KWayVolRefine(ctrl, graph, nparts, tpwgts, ubfactor, 10, 0);
        break;
      case RTYPE_KWAYRANDOM_MCONN:
        Greedy_KWayVolBalanceMConn(ctrl, graph, nparts, tpwgts, ubfactor, 8);
        Random_KWayVolRefineMConn(ctrl, graph, nparts, tpwgts, ubfactor, 10, 0);
        break;
    }
  }

  EliminateVolComponents(ctrl, graph, nparts, tpwgts, ubfactor);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}



/*************************************************************************
* This function allocates memory for k-way edge refinement
**************************************************************************/
void AllocateVolKWayPartitionMemory(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int nvtxs, pad64;

  nvtxs = graph->nvtxs;

  pad64 = (3*nvtxs+nparts)%2;

  graph->rdata = idxmalloc(3*nvtxs+nparts+(sizeof(VRInfoType)/sizeof(idxtype))*nvtxs+pad64, "AllocateVolKWayPartitionMemory: rdata");
  graph->pwgts          = graph->rdata;
  graph->where          = graph->rdata + nparts;
  graph->bndptr         = graph->rdata + nvtxs + nparts;
  graph->bndind         = graph->rdata + 2*nvtxs + nparts;
  graph->vrinfo          = (VRInfoType *)(graph->rdata + 3*nvtxs+nparts + pad64);

}



/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void ComputeVolKWayPartitionParams(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, j, k, nvtxs, mincut, me, other;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *pwgts, *where;
  VRInfoType *rinfo, *myrinfo;
  VEDegreeType *myedegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = idxset(nparts, 0, graph->pwgts);
  rinfo = graph->vrinfo;


  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  ctrl->wspace.cdegree = 0;
  mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];

    myrinfo = rinfo+i;
    myrinfo->id = myrinfo->ed = myrinfo->nid = myrinfo->ndegrees = 0;
    myrinfo->edegrees = NULL;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]]) {
        myrinfo->id += adjwgt[j];
        myrinfo->nid++;
      }
    }
    myrinfo->ed = graph->adjwgtsum[i] - myrinfo->id;

    mincut += myrinfo->ed;

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      myedegrees = myrinfo->edegrees = ctrl->wspace.vedegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += xadj[i+1]-xadj[i];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == other) {
              myedegrees[k].ed += adjwgt[j];
              myedegrees[k].ned++;
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].gv = 0;
            myedegrees[myrinfo->ndegrees].pid = other;
            myedegrees[myrinfo->ndegrees].ed = adjwgt[j];
            myedegrees[myrinfo->ndegrees++].ned = 1;
          }
        }
      }

      ASSERT(myrinfo->ndegrees <= xadj[i+1]-xadj[i]);
    }
  }
  graph->mincut = mincut/2;


  ComputeKWayVolGains(ctrl, graph, nparts);

}



/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void ComputeKWayVolGains(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, ii, j, k, kk, nvtxs, me, other, myndegrees;
  idxtype *xadj, *vsize, *adjncy, *adjwgt, *where, *bndind, *bndptr, *ophtable;
  VRInfoType *rinfo, *myrinfo, *orinfo;
  VEDegreeType *myedegrees, *oedegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  rinfo = graph->vrinfo;

  ophtable = idxset(nparts, -1, idxwspacemalloc(ctrl, nparts));

  /*------------------------------------------------------------
  / Compute now the iv/ev degrees
  /------------------------------------------------------------*/
  graph->minvol = graph->nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    myrinfo = rinfo+i;
    myrinfo->gv = -MAXIDX;

    if (myrinfo->ndegrees > 0) {
      me = where[i];
      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      graph->minvol += myndegrees*vsize[i];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        other = where[ii];
        orinfo = rinfo+ii;
        oedegrees = orinfo->edegrees;

        for (k=0; k<orinfo->ndegrees; k++)
          ophtable[oedegrees[k].pid] = k;
        ophtable[other] = 1;  /* this is to simplify coding */

        if (me == other) {
          /* Find which domains 'i' is connected and 'ii' is not and update their gain */
          for (k=0; k<myndegrees; k++) {
            if (ophtable[myedegrees[k].pid] == -1)
              myedegrees[k].gv -= vsize[ii];
          }
        }
        else {
          ASSERT(ophtable[me] != -1);

          if (oedegrees[ophtable[me]].ned == 1) { /* I'm the only connection of 'ii' in 'me' */
            /* Increase the gains for all the common domains between 'i' and 'ii' */
            for (k=0; k<myndegrees; k++) {
              if (ophtable[myedegrees[k].pid] != -1)
                myedegrees[k].gv += vsize[ii];
            }
          }
          else {
            /* Find which domains 'i' is connected and 'ii' is not and update their gain */
            for (k=0; k<myndegrees; k++) {
              if (ophtable[myedegrees[k].pid] == -1)
                myedegrees[k].gv -= vsize[ii];
            }
          }
        }

        for (kk=0; kk<orinfo->ndegrees; kk++)
          ophtable[oedegrees[kk].pid] = -1;
        ophtable[other] = -1;
      }

      /* Compute the max vgain */
      for (k=0; k<myndegrees; k++) {
        if (myedegrees[k].gv > myrinfo->gv)
          myrinfo->gv = myedegrees[k].gv;
      }
    }

    if (myrinfo->ed > 0 && myrinfo->id == 0)
      myrinfo->gv += vsize[i];

    if (myrinfo->gv >= 0 || myrinfo->ed-myrinfo->id >= 0)
      BNDInsert(graph->nbnd, bndind, bndptr, i);
  }

  idxwspacefree(ctrl, nparts);

}



/*************************************************************************
* This function projects a partition, and at the same time computes the
* parameters for refinement.
**************************************************************************/
void ProjectVolKWayPartition(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, j, k, nvtxs, me, other, istart, iend, ndegrees;
  idxtype *xadj, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *cmap, *where;
  idxtype *cwhere;
  GraphType *cgraph;
  VRInfoType *crinfo, *rinfo, *myrinfo;
  VEDegreeType *myedegrees;
  idxtype *htable;

  cgraph = graph->coarser;
  cwhere = cgraph->where;
  crinfo = cgraph->vrinfo;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  AllocateVolKWayPartitionMemory(ctrl, graph, nparts);
  where = graph->where;
  rinfo = graph->vrinfo;

  /* Go through and project partition and compute id/ed for the nodes */
  for (i=0; i<nvtxs; i++) {
    k = cmap[i];
    where[i] = cwhere[k];
    cmap[i] = crinfo[k].ed;  /* For optimization */
  }

  htable = idxset(nparts, -1, idxwspacemalloc(ctrl, nparts));

  ctrl->wspace.cdegree = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];

    myrinfo = rinfo+i;
    myrinfo->id = myrinfo->ed = myrinfo->nid = myrinfo->ndegrees = 0;
    myrinfo->edegrees = NULL;

    myrinfo->id = adjwgtsum[i];
    myrinfo->nid = xadj[i+1]-xadj[i];

    if (cmap[i] > 0) { /* If it is an interface node. Note cmap[i] = crinfo[cmap[i]].ed */
      istart = xadj[i];
      iend = xadj[i+1];

      myedegrees = myrinfo->edegrees = ctrl->wspace.vedegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += iend-istart;

      ndegrees = 0;
      for (j=istart; j<iend; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          myrinfo->ed += adjwgt[j];
          myrinfo->nid--;
          if ((k = htable[other]) == -1) {
            htable[other] = ndegrees;
            myedegrees[ndegrees].gv = 0;
            myedegrees[ndegrees].pid = other;
            myedegrees[ndegrees].ed = adjwgt[j];
            myedegrees[ndegrees++].ned = 1;
          }
          else {
            myedegrees[k].ed += adjwgt[j];
            myedegrees[k].ned++;
          }
        }
      }
      myrinfo->id -= myrinfo->ed;

      /* Remove space for edegrees if it was interior */
      if (myrinfo->ed == 0) {
        myrinfo->edegrees = NULL;
        ctrl->wspace.cdegree -= iend-istart;
      }
      else {
        myrinfo->ndegrees = ndegrees;

        for (j=0; j<ndegrees; j++)
          htable[myedegrees[j].pid] = -1;
      }
    }
  }

  ComputeKWayVolGains(ctrl, graph, nparts);

  idxcopy(nparts, cgraph->pwgts, graph->pwgts);
  graph->mincut = cgraph->mincut;

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

  idxwspacefree(ctrl, nparts);

}



/*************************************************************************
* This function computes the boundary definition for balancing
**************************************************************************/
void ComputeVolKWayBoundary(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, nvtxs, nbnd;
  idxtype *bndind, *bndptr;

  nvtxs = graph->nvtxs;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);


  /*------------------------------------------------------------
  / Compute the new boundary
  /------------------------------------------------------------*/
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    if (graph->vrinfo[i].gv >=0 || graph->vrinfo[i].ed-graph->vrinfo[i].id >= 0)
      BNDInsert(nbnd, bndind, bndptr, i);
  }

  graph->nbnd = nbnd;
}

/*************************************************************************
* This function computes the boundary definition for balancing
**************************************************************************/
void ComputeVolKWayBalanceBoundary(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, nvtxs, nbnd;
  idxtype *bndind, *bndptr;

  nvtxs = graph->nvtxs;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);


  /*------------------------------------------------------------
  / Compute the new boundary
  /------------------------------------------------------------*/
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    if (graph->vrinfo[i].ed > 0)
      BNDInsert(nbnd, bndind, bndptr, i);
  }

  graph->nbnd = nbnd;
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * match.c
 *
 * This file contains the code that computes matchings and creates the next
 * level coarse graph.
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_RM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED && vwgt[i]+vwgt[adjncy[j]] <= ctrl->maxvwgt) {
          maxidx = adjncy[j];
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_RM_NVW(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *adjncy;
  idxtype *match, *cmap, *perm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED) {
          maxidx = adjncy[j];
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph_NVW(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_HEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx, maxwgt;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (match[k] == UNMATCHED && maxwgt < adjwgt[j] && vwgt[i]+vwgt[k] <= ctrl->maxvwgt) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_SHEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx, maxwgt, avgdegree;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *degrees, *perm, *tperm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  tperm = idxwspacemalloc(ctrl, nvtxs);
  degrees = idxwspacemalloc(ctrl, nvtxs);

  RandomPermute(nvtxs, tperm, 1);
  avgdegree = 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++)
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

  cnvtxs = 0;

  /* Take care any islands. Islands are matched with non-islands due to coarsening */
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      if (xadj[i] < xadj[i+1])
        break;

      maxidx = i;
      for (j=nvtxs-1; j>ii; j--) {
        k = perm[j];
        if (match[k] == UNMATCHED && xadj[k] < xadj[k+1]) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  /* Continue with normal matching */
  for (; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED && maxwgt < adjwgt[j] && vwgt[i]+vwgt[adjncy[j]] <= ctrl->maxvwgt) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  idxwspacefree(ctrl, nvtxs);  /* degrees */
  idxwspacefree(ctrl, nvtxs);  /* tperm */

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mbalance2.c
 *
 * This file contains code that is used to forcefully balance either
 * bisections or k-sections
 *
 * Started 7/29/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point of the bisection balancing algorithms.
**************************************************************************/
void MocBalance2Way2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  float tvec[MAXNCON];

  Compute2WayHLoadImbalanceVec(graph->ncon, graph->npwgts, tpwgts, tvec);
  if (!AreAllBelow(graph->ncon, tvec, ubvec))
    MocGeneral2WayBalance2(ctrl, graph, tpwgts, ubvec);
}



/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void MocGeneral2WayBalance2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  int i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, limit, tmp, cnum;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *moved, *swaps, *perm, *qnum;
  float *nvwgt, *npwgts, origbal[MAXNCON], minbal[MAXNCON], newbal[MAXNCON];
  PQueueType parts[MAXNCON][2];
  int higain, oldgain, mincut, newcut, mincutorder;
  float *maxwgt, *minwgt, tvec[MAXNCON];


  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  limit = amin(amax(0.01*nvtxs, 15), 100);

  /* Setup the weight intervals of the two subdomains */
  minwgt = fwspacemalloc(ctrl, 2*ncon);
  maxwgt = fwspacemalloc(ctrl, 2*ncon);

  for (i=0; i<2; i++) {
    for (j=0; j<ncon; j++) {
      maxwgt[i*ncon+j] = tpwgts[i]*ubvec[j];
      minwgt[i*ncon+j] = tpwgts[i]*(1.0/ubvec[j]);
    }
  }


  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
  }
  for (i=0; i<nvtxs; i++)
    qnum[i] = samax(ncon, nvwgt+i*ncon);

  Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, origbal);
  for (i=0; i<ncon; i++)
    minbal[i] = origbal[i];

  newcut = mincut = graph->mincut;
  mincutorder = -1;

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Parts: [");
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("] T[%.3f %.3f], Nv-Nb[%5d, %5d]. ICut: %6d, LB: ", tpwgts[0], tpwgts[1],
            graph->nvtxs, graph->nbnd, graph->mincut);
    for (i=0; i<ncon; i++)
      printf("%.3f ", origbal[i]);
    printf("[B]\n");
  }

  idxset(nvtxs, -1, moved);

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  /* Insert all nodes in the priority queues */
  nbnd = graph->nbnd;
  RandomPermute(nvtxs, perm, 1);
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    PQueueInsert(&parts[qnum[i]][where[i]], i, ed[i]-id[i]);
  }


  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if (AreAllBelow(ncon, minbal, ubvec))
      break;

    SelectQueue3(ncon, npwgts, tpwgts, &from, &cnum, parts, maxwgt);
    to = (from+1)%2;

    if (from == -1 || (higain = PQueueGetMax(&parts[cnum][from])) == -1)
      break;

    saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
    newcut -= (ed[higain]-id[higain]);
    Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, newbal);

    if (IsBetter2wayBalance(ncon, newbal, minbal, ubvec) ||
        (IsBetter2wayBalance(ncon, newbal, origbal, ubvec) && newcut < mincut)) {
      mincut = newcut;
      for (i=0; i<ncon; i++)
        minbal[i] = newbal[i];
      mincutorder = nswaps;
    }
    else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
      newcut += (ed[higain]-id[higain]);
      saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
      saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      break;
    }

    where[higain] = to;
    moved[higain] = nswaps;
    swaps[nswaps] = higain;

    if (ctrl->dbglvl&DBG_MOVEINFO) {
      printf("Moved %6d from %d(%d). Gain: %5d, Cut: %5d, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], newcut);
      for (i=0; i<ncon; i++)
        printf("(%.3f, %.3f) ", npwgts[i], npwgts[ncon+i]);

      Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, tvec);
      printf(", LB: ");
      for (i=0; i<ncon; i++)
        printf("%.3f ", tvec[i]);
      if (mincutorder == nswaps)
        printf(" *\n");
      else
        printf("\n");
    }


    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (moved[k] == -1)
        PQueueUpdate(&parts[qnum[k]][where[k]], k, oldgain, ed[k]-id[k]);

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }

  }



  /****************************************************************
  * Roll back computations
  *****************************************************************/
  for (i=0; i<nswaps; i++)
    moved[swaps[i]] = -1;  /* reset moved array */
  for (nswaps--; nswaps>mincutorder; nswaps--) {
    higain = swaps[nswaps];

    to = where[higain] = (where[higain]+1)%2;
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    else if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      if (bndptr[k] != -1 && ed[k] == 0)
        BNDDelete(nbnd, bndind, bndptr, k);
      if (bndptr[k] == -1 && ed[k] > 0)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("\tMincut: %6d at %5d, NBND: %6d, NPwgts: [", mincut, mincutorder, nbnd);
    for (i=0; i<ncon; i++)
      printf("(%.3f, %.3f) ", npwgts[i], npwgts[ncon+i]);
    printf("], LB: ");
    Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, tvec);
    for (i=0; i<ncon; i++)
      printf("%.3f ", tvec[i]);
    printf("\n");
  }

  graph->mincut = mincut;
  graph->nbnd = nbnd;


  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  fwspacefree(ctrl, 2*ncon);
  fwspacefree(ctrl, 2*ncon);

}




/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
void SelectQueue3(int ncon, float *npwgts, float *tpwgts, int *from, int *cnum,
       PQueueType queues[MAXNCON][2], float *maxwgt)
{
  int i, j, maxgain=0;
  float maxdiff=0.0, diff;

  *from = -1;
  *cnum = -1;

  /* First determine the side and the queue, irrespective of the presence of nodes */
  for (j=0; j<2; j++) {
    for (i=0; i<ncon; i++) {
      diff = npwgts[j*ncon+i]-maxwgt[j*ncon+i];
      if (diff >= maxdiff) {
        maxdiff = diff;
        *from = j;
        *cnum = i;
      }
    }
  }

/* DELETE
j = *from;
for (i=0; i<ncon; i++)
  printf("[%5d %5d %.4f %.4f] ", i, PQueueGetSize(&queues[i][j]), npwgts[j*ncon+i], maxwgt[j*ncon+i]);
printf("***[%5d %5d]\n", *cnum, *from);
*/

  /* If the desired queue is empty, select a node from that side anyway */
  if (*from != -1 && PQueueGetSize(&queues[*cnum][*from]) == 0) {
    for (i=0; i<ncon; i++) {
      if (PQueueGetSize(&queues[i][*from]) > 0) {
        maxdiff = (npwgts[(*from)*ncon+i] - maxwgt[(*from)*ncon+i]);
        *cnum = i;
        break;
      }
    }

    for (i++; i<ncon; i++) {
      diff = npwgts[(*from)*ncon+i] - maxwgt[(*from)*ncon+i];
      if (diff > maxdiff && PQueueGetSize(&queues[i][*from]) > 0) {
        maxdiff = diff;
        *cnum = i;
      }
    }
  }

  /* If the constraints ar OK, select a high gain vertex */
  if (*from == -1) {
    maxgain = -100000;
    for (j=0; j<2; j++) {
      for (i=0; i<ncon; i++) {
        if (PQueueGetSize(&queues[i][j]) > 0 && PQueueGetKey(&queues[i][j]) > maxgain) {
          maxgain = PQueueGetKey(&queues[i][0]);
          *from = j;
          *cnum = i;
        }
      }
    }

    /* printf("(%2d %2d) %3d\n", *from, *cnum, maxgain); */
  }
}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mbalance.c
 *
 * This file contains code that is used to forcefully balance either
 * bisections or k-sections
 *
 * Started 7/29/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point of the bisection balancing algorithms.
**************************************************************************/
void MocBalance2Way(CtrlType *ctrl, GraphType *graph, float *tpwgts, float lbfactor)
{

  if (Compute2WayHLoadImbalance(graph->ncon, graph->npwgts, tpwgts) < lbfactor)
    return;

  MocGeneral2WayBalance(ctrl, graph, tpwgts, lbfactor);

}


/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void MocGeneral2WayBalance(CtrlType *ctrl, GraphType *graph, float *tpwgts, float lbfactor)
{
  int i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, limit, tmp, cnum;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *moved, *swaps, *perm, *qnum;
  float *nvwgt, *npwgts, mindiff[MAXNCON], origbal, minbal, newbal;
  PQueueType parts[MAXNCON][2];
  int higain, oldgain, mincut, newcut, mincutorder;
  int qsizes[MAXNCON][2];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  limit = amin(amax(0.01*nvtxs, 15), 100);

  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
    qsizes[i][0] = qsizes[i][1] = 0;
  }

  for (i=0; i<nvtxs; i++) {
    qnum[i] = samax(ncon, nvwgt+i*ncon);
    qsizes[qnum[i]][where[i]]++;
  }

/*
  printf("Weight Distribution:    \t");
  for (i=0; i<ncon; i++)
    printf(" [%d %d]", qsizes[i][0], qsizes[i][1]);
  printf("\n");
*/

  for (from=0; from<2; from++) {
    for (j=0; j<ncon; j++) {
      if (qsizes[j][from] == 0) {
        for (i=0; i<nvtxs; i++) {
          if (where[i] != from)
            continue;

          k = samax2(ncon, nvwgt+i*ncon);
          if (k == j && qsizes[qnum[i]][from] > qsizes[j][from] && nvwgt[i*ncon+qnum[i]] < 1.3*nvwgt[i*ncon+j]) {
            qsizes[qnum[i]][from]--;
            qsizes[j][from]++;
            qnum[i] = j;
          }
        }
      }
    }
  }

/*
  printf("Weight Distribution (after):\t ");
  for (i=0; i<ncon; i++)
    printf(" [%d %d]", qsizes[i][0], qsizes[i][1]);
  printf("\n");
*/



  for (i=0; i<ncon; i++)
    mindiff[i] = fabs(tpwgts[0]-npwgts[i]);
  minbal = origbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);
  newcut = mincut = graph->mincut;
  mincutorder = -1;

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Parts: [");
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("] T[%.3f %.3f], Nv-Nb[%5d, %5d]. ICut: %6d, LB: %.3f [B]\n", tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut, origbal);
  }

  idxset(nvtxs, -1, moved);

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  /* Insert all nodes in the priority queues */
  nbnd = graph->nbnd;
  RandomPermute(nvtxs, perm, 1);
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    PQueueInsert(&parts[qnum[i]][where[i]], i, ed[i]-id[i]);
  }

  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if (minbal < lbfactor)
      break;

    SelectQueue(ncon, npwgts, tpwgts, &from, &cnum, parts);
    to = (from+1)%2;

    if (from == -1 || (higain = PQueueGetMax(&parts[cnum][from])) == -1)
      break;

    saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
    newcut -= (ed[higain]-id[higain]);
    newbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

    if (newbal < minbal || (newbal == minbal &&
        (newcut < mincut || (newcut == mincut && BetterBalance(ncon, npwgts, tpwgts, mindiff))))) {
      mincut = newcut;
      minbal = newbal;
      mincutorder = nswaps;
      for (i=0; i<ncon; i++)
        mindiff[i] = fabs(tpwgts[0]-npwgts[i]);
    }
    else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
      newcut += (ed[higain]-id[higain]);
      saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
      saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      break;
    }

    where[higain] = to;
    moved[higain] = nswaps;
    swaps[nswaps] = higain;

    if (ctrl->dbglvl&DBG_MOVEINFO) {
      printf("Moved %6d from %d(%d). Gain: %5d, Cut: %5d, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], newcut);
      for (l=0; l<ncon; l++)
        printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
      printf(", %.3f LB: %.3f\n", minbal, newbal);
    }


    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (moved[k] == -1)
        PQueueUpdate(&parts[qnum[k]][where[k]], k, oldgain, ed[k]-id[k]);

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }



  /****************************************************************
  * Roll back computations
  *****************************************************************/
  for (nswaps--; nswaps>mincutorder; nswaps--) {
    higain = swaps[nswaps];

    to = where[higain] = (where[higain]+1)%2;
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    else if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      if (bndptr[k] != -1 && ed[k] == 0)
        BNDDelete(nbnd, bndind, bndptr, k);
      if (bndptr[k] == -1 && ed[k] > 0)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("\tMincut: %6d at %5d, NBND: %6d, NPwgts: [", mincut, mincutorder, nbnd);
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("], LB: %.3f\n", Compute2WayHLoadImbalance(ncon, npwgts, tpwgts));
  }

  graph->mincut = mincut;
  graph->nbnd = nbnd;


  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}

/*
 * mcoarsen.c
 *
 * This file contains the driving routines for the coarsening process
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *MCCoarsen2Way(CtrlType *ctrl, GraphType *graph)
{
  int i, clevel;
  GraphType *cgraph;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->CoarsenTmr));

  cgraph = graph;

  clevel = 0;
  do {
    if (ctrl->dbglvl&DBG_COARSEN) {
      printf("%6d %7d %10d [%d] [%6.4f", cgraph->nvtxs, cgraph->nedges,
              idxsum(cgraph->nvtxs, cgraph->adjwgtsum), ctrl->CoarsenTo, ctrl->nmaxvwgt);
      for (i=0; i<graph->ncon; i++)
        printf(" %5.3f", ssum_strd(cgraph->nvtxs, cgraph->nvwgt+i, cgraph->ncon));
      printf("]\n");
    }

    switch (ctrl->CType) {
      case MATCH_RM:
        MCMatch_RM(ctrl, cgraph);
        break;
      case MATCH_HEM:
        if (clevel < 1)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_HEM(ctrl, cgraph);
        break;
      case MATCH_SHEM:
        if (clevel < 1)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SHEM(ctrl, cgraph);
        break;
      case MATCH_SHEMKWAY:
        MCMatch_SHEM(ctrl, cgraph);
        break;
      case MATCH_SHEBM_ONENORM:
        MCMatch_SHEBM(ctrl, cgraph, 1);
        break;
      case MATCH_SHEBM_INFNORM:
        MCMatch_SHEBM(ctrl, cgraph, -1);
        break;
      case MATCH_SBHEM_ONENORM:
        MCMatch_SBHEM(ctrl, cgraph, 1);
        break;
      case MATCH_SBHEM_INFNORM:
        MCMatch_SBHEM(ctrl, cgraph, -1);
        break;
      default:
        errexit("Unknown CType: %d\n", ctrl->CType);
    }

    cgraph = cgraph->coarser;
    clevel++;

  } while (cgraph->nvtxs > ctrl->CoarsenTo && cgraph->nvtxs < COARSEN_FRACTION2*cgraph->finer->nvtxs && cgraph->nedges > cgraph->nvtxs/2);

  if (ctrl->dbglvl&DBG_COARSEN) {
    printf("%6d %7d %10d [%d] [%6.4f", cgraph->nvtxs, cgraph->nedges,
            idxsum(cgraph->nvtxs, cgraph->adjwgtsum), ctrl->CoarsenTo, ctrl->nmaxvwgt);
    for (i=0; i<graph->ncon; i++)
      printf(" %5.3f", ssum_strd(cgraph->nvtxs, cgraph->nvwgt+i, cgraph->ncon));
    printf("]\n");
  }


  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->CoarsenTmr));

  return cgraph;
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 2/24/96
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void AllocateWorkSpace(CtrlType *ctrl, GraphType *graph, int nparts)
{
  ctrl->wspace.pmat = NULL;

  if (ctrl->optype == OP_KMETIS) {
    ctrl->wspace.edegrees = (EDegreeType *)GKmalloc(graph->nedges*sizeof(EDegreeType), "AllocateWorkSpace: edegrees");
    ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;

    ctrl->wspace.pmat = idxmalloc(nparts*nparts, "AllocateWorkSpace: pmat");

    /* Memory requirements for different phases
          Coarsening
                    Matching: 4*nvtxs vectors
                 Contraction: 2*nvtxs vectors (from the above 4), 1*nparts, 1*Nedges
            Total = MAX(4*nvtxs, 2*nvtxs+nparts+nedges)

          Refinement
                Random Refinement/Balance: 5*nparts + 1*nvtxs + 2*nedges
                Greedy Refinement/Balance: 5*nparts + 2*nvtxs + 2*nedges + 1*PQueue(==Nvtxs)
            Total = 5*nparts + 3*nvtxs + 2*nedges

         Total = 5*nparts + 3*nvtxs + 2*nedges
    */
    ctrl->wspace.maxcore = 3*(graph->nvtxs+1) +                 /* Match/Refinement vectors */
                           5*(nparts+1) +                       /* Partition weights etc */
                           graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* Greedy k-way balance/refine */
                           20  /* padding for 64 bit machines */
                           ;
  }
  else if (ctrl->optype == OP_KVMETIS) {
    ctrl->wspace.edegrees = NULL;
    ctrl->wspace.vedegrees = (VEDegreeType *)GKmalloc(graph->nedges*sizeof(VEDegreeType), "AllocateWorkSpace: vedegrees");
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.vedegrees;

    ctrl->wspace.pmat = idxmalloc(nparts*nparts, "AllocateWorkSpace: pmat");

    /* Memory requirements for different phases are identical to KMETIS */
    ctrl->wspace.maxcore = 3*(graph->nvtxs+1) +                 /* Match/Refinement vectors */
                           3*(nparts+1) +                       /* Partition weights etc */
                           graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* Greedy k-way balance/refine */
                           20  /* padding for 64 bit machines */
                           ;
  }
  else {
    ctrl->wspace.edegrees = (EDegreeType *)idxmalloc(graph->nedges, "AllocateWorkSpace: edegrees");
    ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;

    ctrl->wspace.maxcore = 5*(graph->nvtxs+1) +                 /* Refinement vectors */
                           4*(nparts+1) +                       /* Partition weights etc */
                           2*graph->ncon*graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* 2-way refinement */
                           2*graph->ncon*(NEG_GAINSPAN+PLUS_GAINSPAN+1)*(sizeof(ListNodeType *)/sizeof(idxtype)) + /* 2-way refinement */
                           20  /* padding for 64 bit machines */
                           ;
  }

  ctrl->wspace.maxcore += HTLENGTH;
  ctrl->wspace.core = idxmalloc(ctrl->wspace.maxcore, "AllocateWorkSpace: maxcore");
  ctrl->wspace.ccore = 0;
}


/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void FreeWorkSpace(CtrlType *ctrl, GraphType *graph)
{
  GKfree((void **) &ctrl->wspace.edegrees, (void **) &ctrl->wspace.vedegrees,
         (void **) &ctrl->wspace.core,     (void **) &ctrl->wspace.pmat, LTERM);
}

/*************************************************************************
* This function returns how may words are left in the workspace
**************************************************************************/
int WspaceAvail(CtrlType *ctrl)
{
  return ctrl->wspace.maxcore - ctrl->wspace.ccore;
}


/*************************************************************************
* This function allocate space from the core
**************************************************************************/
idxtype *idxwspacemalloc(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return ctrl->wspace.core + ctrl->wspace.ccore - n;
}

/*************************************************************************
* This function frees space from the core
**************************************************************************/
void idxwspacefree(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}


/*************************************************************************
* This function allocate space from the core
**************************************************************************/
float *fwspacemalloc(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return (float *) (ctrl->wspace.core + ctrl->wspace.ccore - n);
}

/*************************************************************************
* This function frees space from the core
**************************************************************************/
void fwspacefree(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}



/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
GraphType *CreateGraph(void)
{
  GraphType *graph;

  graph = (GraphType *)GKmalloc(sizeof(GraphType), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
void InitGraph(GraphType *graph)
{
  graph->gdata = graph->rdata = NULL;

  graph->nvtxs = graph->nedges = -1;
  graph->mincut = graph->minvol = -1;

  graph->xadj = graph->vwgt = graph->adjncy = graph->adjwgt = NULL;
  graph->adjwgtsum = NULL;
  graph->label = NULL;
  graph->cmap = NULL;

  graph->where = graph->pwgts = NULL;
  graph->id = graph->ed = NULL;
  graph->bndptr = graph->bndind = NULL;
  graph->rinfo = NULL;
  graph->vrinfo = NULL;
  graph->nrinfo = NULL;

  graph->ncon = -1;
  graph->nvwgt = NULL;
  graph->npwgts = NULL;

  graph->vsize = NULL;

  graph->coarser = graph->finer = NULL;

}

/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeGraph(GraphType *graph)
{

  GKfree((void **) &graph->gdata, (void **) &graph->nvwgt,
         (void **) &graph->rdata, (void **) &graph->npwgts, LTERM);
  free(graph);
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh.c
 *
 * This file contains routines for converting 3D and 4D finite element
 * meshes into dual or nodal graphs
 *
 * Started 8/18/97
 * George
 *
 * $Id$
 *
 */



/*****************************************************************************
* This function creates a graph corresponding to the dual of a finite element
* mesh. At this point the supported elements are triangles, tetrahedrons, and
* bricks.
******************************************************************************/
void METIS_MeshToDual(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,
                      idxtype *dxadj, idxtype *dadjncy)
{
  int esizes[] = {-1, 3, 4, 8, 4};

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esizes[*etype], elmnts);

  GENDUALMETIS(*ne, *nn, *etype, elmnts, dxadj, dadjncy);

  if (*numflag == 1)
    ChangeMesh2FNumbering((*ne)*esizes[*etype], elmnts, *ne, dxadj, dadjncy);
}


/*****************************************************************************
* This function creates a graph corresponding to the finite element mesh.
* At this point the supported elements are triangles, tetrahedrons.
******************************************************************************/
void METIS_MeshToNodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,
                       idxtype *dxadj, idxtype *dadjncy)
{
  int esizes[] = {-1, 3, 4, 8, 4};

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esizes[*etype], elmnts);

  switch (*etype) {
    case 1:
      TRINODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 2:
      TETNODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 3:
      HEXNODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 4:
      QUADNODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering((*ne)*esizes[*etype], elmnts, *nn, dxadj, dadjncy);
}



/*****************************************************************************
* This function creates the dual of a finite element mesh
******************************************************************************/
void GENDUALMETIS(int nelmnts, int nvtxs, int etype, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   int i, j, jj, k, kk, kkk, l, m, n, mask;
   idxtype *nptr, *nind;
   idxtype *mark, ind[200], wgt[200];
   int esize, esizes[] = {-1, 3, 4, 8, 4},
       mgcnum, mgcnums[] = {-1, 2, 3, 4, 2};

   mask = (1<<11)-1;
   mark = idxsmalloc(mask+1, -1, "GENDUALMETIS: mark");

   /* Get the element size and magic number for the particular element */
   esize = esizes[etype];
   mgcnum = mgcnums[etype];

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "GENDUALMETIS: nptr");
   for (j=esize*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "GENDUALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<esize; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;

   for (i=0; i<nelmnts; i++)
     dxadj[i] = esize*i;

   for (i=0; i<nelmnts; i++) {
     for (m=j=0; j<esize; j++) {
       n = elmnts[esize*i+j];
       for (k=nptr[n+1]-1; k>=nptr[n]; k--) {
         if ((kk = nind[k]) <= i)
           break;

         kkk = kk&mask;
         if ((l = mark[kkk]) == -1) {
           ind[m] = kk;
           wgt[m] = 1;
           mark[kkk] = m++;
         }
         else if (ind[l] == kk) {
           wgt[l]++;
         }
         else {
           for (jj=0; jj<m; jj++) {
             if (ind[jj] == kk) {
               wgt[jj]++;
               break;
             }
           }
           if (jj == m) {
             ind[m] = kk;
             wgt[m++] = 1;
           }
         }
       }
     }
     for (j=0; j<m; j++) {
       if (wgt[j] == mgcnum) {
         k = ind[j];
         dadjncy[dxadj[i]++] = k;
         dadjncy[dxadj[k]++] = i;
       }
       mark[ind[j]&mask] = -1;
     }
   }

   /* Go and consolidate the dxadj and dadjncy */
   for (j=i=0; i<nelmnts; i++) {
     for (k=esize*i; k<dxadj[i]; k++, j++)
       dadjncy[j] = dadjncy[k];
     dxadj[i] = j;
   }
   for (i=nelmnts; i>0; i--)
     dxadj[i] = dxadj[i-1];
   dxadj[0] = 0;

   free(mark);
   free(nptr);
   free(nind);

}




/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void TRINODALMETIS(int nelmnts, int nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   int i, j, jj, k, kk, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "TRINODALMETIS: nptr");
   for (j=3*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "TRINODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<3; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "TRINODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       for (jj=3*nind[j], k=0; k<3; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
     }
     dxadj[i+1] = nedges;
   }

   free(mark);
   free(nptr);
   free(nind);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void TETNODALMETIS(int nelmnts, int nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   int i, j, jj, k, kk, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "TETNODALMETIS: nptr");
   for (j=4*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "TETNODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<4; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "TETNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       for (jj=4*nind[j], k=0; k<4; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
     }
     dxadj[i+1] = nedges;
   }

   free(mark);
   free(nptr);
   free(nind);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void HEXNODALMETIS(int nelmnts, int nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   int i, j, jj, k, kk, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;
   int table[8][3] = {{1, 3, 4},
                      {0, 2, 5},
                      {1, 3, 6},
                      {0, 2, 7},
                      {0, 5, 7},
                      {1, 4, 6},
                      {2, 5, 7},
                      {3, 4, 6}};

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "HEXNODALMETIS: nptr");
   for (j=8*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "HEXNODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<8; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "HEXNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       jj=8*nind[j];
       for (k=0; k<8; k++) {
         if (elmnts[jj+k] == i)
           break;
       }
       ASSERT(k != 8);

       /* You found the index, now go and put the 3 neighbors */
       kk = elmnts[jj+table[k][0]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+table[k][1]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+table[k][2]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
     }
     dxadj[i+1] = nedges;
   }

   free(mark);
   free(nptr);
   free(nind);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void QUADNODALMETIS(int nelmnts, int nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   int i, j, jj, k, kk, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;
   int table[4][2] = {{1, 3},
                      {0, 2},
                      {1, 3},
                      {0, 2}};

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "QUADNODALMETIS: nptr");
   for (j=4*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "QUADNODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<4; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "QUADNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       jj=4*nind[j];
       for (k=0; k<4; k++) {
         if (elmnts[jj+k] == i)
           break;
       }
       ASSERT(k != 4);

       /* You found the index, now go and put the 2 neighbors */
       kk = elmnts[jj+table[k][0]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+table[k][1]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
     }
     dxadj[i+1] = nedges;
   }

   free(mark);
   free(nptr);
   free(nind);

}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * meshpart.c
 *
 * This file contains routines for partitioning finite element meshes.
 *
 * Started 9/29/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function partitions a finite element mesh by partitioning its nodal
* graph using KMETIS and then assigning elements in a load balanced fashion.
**************************************************************************/
void METIS_PartMeshNodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,
                         int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  int i, j, k, me;
  idxtype *xadj, *adjncy, *pwgts;
  int options[10], pnumflag=0, wgtflag=0;
  int nnbrs, nbrind[200], nbrwgt[200], maxpwgt;
  int esize, esizes[] = {-1, 3, 4, 8, 4};

  esize = esizes[*etype];

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esize, elmnts);

  xadj = idxmalloc(*nn+1, "METIS_MESHPARTNODAL: xadj");
  adjncy = idxmalloc(20*(*nn), "METIS_MESHPARTNODAL: adjncy");

  METIS_MeshToNodal(ne, nn, elmnts, etype, &pnumflag, xadj, adjncy);

  adjncy = realloc(adjncy, xadj[*nn]*sizeof(idxtype));

  options[0] = 0;
  METIS_PartGraphKway(nn, xadj, adjncy, NULL, NULL, &wgtflag, &pnumflag, nparts, options, edgecut, npart);

  /* OK, now compute an element partition based on the nodal partition npart */
  idxset(*ne, -1, epart);
  pwgts = idxsmalloc(*nparts, 0, "METIS_MESHPARTNODAL: pwgts");
  for (i=0; i<*ne; i++) {
    me = npart[elmnts[i*esize]];
    for (j=1; j<esize; j++) {
      if (npart[elmnts[i*esize+j]] != me)
        break;
    }
    if (j == esize) {
      epart[i] = me;
      pwgts[me]++;
    }
  }

  maxpwgt = 1.03*(*ne)/(*nparts);
  for (i=0; i<*ne; i++) {
    if (epart[i] == -1) { /* Assign the boundary element */
      nnbrs = 0;
      for (j=0; j<esize; j++) {
        me = npart[elmnts[i*esize+j]];
        for (k=0; k<nnbrs; k++) {
          if (nbrind[k] == me) {
            nbrwgt[k]++;
            break;
          }
        }
        if (k == nnbrs) {
          nbrind[nnbrs] = me;
          nbrwgt[nnbrs++] = 1;
        }
      }
      /* Try to assign it first to the domain with most things in common */
      j = iamax(nnbrs, nbrwgt);
      if (pwgts[nbrind[j]] < maxpwgt) {
        epart[i] = nbrind[j];
      }
      else {
        /* If that fails, assign it to a light domain */
        for (j=0; j<nnbrs; j++) {
          if (pwgts[nbrind[j]] < maxpwgt) {
            epart[i] = nbrind[j];
            break;
          }
        }
        if (j == nnbrs)
          epart[i] = nbrind[iamax(nnbrs, nbrwgt)];
      }
      pwgts[epart[i]]++;
    }
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering2((*ne)*esize, elmnts, *ne, *nn, epart, npart);

  GKfree((void **) &xadj, (void **) &adjncy, (void **) &pwgts, LTERM);

}


/*************************************************************************
* This function partitions a finite element mesh by partitioning its dual
* graph using KMETIS and then assigning nodes in a load balanced fashion.
**************************************************************************/
void METIS_PartMeshDual(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,
                        int *nparts, int *edgecut, idxtype *epart, idxtype *npart)
{
  int i, j, k, me;
  idxtype *xadj, *adjncy, *pwgts, *nptr, *nind;
  int options[10], pnumflag=0, wgtflag=0;
  int nnbrs, nbrind[200], nbrwgt[200], maxpwgt;
  int esize, esizes[] = {-1, 3, 4, 8, 4};

  esize = esizes[*etype];

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esize, elmnts);

  xadj = idxmalloc(*ne+1, "METIS_MESHPARTNODAL: xadj");
  adjncy = idxmalloc(esize*(*ne), "METIS_MESHPARTNODAL: adjncy");

  METIS_MeshToDual(ne, nn, elmnts, etype, &pnumflag, xadj, adjncy);

  options[0] = 0;
  METIS_PartGraphKway(ne, xadj, adjncy, NULL, NULL, &wgtflag, &pnumflag, nparts, options, edgecut, epart);

  /* Construct the node-element list */
  nptr = idxsmalloc(*nn+1, 0, "METIS_MESHPARTDUAL: nptr");
  for (j=esize*(*ne), i=0; i<j; i++)
    nptr[elmnts[i]]++;
  MAKECSR(i, *nn, nptr);

  nind = idxmalloc(nptr[*nn], "METIS_MESHPARTDUAL: nind");
  for (k=i=0; i<(*ne); i++) {
    for (j=0; j<esize; j++, k++)
      nind[nptr[elmnts[k]]++] = i;
  }
  for (i=(*nn); i>0; i--)
    nptr[i] = nptr[i-1];
  nptr[0] = 0;


  /* OK, now compute a nodal partition based on the element partition npart */
  idxset(*nn, -1, npart);
  pwgts = idxsmalloc(*nparts, 0, "METIS_MESHPARTDUAL: pwgts");
  for (i=0; i<*nn; i++) {
    me = epart[nind[nptr[i]]];
    for (j=nptr[i]+1; j<nptr[i+1]; j++) {
      if (epart[nind[j]] != me)
        break;
    }
    if (j == nptr[i+1]) {
      npart[i] = me;
      pwgts[me]++;
    }
  }

  maxpwgt = 1.03*(*nn)/(*nparts);
  for (i=0; i<*nn; i++) {
    if (npart[i] == -1) { /* Assign the boundary element */
      nnbrs = 0;
      for (j=nptr[i]; j<nptr[i+1]; j++) {
        me = epart[nind[j]];
        for (k=0; k<nnbrs; k++) {
          if (nbrind[k] == me) {
            nbrwgt[k]++;
            break;
          }
        }
        if (k == nnbrs) {
          nbrind[nnbrs] = me;
          nbrwgt[nnbrs++] = 1;
        }
      }
      /* Try to assign it first to the domain with most things in common */
      j = iamax(nnbrs, nbrwgt);
      if (pwgts[nbrind[j]] < maxpwgt) {
        npart[i] = nbrind[j];
      }
      else {
        /* If that fails, assign it to a light domain */
        npart[i] = nbrind[0];
        for (j=0; j<nnbrs; j++) {
          if (pwgts[nbrind[j]] < maxpwgt) {
            npart[i] = nbrind[j];
            break;
          }
        }
      }
      pwgts[npart[i]]++;
    }
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering2((*ne)*esize, elmnts, *ne, *nn, epart, npart);

  GKfree((void **) &xadj, (void **) &adjncy, (void **) &pwgts,
         (void **) &nptr, (void **) &nind, LTERM);

}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mfm2.c
 *
 * This file contains code that implements the edge-based FM refinement
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void MocFM_2WayEdgeRefine2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *orgubvec,
       int npasses)
{
  int i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, pass, limit, tmp, cnum;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *moved, *swaps, *perm, *qnum;
  float *nvwgt, *npwgts, origdiff[MAXNCON], origbal[MAXNCON], minbal[MAXNCON];
  PQueueType parts[MAXNCON][2];
  int higain, oldgain, mincut, initcut, newcut, mincutorder;
  float *maxwgt, *minwgt, ubvec[MAXNCON], tvec[MAXNCON];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  limit = amin(amax(0.01*nvtxs, 15), 100);

  Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, origbal);
  for (i=0; i<ncon; i++) {
    origdiff[i] = fabs(tpwgts[0]-npwgts[i]);
    ubvec[i] = amax(origbal[i], orgubvec[i]);
  }

  /* Setup the weight intervals of the two subdomains */
  minwgt = fwspacemalloc(ctrl, 2*ncon);
  maxwgt = fwspacemalloc(ctrl, 2*ncon);

  for (i=0; i<2; i++) {
    for (j=0; j<ncon; j++) {
      maxwgt[i*ncon+j] = tpwgts[i]*ubvec[j];
      minwgt[i*ncon+j] = tpwgts[i]*(1.0/ubvec[j]);
    }
  }

  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
  }
  for (i=0; i<nvtxs; i++)
    qnum[i] = samax(ncon, nvwgt+i*ncon);


  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Parts: [");
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("] T[%.3f %.3f], Nv-Nb[%5d, %5d]. ICut: %6d, LB: ", tpwgts[0], tpwgts[1],
            graph->nvtxs, graph->nbnd, graph->mincut);
    for (i=0; i<ncon; i++)
      printf("%.3f ", origbal[i]);
    printf("\n");
  }

  idxset(nvtxs, -1, moved);
  for (pass=0; pass<npasses; pass++) { /* Do a number of passes */
    for (i=0; i<ncon; i++) {
      PQueueReset(&parts[i][0]);
      PQueueReset(&parts[i][1]);
    }

    mincutorder = -1;
    newcut = mincut = initcut = graph->mincut;
    Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, minbal);

    ASSERT(ComputeCut(graph, where) == graph->mincut);
    ASSERT(CheckBnd(graph));

    /* Insert boundary nodes in the priority queues */
    nbnd = graph->nbnd;
    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(ed[i] > 0 || id[i] == 0);
      ASSERT(bndptr[i] != -1);
      PQueueInsert(&parts[qnum[i]][where[i]], i, ed[i]-id[i]);
    }

    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      SelectQueue2(ncon, npwgts, tpwgts, &from, &cnum, parts, maxwgt);
      to = (from+1)%2;

      if (from == -1 || (higain = PQueueGetMax(&parts[cnum][from])) == -1)
        break;
      ASSERT(bndptr[higain] != -1);

      newcut -= (ed[higain]-id[higain]);
      saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

      Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, tvec);
      if ((newcut < mincut && AreAllBelow(ncon, tvec, ubvec)) ||
          (newcut == mincut && IsBetter2wayBalance(ncon, tvec, minbal, ubvec))) {
        mincut = newcut;
        for (i=0; i<ncon; i++)
          minbal[i] = tvec[i];
        mincutorder = nswaps;
      }
      else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
        newcut += (ed[higain]-id[higain]);
        saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
        saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
        break;
      }

      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;

      if (ctrl->dbglvl&DBG_MOVEINFO) {
        printf("Moved %6d from %d(%d). Gain: %5d, Cut: %5d, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], newcut);
        for (l=0; l<ncon; l++)
          printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);

        printf(", LB: ");
        for (i=0; i<ncon; i++)
          printf("%.3f ", tvec[i]);
        if (mincutorder == nswaps)
          printf(" *\n");
        else
          printf("\n");
      }


      /**************************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***************************************************************/
      SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);

      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        oldgain = ed[k]-id[k];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        /* Update its boundary information and queue position */
        if (bndptr[k] != -1) { /* If k was a boundary vertex */
          if (ed[k] == 0) { /* Not a boundary vertex any more */
            BNDDelete(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)  /* Remove it if in the queues */
              PQueueDelete(&parts[qnum[k]][where[k]], k, oldgain);
          }
          else { /* If it has not been moved, update its position in the queue */
            if (moved[k] == -1)
              PQueueUpdate(&parts[qnum[k]][where[k]], k, oldgain, ed[k]-id[k]);
          }
        }
        else {
          if (ed[k] > 0) {  /* It will now become a boundary vertex */
            BNDInsert(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)
              PQueueInsert(&parts[qnum[k]][where[k]], k, ed[k]-id[k]);
          }
        }
      }

    }


    /****************************************************************
    * Roll back computations
    *****************************************************************/
    for (i=0; i<nswaps; i++)
      moved[swaps[i]] = -1;  /* reset moved array */
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      to = where[higain] = (where[higain]+1)%2;
      SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);
      else if (ed[higain] > 0 && bndptr[higain] == -1)
        BNDInsert(nbnd, bndind,  bndptr, higain);

      saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        if (bndptr[k] != -1 && ed[k] == 0)
          BNDDelete(nbnd, bndind, bndptr, k);
        if (bndptr[k] == -1 && ed[k] > 0)
          BNDInsert(nbnd, bndind, bndptr, k);
      }
    }

    if (ctrl->dbglvl&DBG_REFINE) {
      printf("\tMincut: %6d at %5d, NBND: %6d, NPwgts: [", mincut, mincutorder, nbnd);
      for (l=0; l<ncon; l++)
        printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
      printf("], LB: ");
      Compute2WayHLoadImbalanceVec(ncon, npwgts, tpwgts, tvec);
      for (i=0; i<ncon; i++)
        printf("%.3f ", tvec[i]);
      printf("\n");
    }

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut == initcut)
      break;
  }

  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  fwspacefree(ctrl, 2*ncon);
  fwspacefree(ctrl, 2*ncon);

}


/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
void SelectQueue2(int ncon, float *npwgts, float *tpwgts, int *from, int *cnum,
       PQueueType queues[MAXNCON][2], float *maxwgt)
{
  int i, j, maxgain=0;
  float diff, max, maxdiff=0.0;

  *from = -1;
  *cnum = -1;

  /* First determine the side and the queue, irrespective of the presence of nodes */
  for (j=0; j<2; j++) {
    for (i=0; i<ncon; i++) {
      diff = npwgts[j*ncon+i]-maxwgt[j*ncon+i];
      if (diff >= maxdiff) {
        maxdiff = diff;
        *from = j;
        *cnum = i;
      }
    }
  }

  if (*from != -1 && PQueueGetSize(&queues[*cnum][*from]) == 0) {
    /* The desired queue is empty, select a node from that side anyway */
    for (i=0; i<ncon; i++) {
      if (PQueueGetSize(&queues[i][*from]) > 0) {
        max = (npwgts[(*from)*ncon+i] - maxwgt[(*from)*ncon+i]);
        *cnum = i;
        break;
      }
    }

    for (i++; i<ncon; i++) {
      diff = npwgts[(*from)*ncon+i] - maxwgt[(*from)*ncon+i];
      if (diff > max && PQueueGetSize(&queues[i][*from]) > 0) {
        max = diff;
        *cnum = i;
      }
    }
  }

  /* Check to see if you can focus on the cut */
  if (maxdiff <= 0.0 || *from == -1) {
    maxgain = -100000;

    for (j=0; j<2; j++) {
      for (i=0; i<ncon; i++) {
        if (PQueueGetSize(&queues[i][j]) > 0 && PQueueGetKey(&queues[i][j]) > maxgain) {
          maxgain = PQueueGetKey(&queues[i][j]);
          *from = j;
          *cnum = i;
        }
      }
    }

    /* printf("(%2d %2d) %3d\n", *from, *cnum, maxgain); */
  }
}


/*************************************************************************
* This function checks if the newbal is better than oldbal given the
* ubvector ubvec
**************************************************************************/
int IsBetter2wayBalance(int ncon, float *newbal, float *oldbal, float *ubvec)
{
  int i;
  float max1=0.0, max2=0.0, sum1=0.0, sum2=0.0, tmp;

  for (i=0; i<ncon; i++) {
    tmp = (newbal[i]-1)/(ubvec[i]-1);
    max1 = (max1 < tmp ? tmp : max1);
    sum1 += tmp;

    tmp = (oldbal[i]-1)/(ubvec[i]-1);
    max2 = (max2 < tmp ? tmp : max2);
    sum2 += tmp;
  }

  if (max1 < max2)
    return 1;
  else if (max1 > max2)
    return 0;
  else
    return sum1 <= sum2;
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mfm.c
 *
 * This file contains code that implements the edge-based FM refinement
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void MocFM_2WayEdgeRefine(CtrlType *ctrl, GraphType *graph, float *tpwgts, int npasses)
{
  int i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, pass, limit, tmp, cnum;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *moved, *swaps, *perm, *qnum;
  float *nvwgt, *npwgts, mindiff[MAXNCON], origbal, minbal, newbal;
  PQueueType parts[MAXNCON][2];
  int higain, oldgain, mincut, initcut, newcut, mincutorder;
  float rtpwgts[2];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  limit = amin(amax(0.01*nvtxs, 25), 150);

  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
  }
  for (i=0; i<nvtxs; i++)
    qnum[i] = samax(ncon, nvwgt+i*ncon);

  origbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

  rtpwgts[0] = origbal*tpwgts[0];
  rtpwgts[1] = origbal*tpwgts[1];


  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Parts: [");
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("] T[%.3f %.3f], Nv-Nb[%5d, %5d]. ICut: %6d, LB: %.3f\n", tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut, origbal);
  }

  idxset(nvtxs, -1, moved);
  for (pass=0; pass<npasses; pass++) { /* Do a number of passes */
    for (i=0; i<ncon; i++) {
      PQueueReset(&parts[i][0]);
      PQueueReset(&parts[i][1]);
    }

    mincutorder = -1;
    newcut = mincut = initcut = graph->mincut;
    for (i=0; i<ncon; i++)
      mindiff[i] = fabs(tpwgts[0]-npwgts[i]);
    minbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

    ASSERT(ComputeCut(graph, where) == graph->mincut);
    ASSERT(CheckBnd(graph));

    /* Insert boundary nodes in the priority queues */
    nbnd = graph->nbnd;
    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(ed[i] > 0 || id[i] == 0);
      ASSERT(bndptr[i] != -1);
      PQueueInsert(&parts[qnum[i]][where[i]], i, ed[i]-id[i]);
    }

    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      SelectQueue(ncon, npwgts, rtpwgts, &from, &cnum, parts);
      to = (from+1)%2;

      if (from == -1 || (higain = PQueueGetMax(&parts[cnum][from])) == -1)
        break;
      ASSERT(bndptr[higain] != -1);

      saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

      newcut -= (ed[higain]-id[higain]);
      newbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

      if ((newcut < mincut && newbal-origbal <= .00001) ||
          (newcut == mincut && (newbal < minbal ||
                                (newbal == minbal && BetterBalance(ncon, npwgts, tpwgts, mindiff))))) {
        mincut = newcut;
        minbal = newbal;
        mincutorder = nswaps;
        for (i=0; i<ncon; i++)
          mindiff[i] = fabs(tpwgts[0]-npwgts[i]);
      }
      else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
        newcut += (ed[higain]-id[higain]);
        saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
        saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
        break;
      }

      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;

      if (ctrl->dbglvl&DBG_MOVEINFO) {
        printf("Moved %6d from %d(%d). Gain: %5d, Cut: %5d, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], newcut);
        for (l=0; l<ncon; l++)
          printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
        printf(", %.3f LB: %.3f\n", minbal, newbal);
      }


      /**************************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***************************************************************/
      SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);

      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        oldgain = ed[k]-id[k];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        /* Update its boundary information and queue position */
        if (bndptr[k] != -1) { /* If k was a boundary vertex */
          if (ed[k] == 0) { /* Not a boundary vertex any more */
            BNDDelete(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)  /* Remove it if in the queues */
              PQueueDelete(&parts[qnum[k]][where[k]], k, oldgain);
          }
          else { /* If it has not been moved, update its position in the queue */
            if (moved[k] == -1)
              PQueueUpdate(&parts[qnum[k]][where[k]], k, oldgain, ed[k]-id[k]);
          }
        }
        else {
          if (ed[k] > 0) {  /* It will now become a boundary vertex */
            BNDInsert(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)
              PQueueInsert(&parts[qnum[k]][where[k]], k, ed[k]-id[k]);
          }
        }
      }

    }


    /****************************************************************
    * Roll back computations
    *****************************************************************/
    for (i=0; i<nswaps; i++)
      moved[swaps[i]] = -1;  /* reset moved array */
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      to = where[higain] = (where[higain]+1)%2;
      SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);
      else if (ed[higain] > 0 && bndptr[higain] == -1)
        BNDInsert(nbnd, bndind,  bndptr, higain);

      saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        if (bndptr[k] != -1 && ed[k] == 0)
          BNDDelete(nbnd, bndind, bndptr, k);
        if (bndptr[k] == -1 && ed[k] > 0)
          BNDInsert(nbnd, bndind, bndptr, k);
      }
    }

    if (ctrl->dbglvl&DBG_REFINE) {
      printf("\tMincut: %6d at %5d, NBND: %6d, NPwgts: [", mincut, mincutorder, nbnd);
      for (l=0; l<ncon; l++)
        printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
      printf("], LB: %.3f\n", Compute2WayHLoadImbalance(ncon, npwgts, tpwgts));
    }

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut == initcut)
      break;
  }

  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}


/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
void SelectQueue(int ncon, float *npwgts, float *tpwgts, int *from, int *cnum, PQueueType queues[MAXNCON][2])
{
  int i, part, maxgain=0;
  float max, maxdiff=0.0;

  *from = -1;
  *cnum = -1;

  /* First determine the side and the queue, irrespective of the presence of nodes */
  for (part=0; part<2; part++) {
    for (i=0; i<ncon; i++) {
      if (npwgts[part*ncon+i]-tpwgts[part] >= maxdiff) {
        maxdiff = npwgts[part*ncon+i]-tpwgts[part];
        *from = part;
        *cnum = i;
      }
    }
  }

  /* printf("Selected1 %d(%d) -> %d [%5f]\n", *from, *cnum, PQueueGetSize(&queues[*cnum][*from]), maxdiff); */

  if (*from != -1 && PQueueGetSize(&queues[*cnum][*from]) == 0) {
    /* The desired queue is empty, select a node from that side anyway */
    for (i=0; i<ncon; i++) {
      if (PQueueGetSize(&queues[i][*from]) > 0) {
        max = npwgts[(*from)*ncon + i];
        *cnum = i;
        break;
      }
    }

    for (i++; i<ncon; i++) {
      if (npwgts[(*from)*ncon + i] > max && PQueueGetSize(&queues[i][*from]) > 0) {
        max = npwgts[(*from)*ncon + i];
        *cnum = i;
      }
    }
  }

  /* Check to see if you can focus on the cut */
  if (maxdiff <= 0.0 || *from == -1) {
    maxgain = -100000;

    for (part=0; part<2; part++) {
      for (i=0; i<ncon; i++) {
        if (PQueueGetSize(&queues[i][part]) > 0 && PQueueGetKey(&queues[i][part]) > maxgain) {
          maxgain = PQueueGetKey(&queues[i][part]);
          *from = part;
          *cnum = i;
        }
      }
    }
  }

  /* printf("Selected2 %d(%d) -> %d\n", *from, *cnum, PQueueGetSize(&queues[*cnum][*from])); */
}





/*************************************************************************
* This function checks if the balance achieved is better than the diff
* For now, it uses a 2-norm measure
**************************************************************************/
int BetterBalance(int ncon, float *npwgts, float *tpwgts, float *diff)
{
  int i;
  float ndiff[MAXNCON];

  for (i=0; i<ncon; i++)
    ndiff[i] = fabs(tpwgts[0]-npwgts[i]);

  return snorm2(ncon, ndiff) < snorm2(ncon, diff);
}



/*************************************************************************
* This function computes the load imbalance over all the constrains
**************************************************************************/
float Compute2WayHLoadImbalance(int ncon, float *npwgts, float *tpwgts)
{
  int i;
  float max=0.0, temp;

  for (i=0; i<ncon; i++) {
    /* temp = amax(npwgts[i]/tpwgts[0], npwgts[ncon+i]/tpwgts[1]); */
    temp = fabs(tpwgts[0]-npwgts[i])/tpwgts[0];
    max = (max < temp ? temp : max);
  }
  return 1.0+max;
}


/*************************************************************************
* This function computes the load imbalance over all the constrains
* For now assume that we just want balanced partitionings
**************************************************************************/
void Compute2WayHLoadImbalanceVec(int ncon, float *npwgts, float *tpwgts, float *lbvec)
{
  int i;

  for (i=0; i<ncon; i++)
    lbvec[i] = 1.0 + fabs(tpwgts[0]-npwgts[i])/tpwgts[0];
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mincover.c
 *
 * This file implements the minimum cover algorithm
 *
 * Started 8/1/97
 * George
 *
 * $Id$
 */



/*************************************************************************
* Constants used by mincover algorithm
**************************************************************************/
#define INCOL 10
#define INROW 20
#define VC 1
#define SC 2
#define HC 3
#define VR 4
#define SR 5
#define HR 6


/*************************************************************************
* This function returns the min-cover of a bipartite graph.
* The algorithm used is due to Hopcroft and Karp as modified by Duff etal
* adj: the adjacency list of the bipartite graph
*       asize: the number of vertices in the first part of the bipartite graph
* bsize-asize: the number of vertices in the second part
*        0..(asize-1) > A vertices
*        asize..bsize > B vertices
*
* Returns:
*  cover : the actual cover (array)
*  csize : the size of the cover
**************************************************************************/
void MinCover(idxtype *xadj, idxtype *adjncy, int asize, int bsize, idxtype *cover, int *csize)
{
  int i, j;
  idxtype *mate, *queue, *flag, *level, *lst;
  int fptr, rptr, lstptr;
  int row, maxlevel, col;

  mate = idxsmalloc(bsize, -1, "MinCover: mate");
  flag = idxmalloc(bsize, "MinCover: flag");
  level = idxmalloc(bsize, "MinCover: level");
  queue = idxmalloc(bsize, "MinCover: queue");
  lst = idxmalloc(bsize, "MinCover: lst");

  /* Get a cheap matching */
  for (i=0; i<asize; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (mate[adjncy[j]] == -1) {
        mate[i] = adjncy[j];
        mate[adjncy[j]] = i;
        break;
      }
    }
  }

  /* Get into the main loop */
  while (1) {
    /* Initialization */
    fptr = rptr = 0;   /* Empty Queue */
    lstptr = 0;        /* Empty List */
    for (i=0; i<bsize; i++) {
      level[i] = -1;
      flag[i] = 0;
    }
    maxlevel = bsize;

    /* Insert free nodes into the queue */
    for (i=0; i<asize; i++)
      if (mate[i] == -1) {
        queue[rptr++] = i;
        level[i] = 0;
      }

    /* Perform the BFS */
    while (fptr != rptr) {
      row = queue[fptr++];
      if (level[row] < maxlevel) {
        flag[row] = 1;
        for (j=xadj[row]; j<xadj[row+1]; j++) {
          col = adjncy[j];
          if (!flag[col]) {  /* If this column has not been accessed yet */
            flag[col] = 1;
            if (mate[col] == -1) { /* Free column node was found */
              maxlevel = level[row];
              lst[lstptr++] = col;
            }
            else { /* This column node is matched */
              if (flag[mate[col]])
                printf("\nSomething wrong, flag[%d] is 1",mate[col]);
              queue[rptr++] = mate[col];
              level[mate[col]] = level[row] + 1;
            }
          }
        }
      }
    }

    if (lstptr == 0)
      break;   /* No free columns can be reached */

    /* Perform restricted DFS from the free column nodes */
    for (i=0; i<lstptr; i++)
      MinCover_Augment(xadj, adjncy, lst[i], mate, flag, level, maxlevel);
  }

  MinCover_Decompose(xadj, adjncy, asize, bsize, mate, cover, csize);

  GKfree((void **) &mate,  (void **) &flag, (void **) &level,
         (void **) &queue, (void **) &lst, LTERM);

}


/*************************************************************************
* This function perfoms a restricted DFS and augments matchings
**************************************************************************/
int MinCover_Augment(idxtype *xadj, idxtype *adjncy, int col, idxtype *mate, idxtype *flag, idxtype *level, int maxlevel)
{
  int i;
  int row = -1;
  int status;

  flag[col] = 2;
  for (i=xadj[col]; i<xadj[col+1]; i++) {
    row = adjncy[i];

    if (flag[row] == 1) { /* First time through this row node */
      if (level[row] == maxlevel) {  /* (col, row) is an edge of the G^T */
        flag[row] = 2;  /* Mark this node as being visited */
        if (maxlevel != 0)
          status = MinCover_Augment(xadj, adjncy, mate[row], mate, flag, level, maxlevel-1);
        else
          status = 1;

        if (status) {
          mate[col] = row;
          mate[row] = col;
          return 1;
        }
      }
    }
  }

  return 0;
}



/*************************************************************************
* This function performs a coarse decomposition and determines the
* min-cover.
* REF: Pothen ACMTrans. on Amth Software
**************************************************************************/
void MinCover_Decompose(idxtype *xadj, idxtype *adjncy, int asize, int bsize, idxtype *mate, idxtype *cover, int *csize)
{
  int i, k;
  idxtype *where;
  int card[10];

  where = idxmalloc(bsize, "MinCover_Decompose: where");
  for (i=0; i<10; i++)
    card[i] = 0;

  for (i=0; i<asize; i++)
    where[i] = SC;
  for (; i<bsize; i++)
    where[i] = SR;

  for (i=0; i<asize; i++)
    if (mate[i] == -1)
      MinCover_ColDFS(xadj, adjncy, i, mate, where, INCOL);
  for (; i<bsize; i++)
    if (mate[i] == -1)
      MinCover_RowDFS(xadj, adjncy, i, mate, where, INROW);

  for (i=0; i<bsize; i++)
    card[where[i]]++;

  k = 0;
  if (abs(card[VC]+card[SC]-card[HR]) < abs(card[VC]-card[SR]-card[HR])) {  /* S = VC+SC+HR */
    /* printf("%d %d ",vc+sc, hr); */
    for (i=0; i<bsize; i++)
      if (where[i] == VC || where[i] == SC || where[i] == HR)
        cover[k++] = i;
  }
  else {  /* S = VC+SR+HR */
    /* printf("%d %d ",vc, hr+sr); */
    for (i=0; i<bsize; i++)
      if (where[i] == VC || where[i] == SR || where[i] == HR)
        cover[k++] = i;
  }

  *csize = k;
  free(where);

}


/*************************************************************************
* This function perfoms a dfs starting from an unmatched col node
* forming alternate paths
**************************************************************************/
void MinCover_ColDFS(idxtype *xadj, idxtype *adjncy, int root, idxtype *mate, idxtype *where, int flag)
{
  int i;

  if (flag == INCOL) {
    if (where[root] == HC)
      return;
    where[root] = HC;
    for (i=xadj[root]; i<xadj[root+1]; i++)
      MinCover_ColDFS(xadj, adjncy, adjncy[i], mate, where, INROW);
  }
  else {
    if (where[root] == HR)
      return;
    where[root] = HR;
    if (mate[root] != -1)
      MinCover_ColDFS(xadj, adjncy, mate[root], mate, where, INCOL);
  }

}

/*************************************************************************
* This function perfoms a dfs starting from an unmatched col node
* forming alternate paths
**************************************************************************/
void MinCover_RowDFS(idxtype *xadj, idxtype *adjncy, int root, idxtype *mate, idxtype *where, int flag)
{
  int i;

  if (flag == INROW) {
    if (where[root] == VR)
      return;
    where[root] = VR;
    for (i=xadj[root]; i<xadj[root+1]; i++)
      MinCover_RowDFS(xadj, adjncy, adjncy[i], mate, where, INCOL);
  }
  else {
    if (where[root] == VC)
      return;
    where[root] = VC;
    if (mate[root] != -1)
      MinCover_RowDFS(xadj, adjncy, mate[root], mate, where, INROW);
  }

}



/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * minitpart2.c
 *
 * This file contains code that performs the initial partition of the
 * coarsest graph
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function computes the initial bisection of the coarsest graph
**************************************************************************/
void MocInit2WayPartition2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  int dbglvl;

  dbglvl = ctrl->dbglvl;
  IFSET(ctrl->dbglvl, DBG_REFINE, ctrl->dbglvl -= DBG_REFINE);
  IFSET(ctrl->dbglvl, DBG_MOVEINFO, ctrl->dbglvl -= DBG_MOVEINFO);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  switch (ctrl->IType) {
    case IPART_GGPKL:
    case IPART_RANDOM:
      MocGrowBisection2(ctrl, graph, tpwgts, ubvec);
      break;
    case 3:
      MocGrowBisectionNew2(ctrl, graph, tpwgts, ubvec);
      break;
    default:
      errexit("Unknown initial partition type: %d\n", ctrl->IType);
  }

  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial Cut: %d\n", graph->mincut));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  ctrl->dbglvl = dbglvl;

}




/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void MocGrowBisection2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  int nvtxs, bestcut, nbfs;
  idxtype *bestwhere, *where;

  nvtxs = graph->nvtxs;

  MocAllocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  nbfs = 2*(nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = idxsum(graph->nedges, graph->adjwgt);

  for (; nbfs>0; nbfs--) {
    idxset(nvtxs, 1, where);
    where[RandomInRange(nvtxs)] = 0;

    MocCompute2WayPartitionParams(ctrl, graph);

    MocBalance2Way2(ctrl, graph, tpwgts, ubvec);

    MocFM_2WayEdgeRefine2(ctrl, graph, tpwgts, ubvec, 4);

    MocBalance2Way2(ctrl, graph, tpwgts, ubvec);
    MocFM_2WayEdgeRefine2(ctrl, graph, tpwgts, ubvec, 4);

    if (bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  GKfree((void **) &bestwhere, LTERM);
}






/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void MocGrowBisectionNew2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  int nvtxs, bestcut, nbfs;
  idxtype *bestwhere, *where;

  nvtxs = graph->nvtxs;

  MocAllocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  nbfs = 2*(nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = idxsum(graph->nedges, graph->adjwgt);

  for (; nbfs>0; nbfs--) {
    idxset(nvtxs, 1, where);
    where[RandomInRange(nvtxs)] = 0;

    MocCompute2WayPartitionParams(ctrl, graph);

    MocInit2WayBalance2(ctrl, graph, tpwgts, ubvec);

    MocFM_2WayEdgeRefine2(ctrl, graph, tpwgts, ubvec, 4);

    if (bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  GKfree((void **) &bestwhere, LTERM);
}



/*************************************************************************
* This function balances two partitions by moving the highest gain
* (including negative gain) vertices to the other domain.
* It is used only when tha unbalance is due to non contigous
* subdomains. That is, the are no boundary vertices.
* It moves vertices from the domain that is overweight to the one that
* is underweight.
**************************************************************************/
void MocInit2WayBalance2(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  int i, ii, j, k, l, kwgt, nvtxs, nbnd, ncon, nswaps, from, to, cnum, tmp, imin;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *moved, *perm, *qnum;
  float *nvwgt, *npwgts, minwgt;
  PQueueType parts[MAXNCON][2];
  int higain, oldgain, mincut;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  nvwgt = graph->nvwgt;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  /* This is called for initial partitioning so we know from where to pick nodes */
  from = 1;
  to = (from+1)%2;

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Parts: [");
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("] T[%.3f %.3f], Nv-Nb[%5d, %5d]. ICut: %6d, LB: %.3f [B]\n", tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut, ComputeLoadImbalance(ncon, 2, npwgts, tpwgts));
  }

  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
  }

  idxset(nvtxs, -1, moved);

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));
  ASSERT(CheckGraph(graph));

  /* Compute the queues in which each vertex will be assigned to */
  for (i=0; i<nvtxs; i++)
    qnum[i] = samax(ncon, nvwgt+i*ncon);

  /* Insert the nodes of the proper partition in the appropriate priority queue */
  RandomPermute(nvtxs, perm, 1);
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    if (where[i] == from) {
      if (ed[i] > 0)
        PQueueInsert(&parts[qnum[i]][0], i, ed[i]-id[i]);
      else
        PQueueInsert(&parts[qnum[i]][1], i, ed[i]-id[i]);
    }
  }

/*
  for (i=0; i<ncon; i++)
    printf("Queue #%d has %d %d\n", i, parts[i][0].nnodes, parts[i][1].nnodes);
*/

  /* Determine the termination criterion */
  imin = 0;
  for (i=1; i<ncon; i++)
    imin = (ubvec[i] < ubvec[imin] ? i : imin);
  minwgt = .5/ubvec[imin];

  mincut = graph->mincut;
  nbnd = graph->nbnd;
  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    /* Exit as soon as the minimum weight crossed over */
    if (npwgts[to*ncon+imin] > minwgt)
      break;

    if ((cnum = SelectQueueOneWay2(ncon, npwgts+to*ncon, parts, ubvec)) == -1)
      break;

    if ((higain = PQueueGetMax(&parts[cnum][0])) == -1)
      higain = PQueueGetMax(&parts[cnum][1]);

    mincut -= (ed[higain]-id[higain]);
    saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

    where[higain] = to;
    moved[higain] = nswaps;

    if (ctrl->dbglvl&DBG_MOVEINFO) {
      printf("Moved %6d from %d(%d). [%5d] %5d, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], mincut);
      for (l=0; l<ncon; l++)
        printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
      printf(", LB: %.3f\n", ComputeLoadImbalance(ncon, 2, npwgts, tpwgts));
      if (ed[higain] == 0 && id[higain] > 0)
        printf("\t Pulled from the interior!\n");
    }


    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (moved[k] == -1 && where[k] == from) {
        if (ed[k] > 0 && bndptr[k] == -1) {  /* It moves in boundary */
          PQueueDelete(&parts[qnum[k]][1], k, oldgain);
          PQueueInsert(&parts[qnum[k]][0], k, ed[k]-id[k]);
        }
        else { /* It must be in the boundary already */
          if (bndptr[k] == -1)
            printf("What you thought was wrong!\n");
          PQueueUpdate(&parts[qnum[k]][0], k, oldgain, ed[k]-id[k]);
        }
      }

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }

    ASSERTP(ComputeCut(graph, where) == mincut, ("%d != %d\n", ComputeCut(graph, where), mincut));

  }

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("\tMincut: %6d, NBND: %6d, NPwgts: ", mincut, nbnd);
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf(", LB: %.3f\n", ComputeLoadImbalance(ncon, 2, npwgts, tpwgts));
  }

  graph->mincut = mincut;
  graph->nbnd = nbnd;

  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
int SelectQueueOneWay2(int ncon, float *pto, PQueueType queues[MAXNCON][2], float *ubvec)
{
  int i, cnum=-1, imax, maxgain;
  float max=0.0;
  float twgt[MAXNCON];

  for (i=0; i<ncon; i++) {
    if (max < pto[i]) {
      imax = i;
      max = pto[i];
    }
  }
  for (i=0; i<ncon; i++)
    twgt[i] = (max/(ubvec[imax]*ubvec[i]))/pto[i];
  twgt[imax] = 0.0;

  max = 0.0;
  for (i=0; i<ncon; i++) {
    if (max < twgt[i] && (PQueueGetSize(&queues[i][0]) > 0 || PQueueGetSize(&queues[i][1]) > 0)) {
      max = twgt[i];
      cnum = i;
    }
  }
  if (max > 1)
    return cnum;

  /* optimize of cut */
  maxgain = -10000000;
  for (i=0; i<ncon; i++) {
    if (PQueueGetSize(&queues[i][0]) > 0 && PQueueGetKey(&queues[i][0]) > maxgain) {
      maxgain = PQueueGetKey(&queues[i][0]);
      cnum = i;
    }
  }

  return cnum;

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * minitpart.c
 *
 * This file contains code that performs the initial partition of the
 * coarsest graph
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function computes the initial bisection of the coarsest graph
**************************************************************************/
void MocInit2WayPartition(CtrlType *ctrl, GraphType *graph, float *tpwgts, float ubfactor)
{
  int dbglvl;

  dbglvl = ctrl->dbglvl;
  IFSET(ctrl->dbglvl, DBG_REFINE, ctrl->dbglvl -= DBG_REFINE);
  IFSET(ctrl->dbglvl, DBG_MOVEINFO, ctrl->dbglvl -= DBG_MOVEINFO);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  switch (ctrl->IType) {
    case IPART_GGPKL:
      MocGrowBisection(ctrl, graph, tpwgts, ubfactor);
      break;
    case IPART_RANDOM:
      MocRandomBisection(ctrl, graph, tpwgts, ubfactor);
      break;
    default:
      errexit("Unknown initial partition type: %d\n", ctrl->IType);
  }

  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial Cut: %d [%d]\n", graph->mincut, graph->where[0]));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  ctrl->dbglvl = dbglvl;

}





/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void MocGrowBisection(CtrlType *ctrl, GraphType *graph, float *tpwgts, float ubfactor)
{
  int nvtxs, bestcut, nbfs;
  idxtype *bestwhere, *where;

  nvtxs = graph->nvtxs;

  MocAllocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  nbfs = 2*(nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = idxsum(graph->nedges, graph->adjwgt);

  for (; nbfs>0; nbfs--) {
    idxset(nvtxs, 1, where);
    where[RandomInRange(nvtxs)] = 0;

    MocCompute2WayPartitionParams(ctrl, graph);

    MocInit2WayBalance(ctrl, graph, tpwgts);

    MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);

    MocBalance2Way(ctrl, graph, tpwgts, 1.02);
    MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);

    if (bestcut >= graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  GKfree((void **) &bestwhere, LTERM);
}



/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void MocRandomBisection(CtrlType *ctrl, GraphType *graph, float *tpwgts, float ubfactor)
{
  int i, ii, nvtxs, ncon, bestcut, nbfs, qnum;
  idxtype *bestwhere, *where, *perm;
  int counts[MAXNCON];
  float *nvwgt;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  nvwgt = graph->nvwgt;

  MocAllocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  nbfs = 2*(nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  bestcut = idxsum(graph->nedges, graph->adjwgt);
  perm = idxmalloc(nvtxs, "BisectGraph: perm");

  for (; nbfs>0; nbfs--) {
    for (i=0; i<ncon; i++)
      counts[i] = 0;

    RandomPermute(nvtxs, perm, 1);

    /* Partition by spliting the queues randomly */
    for (ii=0; ii<nvtxs; ii++) {
      i = perm[ii];
      qnum = samax(ncon, nvwgt+i*ncon);
      where[i] = counts[qnum];
      counts[qnum] = (counts[qnum]+1)%2;
    }

    MocCompute2WayPartitionParams(ctrl, graph);

    MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 6);
    MocBalance2Way(ctrl, graph, tpwgts, 1.02);
    MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 6);
    MocBalance2Way(ctrl, graph, tpwgts, 1.02);
    MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 6);

    /*
    printf("Edgecut: %6d, NPwgts: [", graph->mincut);
    for (i=0; i<graph->ncon; i++)
      printf("(%.3f %.3f) ", graph->npwgts[i], graph->npwgts[graph->ncon+i]);
    printf("]\n");
    */

    if (bestcut >= graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  GKfree((void **) &bestwhere, &perm, LTERM);
}




/*************************************************************************
* This function balances two partitions by moving the highest gain
* (including negative gain) vertices to the other domain.
* It is used only when tha unbalance is due to non contigous
* subdomains. That is, the are no boundary vertices.
* It moves vertices from the domain that is overweight to the one that
* is underweight.
**************************************************************************/
void MocInit2WayBalance(CtrlType *ctrl, GraphType *graph, float *tpwgts)
{
  int i, ii, j, k, l, kwgt, nvtxs, nbnd, ncon, nswaps, from, to, cnum, tmp;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *perm, *qnum;
  float *nvwgt, *npwgts;
  PQueueType parts[MAXNCON][2];
  int higain, oldgain, mincut;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  nvwgt = graph->nvwgt;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  /* This is called for initial partitioning so we know from where to pick nodes */
  from = 1;
  to = (from+1)%2;

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Parts: [");
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf("] T[%.3f %.3f], Nv-Nb[%5d, %5d]. ICut: %6d, LB: %.3f [B]\n", tpwgts[0], tpwgts[1],
           graph->nvtxs, graph->nbnd, graph->mincut,
           Compute2WayHLoadImbalance(ncon, npwgts, tpwgts));
  }

  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
  }

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));
  ASSERT(CheckGraph(graph));

  /* Compute the queues in which each vertex will be assigned to */
  for (i=0; i<nvtxs; i++)
    qnum[i] = samax(ncon, nvwgt+i*ncon);

  /* Insert the nodes of the proper partition in the appropriate priority queue */
  RandomPermute(nvtxs, perm, 1);
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    if (where[i] == from) {
      if (ed[i] > 0)
        PQueueInsert(&parts[qnum[i]][0], i, ed[i]-id[i]);
      else
        PQueueInsert(&parts[qnum[i]][1], i, ed[i]-id[i]);
    }
  }


  mincut = graph->mincut;
  nbnd = graph->nbnd;
  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if (AreAnyVwgtsBelow(ncon, 1.0, npwgts+from*ncon, 0.0, nvwgt, tpwgts[from]))
      break;

    if ((cnum = SelectQueueOneWay(ncon, npwgts, tpwgts, from, parts)) == -1)
      break;

    if ((higain = PQueueGetMax(&parts[cnum][0])) == -1)
      higain = PQueueGetMax(&parts[cnum][1]);

    mincut -= (ed[higain]-id[higain]);
    saxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    saxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

    where[higain] = to;

    if (ctrl->dbglvl&DBG_MOVEINFO) {
      printf("Moved %6d from %d(%d). [%5d] %5d, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], mincut);
      for (l=0; l<ncon; l++)
        printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
      printf(", LB: %.3f\n", Compute2WayHLoadImbalance(ncon, npwgts, tpwgts));
      if (ed[higain] == 0 && id[higain] > 0)
        printf("\t Pulled from the interior!\n");
    }


    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (where[k] == from) {
        if (ed[k] > 0 && bndptr[k] == -1) {  /* It moves in boundary */
          PQueueDelete(&parts[qnum[k]][1], k, oldgain);
          PQueueInsert(&parts[qnum[k]][0], k, ed[k]-id[k]);
        }
        else { /* It must be in the boundary already */
          if (bndptr[k] == -1)
            printf("What you thought was wrong!\n");
          PQueueUpdate(&parts[qnum[k]][0], k, oldgain, ed[k]-id[k]);
        }
      }

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }

    ASSERTP(ComputeCut(graph, where) == mincut, ("%d != %d\n", ComputeCut(graph, where), mincut));

  }

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("\tMincut: %6d, NBND: %6d, NPwgts: ", mincut, nbnd);
    for (l=0; l<ncon; l++)
      printf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    printf(", LB: %.3f\n", Compute2WayHLoadImbalance(ncon, npwgts, tpwgts));
  }

  graph->mincut = mincut;
  graph->nbnd = nbnd;

  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}




/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
int SelectQueueOneWay(int ncon, float *npwgts, float *tpwgts, int from, PQueueType queues[MAXNCON][2])
{
  int i, cnum=-1;
  float max=0.0;

  for (i=0; i<ncon; i++) {
    if (npwgts[from*ncon+i]-tpwgts[from] >= max &&
        PQueueGetSize(&queues[i][0]) + PQueueGetSize(&queues[i][1]) > 0) {
      max = npwgts[from*ncon+i]-tpwgts[0];
      cnum = i;
    }
  }

  return cnum;
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mkmetis.c
 *
 * This file contains the top level routines for the multilevel k-way partitioning
 * algorithm KMETIS.
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 *
 */





/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_mCPartGraphKway(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
                          idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag,
                          int *nparts, float *rubvec, int *options, int *edgecut,
                          idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, *ncon, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType  = McKMETIS_CTYPE;
    ctrl.IType  = McKMETIS_ITYPE;
    ctrl.RType  = McKMETIS_RTYPE;
    ctrl.dbglvl = McKMETIS_DBGLVL;
  }
  else {
    ctrl.CType  = options[OPTION_CTYPE];
    ctrl.IType  = options[OPTION_ITYPE];
    ctrl.RType  = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = amax((*nvtxs)/(20*log2(*nparts)), 30*(*nparts));

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MCMlevelKWayPartitioning(&ctrl, &graph, *nparts, part, rubvec);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MCMlevelKWayPartitioning(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part,
      float *rubvec)
{
  int i;
  GraphType *cgraph;
  int options[10], edgecut;

  cgraph = MCCoarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
  MocAllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  options[0] = 1;
  options[OPTION_CTYPE] = MATCH_SBHEM_INFNORM;
  options[OPTION_ITYPE] = IPART_RANDOM;
  options[OPTION_RTYPE] = RTYPE_FM;
  options[OPTION_DBGLVL] = 0;

  /* Determine what you will use as the initial partitioner, based on tolerances */
  for (i=0; i<graph->ncon; i++) {
    if (rubvec[i] > 1.2)
      break;
  }
  if (i == graph->ncon)
    METIS_mCPartGraphRecursiveInternal(&cgraph->nvtxs, &cgraph->ncon,
          cgraph->xadj, cgraph->adjncy, cgraph->nvwgt, cgraph->adjwgt, &nparts,
          options, &edgecut, cgraph->where);
  else
    METIS_mCHPartGraphRecursiveInternal(&cgraph->nvtxs, &cgraph->ncon,
          cgraph->xadj, cgraph->adjncy, cgraph->nvwgt, cgraph->adjwgt, &nparts,
          rubvec, options, &edgecut, cgraph->where);


  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial %d-way partitioning cut: %d\n", nparts, edgecut));

  IFSET(ctrl->dbglvl, DBG_KWAYPINFO, ComputePartitionInfo(cgraph, nparts, cgraph->where));

  MocRefineKWayHorizontal(ctrl, graph, cgraph, nparts, rubvec);

  idxcopy(graph->nvtxs, graph->where, part);

  GKfree((void **) &graph->nvwgt, (void **) &graph->npwgts,
         (void **) &graph->gdata, (void **) &graph->rdata, LTERM);

  return graph->mincut;

}

/*
 * mkwayfmh.c
 *
 * This file contains code that implements the multilevel k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 *
 */





/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void MCRandom_KWayEdgeRefineHorizontal(CtrlType *ctrl, GraphType *graph, int nparts,
       float *orgubvec, int npasses)
{
  int i, ii, iii, j, k, pass, nvtxs, ncon, nmoves, nbnd, myndegrees, same;
  int from, me, to, oldcut, gain;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *perm, *bndptr, *bndind;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;
  float *npwgts, *nvwgt, *minwgt, *maxwgt, maxlb, minlb, ubvec[MAXNCON], tvec[MAXNCON];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  npwgts = graph->npwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  fwspacemalloc(ctrl, nparts*ncon);
  maxwgt = fwspacemalloc(ctrl, nparts*ncon);

  /* See if the orgubvec consists of identical constraints */
  maxlb = minlb = orgubvec[0];
  for (i=1; i<ncon; i++) {
    minlb = (orgubvec[i] < minlb ? orgubvec[i] : minlb);
    maxlb = (orgubvec[i] > maxlb ? orgubvec[i] : maxlb);
  }
  same = (fabs(maxlb-minlb) < .01 ? 1 : 0);


  /* Let's not get very optimistic. Let Balancing do the work */
  ComputeHKWayLoadImbalance(ncon, nparts, npwgts, ubvec);
  for (i=0; i<ncon; i++)
    ubvec[i] = amax(ubvec[i], orgubvec[i]);

  if (!same) {
    for (i=0; i<nparts; i++) {
      for (j=0; j<ncon; j++) {
        maxwgt[i*ncon+j] = ubvec[j]/nparts;
        minwgt[i*ncon+j] = 1.0/(ubvec[j]*nparts);
      }
    }
  }
  else {
    maxlb = ubvec[0];
    for (i=1; i<ncon; i++)
      maxlb = (ubvec[i] > maxlb ? ubvec[i] : maxlb);

    for (i=0; i<nparts; i++) {
      for (j=0; j<ncon; j++) {
        maxwgt[i*ncon+j] = maxlb/nparts;
        minwgt[i*ncon+j] = 1.0/(maxlb*nparts);
      }
    }
  }


  perm = idxwspacemalloc(ctrl, nvtxs);

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Partitions: [%5.4f %5.4f], Nv-Nb[%6d %6d]. Cut: %6d, LB: ",
            npwgts[samin(ncon*nparts, npwgts)], npwgts[samax(ncon*nparts, npwgts)],
            graph->nvtxs, graph->nbnd, graph->mincut);
    ComputeHKWayLoadImbalance(ncon, nparts, npwgts, tvec);
    for (i=0; i<ncon; i++)
      printf("%.3f ", tvec[i]);
    printf("\n");
  }

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (nmoves=iii=0; iii<graph->nbnd; iii++) {
      ii = perm[iii];
      if (ii >= nbnd)
        continue;
      i = bndind[ii];

      myrinfo = graph->rinfo+i;

      if (myrinfo->ed >= myrinfo->id) { /* Total ED is too high */
        from = where[i];
        nvwgt = graph->nvwgt+i*ncon;

        if (myrinfo->id > 0 && AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, -1.0, nvwgt, minwgt+from*ncon))
          continue;   /* This cannot be moved! */

        myedegrees = myrinfo->edegrees;
        myndegrees = myrinfo->ndegrees;

        for (k=0; k<myndegrees; k++) {
          to = myedegrees[k].pid;
          gain = myedegrees[k].ed - myrinfo->id;
          if (gain >= 0 &&
              (AreAllHVwgtsBelow(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, maxwgt+to*ncon) ||
               IsHBalanceBetterFT(ncon, nparts, npwgts+from*ncon, npwgts+to*ncon, nvwgt, ubvec)))
            break;
        }
        if (k == myndegrees)
          continue;  /* break out if you did not find a candidate */

        for (j=k+1; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          if ((myedegrees[j].ed > myedegrees[k].ed &&
               (AreAllHVwgtsBelow(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, maxwgt+to*ncon) ||
               IsHBalanceBetterFT(ncon, nparts, npwgts+from*ncon, npwgts+to*ncon, nvwgt, ubvec))) ||
              (myedegrees[j].ed == myedegrees[k].ed &&
               IsHBalanceBetterTT(ncon, nparts, npwgts+myedegrees[k].pid*ncon, npwgts+to*ncon, nvwgt, ubvec)))
            k = j;
        }

        to = myedegrees[k].pid;

        if (myedegrees[k].ed-myrinfo->id == 0
            && !IsHBalanceBetterFT(ncon, nparts, npwgts+from*ncon, npwgts+to*ncon, nvwgt, ubvec)
            && AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, 0.0, npwgts+from*ncon, maxwgt+from*ncon))
          continue;

        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to'
        *======================================================================*/
        graph->mincut -= myedegrees[k].ed-myrinfo->id;

        IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

        /* Update where, weight, and ID/ED information of the vertex you moved */
        saxpy(ncon, 1.0, nvwgt, 1, npwgts+to*ncon, 1);
        saxpy(ncon, -1.0, nvwgt, 1, npwgts+from*ncon, 1);
        where[i] = to;
        myrinfo->ed += myrinfo->id-myedegrees[k].ed;
        SWAP(myrinfo->id, myedegrees[k].ed, j);
        if (myedegrees[k].ed == 0)
          myedegrees[k] = myedegrees[--myrinfo->ndegrees];
        else
          myedegrees[k].pid = from;

        if (myrinfo->ed-myrinfo->id < 0)
          BNDDelete(nbnd, bndind, bndptr, i);

        /* Update the degrees of adjacent vertices */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          ii = adjncy[j];
          me = where[ii];

          myrinfo = graph->rinfo+ii;
          if (myrinfo->edegrees == NULL) {
            myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
            ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
          }
          myedegrees = myrinfo->edegrees;

          ASSERT(CheckRInfo(myrinfo));

          if (me == from) {
            INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

            if (myrinfo->ed-myrinfo->id >= 0 && bndptr[ii] == -1)
              BNDInsert(nbnd, bndind, bndptr, ii);
          }
          else if (me == to) {
            INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

            if (myrinfo->ed-myrinfo->id < 0 && bndptr[ii] != -1)
              BNDDelete(nbnd, bndind, bndptr, ii);
          }

          /* Remove contribution from the .ed of 'from' */
          if (me != from) {
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (myedegrees[k].pid == from) {
                if (myedegrees[k].ed == adjwgt[j])
                  myedegrees[k] = myedegrees[--myrinfo->ndegrees];
                else
                  myedegrees[k].ed -= adjwgt[j];
                break;
              }
            }
          }

          /* Add contribution to the .ed of 'to' */
          if (me != to) {
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (myedegrees[k].pid == to) {
                myedegrees[k].ed += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->ndegrees) {
              myedegrees[myrinfo->ndegrees].pid = to;
              myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
            }
          }

          ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
          ASSERT(CheckRInfo(myrinfo));

        }
        nmoves++;
      }
    }

    graph->nbnd = nbnd;

    if (ctrl->dbglvl&DBG_REFINE) {
      printf("\t [%5.4f %5.4f], Nb: %6d, Nmoves: %5d, Cut: %6d, LB: ",
              npwgts[samin(ncon*nparts, npwgts)], npwgts[samax(ncon*nparts, npwgts)],
              nbnd, nmoves, graph->mincut);
      ComputeHKWayLoadImbalance(ncon, nparts, npwgts, tvec);
      for (i=0; i<ncon; i++)
        printf("%.3f ", tvec[i]);
      printf("\n");
    }

    if (graph->mincut == oldcut)
      break;
  }

  fwspacefree(ctrl, ncon*nparts);
  fwspacefree(ctrl, ncon*nparts);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void MCGreedy_KWayEdgeBalanceHorizontal(CtrlType *ctrl, GraphType *graph, int nparts,
       float *ubvec, int npasses)
{
  int i, ii, j, k, pass, nvtxs, ncon, nbnd, myndegrees, oldgain, gain, nmoves;
  int from, me, to, oldcut;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *perm, *bndptr, *bndind, *moved;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;
  PQueueType queue;
  float *npwgts, *nvwgt, *minwgt, *maxwgt, tvec[MAXNCON];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;

  where = graph->where;
  npwgts = graph->npwgts;

  /* Setup the weight intervals of the various subdomains */
  minwgt =  fwspacemalloc(ctrl, ncon*nparts);
  maxwgt = fwspacemalloc(ctrl, ncon*nparts);

  for (i=0; i<nparts; i++) {
    for (j=0; j<ncon; j++) {
      maxwgt[i*ncon+j] = ubvec[j]/nparts;
      minwgt[i*ncon+j] = 1.0/(ubvec[j]*nparts);
    }
  }

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxwspacemalloc(ctrl, nvtxs);

  PQueueInit(ctrl, &queue, nvtxs, graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)]);

  if (ctrl->dbglvl&DBG_REFINE) {
    printf("Partitions: [%5.4f %5.4f], Nv-Nb[%6d %6d]. Cut: %6d, LB: ",
            npwgts[samin(ncon*nparts, npwgts)], npwgts[samax(ncon*nparts, npwgts)],
            graph->nvtxs, graph->nbnd, graph->mincut);
    ComputeHKWayLoadImbalance(ncon, nparts, npwgts, tvec);
    for (i=0; i<ncon; i++)
      printf("%.3f ", tvec[i]);
    printf("[B]\n");
  }


  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    /* Check to see if things are out of balance, given the tolerance */
    if (MocIsHBalanced(ncon, nparts, npwgts, ubvec))
      break;

    PQueueReset(&queue);
    idxset(nvtxs, -1, moved);

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      PQueueInsert(&queue, i, graph->rinfo[i].ed - graph->rinfo[i].id);
      moved[i] = 2;
    }

    nmoves = 0;
    for (;;) {
      if ((i = PQueueGetMax(&queue)) == -1)
        break;
      moved[i] = 1;

      myrinfo = graph->rinfo+i;
      from = where[i];
      nvwgt = graph->nvwgt+i*ncon;

      if (AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, -1.0, nvwgt, minwgt+from*ncon))
        continue;   /* This cannot be moved! */

      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      for (k=0; k<myndegrees; k++) {
        to = myedegrees[k].pid;
        if (IsHBalanceBetterFT(ncon, nparts, npwgts+from*ncon, npwgts+to*ncon, nvwgt, ubvec))
          break;
      }
      if (k == myndegrees)
        continue;  /* break out if you did not find a candidate */

      for (j=k+1; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        if (IsHBalanceBetterTT(ncon, nparts, npwgts+myedegrees[k].pid*ncon, npwgts+to*ncon, nvwgt, ubvec))
          k = j;
      }

      to = myedegrees[k].pid;

      j = 0;
      if (!AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, 0.0, nvwgt, maxwgt+from*ncon))
        j++;
      if (myedegrees[k].ed-myrinfo->id >= 0)
        j++;
      if (!AreAllHVwgtsAbove(ncon, 1.0, npwgts+to*ncon, 0.0, nvwgt, minwgt+to*ncon) &&
          AreAllHVwgtsBelow(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, maxwgt+to*ncon))
        j++;
      if (j == 0)
        continue;

/* DELETE
      if (myedegrees[k].ed-myrinfo->id < 0 &&
          AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, 0.0, nvwgt, maxwgt+from*ncon) &&
          AreAllHVwgtsAbove(ncon, 1.0, npwgts+to*ncon, 0.0, nvwgt, minwgt+to*ncon) &&
          AreAllHVwgtsBelow(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, maxwgt+to*ncon))
        continue;
*/
      /*=====================================================================
      * If we got here, we can now move the vertex from 'from' to 'to'
      *======================================================================*/
      graph->mincut -= myedegrees[k].ed-myrinfo->id;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

      /* Update where, weight, and ID/ED information of the vertex you moved */
      saxpy(ncon, 1.0, nvwgt, 1, npwgts+to*ncon, 1);
      saxpy(ncon, -1.0, nvwgt, 1, npwgts+from*ncon, 1);
      where[i] = to;
      myrinfo->ed += myrinfo->id-myedegrees[k].ed;
      SWAP(myrinfo->id, myedegrees[k].ed, j);
      if (myedegrees[k].ed == 0)
        myedegrees[k] = myedegrees[--myrinfo->ndegrees];
      else
        myedegrees[k].pid = from;

      if (myrinfo->ed == 0)
        BNDDelete(nbnd, bndind, bndptr, i);

      /* Update the degrees of adjacent vertices */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        me = where[ii];

        myrinfo = graph->rinfo+ii;
        if (myrinfo->edegrees == NULL) {
          myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
          ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
        }
        myedegrees = myrinfo->edegrees;

        ASSERT(CheckRInfo(myrinfo));

        oldgain = (myrinfo->ed-myrinfo->id);

        if (me == from) {
          INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

          if (myrinfo->ed > 0 && bndptr[ii] == -1)
            BNDInsert(nbnd, bndind, bndptr, ii);
        }
        else if (me == to) {
          INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

          if (myrinfo->ed == 0 && bndptr[ii] != -1)
            BNDDelete(nbnd, bndind, bndptr, ii);
        }

        /* Remove contribution from the .ed of 'from' */
        if (me != from) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == from) {
              if (myedegrees[k].ed == adjwgt[j])
                myedegrees[k] = myedegrees[--myrinfo->ndegrees];
              else
                myedegrees[k].ed -= adjwgt[j];
              break;
            }
          }
        }

        /* Add contribution to the .ed of 'to' */
        if (me != to) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == to) {
              myedegrees[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].pid = to;
            myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
          }
        }


        /* Update the queue */
        if (me == to || me == from) {
          gain = myrinfo->ed-myrinfo->id;
          if (moved[ii] == 2) {
            if (myrinfo->ed > 0)
              PQueueUpdate(&queue, ii, oldgain, gain);
            else {
              PQueueDelete(&queue, ii, oldgain);
              moved[ii] = -1;
            }
          }
          else if (moved[ii] == -1 && myrinfo->ed > 0) {
            PQueueInsert(&queue, ii, gain);
            moved[ii] = 2;
          }
        }

        ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
        ASSERT(CheckRInfo(myrinfo));
      }
      nmoves++;
    }

    graph->nbnd = nbnd;

    if (ctrl->dbglvl&DBG_REFINE) {
      printf("\t [%5.4f %5.4f], Nb: %6d, Nmoves: %5d, Cut: %6d, LB: ",
              npwgts[samin(ncon*nparts, npwgts)], npwgts[samax(ncon*nparts, npwgts)],
              nbnd, nmoves, graph->mincut);
      ComputeHKWayLoadImbalance(ncon, nparts, npwgts, tvec);
      for (i=0; i<ncon; i++)
        printf("%.3f ", tvec[i]);
      printf("\n");
    }

    if (nmoves == 0)
      break;
  }

  PQueueFree(ctrl, &queue);

  fwspacefree(ctrl, ncon*nparts);
  fwspacefree(ctrl, ncon*nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}





/*************************************************************************
* This function checks if the vertex weights of two vertices are below
* a given set of values
**************************************************************************/
int AreAllHVwgtsBelow(int ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float *limit)
{
  int i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] > limit[i])
      return 0;

  return 1;
}



/*************************************************************************
* This function checks if the vertex weights of two vertices are above
* a given set of values
**************************************************************************/
int AreAllHVwgtsAbove(int ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float *limit)
{
  int i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] < limit[i])
      return 0;

  return 1;
}


/*************************************************************************
* This function computes the load imbalance over all the constrains
* For now assume that we just want balanced partitionings
**************************************************************************/
void ComputeHKWayLoadImbalance(int ncon, int nparts, float *npwgts, float *lbvec)
{
  int i, j;
  float max;

  for (i=0; i<ncon; i++) {
    max = 0.0;
    for (j=0; j<nparts; j++) {
      if (npwgts[j*ncon+i] > max)
        max = npwgts[j*ncon+i];
    }

    lbvec[i] = max*nparts;
  }
}


/*************************************************************************
* This function determines if a partitioning is horizontally balanced
**************************************************************************/
int MocIsHBalanced(int ncon, int nparts, float *npwgts, float *ubvec)
{
  int i, j;
  float max;

  for (i=0; i<ncon; i++) {
    max = 0.0;
    for (j=0; j<nparts; j++) {
      if (npwgts[j*ncon+i] > max)
        max = npwgts[j*ncon+i];
    }

    if (ubvec[i] < max*nparts)
      return 0;
  }

  return 1;
}





/*************************************************************************
* This function checks if the pairwise balance of the between the two
* partitions will improve by moving the vertex v from pfrom to pto,
* subject to the target partition weights of tfrom, and tto respectively
**************************************************************************/
int IsHBalanceBetterFT(int ncon, int nparts, float *pfrom, float *pto, float *vwgt, float *ubvec)
{
  int i;
  float blb1=0.0, alb1=0.0, sblb=0.0, salb=0.0;
  float blb2=0.0, alb2=0.0;
  float temp;

  for (i=0; i<ncon; i++) {
    temp = amax(pfrom[i], pto[i])*nparts/ubvec[i];
    if (blb1 < temp) {
      blb2 = blb1;
      blb1 = temp;
    }
    else if (blb2 < temp)
      blb2 = temp;
    sblb += temp;

    temp = amax(pfrom[i]-vwgt[i], pto[i]+vwgt[i])*nparts/ubvec[i];
    if (alb1 < temp) {
      alb2 = alb1;
      alb1 = temp;
    }
    else if (alb2 < temp)
      alb2 = temp;
    salb += temp;
  }

  if (alb1 < blb1)
    return 1;
  if (blb1 < alb1)
    return 0;
  if (alb2 < blb2)
    return 1;
  if (blb2 < alb2)
    return 0;

  return salb < sblb;

}




/*************************************************************************
* This function checks if it will be better to move a vertex to pt2 than
* to pt1 subject to their target weights of tt1 and tt2, respectively
* This routine takes into account the weight of the vertex in question
**************************************************************************/
int IsHBalanceBetterTT(int ncon, int nparts, float *pt1, float *pt2, float *vwgt, float *ubvec)
{
  int i;
  float m11=0.0, m12=0.0, m21=0.0, m22=0.0, sm1=0.0, sm2=0.0, temp;

  for (i=0; i<ncon; i++) {
    temp = (pt1[i]+vwgt[i])*nparts/ubvec[i];
    if (m11 < temp) {
      m12 = m11;
      m11 = temp;
    }
    else if (m12 < temp)
      m12 = temp;
    sm1 += temp;

    temp = (pt2[i]+vwgt[i])*nparts/ubvec[i];
    if (m21 < temp) {
      m22 = m21;
      m21 = temp;
    }
    else if (m22 < temp)
      m22 = temp;
    sm2 += temp;
  }

  if (m21 < m11)
    return 1;
  if (m21 > m11)
    return 0;
  if (m22 < m12)
    return 1;
  if (m22 > m12)
    return 0;

  return sm2 < sm1;
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mkwayrefine.c
 *
 * This file contains the driving routines for multilevel k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void MocRefineKWayHorizontal(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int nparts,
       float *ubvec)
{

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  MocComputeKWayPartitionParams(ctrl, graph, nparts);

  for (;;) {
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));

    if (!MocIsHBalanced(graph->ncon, nparts, graph->npwgts, ubvec)) {
      MocComputeKWayBalanceBoundary(ctrl, graph, nparts);
      MCGreedy_KWayEdgeBalanceHorizontal(ctrl, graph, nparts, ubvec, 4);
      ComputeKWayBoundary(ctrl, graph, nparts);
    }

    MCRandom_KWayEdgeRefineHorizontal(ctrl, graph, nparts, ubvec, 10);

    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    MocProjectKWayPartition(ctrl, graph, nparts);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  if (!MocIsHBalanced(graph->ncon, nparts, graph->npwgts, ubvec)) {
    MocComputeKWayBalanceBoundary(ctrl, graph, nparts);
    MCGreedy_KWayEdgeBalanceHorizontal(ctrl, graph, nparts, ubvec, 4);
    ComputeKWayBoundary(ctrl, graph, nparts);
    MCRandom_KWayEdgeRefineHorizontal(ctrl, graph, nparts, ubvec, 10);
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}




/*************************************************************************
* This function allocates memory for k-way edge refinement
**************************************************************************/
void MocAllocateKWayPartitionMemory(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int nvtxs, ncon, pad64;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;

  pad64 = (3*nvtxs)%2;

  graph->rdata = idxmalloc(3*nvtxs+(sizeof(RInfoType)/sizeof(idxtype))*nvtxs+pad64, "AllocateKWayPartitionMemory: rdata");
  graph->where          = graph->rdata;
  graph->bndptr         = graph->rdata + nvtxs;
  graph->bndind         = graph->rdata + 2*nvtxs;
  graph->rinfo          = (RInfoType *)(graph->rdata + 3*nvtxs + pad64);

  graph->npwgts         = fmalloc(ncon*nparts, "MocAllocateKWayPartitionMemory: npwgts");
}


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void MocComputeKWayPartitionParams(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, j, k, nvtxs, ncon, nbnd, mincut, me, other;
  idxtype *xadj, *adjncy, *adjwgt, *where, *bndind, *bndptr;
  RInfoType *rinfo, *myrinfo;
  EDegreeType *myedegrees;
  float *nvwgt, *npwgts;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  npwgts = sset(ncon*nparts, 0.0, graph->npwgts);
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  rinfo = graph->rinfo;


  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  ctrl->wspace.cdegree = 0;
  nbnd = mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    saxpy(ncon, 1.0, nvwgt+i*ncon, 1, npwgts+me*ncon, 1);

    myrinfo = rinfo+i;
    myrinfo->id = myrinfo->ed = myrinfo->ndegrees = 0;
    myrinfo->edegrees = NULL;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me != where[adjncy[j]])
        myrinfo->ed += adjwgt[j];
    }
    myrinfo->id = graph->adjwgtsum[i] - myrinfo->ed;

    if (myrinfo->ed > 0)
      mincut += myrinfo->ed;

    if (myrinfo->ed-myrinfo->id >= 0)
      BNDInsert(nbnd, bndind, bndptr, i);

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      myedegrees = myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += xadj[i+1]-xadj[i];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == other) {
              myedegrees[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].pid = other;
            myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
          }
        }
      }

      ASSERT(myrinfo->ndegrees <= xadj[i+1]-xadj[i]);
    }
  }

  graph->mincut = mincut/2;
  graph->nbnd = nbnd;

}



/*************************************************************************
* This function projects a partition, and at the same time computes the
* parameters for refinement.
**************************************************************************/
void MocProjectKWayPartition(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, j, k, nvtxs, nbnd, me, other, istart, iend, ndegrees;
  idxtype *xadj, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *cmap, *where, *bndptr, *bndind;
  idxtype *cwhere;
  GraphType *cgraph;
  RInfoType *crinfo, *rinfo, *myrinfo;
  EDegreeType *myedegrees;
  idxtype *htable;

  cgraph = graph->coarser;
  cwhere = cgraph->where;
  crinfo = cgraph->rinfo;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  MocAllocateKWayPartitionMemory(ctrl, graph, nparts);
  where = graph->where;
  rinfo = graph->rinfo;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);

  /* Go through and project partition and compute id/ed for the nodes */
  for (i=0; i<nvtxs; i++) {
    k = cmap[i];
    where[i] = cwhere[k];
    cmap[i] = crinfo[k].ed;  /* For optimization */
  }

  htable = idxset(nparts, -1, idxwspacemalloc(ctrl, nparts));

  ctrl->wspace.cdegree = 0;
  for (nbnd=0, i=0; i<nvtxs; i++) {
    me = where[i];

    myrinfo = rinfo+i;
    myrinfo->id = myrinfo->ed = myrinfo->ndegrees = 0;
    myrinfo->edegrees = NULL;

    myrinfo->id = adjwgtsum[i];

    if (cmap[i] > 0) { /* If it is an interface node. Note cmap[i] = crinfo[cmap[i]].ed */
      istart = xadj[i];
      iend = xadj[i+1];

      myedegrees = myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += iend-istart;

      ndegrees = 0;
      for (j=istart; j<iend; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          myrinfo->ed += adjwgt[j];
          if ((k = htable[other]) == -1) {
            htable[other] = ndegrees;
            myedegrees[ndegrees].pid = other;
            myedegrees[ndegrees++].ed = adjwgt[j];
          }
          else {
            myedegrees[k].ed += adjwgt[j];
          }
        }
      }
      myrinfo->id -= myrinfo->ed;

      /* Remove space for edegrees if it was interior */
      if (myrinfo->ed == 0) {
        myrinfo->edegrees = NULL;
        ctrl->wspace.cdegree -= iend-istart;
      }
      else {
        if (myrinfo->ed-myrinfo->id >= 0)
          BNDInsert(nbnd, bndind, bndptr, i);

        myrinfo->ndegrees = ndegrees;

        for (j=0; j<ndegrees; j++)
          htable[myedegrees[j].pid] = -1;
      }
    }
  }

  scopy(graph->ncon*nparts, cgraph->npwgts, graph->npwgts);
  graph->mincut = cgraph->mincut;
  graph->nbnd = nbnd;

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

  idxwspacefree(ctrl, nparts);

  ASSERT(CheckBnd2(graph));

}



/*************************************************************************
* This function computes the boundary definition for balancing
**************************************************************************/
void MocComputeKWayBalanceBoundary(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int i, nvtxs, nbnd;
  idxtype *bndind, *bndptr;

  nvtxs = graph->nvtxs;
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);


  /* Compute the new boundary */
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    if (graph->rinfo[i].ed > 0)
      BNDInsert(nbnd, bndind, bndptr, i);
  }

  graph->nbnd = nbnd;
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mmatch.c
 *
 * This file contains the code that computes matchings and creates the next
 * level coarse graph.
 *
 * Started 7/23/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void MCMatch_RM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, ncon, cnvtxs, maxidx;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;
  float *nvwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (match[k] == UNMATCHED && AreAllVwgtsBelowFast(ncon, nvwgt+i*ncon, nvwgt+k*ncon, ctrl->nmaxvwgt)) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void MCMatch_HEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, ncon, maxidx, maxwgt;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;
  float *nvwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (match[k] == UNMATCHED && maxwgt <= adjwgt[j] &&
               AreAllVwgtsBelowFast(ncon, nvwgt+i*ncon, nvwgt+k*ncon, ctrl->nmaxvwgt)) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void MCMatch_SHEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, ncon, maxidx, maxwgt, avgdegree;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *match, *cmap, *degrees, *perm, *tperm;
  float *nvwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  tperm = idxwspacemalloc(ctrl, nvtxs);
  degrees = idxwspacemalloc(ctrl, nvtxs);

  RandomPermute(nvtxs, tperm, 1);
  avgdegree = 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++)
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

  cnvtxs = 0;

  /* Take care any islands. Islands are matched with non-islands due to coarsening */
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      if (xadj[i] < xadj[i+1])
        break;

      maxidx = i;
      for (j=nvtxs-1; j>ii; j--) {
        k = perm[j];
        if (match[k] == UNMATCHED && xadj[k] < xadj[k+1]) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  /* Continue with normal matching */
  for (; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (match[k] == UNMATCHED && maxwgt <= adjwgt[j] &&
               AreAllVwgtsBelowFast(ncon, nvwgt+i*ncon, nvwgt+k*ncon, ctrl->nmaxvwgt)) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  idxwspacefree(ctrl, nvtxs);  /* degrees */
  idxwspacefree(ctrl, nvtxs);  /* tperm */

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void MCMatch_SHEBM(CtrlType *ctrl, GraphType *graph, int norm)
{
  int i, ii, j, k, nvtxs, cnvtxs, ncon, maxidx, maxwgt, avgdegree;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *match, *cmap, *degrees, *perm, *tperm;
  float *nvwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  tperm = idxwspacemalloc(ctrl, nvtxs);
  degrees = idxwspacemalloc(ctrl, nvtxs);

  RandomPermute(nvtxs, tperm, 1);
  avgdegree = 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++)
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

  cnvtxs = 0;

  /* Take care any islands. Islands are matched with non-islands due to coarsening */
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      if (xadj[i] < xadj[i+1])
        break;

      maxidx = i;
      for (j=nvtxs-1; j>ii; j--) {
        k = perm[j];
        if (match[k] == UNMATCHED && xadj[k] < xadj[k+1]) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  /* Continue with normal matching */
  for (; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = -1;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];

        if (match[k] == UNMATCHED &&
            AreAllVwgtsBelowFast(ncon, nvwgt+i*ncon, nvwgt+k*ncon, ctrl->nmaxvwgt) &&
            (maxwgt < adjwgt[j] ||
              (maxwgt == adjwgt[j] &&
               BetterVBalance(ncon, norm, nvwgt+i*ncon, nvwgt+maxidx*ncon, nvwgt+k*ncon) >= 0
              )
            )
           ) {
          maxwgt = adjwgt[j];
          maxidx = k;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  idxwspacefree(ctrl, nvtxs);  /* degrees */
  idxwspacefree(ctrl, nvtxs);  /* tperm */

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void MCMatch_SBHEM(CtrlType *ctrl, GraphType *graph, int norm)
{
  int i, ii, j, k, nvtxs, cnvtxs, ncon, maxidx, maxwgt, avgdegree;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *match, *cmap, *degrees, *perm, *tperm;
  float *nvwgt, vbal;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  tperm = idxwspacemalloc(ctrl, nvtxs);
  degrees = idxwspacemalloc(ctrl, nvtxs);

  RandomPermute(nvtxs, tperm, 1);
  avgdegree = 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++)
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

  cnvtxs = 0;

  /* Take care any islands. Islands are matched with non-islands due to coarsening */
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      if (xadj[i] < xadj[i+1])
        break;

      maxidx = i;
      for (j=nvtxs-1; j>ii; j--) {
        k = perm[j];
        if (match[k] == UNMATCHED && xadj[k] < xadj[k+1]) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  /* Continue with normal matching */
  for (; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = -1;
      vbal = 0.0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (match[k] == UNMATCHED && AreAllVwgtsBelowFast(ncon, nvwgt+i*ncon, nvwgt+k*ncon, ctrl->nmaxvwgt)) {
          if (maxidx != i)
            vbal = BetterVBalance(ncon, norm, nvwgt+i*ncon, nvwgt+maxidx*ncon, nvwgt+k*ncon);

          if (vbal > 0 || (vbal > -.01 && maxwgt < adjwgt[j])) {
            maxwgt = adjwgt[j];
            maxidx = k;
          }
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  idxwspacefree(ctrl, nvtxs);  /* degrees */
  idxwspacefree(ctrl, nvtxs);  /* tperm */

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}





/*************************************************************************
* This function checks if v+u2 provides a better balance in the weight
* vector that v+u1
**************************************************************************/
float BetterVBalance(int ncon, int norm, float *vwgt, float *u1wgt, float *u2wgt)
{
  int i;
  float sum1, sum2, max1, max2, min1, min2, diff1, diff2;

  if (norm == -1) {
    max1 = min1 = vwgt[0]+u1wgt[0];
    max2 = min2 = vwgt[0]+u2wgt[0];
    sum1 = vwgt[0]+u1wgt[0];
    sum2 = vwgt[0]+u2wgt[0];

    for (i=1; i<ncon; i++) {
      if (max1 < vwgt[i]+u1wgt[i])
        max1 = vwgt[i]+u1wgt[i];
      if (min1 > vwgt[i]+u1wgt[i])
        min1 = vwgt[i]+u1wgt[i];

      if (max2 < vwgt[i]+u2wgt[i])
        max2 = vwgt[i]+u2wgt[i];
      if (min2 > vwgt[i]+u2wgt[i])
        min2 = vwgt[i]+u2wgt[i];

      sum1 += vwgt[i]+u1wgt[i];
      sum2 += vwgt[i]+u2wgt[i];
    }

    if (sum1 == 0.0)
      return 1;
    else if (sum2 == 0.0)
      return -1;
    else
      return ((max1-min1)/sum1) - ((max2-min2)/sum2);
  }
  else if (norm == 1) {
    sum1 = sum2 = 0.0;
    for (i=0; i<ncon; i++) {
      sum1 += vwgt[i]+u1wgt[i];
      sum2 += vwgt[i]+u2wgt[i];
    }
    sum1 = sum1/(1.0*ncon);
    sum2 = sum2/(1.0*ncon);

    diff1 = diff2 = 0.0;
    for (i=0; i<ncon; i++) {
      diff1 += fabs(sum1 - (vwgt[i]+u1wgt[i]));
      diff2 += fabs(sum2 - (vwgt[i]+u2wgt[i]));
    }

    return diff1 - diff2;
  }
  else {
    errexit("Unknown norm: %d\n", norm);
  }
  return 0.0;
}


/*************************************************************************
* This function checks if the vertex weights of two vertices are below
* a given set of values
**************************************************************************/
int AreAllVwgtsBelowFast(int ncon, float *vwgt1, float *vwgt2, float limit)
{
  int i;

  for (i=0; i<ncon; i++)
    if (vwgt1[i] + vwgt2[i] > limit)
      return 0;

  return 1;
}

/*
 * mmd.c
 *
 * **************************************************************
 * The following C function was developed from a FORTRAN subroutine
 * in SPARSPAK written by Eleanor Chu, Alan George, Joseph Liu
 * and Esmond Ng.
 *
 * The FORTRAN-to-C transformation and modifications such as dynamic
 * memory allocation and deallocation were performed by Chunguang
 * Sun.
 * **************************************************************
 *
 * Taken from SMMS, George 12/13/94
 *
 * The meaning of invperm, and perm vectors is different from that
 * in genqmd_ of SparsPak
 *
 * $Id$
 */




/*************************************************************************
*  genmmd  -- multiple minimum external degree
*  purpose -- this routine implements the minimum degree
*     algorithm. it makes use of the implicit representation
*     of elimination graphs by quotient graphs, and the notion
*     of indistinguishable nodes. It also implements the modifications
*     by multiple elimination and minimum external degree.
*     Caution -- the adjacency vector adjncy will be destroyed.
*  Input parameters --
*     neqns -- number of equations.
*     (xadj, adjncy) -- the adjacency structure.
*     delta  -- tolerance value for multiple elimination.
*     maxint -- maximum machine representable (short) integer
*               (any smaller estimate will do) for marking nodes.
*  Output parameters --
*     perm -- the minimum degree ordering.
*     invp -- the inverse of perm.
*     *ncsub -- an upper bound on the number of nonzero subscripts
*               for the compressed storage scheme.
*  Working parameters --
*     head -- vector for head of degree lists.
*     invp  -- used temporarily for degree forward link.
*     perm  -- used temporarily for degree backward link.
*     qsize -- vector for size of supernodes.
*     list -- vector for temporary linked lists.
*     marker -- a temporary marker vector.
*  Subroutines used -- mmdelm, mmdint, mmdnum, mmdupd.
**************************************************************************/
void genmmd(int neqns, idxtype *xadj, idxtype *adjncy, idxtype *invp, idxtype *perm,
     int delta, idxtype *head, idxtype *qsize, idxtype *list, idxtype *marker,
     int maxint, int *ncsub)
{
    int  ehead, i, mdeg, mdlmt, mdeg_node, nextmd, num, tag;

    if (neqns <= 0)
      return;

    /* Adjust from C to Fortran */
    xadj--; adjncy--; invp--; perm--; head--; qsize--; list--; marker--;

    /* initialization for the minimum degree algorithm. */
    *ncsub = 0;
    mmdint(neqns, xadj, adjncy, head, invp, perm, qsize, list, marker);

    /*  'num' counts the number of ordered nodes plus 1. */
    num = 1;

    /* eliminate all isolated nodes. */
    nextmd = head[1];
    while (nextmd > 0) {
      mdeg_node = nextmd;
      nextmd = invp[mdeg_node];
      marker[mdeg_node] = maxint;
      invp[mdeg_node] = -num;
      num = num + 1;
    }

    /* search for node of the minimum degree. 'mdeg' is the current */
    /* minimum degree; 'tag' is used to facilitate marking nodes.   */
    if (num > neqns)
      goto n1000;
    tag = 1;
    head[1] = 0;
    mdeg = 2;

    /* infinite loop here ! */
    while (1) {
      while (head[mdeg] <= 0)
        mdeg++;

      /* use value of 'delta' to set up 'mdlmt', which governs */
      /* when a degree update is to be performed.              */
      mdlmt = mdeg + delta;
      ehead = 0;

n500:
      mdeg_node = head[mdeg];
      while (mdeg_node <= 0) {
        mdeg++;

        if (mdeg > mdlmt)
          goto n900;
        mdeg_node = head[mdeg];
      };

      /*  remove 'mdeg_node' from the degree structure. */
      nextmd = invp[mdeg_node];
      head[mdeg] = nextmd;
      if (nextmd > 0)
        perm[nextmd] = -mdeg;
      invp[mdeg_node] = -num;
      *ncsub += mdeg + qsize[mdeg_node] - 2;
      if ((num+qsize[mdeg_node]) > neqns)
        goto n1000;

      /*  eliminate 'mdeg_node' and perform quotient graph */
      /*  transformation. reset 'tag' value if necessary.    */
      tag++;
      if (tag >= maxint) {
        tag = 1;
        for (i = 1; i <= neqns; i++)
          if (marker[i] < maxint)
            marker[i] = 0;
      };

      mmdelm(mdeg_node, xadj, adjncy, head, invp, perm, qsize, list, marker, maxint, tag);

      num += qsize[mdeg_node];
      list[mdeg_node] = ehead;
      ehead = mdeg_node;
      if (delta >= 0)
        goto n500;

 n900:
      /* update degrees of the nodes involved in the  */
      /* minimum degree nodes elimination.            */
      if (num > neqns)
        goto n1000;
      mmdupd( ehead, neqns, xadj, adjncy, delta, &mdeg, head, invp, perm, qsize, list, marker, maxint, &tag);
    }; /* end of -- while ( 1 ) -- */

n1000:
    mmdnum( neqns, perm, invp, qsize );

    /* Adjust from Fortran back to C*/
    xadj++; adjncy++; invp++; perm++; head++; qsize++; list++; marker++;
}


/**************************************************************************
*           mmdelm ...... multiple minimum degree elimination
* Purpose -- This routine eliminates the node mdeg_node of minimum degree
*     from the adjacency structure, which is stored in the quotient
*     graph format. It also transforms the quotient graph representation
*     of the elimination graph.
* Input parameters --
*     mdeg_node -- node of minimum degree.
*     maxint -- estimate of maximum representable (short) integer.
*     tag    -- tag value.
* Updated parameters --
*     (xadj, adjncy) -- updated adjacency structure.
*     (head, forward, backward) -- degree doubly linked structure.
*     qsize -- size of supernode.
*     marker -- marker vector.
*     list -- temporary linked list of eliminated nabors.
***************************************************************************/
void mmdelm(int mdeg_node, idxtype *xadj, idxtype *adjncy, idxtype *head, idxtype *forward,
     idxtype *backward, idxtype *qsize, idxtype *list, idxtype *marker, int maxint,int tag)
{
    int   element, i,   istop, istart, j,
          jstop, jstart, link,
          nabor, node, npv, nqnbrs, nxnode,
          pvnode, rlmt, rloc, rnode, xqnbr;

    /* find the reachable set of 'mdeg_node' and */
    /* place it in the data structure.           */
    marker[mdeg_node] = tag;
    istart = xadj[mdeg_node];
    istop = xadj[mdeg_node+1] - 1;

    /* 'element' points to the beginning of the list of  */
    /* eliminated nabors of 'mdeg_node', and 'rloc' gives the */
    /* storage location for the next reachable node.   */
    element = 0;
    rloc = istart;
    rlmt = istop;
    for ( i = istart; i <= istop; i++ ) {
        nabor = adjncy[i];
        if ( nabor == 0 ) break;
        if ( marker[nabor] < tag ) {
           marker[nabor] = tag;
           if ( forward[nabor] < 0 )  {
              list[nabor] = element;
              element = nabor;
           } else {
              adjncy[rloc] = nabor;
              rloc++;
           };
        }; /* end of -- if -- */
    }; /* end of -- for -- */

  /* merge with reachable nodes from generalized elements. */
  while ( element > 0 ) {
      adjncy[rlmt] = -element;
      link = element;

n400:
      jstart = xadj[link];
      jstop = xadj[link+1] - 1;
      for ( j = jstart; j <= jstop; j++ ) {
          node = adjncy[j];
          link = -node;
          if ( node < 0 )  goto n400;
          if ( node == 0 ) break;
          if ((marker[node]<tag)&&(forward[node]>=0)) {
             marker[node] = tag;
             /*use storage from eliminated nodes if necessary.*/
             while ( rloc >= rlmt ) {
                   link = -adjncy[rlmt];
                   rloc = xadj[link];
                   rlmt = xadj[link+1] - 1;
             };
             adjncy[rloc] = node;
             rloc++;
          };
      }; /* end of -- for ( j = jstart; -- */
      element = list[element];
    };  /* end of -- while ( element > 0 ) -- */
    if ( rloc <= rlmt ) adjncy[rloc] = 0;
    /* for each node in the reachable set, do the following. */
    link = mdeg_node;

n1100:
    istart = xadj[link];
    istop = xadj[link+1] - 1;
    for ( i = istart; i <= istop; i++ ) {
        rnode = adjncy[i];
        link = -rnode;
        if ( rnode < 0 ) goto n1100;
        if ( rnode == 0 ) return;

        /* 'rnode' is in the degree list structure. */
        pvnode = backward[rnode];
        if (( pvnode != 0 ) && ( pvnode != (-maxint) )) {
           /* then remove 'rnode' from the structure. */
           nxnode = forward[rnode];
           if ( nxnode > 0 ) backward[nxnode] = pvnode;
           if ( pvnode > 0 ) forward[pvnode] = nxnode;
           npv = -pvnode;
           if ( pvnode < 0 ) head[npv] = nxnode;
        };

        /* purge inactive quotient nabors of 'rnode'. */
        jstart = xadj[rnode];
        jstop = xadj[rnode+1] - 1;
        xqnbr = jstart;
        for ( j = jstart; j <= jstop; j++ ) {
            nabor = adjncy[j];
            if ( nabor == 0 ) break;
            if ( marker[nabor] < tag ) {
                adjncy[xqnbr] = nabor;
                xqnbr++;
            };
        };

        /* no active nabor after the purging. */
        nqnbrs = xqnbr - jstart;
        if ( nqnbrs <= 0 ) {
           /* merge 'rnode' with 'mdeg_node'. */
           qsize[mdeg_node] += qsize[rnode];
           qsize[rnode] = 0;
           marker[rnode] = maxint;
           forward[rnode] = -mdeg_node;
           backward[rnode] = -maxint;
        } else {
           /* flag 'rnode' for degree update, and  */
           /* add 'mdeg_node' as a nabor of 'rnode'.      */
           forward[rnode] = nqnbrs + 1;
           backward[rnode] = 0;
           adjncy[xqnbr] = mdeg_node;
           xqnbr++;
           if ( xqnbr <= jstop )  adjncy[xqnbr] = 0;
        };
      }; /* end of -- for ( i = istart; -- */
      return;
 }

/***************************************************************************
*    mmdint ---- mult minimum degree initialization
*    purpose -- this routine performs initialization for the
*       multiple elimination version of the minimum degree algorithm.
*    input parameters --
*       neqns  -- number of equations.
*       (xadj, adjncy) -- adjacency structure.
*    output parameters --
*       (head, dfrow, backward) -- degree doubly linked structure.
*       qsize -- size of supernode ( initialized to one).
*       list -- linked list.
*       marker -- marker vector.
****************************************************************************/
int  mmdint(int neqns, idxtype *xadj, idxtype *adjncy, idxtype *head, idxtype *forward,
     idxtype *backward, idxtype *qsize, idxtype *list, idxtype *marker)
{
    int  fnode, ndeg, node;

    for ( node = 1; node <= neqns; node++ ) {
        head[node] = 0;
        qsize[node] = 1;
        marker[node] = 0;
        list[node] = 0;
    };

    /* initialize the degree doubly linked lists. */
    for ( node = 1; node <= neqns; node++ ) {
        ndeg = xadj[node+1] - xadj[node]/* + 1*/;   /* george */
        if (ndeg == 0)
          ndeg = 1;
        fnode = head[ndeg];
        forward[node] = fnode;
        head[ndeg] = node;
        if ( fnode > 0 ) backward[fnode] = node;
        backward[node] = -ndeg;
    };
    return 0;
}

/****************************************************************************
* mmdnum --- multi minimum degree numbering
* purpose -- this routine performs the final step in producing
*    the permutation and inverse permutation vectors in the
*    multiple elimination version of the minimum degree
*    ordering algorithm.
* input parameters --
*     neqns -- number of equations.
*     qsize -- size of supernodes at elimination.
* updated parameters --
*     invp -- inverse permutation vector. on input,
*             if qsize[node] = 0, then node has been merged
*             into the node -invp[node]; otherwise,
*            -invp[node] is its inverse labelling.
* output parameters --
*     perm -- the permutation vector.
****************************************************************************/
void mmdnum(int neqns, idxtype *perm, idxtype *invp, idxtype *qsize)
{
  int father, nextf, node, nqsize, num, root;

  for ( node = 1; node <= neqns; node++ ) {
      nqsize = qsize[node];
      if ( nqsize <= 0 ) perm[node] = invp[node];
      if ( nqsize > 0 )  perm[node] = -invp[node];
  };

  /* for each node which has been merged, do the following. */
  for ( node = 1; node <= neqns; node++ ) {
      if ( perm[node] <= 0 )  {

	 /* trace the merged tree until one which has not */
         /* been merged, call it root.                    */
         father = node;
         while ( perm[father] <= 0 )
            father = - perm[father];

         /* number node after root. */
         root = father;
         num = perm[root] + 1;
         invp[node] = -num;
         perm[root] = num;

         /* shorten the merged tree. */
         father = node;
         nextf = - perm[father];
         while ( nextf > 0 ) {
            perm[father] = -root;
            father = nextf;
            nextf = -perm[father];
         };
      };  /* end of -- if ( perm[node] <= 0 ) -- */
  }; /* end of -- for ( node = 1; -- */

  /* ready to compute perm. */
  for ( node = 1; node <= neqns; node++ ) {
        num = -invp[node];
        invp[node] = num;
        perm[num] = node;
  };
  return;
}

/****************************************************************************
* mmdupd ---- multiple minimum degree update
* purpose -- this routine updates the degrees of nodes after a
*            multiple elimination step.
* input parameters --
*    ehead -- the beginning of the list of eliminated nodes
*             (i.e., newly formed elements).
*    neqns -- number of equations.
*    (xadj, adjncy) -- adjacency structure.
*    delta -- tolerance value for multiple elimination.
*    maxint -- maximum machine representable (short) integer.
* updated parameters --
*    mdeg -- new minimum degree after degree update.
*    (head, forward, backward) -- degree doubly linked structure.
*    qsize -- size of supernode.
*    list -- marker vector for degree update.
*    *tag   -- tag value.
****************************************************************************/
void mmdupd(int ehead, int neqns, idxtype *xadj, idxtype *adjncy, int delta, int *mdeg,
     idxtype *head, idxtype *forward, idxtype *backward, idxtype *qsize, idxtype *list,
     idxtype *marker, int maxint,int *tag)
{
 int  deg, deg0, element, enode, fnode, i, iq2, istop,
      istart, j, jstop, jstart, link, mdeg0, mtag, nabor,
      node, q2head, qxhead;

      mdeg0 = *mdeg + delta;
      element = ehead;

n100:
      if ( element <= 0 ) return;

      /* for each of the newly formed element, do the following. */
      /* reset tag value if necessary.                           */
      mtag = *tag + mdeg0;
      if ( mtag >= maxint ) {
         *tag = 1;
         for ( i = 1; i <= neqns; i++ )
             if ( marker[i] < maxint ) marker[i] = 0;
         mtag = *tag + mdeg0;
      };

      /* create two linked lists from nodes associated with 'element': */
      /* one with two nabors (q2head) in the adjacency structure, and the*/
      /* other with more than two nabors (qxhead). also compute 'deg0',*/
      /* number of nodes in this element.                              */
      q2head = 0;
      qxhead = 0;
      deg0 = 0;
      link =element;

n400:
      istart = xadj[link];
      istop = xadj[link+1] - 1;
      for ( i = istart; i <= istop; i++ ) {
          enode = adjncy[i];
          link = -enode;
          if ( enode < 0 )  goto n400;
          if ( enode == 0 ) break;
          if ( qsize[enode] != 0 ) {
             deg0 += qsize[enode];
             marker[enode] = mtag;

             /*'enode' requires a degree update*/
             if ( backward[enode] == 0 ) {
                /* place either in qxhead or q2head list. */
                if ( forward[enode] != 2 ) {
                     list[enode] = qxhead;
                     qxhead = enode;
                } else {
                     list[enode] = q2head;
                     q2head = enode;
                };
             };
          }; /* enf of -- if ( qsize[enode] != 0 ) -- */
      }; /* end of -- for ( i = istart; -- */

      /* for each node in q2 list, do the following. */
      enode = q2head;
      iq2 = 1;

n900:
      if ( enode <= 0 ) goto n1500;
      if ( backward[enode] != 0 ) goto n2200;
      (*tag)++;
      deg = deg0;

      /* identify the other adjacent element nabor. */
      istart = xadj[enode];
      nabor = adjncy[istart];
      if ( nabor == element ) nabor = adjncy[istart+1];
      link = nabor;
      if ( forward[nabor] >= 0 ) {
           /* nabor is uneliminated, increase degree count. */
           deg += qsize[nabor];
           goto n2100;
      };

       /* the nabor is eliminated. for each node in the 2nd element */
       /* do the following.                                         */
n1000:
       istart = xadj[link];
       istop = xadj[link+1] - 1;
       for ( i = istart; i <= istop; i++ ) {
           node = adjncy[i];
           link = -node;
           if ( node != enode ) {
                if ( node < 0 ) goto n1000;
                if ( node == 0 )  goto n2100;
                if ( qsize[node] != 0 ) {
                     if ( marker[node] < *tag ) {
                        /* 'node' is not yet considered. */
                        marker[node] = *tag;
                        deg += qsize[node];
                     } else {
                        if ( backward[node] == 0 ) {
                             if ( forward[node] == 2 ) {
                                /* 'node' is indistinguishable from 'enode'.*/
                                /* merge them into a new supernode.         */
                                qsize[enode] += qsize[node];
                                qsize[node] = 0;
                                marker[node] = maxint;
                                forward[node] = -enode;
                                backward[node] = -maxint;
                             } else {
                                /* 'node' is outmacthed by 'enode' */
				if (backward[node]==0) backward[node] = -maxint;
                             };
                        }; /* end of -- if ( backward[node] == 0 ) -- */
                    }; /* end of -- if ( marker[node] < *tag ) -- */
                }; /* end of -- if ( qsize[node] != 0 ) -- */
              }; /* end of -- if ( node != enode ) -- */
          }; /* end of -- for ( i = istart; -- */
          goto n2100;

n1500:
          /* for each 'enode' in the 'qx' list, do the following. */
          enode = qxhead;
          iq2 = 0;

n1600:    if ( enode <= 0 )  goto n2300;
          if ( backward[enode] != 0 )  goto n2200;
          (*tag)++;
          deg = deg0;

          /*for each unmarked nabor of 'enode', do the following.*/
          istart = xadj[enode];
          istop = xadj[enode+1] - 1;
          for ( i = istart; i <= istop; i++ ) {
                nabor = adjncy[i];
                if ( nabor == 0 ) break;
                if ( marker[nabor] < *tag ) {
                     marker[nabor] = *tag;
                     link = nabor;
                     if ( forward[nabor] >= 0 )
                          /*if uneliminated, include it in deg count.*/
                          deg += qsize[nabor];
                     else {
n1700:
                          /* if eliminated, include unmarked nodes in this*/
                          /* element into the degree count.             */
                          jstart = xadj[link];
                          jstop = xadj[link+1] - 1;
                          for ( j = jstart; j <= jstop; j++ ) {
                                node = adjncy[j];
                                link = -node;
                                if ( node < 0 ) goto n1700;
                                if ( node == 0 ) break;
                                if ( marker[node] < *tag ) {
                                    marker[node] = *tag;
                                    deg += qsize[node];
                                };
                          }; /* end of -- for ( j = jstart; -- */
                     }; /* end of -- if ( forward[nabor] >= 0 ) -- */
                  }; /* end of -- if ( marker[nabor] < *tag ) -- */
          }; /* end of -- for ( i = istart; -- */

n2100:
          /* update external degree of 'enode' in degree structure, */
          /* and '*mdeg' if necessary.                     */
          deg = deg - qsize[enode] + 1;
          fnode = head[deg];
          forward[enode] = fnode;
          backward[enode] = -deg;
          if ( fnode > 0 ) backward[fnode] = enode;
          head[deg] = enode;
          if ( deg < *mdeg ) *mdeg = deg;

n2200:
          /* get next enode in current element. */
          enode = list[enode];
          if ( iq2 == 1 ) goto n900;
          goto n1600;

n2300:
          /* get next element in the list. */
          *tag = mtag;
          element = list[element];
          goto n100;
    }
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mpmetis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 *
 */





/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_mCPartGraphRecursive(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
       idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
       int *options, int *edgecut, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_PMETIS, *nvtxs, *ncon, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType  = McPMETIS_CTYPE;
    ctrl.IType  = McPMETIS_ITYPE;
    ctrl.RType  = McPMETIS_RTYPE;
    ctrl.dbglvl = McPMETIS_DBGLVL;
  }
  else {
    ctrl.CType  = options[OPTION_CTYPE];
    ctrl.IType  = options[OPTION_ITYPE];
    ctrl.RType  = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 100;

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MCMlevelRecursiveBisection(&ctrl, &graph, *nparts, part, 1.000, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}



/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_mCHPartGraphRecursive(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
       idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
       float *ubvec, int *options, int *edgecut, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;
  float *myubvec;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_PMETIS, *nvtxs, *ncon, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = PMETIS_CTYPE;
    ctrl.IType = PMETIS_ITYPE;
    ctrl.RType = PMETIS_RTYPE;
    ctrl.dbglvl = PMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 100;

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  myubvec = fmalloc(*ncon, "PWMETIS: mytpwgts");
  scopy(*ncon, ubvec, myubvec);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MCHMlevelRecursiveBisection(&ctrl, &graph, *nparts, part, myubvec, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);
  GKfree((void **) &myubvec, LTERM);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}



/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_mCPartGraphRecursiveInternal(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
       float *nvwgt, idxtype *adjwgt, int *nparts, int *options, int *edgecut, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph2(&graph, *nvtxs, *ncon, xadj, adjncy, nvwgt, adjwgt);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = PMETIS_CTYPE;
    ctrl.IType = PMETIS_ITYPE;
    ctrl.RType = PMETIS_RTYPE;
    ctrl.dbglvl = PMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 100;

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MCMlevelRecursiveBisection(&ctrl, &graph, *nparts, part, 1.000, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

}


/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_mCHPartGraphRecursiveInternal(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
       float *nvwgt, idxtype *adjwgt, int *nparts, float *ubvec, int *options, int *edgecut,
       idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;
  float *myubvec;

  SetUpGraph2(&graph, *nvtxs, *ncon, xadj, adjncy, nvwgt, adjwgt);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = PMETIS_CTYPE;
    ctrl.IType = PMETIS_ITYPE;
    ctrl.RType = PMETIS_RTYPE;
    ctrl.dbglvl = PMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 100;

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  myubvec = fmalloc(*ncon, "PWMETIS: mytpwgts");
  scopy(*ncon, ubvec, myubvec);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MCHMlevelRecursiveBisection(&ctrl, &graph, *nparts, part, myubvec, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);
  GKfree((void **) &myubvec, LTERM);

}




/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MCMlevelRecursiveBisection(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part,
       float ubfactor, int fpart)
{
  int i, nvtxs, cut;
  idxtype *label, *where;
  GraphType lgraph, rgraph;
  float tpwgts[2];

  nvtxs = graph->nvtxs;
  if (nvtxs == 0) {
    printf("\t***Cannot bisect a graph with 0 vertices!\n\t***You are trying to partition a graph into too many parts!\n");
    return 0;
  }

  /* Determine the weights of the partitions */
  tpwgts[0] = 1.0*(nparts>>1)/(1.0*nparts);
  tpwgts[1] = 1.0 - tpwgts[0];

  MCMlevelEdgeBisection(ctrl, graph, tpwgts, ubfactor);
  cut = graph->mincut;

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

  if (nparts > 2)
    SplitGraphPart(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  GKfree((void **) &graph->gdata,  (void **) &graph->nvwgt, (void **) &graph->rdata,
         (void **) &graph->npwgts, (void **) &graph->label, LTERM);


  /* Do the recursive call */
  if (nparts > 3) {
    cut += MCMlevelRecursiveBisection(ctrl, &lgraph, nparts/2, part, ubfactor, fpart);
    cut += MCMlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, ubfactor, fpart+nparts/2);
  }
  else if (nparts == 3) {
    cut += MCMlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, ubfactor, fpart+nparts/2);
    GKfree((void **) &lgraph.gdata, (void **) &lgraph.nvwgt, (void **) &lgraph.label, LTERM);
  }

  return cut;

}



/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MCHMlevelRecursiveBisection(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part,
      float *ubvec, int fpart)
{
  int i, nvtxs, ncon, cut;
  idxtype *label, *where;
  GraphType lgraph, rgraph;
  float tpwgts[2], *npwgts, *lubvec, *rubvec;

  lubvec = rubvec = NULL;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  if (nvtxs == 0) {
    printf("\t***Cannot bisect a graph with 0 vertices!\n\t***You are trying to partition a graph into too many parts!\n");
    return 0;
  }

  /* Determine the weights of the partitions */
  tpwgts[0] = 1.0*(nparts>>1)/(1.0*nparts);
  tpwgts[1] = 1.0 - tpwgts[0];

  /* For now, relax at the coarsest level only */
  if (nparts == 2)
    MCHMlevelEdgeBisection(ctrl, graph, tpwgts, ubvec);
  else
    MCMlevelEdgeBisection(ctrl, graph, tpwgts, 1.000);
  cut = graph->mincut;

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

  if (nparts > 2) {
    /* Adjust the ubvecs before the split */
    npwgts = graph->npwgts;
    lubvec = fmalloc(ncon, "MCHMlevelRecursiveBisection");
    rubvec = fmalloc(ncon, "MCHMlevelRecursiveBisection");

    for (i=0; i<ncon; i++) {
      lubvec[i] = ubvec[i]*tpwgts[0]/npwgts[i];
      lubvec[i] = amax(lubvec[i], 1.01);

      rubvec[i] = ubvec[i]*tpwgts[1]/npwgts[ncon+i];
      rubvec[i] = amax(rubvec[i], 1.01);
    }

    SplitGraphPart(ctrl, graph, &lgraph, &rgraph);
  }

  /* Free the memory of the top level graph */
  GKfree((void **) &graph->gdata,  (void **) &graph->nvwgt, (void **) &graph->rdata,
         (void **) &graph->npwgts, (void **) &graph->label, LTERM);


  /* Do the recursive call */
  if (nparts > 3) {
    cut += MCHMlevelRecursiveBisection(ctrl, &lgraph, nparts/2, part, lubvec, fpart);
    cut += MCHMlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, rubvec, fpart+nparts/2);
  }
  else if (nparts == 3) {
    cut += MCHMlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, rubvec, fpart+nparts/2);
    GKfree((void **) &lgraph.gdata, (void **) &lgraph.nvwgt, (void **) &lgraph.label, LTERM);
  }

  GKfree((void **) &lubvec, (void **) &rubvec, LTERM);

  return cut;

}




/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void MCMlevelEdgeBisection(CtrlType *ctrl, GraphType *graph, float *tpwgts, float ubfactor)
{
  GraphType *cgraph;

  cgraph = MCCoarsen2Way(ctrl, graph);

  MocInit2WayPartition(ctrl, cgraph, tpwgts, ubfactor);

  MocRefine2Way(ctrl, graph, cgraph, tpwgts, ubfactor);

}



/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void MCHMlevelEdgeBisection(CtrlType *ctrl, GraphType *graph, float *tpwgts, float *ubvec)
{
  GraphType *cgraph;

  cgraph = MCCoarsen2Way(ctrl, graph);

  MocInit2WayPartition2(ctrl, cgraph, tpwgts, ubvec);

  MocRefine2Way2(ctrl, graph, cgraph, tpwgts, ubvec);

}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mrefine2.c
 *
 * This file contains the driving routines for multilevel refinement
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void MocRefine2Way2(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, float *tpwgts,
       float *ubvec)
{

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  MocCompute2WayPartitionParams(ctrl, graph);

  for (;;) {
    ASSERT(CheckBnd(graph));

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));
    switch (ctrl->RType) {
      case RTYPE_FM:
        MocBalance2Way2(ctrl, graph, tpwgts, ubvec);
        MocFM_2WayEdgeRefine2(ctrl, graph, tpwgts, ubvec, 8);
        break;
      default:
        errexit("Unknown refinement type: %d\n", ctrl->RType);
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    MocProject2WayPartition(ctrl, graph);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * refine.c
 *
 * This file contains the driving routines for multilevel refinement
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void MocRefine2Way(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, float *tpwgts, float ubfactor)
{
  int i;
  float tubvec[MAXNCON];

  for (i=0; i<graph->ncon; i++)
    tubvec[i] = 1.0;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  MocCompute2WayPartitionParams(ctrl, graph);

  for (;;) {
    ASSERT(CheckBnd(graph));

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));
    switch (ctrl->RType) {
      case RTYPE_FM:
        MocBalance2Way(ctrl, graph, tpwgts, 1.03);
        MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 8);
        break;
      case 2:
        MocBalance2Way(ctrl, graph, tpwgts, 1.03);
        MocFM_2WayEdgeRefine2(ctrl, graph, tpwgts, tubvec, 8);
        break;
      default:
        errexit("Unknown refinement type: %d\n", ctrl->RType);
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    MocProject2WayPartition(ctrl, graph);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  MocBalance2Way(ctrl, graph, tpwgts, 1.01);
  MocFM_2WayEdgeRefine(ctrl, graph, tpwgts, 8);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}


/*************************************************************************
* This function allocates memory for 2-way edge refinement
**************************************************************************/
void MocAllocate2WayPartitionMemory(CtrlType *ctrl, GraphType *graph)
{
  int nvtxs, ncon;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;

  graph->rdata = idxmalloc(5*nvtxs, "Allocate2WayPartitionMemory: rdata");
  graph->where		= graph->rdata;
  graph->id		= graph->rdata + nvtxs;
  graph->ed		= graph->rdata + 2*nvtxs;
  graph->bndptr		= graph->rdata + 3*nvtxs;
  graph->bndind		= graph->rdata + 4*nvtxs;

  graph->npwgts 	= fmalloc(2*ncon, "npwgts");
}


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void MocCompute2WayPartitionParams(CtrlType *ctrl, GraphType *graph)
{
  int i, j, nvtxs, ncon, nbnd, mincut;
  idxtype *xadj, *adjncy, *adjwgt;
  float *nvwgt, *npwgts;
  idxtype *id, *ed, *where;
  idxtype *bndptr, *bndind;
  int me;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  npwgts = sset(2*ncon, 0.0, graph->npwgts);
  id = idxset(nvtxs, 0, graph->id);
  ed = idxset(nvtxs, 0, graph->ed);
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;


  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  nbnd = mincut = 0;
  for (i=0; i<nvtxs; i++) {
    ASSERT(where[i] >= 0 && where[i] <= 1);
    me = where[i];
    /*printf("saxpy: %d",i);*/
    saxpy(ncon, 1.0, nvwgt+i*ncon, 1, npwgts+me*ncon, 1);

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        id[i] += adjwgt[j];
      else
        ed[i] += adjwgt[j];
    }

    if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
      mincut += ed[i];
      bndptr[i] = nbnd;
      bndind[nbnd++] = i;
    }
  }

  graph->mincut = mincut/2;
  graph->nbnd = nbnd;

}



/*************************************************************************
* This function projects a partition, and at the same time computes the
* parameters for refinement.
**************************************************************************/
void MocProject2WayPartition(CtrlType *ctrl, GraphType *graph)
{
  int i, j, k, nvtxs, nbnd, me;
  idxtype *xadj, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *cmap, *where, *id, *ed, *bndptr, *bndind;
  idxtype *cwhere, *cid, *ced, *cbndptr;
  GraphType *cgraph;

  cgraph = graph->coarser;
  cwhere = cgraph->where;
  cid = cgraph->id;
  ced = cgraph->ed;
  cbndptr = cgraph->bndptr;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  MocAllocate2WayPartitionMemory(ctrl, graph);

  where = graph->where;
  id = idxset(nvtxs, 0, graph->id);
  ed = idxset(nvtxs, 0, graph->ed);
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;


  /* Go through and project partition and compute id/ed for the nodes */
  for (i=0; i<nvtxs; i++) {
    k = cmap[i];
    where[i] = cwhere[k];
    cmap[i] = cbndptr[k];
  }

  for (nbnd=0, i=0; i<nvtxs; i++) {
    me = where[i];

    id[i] = adjwgtsum[i];

    if (xadj[i] == xadj[i+1]) {
      bndptr[i] = nbnd;
      bndind[nbnd++] = i;
    }
    else {
      if (cmap[i] != -1) { /* If it is an interface node. Note that cmap[i] = cbndptr[cmap[i]] */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          if (me != where[adjncy[j]])
            ed[i] += adjwgt[j];
        }
        id[i] -= ed[i];

        if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
          bndptr[i] = nbnd;
          bndind[nbnd++] = i;
        }
      }
    }
  }

  graph->mincut = cgraph->mincut;
  graph->nbnd = nbnd;
  scopy(2*graph->ncon, cgraph->npwgts, graph->npwgts);

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

}

/*
 * mutil.c
 *
 * This file contains various utility functions for the MOC portion of the
 * code
 *
 * Started 2/15/98
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function checks if the vertex weights of two vertices are below
* a given set of values
**************************************************************************/
int AreAllVwgtsBelow(int ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float limit)
{
  int i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] > limit)
      return 0;

  return 1;
}


/*************************************************************************
* This function checks if the vertex weights of two vertices are below
* a given set of values
**************************************************************************/
int AreAnyVwgtsBelow(int ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float limit)
{
  int i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] < limit)
      return 1;

  return 0;
}



/*************************************************************************
* This function checks if the vertex weights of two vertices are above
* a given set of values
**************************************************************************/
int AreAllVwgtsAbove(int ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float limit)
{
  int i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] < limit)
      return 0;

  return 1;
}


/*************************************************************************
* This function computes the load imbalance over all the constrains
* For now assume that we just want balanced partitionings
**************************************************************************/
float ComputeLoadImbalance(int ncon, int nparts, float *npwgts, float *tpwgts)
{
  int i, j;
  float max, lb=0.0;

  for (i=0; i<ncon; i++) {
    max = 0.0;
    for (j=0; j<nparts; j++) {
      if (npwgts[j*ncon+i] > max)
        max = npwgts[j*ncon+i];
    }
    if (max*nparts > lb)
      lb = max*nparts;
  }

  return lb;
}

/*************************************************************************
* This function checks if the vertex weights of two vertices are below
* a given set of values
**************************************************************************/
int AreAllBelow(int ncon, float *v1, float *v2)
{
  int i;

  for (i=0; i<ncon; i++)
    if (v1[i] > v2[i])
      return 0;

  return 1;
}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * myqsort.c
 *
 * This file contains a fast idxtype increasing qsort algorithm.
 * Addopted from TeX
 *
 * Started 10/18/96
 * George
 *
 * $Id$
 */

			/* only for type declarations */

#define		THRESH		1	/* threshold for insertion */
#define		MTHRESH		6	/* threshold for median */




static void siqst(idxtype *, idxtype *);
static void iiqst(int *, int *);
static void keyiqst(KeyValueType *, KeyValueType *);
static void keyvaliqst(KeyValueType *, KeyValueType *);


/*************************************************************************
* Entry point of idxtype increasing sort
**************************************************************************/
void iidxsort(int n, idxtype *base)
{
  register idxtype *i;
  register idxtype *j;
  register idxtype *lo;
  register idxtype *hi;
  register idxtype *min;
  register idxtype c;
  idxtype *max;

  if (n <= 1)
    return;

  max = base + n;

  if (n >= THRESH) {
    siqst(base, max);
    hi = base + THRESH;
  }
  else
    hi = max;

  for (j = lo = base; lo++ < hi;) {
    if (*j > *lo)
      j = lo;
  }
  if (j != base) { /* swap j into place */
    c = *base;
    *base = *j;
    *j = c;
  }

  for (min = base; (hi = min += 1) < max;) {
    while (*(--hi) > *min);
    if ((hi += 1) != min) {
      for (lo = min + 1; --lo >= min;) {
	c = *lo;
	for (i = j = lo; (j -= 1) >= hi; i = j)
	   *i = *j;
	*i = c;
      }
    }
  }
}

static void siqst(idxtype *base, idxtype *max)
{
  register idxtype *i;
  register idxtype *j;
  register idxtype *jj;
  register idxtype *mid;
  register idxtype c;
  idxtype *tmp;
  int lo;
  int hi;

  lo = max - base;		/* number of elements as idxtype */
  do {
    mid = base + ((unsigned) lo>>1);
    if (lo >= MTHRESH) {
      j = (*base > *mid ? base : mid);
      tmp = max - 1;
      if (*j > *tmp) {
        j = (j == base ? mid : base); /* switch to first loser */
        if (*j < *tmp)
          j = tmp;
      }

      if (j != mid) {  /* SWAP */
        c = *mid;
        *mid = *j;
        *j = c;
      }
    }

    /* Semi-standard quicksort partitioning/swapping */
    for (i = base, j = max - 1;;) {
      while (i < mid && *i <= *mid)
        i++;
      while (j > mid) {
        if (*mid <= *j) {
          j--;
          continue;
        }
        tmp = i + 1;	/* value of i after swap */
        if (i == mid) 	/* j <-> mid, new mid is j */
          mid = jj = j;
        else 		/* i <-> j */
          jj = j--;
        goto swap;
      }

      if (i == mid)
	break;
      else {		/* i <-> mid, new mid is i */
        jj = mid;
        tmp = mid = i;	/* value of i after swap */
        j--;
      }
swap:
      c = *i;
      *i = *jj;
      *jj = c;
      i = tmp;
    }

    i = (j = mid) + 1;
    if ((lo = j - base) <= (hi = max - i)) {
      if (lo >= THRESH)
        siqst(base, j);
      base = i;
      lo = hi;
    }
    else {
      if (hi >= THRESH)
        siqst(i, max);
      max = j;
    }
  } while (lo >= THRESH);
}





/*************************************************************************
* Entry point of int increasing sort
**************************************************************************/
void iintsort(int n, int *base)
{
  register int *i;
  register int *j;
  register int *lo;
  register int *hi;
  register int *min;
  register int c;
  int *max;

  if (n <= 1)
    return;

  max = base + n;

  if (n >= THRESH) {
    iiqst(base, max);
    hi = base + THRESH;
  }
  else
    hi = max;

  for (j = lo = base; lo++ < hi;) {
    if (*j > *lo)
      j = lo;
  }
  if (j != base) { /* swap j into place */
    c = *base;
    *base = *j;
    *j = c;
  }

  for (min = base; (hi = min += 1) < max;) {
    while (*(--hi) > *min);
    if ((hi += 1) != min) {
      for (lo = min + 1; --lo >= min;) {
	c = *lo;
	for (i = j = lo; (j -= 1) >= hi; i = j)
	   *i = *j;
	*i = c;
      }
    }
  }
}


static void iiqst(int *base, int *max)
{
  register int *i;
  register int *j;
  register int *jj;
  register int *mid;
  register int c;
  int *tmp;
  int lo;
  int hi;

  lo = max - base;		/* number of elements as ints */
  do {
    mid = base + ((unsigned) lo>>1);
    if (lo >= MTHRESH) {
      j = (*base > *mid ? base : mid);
      tmp = max - 1;
      if (*j > *tmp) {
        j = (j == base ? mid : base); /* switch to first loser */
        if (*j < *tmp)
          j = tmp;
      }

      if (j != mid) {  /* SWAP */
        c = *mid;
        *mid = *j;
        *j = c;
      }
    }

    /* Semi-standard quicksort partitioning/swapping */
    for (i = base, j = max - 1;;) {
      while (i < mid && *i <= *mid)
        i++;
      while (j > mid) {
        if (*mid <= *j) {
          j--;
          continue;
        }
        tmp = i + 1;	/* value of i after swap */
        if (i == mid) 	/* j <-> mid, new mid is j */
          mid = jj = j;
        else 		/* i <-> j */
          jj = j--;
        goto swap;
      }

      if (i == mid)
	break;
      else {		/* i <-> mid, new mid is i */
        jj = mid;
        tmp = mid = i;	/* value of i after swap */
        j--;
      }
swap:
      c = *i;
      *i = *jj;
      *jj = c;
      i = tmp;
    }

    i = (j = mid) + 1;
    if ((lo = j - base) <= (hi = max - i)) {
      if (lo >= THRESH)
        iiqst(base, j);
      base = i;
      lo = hi;
    }
    else {
      if (hi >= THRESH)
        iiqst(i, max);
      max = j;
    }
  } while (lo >= THRESH);
}





/*************************************************************************
* Entry point of KeyVal increasing sort, ONLY key part
**************************************************************************/
void ikeysort(int n, KeyValueType *base)
{
  register KeyValueType *i;
  register KeyValueType *j;
  register KeyValueType *lo;
  register KeyValueType *hi;
  register KeyValueType *min;
  register KeyValueType c;
  KeyValueType *max;

  if (n <= 1)
    return;

  max = base + n;

  if (n >= THRESH) {
    keyiqst(base, max);
    hi = base + THRESH;
  }
  else
    hi = max;

  for (j = lo = base; lo++ < hi;) {
    if (j->key > lo->key)
      j = lo;
  }
  if (j != base) { /* swap j into place */
    c = *base;
    *base = *j;
    *j = c;
  }

  for (min = base; (hi = min += 1) < max;) {
    while ((--hi)->key > min->key);
    if ((hi += 1) != min) {
      for (lo = min + 1; --lo >= min;) {
	c = *lo;
	for (i = j = lo; (j -= 1) >= hi; i = j)
	   *i = *j;
	*i = c;
      }
    }
  }

  /* Sanity check */
  {
    int i;
    for (i=0; i<n-1; i++)
      if (base[i].key > base[i+1].key)
        printf("Something went wrong!\n");
  }
}


static void keyiqst(KeyValueType *base, KeyValueType *max)
{
  register KeyValueType *i;
  register KeyValueType *j;
  register KeyValueType *jj;
  register KeyValueType *mid;
  register KeyValueType c;
  KeyValueType *tmp;
  int lo;
  int hi;

  lo = (max - base)>>1;		/* number of elements as KeyValueType */
  do {
    mid = base + ((unsigned) lo>>1);
    if (lo >= MTHRESH) {
      j = (base->key > mid->key ? base : mid);
      tmp = max - 1;
      if (j->key > tmp->key) {
        j = (j == base ? mid : base); /* switch to first loser */
        if (j->key < tmp->key)
          j = tmp;
      }

      if (j != mid) {  /* SWAP */
        c = *mid;
        *mid = *j;
        *j = c;
      }
    }

    /* Semi-standard quicksort partitioning/swapping */
    for (i = base, j = max - 1;;) {
      while (i < mid && i->key <= mid->key)
        i++;
      while (j > mid) {
        if (mid->key <= j->key) {
          j--;
          continue;
        }
        tmp = i + 1;	/* value of i after swap */
        if (i == mid) 	/* j <-> mid, new mid is j */
          mid = jj = j;
        else 		/* i <-> j */
          jj = j--;
        goto swap;
      }

      if (i == mid)
	break;
      else {		/* i <-> mid, new mid is i */
        jj = mid;
        tmp = mid = i;	/* value of i after swap */
        j--;
      }
swap:
      c = *i;
      *i = *jj;
      *jj = c;
      i = tmp;
    }

    i = (j = mid) + 1;
    if ((lo = (j - base)>>1) <= (hi = (max - i)>>1)) {
      if (lo >= THRESH)
        keyiqst(base, j);
      base = i;
      lo = hi;
    }
    else {
      if (hi >= THRESH)
        keyiqst(i, max);
      max = j;
    }
  } while (lo >= THRESH);
}




/*************************************************************************
* Entry point of KeyVal increasing sort, BOTH key and val part
**************************************************************************/
void ikeyvalsort(int n, KeyValueType *base)
{
  register KeyValueType *i;
  register KeyValueType *j;
  register KeyValueType *lo;
  register KeyValueType *hi;
  register KeyValueType *min;
  register KeyValueType c;
  KeyValueType *max;

  if (n <= 1)
    return;

  max = base + n;

  if (n >= THRESH) {
    keyvaliqst(base, max);
    hi = base + THRESH;
  }
  else
    hi = max;

  for (j = lo = base; lo++ < hi;) {
    if ((j->key > lo->key) || (j->key == lo->key && j->val > lo->val))
      j = lo;
  }
  if (j != base) { /* swap j into place */
    c = *base;
    *base = *j;
    *j = c;
  }

  for (min = base; (hi = min += 1) < max;) {
    while ((--hi)->key > min->key || (hi->key == min->key && hi->val > min->val));
    if ((hi += 1) != min) {
      for (lo = min + 1; --lo >= min;) {
	c = *lo;
	for (i = j = lo; (j -= 1) >= hi; i = j)
	   *i = *j;
	*i = c;
      }
    }
  }
}


static void keyvaliqst(KeyValueType *base, KeyValueType *max)
{
  register KeyValueType *i;
  register KeyValueType *j;
  register KeyValueType *jj;
  register KeyValueType *mid;
  register KeyValueType c;
  KeyValueType *tmp;
  int lo;
  int hi;

  lo = (max - base)>>1;		/* number of elements as KeyValueType */
  do {
    mid = base + ((unsigned) lo>>1);
    if (lo >= MTHRESH) {
      j = (base->key > mid->key || (base->key == mid->key && base->val > mid->val) ? base : mid);
      tmp = max - 1;
      if (j->key > tmp->key || (j->key == tmp->key && j->val > tmp->val)) {
        j = (j == base ? mid : base); /* switch to first loser */
        if (j->key < tmp->key || (j->key == tmp->key && j->val < tmp->val))
          j = tmp;
      }

      if (j != mid) {  /* SWAP */
        c = *mid;
        *mid = *j;
        *j = c;
      }
    }

    /* Semi-standard quicksort partitioning/swapping */
    for (i = base, j = max - 1;;) {
      while (i < mid && (i->key < mid->key || (i->key == mid->key && i->val <= mid->val)))
        i++;
      while (j > mid) {
        if (mid->key < j->key || (mid->key == j->key && mid->val <= j->val)) {
          j--;
          continue;
        }
        tmp = i + 1;	/* value of i after swap */
        if (i == mid) 	/* j <-> mid, new mid is j */
          mid = jj = j;
        else 		/* i <-> j */
          jj = j--;
        goto swap;
      }

      if (i == mid)
	break;
      else {		/* i <-> mid, new mid is i */
        jj = mid;
        tmp = mid = i;	/* value of i after swap */
        j--;
      }
swap:
      c = *i;
      *i = *jj;
      *jj = c;
      i = tmp;
    }

    i = (j = mid) + 1;
    if ((lo = (j - base)>>1) <= (hi = (max - i)>>1)) {
      if (lo >= THRESH)
        keyvaliqst(base, j);
      base = i;
      lo = hi;
    }
    else {
      if (hi >= THRESH)
        keyvaliqst(i, max);
      max = j;
    }
  } while (lo >= THRESH);
}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ometis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point for OEMETIS
**************************************************************************/
void METIS_EdgeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
                  idxtype *perm, idxtype *iperm)
{
  int i;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_OEMETIS, *nvtxs, 1, xadj, adjncy, NULL, NULL, 0);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = OEMETIS_CTYPE;
    ctrl.IType = OEMETIS_ITYPE;
    ctrl.RType = OEMETIS_RTYPE;
    ctrl.dbglvl = OEMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.oflags  = 0;
  ctrl.pfactor = -1;
  ctrl.nseps   = 1;

  ctrl.optype = OP_OEMETIS;
  ctrl.CoarsenTo = 20;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, 2);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, *nvtxs);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
}


/*************************************************************************
* This function is the entry point for ONCMETIS
**************************************************************************/
void METIS_NodeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
                  idxtype *perm, idxtype *iperm)
{
  int i, ii, j, l;
  GraphType graph;
  CtrlType ctrl;
  idxtype *cptr, *cind, *piperm;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
    ctrl.oflags  = ONMETIS_OFLAGS;
    ctrl.pfactor = ONMETIS_PFACTOR;
    ctrl.nseps   = ONMETIS_NSEPS;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
    ctrl.oflags  = options[OPTION_OFLAGS];
    ctrl.pfactor = options[OPTION_PFACTOR];
    ctrl.nseps   = options[OPTION_NSEPS];
  }
  if (ctrl.nseps < 1)
    ctrl.nseps = 1;

  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  InitRandom(-1);

  if (ctrl.pfactor > 0) {
    /*============================================================
    * Prune the dense columns
    ==============================================================*/
    piperm = idxmalloc(*nvtxs, "ONMETIS: piperm");

    PruneGraph(&ctrl, &graph, *nvtxs, xadj, adjncy, piperm, (float)(0.1*ctrl.pfactor));
  }
  else if (ctrl.oflags&OFLAG_COMPRESS) {
    /*============================================================
    * Compress the graph
    ==============================================================*/
    cptr = idxmalloc(*nvtxs+1, "ONMETIS: cptr");
    cind = idxmalloc(*nvtxs, "ONMETIS: cind");

    CompressGraph(&ctrl, &graph, *nvtxs, xadj, adjncy, cptr, cind);

    if (graph.nvtxs >= COMPRESSION_FRACTION*(*nvtxs)) {
      ctrl.oflags--; /* We actually performed no compression */
      GKfree((void **) &cptr, &cind, LTERM);
    }
    else if (2*graph.nvtxs < *nvtxs && ctrl.nseps == 1)
      ctrl.nseps = 2;
  }
  else {
    SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, NULL, NULL, 0);
  }


  /*=============================================================
  * Do the nested dissection ordering
  --=============================================================*/
  ctrl.maxvwgt = 1.5*(idxsum(graph.nvtxs, graph.vwgt)/ctrl.CoarsenTo);
  AllocateWorkSpace(&ctrl, &graph, 2);

  if (ctrl.oflags&OFLAG_CCMP)
    MlevelNestedDissectionCC(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, graph.nvtxs);
  else
    MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, graph.nvtxs);

  FreeWorkSpace(&ctrl, &graph);

  if (ctrl.pfactor > 0) { /* Order any prunned vertices */
    if (graph.nvtxs < *nvtxs) {
      idxcopy(graph.nvtxs, iperm, perm);  /* Use perm as an auxiliary array */
      for (i=0; i<graph.nvtxs; i++)
        iperm[piperm[i]] = perm[i];
      for (i=graph.nvtxs; i<*nvtxs; i++)
        iperm[piperm[i]] = i;
    }

    GKfree((void **) &piperm, LTERM);
  }
  else if (ctrl.oflags&OFLAG_COMPRESS) { /* Uncompress the ordering */
    if (graph.nvtxs < COMPRESSION_FRACTION*(*nvtxs)) {
      /* construct perm from iperm */
      for (i=0; i<graph.nvtxs; i++)
        perm[iperm[i]] = i;
      for (l=ii=0; ii<graph.nvtxs; ii++) {
        i = perm[ii];
        for (j=cptr[i]; j<cptr[i+1]; j++)
          iperm[cind[j]] = l++;
      }
    }

    GKfree((void **) &cptr, &cind, LTERM);
  }


  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);

}


/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_NodeWND(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag,
                   int *options, idxtype *perm, idxtype *iperm)
{
  int i;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, NULL, 2);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = OFLAG_COMPRESS;
  ctrl.pfactor = 0;
  ctrl.nseps = 2;
  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, 2);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, *nvtxs);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
}




/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissection(CtrlType *ctrl, GraphType *graph, idxtype *order, float ubfactor, int lastvtx)
{
  int i, nvtxs, nbnd, tvwgt, tpwgts2[2];
  idxtype *label, *bndind;
  GraphType lgraph, rgraph;

  nvtxs = graph->nvtxs;

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  switch (ctrl->optype) {
    case OP_OEMETIS:
      MlevelEdgeBisection(ctrl, graph, tpwgts2, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SepTmr));
      ConstructMinCoverSeparator(ctrl, graph, ubfactor);
      IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SepTmr));

      break;
    case OP_ONMETIS:
      MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);

      IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

      break;
  }

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++)
    order[label[bndind[i]]] = --lastvtx;

  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  GKfree((void **) &graph->gdata, (void **) &graph->rdata,
         (void **) &graph->label, LTERM);

  if (rgraph.nvtxs > MMDSWITCH)
    MlevelNestedDissection(ctrl, &rgraph, order, ubfactor, lastvtx);
  else {
    MMDOrder(ctrl, &rgraph, order, lastvtx);
    GKfree((void **) &rgraph.gdata, (void **) &rgraph.rdata,
           (void **) &rgraph.label, LTERM);
  }
  if (lgraph.nvtxs > MMDSWITCH)
    MlevelNestedDissection(ctrl, &lgraph, order, ubfactor, lastvtx-rgraph.nvtxs);
  else {
    MMDOrder(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs);
    GKfree((void **) &lgraph.gdata, (void **) &lgraph.rdata,
           (void **) &lgraph.label, LTERM);
  }
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissectionCC(CtrlType *ctrl, GraphType *graph, idxtype *order, float ubfactor, int lastvtx)
{
  int i, nvtxs, nbnd, tvwgt, tpwgts2[2], nsgraphs, ncmps, rnvtxs;
  idxtype *label, *bndind;
  idxtype *cptr, *cind;
  GraphType *sgraphs;

  nvtxs = graph->nvtxs;

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);
  IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++)
    order[label[bndind[i]]] = --lastvtx;

  cptr = idxmalloc(nvtxs, "MlevelNestedDissectionCC: cptr");
  cind = idxmalloc(nvtxs, "MlevelNestedDissectionCC: cind");
  ncmps = FindComponents(ctrl, graph, cptr, cind);

/*
  if (ncmps > 2)
    printf("[%5d] has %3d components\n", nvtxs, ncmps);
*/

  sgraphs = (GraphType *)GKmalloc(ncmps*sizeof(GraphType), "MlevelNestedDissectionCC: sgraphs");

  nsgraphs = SplitGraphOrderCC(ctrl, graph, sgraphs, ncmps, cptr, cind);

  GKfree((void **) &cptr, (void **) &cind, LTERM);

  /* Free the memory of the top level graph */
  GKfree((void **) &graph->gdata, (void **) &graph->rdata,
         (void **) &graph->label, LTERM);

  /* Go and process the subgraphs */
  for (rnvtxs=i=0; i<nsgraphs; i++) {
    if (sgraphs[i].adjwgt == NULL) {
      MMDOrder(ctrl, sgraphs+i, order, lastvtx-rnvtxs);
      GKfree((void **) &sgraphs[i].gdata, (void **) &sgraphs[i].label, LTERM);
    }
    else {
      MlevelNestedDissectionCC(ctrl, sgraphs+i, order, ubfactor, lastvtx-rnvtxs);
    }
    rnvtxs += sgraphs[i].nvtxs;
  }

  free(sgraphs);
}



/*************************************************************************
* This function performs multilevel bisection. It performs multiple
* bisections and selects the best.
**************************************************************************/
void MlevelNodeBisectionMultiple(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  int i, nvtxs, cnvtxs, mincut;
  GraphType *cgraph;
  idxtype *bestwhere;

  if (ctrl->nseps == 1 || graph->nvtxs < (ctrl->oflags&OFLAG_COMPRESS ? 1000 : 2000)) {
    MlevelNodeBisection(ctrl, graph, tpwgts, ubfactor);
    return;
  }

  nvtxs = graph->nvtxs;

  if (ctrl->oflags&OFLAG_COMPRESS) { /* Multiple separators at the original graph */
    bestwhere = idxmalloc(nvtxs, "MlevelNodeBisection2: bestwhere");
    mincut = nvtxs;

    for (i=ctrl->nseps; i>0; i--) {
      MlevelNodeBisection(ctrl, graph, tpwgts, ubfactor);

      /* printf("%5d ", cgraph->mincut); */

      if (graph->mincut < mincut) {
        mincut = graph->mincut;
        idxcopy(nvtxs, graph->where, bestwhere);
      }

      GKfree((void **) &graph->rdata, LTERM);

      if (mincut == 0)
        break;
    }
    /* printf("[%5d]\n", mincut); */

    Allocate2WayNodePartitionMemory(ctrl, graph);
    idxcopy(nvtxs, bestwhere, graph->where);
    free(bestwhere);

    Compute2WayNodePartitionParams(ctrl, graph);
  }
  else {  /* Coarsen it a bit */
    ctrl->CoarsenTo = nvtxs-1;

    cgraph = Coarsen2Way(ctrl, graph);

    cnvtxs = cgraph->nvtxs;

    bestwhere = idxmalloc(cnvtxs, "MlevelNodeBisection2: bestwhere");
    mincut = nvtxs;

    for (i=ctrl->nseps; i>0; i--) {
      ctrl->CType += 20; /* This is a hack. Look at coarsen.c */
      MlevelNodeBisection(ctrl, cgraph, tpwgts, ubfactor);

      /* printf("%5d ", cgraph->mincut); */

      if (cgraph->mincut < mincut) {
        mincut = cgraph->mincut;
        idxcopy(cnvtxs, cgraph->where, bestwhere);
      }

      GKfree((void **) &cgraph->rdata, LTERM);

      if (mincut == 0)
        break;
    }
    /* printf("[%5d]\n", mincut); */

    Allocate2WayNodePartitionMemory(ctrl, cgraph);
    idxcopy(cnvtxs, bestwhere, cgraph->where);
    free(bestwhere);

    Compute2WayNodePartitionParams(ctrl, cgraph);

    Refine2WayNode(ctrl, graph, cgraph, ubfactor);
  }

}

/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void MlevelNodeBisection(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  GraphType *cgraph;

  ctrl->CoarsenTo = graph->nvtxs/8;
  if (ctrl->CoarsenTo > 100)
    ctrl->CoarsenTo = 100;
  else if (ctrl->CoarsenTo < 40)
    ctrl->CoarsenTo = 40;
  ctrl->maxvwgt = 1.5*((tpwgts[0]+tpwgts[1])/ctrl->CoarsenTo);

  cgraph = Coarsen2Way(ctrl, graph);

  switch (ctrl->IType) {
    case IPART_GGPKL:
      Init2WayPartition(ctrl, cgraph, tpwgts, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SepTmr));

      Compute2WayPartitionParams(ctrl, cgraph);
      ConstructSeparator(ctrl, cgraph, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SepTmr));
      break;
    case IPART_GGPKLNODE:
      InitSeparator(ctrl, cgraph, ubfactor);
      break;
  }

  Refine2WayNode(ctrl, graph, cgraph, ubfactor);

}




/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
* This function relies on the fact that adjwgt is all equal to 1.
**************************************************************************/
void SplitGraphOrder(CtrlType *ctrl, GraphType *graph, GraphType *lgraph, GraphType *rgraph)
{
  int i, ii, j, k, l, istart, iend, mypart, nvtxs, snvtxs[3], snedges[3];
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *label, *where, *bndptr, *bndind;
  idxtype *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *sadjwgtsum[2], *slabel[2];
  idxtype *rename;
  idxtype *auxadjncy;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SplitTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  label = graph->label;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;
  ASSERT(bndptr != NULL);

  rename = idxwspacemalloc(ctrl, nvtxs);

  snvtxs[0] = snvtxs[1] = snvtxs[2] = snedges[0] = snedges[1] = snedges[2] = 0;
  for (i=0; i<nvtxs; i++) {
    k = where[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  SetUpSplitGraph(graph, lgraph, snvtxs[0], snedges[0]);
  sxadj[0] = lgraph->xadj;
  svwgt[0] = lgraph->vwgt;
  sadjwgtsum[0] = lgraph->adjwgtsum;
  sadjncy[0] = lgraph->adjncy;
  sadjwgt[0] = lgraph->adjwgt;
  slabel[0] = lgraph->label;

  SetUpSplitGraph(graph, rgraph, snvtxs[1], snedges[1]);
  sxadj[1] = rgraph->xadj;
  svwgt[1] = rgraph->vwgt;
  sadjwgtsum[1] = rgraph->adjwgtsum;
  sadjncy[1] = rgraph->adjncy;
  sadjwgt[1] = rgraph->adjwgt;
  slabel[1] = rgraph->label;

  /* Go and use bndptr to also mark the boundary nodes in the two partitions */
  for (ii=0; ii<graph->nbnd; ii++) {
    i = bndind[ii];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      bndptr[adjncy[j]] = 1;
  }

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    if ((mypart = where[i]) == 2)
      continue;

    istart = xadj[i];
    iend = xadj[i+1];
    if (bndptr[i] == -1) { /* This is an interior vertex */
      auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
      for(j=istart; j<iend; j++)
        auxadjncy[j] = adjncy[j];
      snedges[mypart] += iend-istart;
    }
    else {
      auxadjncy = sadjncy[mypart];
      l = snedges[mypart];
      for (j=istart; j<iend; j++) {
        k = adjncy[j];
        if (where[k] == mypart)
          auxadjncy[l++] = k;
      }
      snedges[mypart] = l;
    }

    svwgt[mypart][snvtxs[mypart]] = vwgt[i];
    sadjwgtsum[mypart][snvtxs[mypart]] = snedges[mypart]-sxadj[mypart][snvtxs[mypart]];
    slabel[mypart][snvtxs[mypart]] = label[i];
    sxadj[mypart][++snvtxs[mypart]] = snedges[mypart];
  }

  for (mypart=0; mypart<2; mypart++) {
    iend = snedges[mypart];
    idxset(iend, 1, sadjwgt[mypart]);

    auxadjncy = sadjncy[mypart];
    for (i=0; i<iend; i++)
      auxadjncy[i] = rename[auxadjncy[i]];
  }

  lgraph->nvtxs = snvtxs[0];
  lgraph->nedges = snedges[0];
  rgraph->nvtxs = snvtxs[1];
  rgraph->nedges = snedges[1];

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SplitTmr));

  idxwspacefree(ctrl, nvtxs);

}

/*************************************************************************
* This function uses MMD to order the graph. The vertices are numbered
* from lastvtx downwards
**************************************************************************/
void MMDOrder(CtrlType *ctrl, GraphType *graph, idxtype *order, int lastvtx)
{
  int i, k, nvtxs, nofsub, firstvtx;
  idxtype *xadj, *adjncy, *label;
  idxtype *perm, *iperm, *head, *qsize, *list, *marker;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  /* Relabel the vertices so that it starts from 1 */
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]++;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]++;

  perm = idxmalloc(6*(nvtxs+5), "MMDOrder: perm");
  iperm = perm + nvtxs + 5;
  head = iperm + nvtxs + 5;
  qsize = head + nvtxs + 5;
  list = qsize + nvtxs + 5;
  marker = list + nvtxs + 5;

  genmmd(nvtxs, xadj, adjncy, iperm, perm, 1, head, qsize, list, marker, MAXIDX, &nofsub);

  label = graph->label;
  firstvtx = lastvtx-nvtxs;
  for (i=0; i<nvtxs; i++)
    order[label[i]] = firstvtx+iperm[i]-1;

  free(perm);

  /* Relabel the vertices so that it starts from 0 */
  for (i=0; i<nvtxs+1; i++)
    xadj[i]--;
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]--;
}


/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
* It relies on the fact that adjwgt is all set to 1.
**************************************************************************/
int SplitGraphOrderCC(CtrlType *ctrl, GraphType *graph, GraphType *sgraphs, int ncmps, idxtype *cptr, idxtype *cind)
{
  int i, ii, iii, j, k, l, istart, iend, nvtxs, snvtxs, snedges;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *label, *where, *bndptr, *bndind;
  idxtype *sxadj, *svwgt, *sadjncy, *sadjwgt, *sadjwgtsum, *slabel;
  idxtype *rename;
  idxtype *auxadjncy, *auxadjwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SplitTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  label = graph->label;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;
  ASSERT(bndptr != NULL);

  /* Go and use bndptr to also mark the boundary nodes in the two partitions */
  for (ii=0; ii<graph->nbnd; ii++) {
    i = bndind[ii];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      bndptr[adjncy[j]] = 1;
  }

  rename = idxwspacemalloc(ctrl, nvtxs);

  /* Go and split the graph a component at a time */
  for (iii=0; iii<ncmps; iii++) {
    RandomPermute(cptr[iii+1]-cptr[iii], cind+cptr[iii], 0);
    snvtxs = snedges = 0;
    for (j=cptr[iii]; j<cptr[iii+1]; j++) {
      i = cind[j];
      rename[i] = snvtxs++;
      snedges += xadj[i+1]-xadj[i];
    }

    SetUpSplitGraph(graph, sgraphs+iii, snvtxs, snedges);
    sxadj = sgraphs[iii].xadj;
    svwgt = sgraphs[iii].vwgt;
    sadjwgtsum = sgraphs[iii].adjwgtsum;
    sadjncy = sgraphs[iii].adjncy;
    sadjwgt = sgraphs[iii].adjwgt;
    slabel = sgraphs[iii].label;

    snvtxs = snedges = sxadj[0] = 0;
    for (ii=cptr[iii]; ii<cptr[iii+1]; ii++) {
      i = cind[ii];

      istart = xadj[i];
      iend = xadj[i+1];
      if (bndptr[i] == -1) { /* This is an interior vertex */
        auxadjncy = sadjncy + snedges - istart;
        auxadjwgt = sadjwgt + snedges - istart;
        for(j=istart; j<iend; j++)
          auxadjncy[j] = adjncy[j];
        snedges += iend-istart;
      }
      else {
        l = snedges;
        for (j=istart; j<iend; j++) {
          k = adjncy[j];
          if (where[k] != 2)
            sadjncy[l++] = k;
        }
        snedges = l;
      }

      svwgt[snvtxs] = vwgt[i];
      sadjwgtsum[snvtxs] = snedges-sxadj[snvtxs];
      slabel[snvtxs] = label[i];
      sxadj[++snvtxs] = snedges;
    }

    idxset(snedges, 1, sadjwgt);
    for (i=0; i<snedges; i++)
      sadjncy[i] = rename[sadjncy[i]];

    sgraphs[iii].nvtxs = snvtxs;
    sgraphs[iii].nedges = snedges;
    sgraphs[iii].ncon = 1;

    if (snvtxs < MMDSWITCH)
      sgraphs[iii].adjwgt = NULL;  /* A marker to call MMD on the driver */
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SplitTmr));

  idxwspacefree(ctrl, nvtxs);

  return ncmps;

}





/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * parmetis.c
 *
 * This file contains top level routines that are used by ParMETIS
 *
 * Started 10/14/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point for KMETIS with seed specification
* in options[7]
**************************************************************************/
void METIS_PartGraphKway2(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                         idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                         int *options, int *edgecut, idxtype *part)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++)
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphKway2(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts,
                       tpwgts, options, edgecut, part);

  free(tpwgts);
}


/*************************************************************************
* This function is the entry point for KWMETIS with seed specification
* in options[7]
**************************************************************************/
void METIS_WPartGraphKway2(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                          idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                          float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = KMETIS_CTYPE;
    ctrl.IType = KMETIS_ITYPE;
    ctrl.RType = KMETIS_RTYPE;
    ctrl.dbglvl = KMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = 20*(*nparts);
  ctrl.maxvwgt = 1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt) : (*nvtxs))/ctrl.CoarsenTo);

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MlevelKWayPartitioning(&ctrl, &graph, *nparts, part, tpwgts, 1.03);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function is the entry point for the node ND code for ParMETIS
**************************************************************************/
void METIS_NodeNDP(int nvtxs, idxtype *xadj, idxtype *adjncy, int npes,
                   int *options, idxtype *perm, idxtype *iperm, idxtype *sizes)
{
  int i, ii, j, l;
  GraphType graph;
  CtrlType ctrl;
  idxtype *cptr, *cind;

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
    ctrl.oflags  = ONMETIS_OFLAGS;
    ctrl.pfactor = ONMETIS_PFACTOR;
    ctrl.nseps   = ONMETIS_NSEPS;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
    ctrl.oflags  = options[OPTION_OFLAGS];
    ctrl.pfactor = options[OPTION_PFACTOR];
    ctrl.nseps   = options[OPTION_NSEPS];
  }
  if (ctrl.nseps < 1)
    ctrl.nseps = 1;

  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  InitRandom(-1);

  if (ctrl.oflags&OFLAG_COMPRESS) {
    /*============================================================
    * Compress the graph
    ==============================================================*/
    cptr = idxmalloc(nvtxs+1, "ONMETIS: cptr");
    cind = idxmalloc(nvtxs, "ONMETIS: cind");

    CompressGraph(&ctrl, &graph, nvtxs, xadj, adjncy, cptr, cind);

    if (graph.nvtxs >= COMPRESSION_FRACTION*(nvtxs)) {
      ctrl.oflags--; /* We actually performed no compression */
      GKfree((void **) &cptr, &cind, LTERM);
    }
    else if (2*graph.nvtxs < nvtxs && ctrl.nseps == 1)
      ctrl.nseps = 2;
  }
  else {
    SetUpGraph(&graph, OP_ONMETIS, nvtxs, 1, xadj, adjncy, NULL, NULL, 0);
  }


  /*=============================================================
  * Do the nested dissection ordering
  --=============================================================*/
  ctrl.maxvwgt = 1.5*(idxsum(graph.nvtxs, graph.vwgt)/ctrl.CoarsenTo);
  AllocateWorkSpace(&ctrl, &graph, 2);

  idxset(2*npes-1, 0, sizes);
  MlevelNestedDissectionP(&ctrl, &graph, iperm, graph.nvtxs, npes, 0, sizes);

  FreeWorkSpace(&ctrl, &graph);

  if (ctrl.oflags&OFLAG_COMPRESS) { /* Uncompress the ordering */
    if (graph.nvtxs < COMPRESSION_FRACTION*(nvtxs)) {
      /* construct perm from iperm */
      for (i=0; i<graph.nvtxs; i++)
        perm[iperm[i]] = i;
      for (l=ii=0; ii<graph.nvtxs; ii++) {
        i = perm[ii];
        for (j=cptr[i]; j<cptr[i+1]; j++)
          iperm[cind[j]] = l++;
      }
    }

    GKfree((void **) &cptr, (void **) &cind, LTERM);
  }


  for (i=0; i<nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

}



/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissectionP(CtrlType *ctrl, GraphType *graph, idxtype *order, int lastvtx,
                             int npes, int cpos, idxtype *sizes)
{
  int i, nvtxs, nbnd, tvwgt, tpwgts2[2];
  idxtype *label, *bndind;
  GraphType lgraph, rgraph;
  float ubfactor;

  nvtxs = graph->nvtxs;

  if (nvtxs == 0) {
    GKfree((void **) &graph->gdata, (void **) &graph->rdata,
           (void **) &graph->label, LTERM);
    return;
  }

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  if (cpos >= npes-1)
    ubfactor = ORDER_UNBALANCE_FRACTION;
  else
    ubfactor = 1.05;


  MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);

  IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

  if (cpos < npes-1) {
    sizes[2*npes-2-cpos] = graph->pwgts[2];
    sizes[2*npes-2-(2*cpos+1)] = graph->pwgts[1];
    sizes[2*npes-2-(2*cpos+2)] = graph->pwgts[0];
  }

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++)
    order[label[bndind[i]]] = --lastvtx;

  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  GKfree((void **) &graph->gdata, (void **) &graph->rdata,
         (void **) &graph->label, LTERM);

  if (rgraph.nvtxs > MMDSWITCH || 2*cpos+1 < npes-1)
    MlevelNestedDissectionP(ctrl, &rgraph, order, lastvtx, npes, 2*cpos+1, sizes);
  else {
    MMDOrder(ctrl, &rgraph, order, lastvtx);
    GKfree((void **) &rgraph.gdata, (void **) &rgraph.rdata,
           (void **) &rgraph.label, LTERM);
  }
  if (lgraph.nvtxs > MMDSWITCH || 2*cpos+2 < npes-1)
    MlevelNestedDissectionP(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs, npes, 2*cpos+2, sizes);
  else {
    MMDOrder(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs);
    GKfree((void **) &lgraph.gdata, (void **) &lgraph.rdata,
           (void **) &lgraph.label, LTERM);
  }
}




/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_NodeComputeSeparator(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
           idxtype *adjwgt, int *options, int *sepsize, idxtype *part)
{
  int tvwgt, tpwgts[2];
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);
  tvwgt = idxsum(*nvtxs, graph.vwgt);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = 0;
  ctrl.pfactor = 0;
  ctrl.nseps = 1;
  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = amin(100, *nvtxs-1);
  ctrl.maxvwgt = 1.5*tvwgt/ctrl.CoarsenTo;

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, 2);

  /*============================================================
   * Perform the bisection
   *============================================================*/
  tpwgts[0] = tvwgt/2;
  tpwgts[1] = tvwgt-tpwgts[0];

  MlevelNodeBisectionMultiple(&ctrl, &graph, tpwgts, 1.05);

  *sepsize = graph.pwgts[2];
  idxcopy(*nvtxs, graph.where, part);

  GKfree((void **) &graph.gdata, (void **) &graph.rdata,
         (void **) &graph.label, LTERM);


  FreeWorkSpace(&ctrl, &graph);

}



/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_EdgeComputeSeparator(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
           idxtype *adjwgt, int *options, int *sepsize, idxtype *part)
{
  int tvwgt, tpwgts[2];
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);
  tvwgt = idxsum(*nvtxs, graph.vwgt);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = 0;
  ctrl.pfactor = 0;
  ctrl.nseps = 1;
  ctrl.optype = OP_OEMETIS;
  ctrl.CoarsenTo = amin(100, *nvtxs-1);
  ctrl.maxvwgt = 1.5*tvwgt/ctrl.CoarsenTo;

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, 2);

  /*============================================================
   * Perform the bisection
   *============================================================*/
  tpwgts[0] = tvwgt/2;
  tpwgts[1] = tvwgt-tpwgts[0];

  MlevelEdgeBisection(&ctrl, &graph, tpwgts, 1.05);
  ConstructMinCoverSeparator(&ctrl, &graph, 1.05);

  *sepsize = graph.pwgts[2];
  idxcopy(*nvtxs, graph.where, part);

  GKfree((void **) &graph.gdata, (void **) &graph.rdata,
         (void **) &graph.label, LTERM);


  FreeWorkSpace(&ctrl, &graph);

}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pmetis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point for PMETIS
**************************************************************************/
void METIS_PartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                              idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                              int *options, int *edgecut, idxtype *part)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++)
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts,
                            tpwgts, options, edgecut, part);

  free(tpwgts);
}



/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_WPartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                               idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                               float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  int i;
  GraphType graph;
  CtrlType ctrl;
  float *mytpwgts;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_PMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = PMETIS_CTYPE;
    ctrl.IType = PMETIS_ITYPE;
    ctrl.RType = PMETIS_RTYPE;
    ctrl.dbglvl = PMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 20;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  mytpwgts = fmalloc(*nparts, "PWMETIS: mytpwgts");
  for (i=0; i<*nparts; i++)
    mytpwgts[i] = tpwgts[i];

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MlevelRecursiveBisection(&ctrl, &graph, *nparts, part, mytpwgts, 1.000, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);
  free(mytpwgts);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}



/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MlevelRecursiveBisection(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part, float *tpwgts, float ubfactor, int fpart)
{
  int i, nvtxs, cut, tvwgt, tpwgts2[2];
  idxtype *label, *where;
  GraphType lgraph, rgraph;
  float wsum;

  nvtxs = graph->nvtxs;
  if (nvtxs == 0) {
    printf("\t***Cannot bisect a graph with 0 vertices!\n\t***You are trying to partition a graph into too many parts!\n");
    return 0;
  }

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt*ssum(nparts/2, tpwgts);
  tpwgts2[1] = tvwgt-tpwgts2[0];

  MlevelEdgeBisection(ctrl, graph, tpwgts2, ubfactor);
  cut = graph->mincut;

  /* printf("%5d %5d %5d [%5d %f]\n", tpwgts2[0], tpwgts2[1], cut, tvwgt, ssum(nparts/2, tpwgts));*/

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

  if (nparts > 2) {
    SplitGraphPart(ctrl, graph, &lgraph, &rgraph);
    /* printf("%d %d\n", lgraph.nvtxs, rgraph.nvtxs); */
  }


  /* Free the memory of the top level graph */
  GKfree((void **) &graph->gdata, (void **) &graph->rdata,
         (void **) &graph->label, LTERM);

  /* Scale the fractions in the tpwgts according to the true weight */
  wsum = ssum(nparts/2, tpwgts);
  sscale(nparts/2, 1.0/wsum, tpwgts);
  sscale(nparts-nparts/2, 1.0/(1.0-wsum), tpwgts+nparts/2);
  /*
  for (i=0; i<nparts; i++)
    printf("%5.3f ", tpwgts[i]);
  printf("[%5.3f]\n", wsum);
  */

  /* Do the recursive call */
  if (nparts > 3) {
    cut += MlevelRecursiveBisection(ctrl, &lgraph, nparts/2, part, tpwgts, ubfactor, fpart);
    cut += MlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, tpwgts+nparts/2, ubfactor, fpart+nparts/2);
  }
  else if (nparts == 3) {
    cut += MlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, tpwgts+nparts/2, ubfactor, fpart+nparts/2);
    GKfree((void **) &lgraph.gdata, (void **) &lgraph.label, LTERM);
  }

  return cut;

}


/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void MlevelEdgeBisection(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  GraphType *cgraph;

  cgraph = Coarsen2Way(ctrl, graph);

  Init2WayPartition(ctrl, cgraph, tpwgts, ubfactor);

  Refine2Way(ctrl, graph, cgraph, tpwgts, ubfactor);

/*
  IsConnectedSubdomain(ctrl, graph, 0);
  IsConnectedSubdomain(ctrl, graph, 1);
*/
}




/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
**************************************************************************/
void SplitGraphPart(CtrlType *ctrl, GraphType *graph, GraphType *lgraph, GraphType *rgraph)
{
  int i, j, k, kk, l, istart, iend, mypart, nvtxs, ncon, snvtxs[2], snedges[2], sum;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *label, *where, *bndptr;
  idxtype *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *sadjwgtsum[2], *slabel[2];
  idxtype *rename;
  idxtype *auxadjncy, *auxadjwgt;
  float *nvwgt, *snvwgt[2], *npwgts;


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SplitTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  label = graph->label;
  where = graph->where;
  bndptr = graph->bndptr;
  npwgts = graph->npwgts;

  ASSERT(bndptr != NULL);

  rename = idxwspacemalloc(ctrl, nvtxs);

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  for (i=0; i<nvtxs; i++) {
    k = where[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  SetUpSplitGraph(graph, lgraph, snvtxs[0], snedges[0]);
  sxadj[0] = lgraph->xadj;
  svwgt[0] = lgraph->vwgt;
  snvwgt[0] = lgraph->nvwgt;
  sadjwgtsum[0] = lgraph->adjwgtsum;
  sadjncy[0] = lgraph->adjncy;
  sadjwgt[0] = lgraph->adjwgt;
  slabel[0] = lgraph->label;

  SetUpSplitGraph(graph, rgraph, snvtxs[1], snedges[1]);
  sxadj[1] = rgraph->xadj;
  svwgt[1] = rgraph->vwgt;
  snvwgt[1] = rgraph->nvwgt;
  sadjwgtsum[1] = rgraph->adjwgtsum;
  sadjncy[1] = rgraph->adjncy;
  sadjwgt[1] = rgraph->adjwgt;
  slabel[1] = rgraph->label;

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    mypart = where[i];
    sum = adjwgtsum[i];

    istart = xadj[i];
    iend = xadj[i+1];
    if (bndptr[i] == -1) { /* This is an interior vertex */
      auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
      auxadjwgt = sadjwgt[mypart] + snedges[mypart] - istart;
      for(j=istart; j<iend; j++) {
        auxadjncy[j] = adjncy[j];
        auxadjwgt[j] = adjwgt[j];
      }
      snedges[mypart] += iend-istart;
    }
    else {
      auxadjncy = sadjncy[mypart];
      auxadjwgt = sadjwgt[mypart];
      l = snedges[mypart];
      for (j=istart; j<iend; j++) {
        k = adjncy[j];
        if (where[k] == mypart) {
          auxadjncy[l] = k;
          auxadjwgt[l++] = adjwgt[j];
        }
        else {
          sum -= adjwgt[j];
        }
      }
      snedges[mypart] = l;
    }

    if (ncon == 1)
      svwgt[mypart][snvtxs[mypart]] = vwgt[i];
    else {
      for (kk=0; kk<ncon; kk++)
        snvwgt[mypart][snvtxs[mypart]*ncon+kk] = nvwgt[i*ncon+kk]/npwgts[mypart*ncon+kk];
    }

    sadjwgtsum[mypart][snvtxs[mypart]] = sum;
    slabel[mypart][snvtxs[mypart]] = label[i];
    sxadj[mypart][++snvtxs[mypart]] = snedges[mypart];
  }

  for (mypart=0; mypart<2; mypart++) {
    iend = sxadj[mypart][snvtxs[mypart]];
    auxadjncy = sadjncy[mypart];
    for (i=0; i<iend; i++)
      auxadjncy[i] = rename[auxadjncy[i]];
  }

  lgraph->nedges = snedges[0];
  rgraph->nedges = snedges[1];

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SplitTmr));

  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* Setup the various arrays for the splitted graph
**************************************************************************/
void SetUpSplitGraph(GraphType *graph, GraphType *sgraph, int snvtxs, int snedges)
{
  InitGraph(sgraph);
  sgraph->nvtxs = snvtxs;
  sgraph->nedges = snedges;
  sgraph->ncon = graph->ncon;

  /* Allocate memory for the splitted graph */
  if (graph->ncon == 1) {
    sgraph->gdata = idxmalloc(4*snvtxs+1 + 2*snedges, "SetUpSplitGraph: gdata");

    sgraph->xadj        = sgraph->gdata;
    sgraph->vwgt        = sgraph->gdata + snvtxs+1;
    sgraph->adjwgtsum   = sgraph->gdata + 2*snvtxs+1;
    sgraph->cmap        = sgraph->gdata + 3*snvtxs+1;
    sgraph->adjncy      = sgraph->gdata + 4*snvtxs+1;
    sgraph->adjwgt      = sgraph->gdata + 4*snvtxs+1 + snedges;
  }
  else {
    sgraph->gdata = idxmalloc(3*snvtxs+1 + 2*snedges, "SetUpSplitGraph: gdata");

    sgraph->xadj        = sgraph->gdata;
    sgraph->adjwgtsum   = sgraph->gdata + snvtxs+1;
    sgraph->cmap        = sgraph->gdata + 2*snvtxs+1;
    sgraph->adjncy      = sgraph->gdata + 3*snvtxs+1;
    sgraph->adjwgt      = sgraph->gdata + 3*snvtxs+1 + snedges;

    sgraph->nvwgt       = fmalloc(graph->ncon*snvtxs, "SetUpSplitGraph: nvwgt");
  }

  sgraph->label	= idxmalloc(snvtxs, "SetUpSplitGraph: sgraph->label");
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pqueue.c
 *
 * This file contains functions for manipulating the bucket list
 * representation of the gains associated with each vertex in a graph.
 * These functions are used by the refinement algorithms
 *
 * Started 9/2/94
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function initializes the data structures of the priority queue
**************************************************************************/
void PQueueInit(CtrlType *ctrl, PQueueType *queue, int maxnodes, int maxgain)
{
  int i, j, ncore;

  queue->nnodes = 0;
  queue->maxnodes = maxnodes;

  queue->buckets = NULL;
  queue->nodes = NULL;
  queue->heap = NULL;
  queue->locator = NULL;

  if (maxgain > PLUS_GAINSPAN || maxnodes < 500)
    queue->type = 2;
  else
    queue->type = 1;

  if (queue->type == 1) {
    queue->pgainspan = amin(PLUS_GAINSPAN, maxgain);
    queue->ngainspan = amin(NEG_GAINSPAN, maxgain);

    j = queue->ngainspan+queue->pgainspan+1;

    ncore = 2 + (sizeof(ListNodeType)/sizeof(idxtype))*maxnodes + (sizeof(ListNodeType *)/sizeof(idxtype))*j;

    if (WspaceAvail(ctrl) > ncore) {
      queue->nodes = (ListNodeType *)idxwspacemalloc(ctrl, (sizeof(ListNodeType)/sizeof(idxtype))*maxnodes);
      queue->buckets = (ListNodeType **)idxwspacemalloc(ctrl, (sizeof(ListNodeType *)/sizeof(idxtype))*j);
      queue->mustfree = 0;
    }
    else { /* Not enough memory in the wspace, allocate it */
      queue->nodes = (ListNodeType *)idxmalloc((sizeof(ListNodeType)/sizeof(idxtype))*maxnodes, "PQueueInit: queue->nodes");
      queue->buckets = (ListNodeType **)idxmalloc((sizeof(ListNodeType *)/sizeof(idxtype))*j, "PQueueInit: queue->buckets");
      queue->mustfree = 1;
    }

    for (i=0; i<maxnodes; i++)
      queue->nodes[i].id = i;

    for (i=0; i<j; i++)
      queue->buckets[i] = NULL;

    queue->buckets += queue->ngainspan;  /* Advance buckets by the ngainspan proper indexing */
    queue->maxgain = -queue->ngainspan;
  }
  else {
    queue->heap = (KeyValueType *)idxwspacemalloc(ctrl, (sizeof(KeyValueType)/sizeof(idxtype))*maxnodes);
    queue->locator = idxwspacemalloc(ctrl, maxnodes);
    idxset(maxnodes, -1, queue->locator);
  }

}


/*************************************************************************
* This function resets the buckets
**************************************************************************/
void PQueueReset(PQueueType *queue)
{
  int i, j;
  queue->nnodes = 0;

  if (queue->type == 1) {
    queue->maxgain = -queue->ngainspan;

    j = queue->ngainspan+queue->pgainspan+1;
    queue->buckets -= queue->ngainspan;
    for (i=0; i<j; i++)
      queue->buckets[i] = NULL;
    queue->buckets += queue->ngainspan;
  }
  else {
    idxset(queue->maxnodes, -1, queue->locator);
  }

}


/*************************************************************************
* This function frees the buckets
**************************************************************************/
void PQueueFree(CtrlType *ctrl, PQueueType *queue)
{

  if (queue->type == 1) {
    if (queue->mustfree) {
      queue->buckets -= queue->ngainspan;
      GKfree((void **) &queue->nodes, (void **) &queue->buckets, LTERM);
    }
    else {
      idxwspacefree(ctrl, sizeof(ListNodeType *)*(queue->ngainspan+queue->pgainspan+1)/sizeof(idxtype));
      idxwspacefree(ctrl, sizeof(ListNodeType)*queue->maxnodes/sizeof(idxtype));
    }
  }
  else {
    idxwspacefree(ctrl, sizeof(KeyValueType)*queue->maxnodes/sizeof(idxtype));
    idxwspacefree(ctrl, queue->maxnodes);
  }

  queue->maxnodes = 0;
}


/*************************************************************************
* This function returns the number of nodes in the queue
**************************************************************************/
int PQueueGetSize(PQueueType *queue)
{
  return queue->nnodes;
}


/*************************************************************************
* This function adds a node of certain gain into a partition
**************************************************************************/
int PQueueInsert(PQueueType *queue, int node, int gain)
{
  int i, j;
  idxtype *locator;
  ListNodeType *newnode;
  KeyValueType *heap;

  if (queue->type == 1) {
    ASSERT(gain >= -queue->ngainspan && gain <= queue->pgainspan);

    /* Allocate and add the node */
    queue->nnodes++;
    newnode = queue->nodes + node;

    /* Attach this node in the doubly-linked list */
    newnode->next = queue->buckets[gain];
    newnode->prev = NULL;
    if (newnode->next != NULL)
      newnode->next->prev = newnode;
    queue->buckets[gain] = newnode;

    if (queue->maxgain < gain)
      queue->maxgain = gain;
  }
  else {
    ASSERT(CheckHeap(queue));

    heap = queue->heap;
    locator = queue->locator;

    ASSERT(locator[node] == -1);

    i = queue->nnodes++;
    while (i > 0) {
      j = (i-1)/2;
      if (heap[j].key < gain) {
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
        break;
    }
    ASSERT(i >= 0);
    heap[i].key = gain;
    heap[i].val = node;
    locator[node] = i;

    ASSERT(CheckHeap(queue));
  }

  return 0;
}


/*************************************************************************
* This function deletes a node from a partition and reinserts it with
* an updated gain
**************************************************************************/
int PQueueDelete(PQueueType *queue, int node, int gain)
{
  int i, j, newgain, oldgain;
  idxtype *locator;
  ListNodeType *newnode, **buckets;
  KeyValueType *heap;

  if (queue->type == 1) {
    ASSERT(gain >= -queue->ngainspan && gain <= queue->pgainspan);
    ASSERT(queue->nnodes > 0);

    buckets = queue->buckets;
    queue->nnodes--;
    newnode = queue->nodes+node;

    /* Remove newnode from the doubly-linked list */
    if (newnode->prev != NULL)
      newnode->prev->next = newnode->next;
    else
      buckets[gain] = newnode->next;
    if (newnode->next != NULL)
      newnode->next->prev = newnode->prev;

    if (buckets[gain] == NULL && gain == queue->maxgain) {
      if (queue->nnodes == 0)
        queue->maxgain = -queue->ngainspan;
      else
        for (; buckets[queue->maxgain]==NULL; queue->maxgain--);
    }
  }
  else { /* Heap Priority Queue */
    heap = queue->heap;
    locator = queue->locator;

    ASSERT(locator[node] != -1);
    ASSERT(heap[locator[node]].val == node);

    ASSERT(CheckHeap(queue));

    i = locator[node];
    locator[node] = -1;

    if (--queue->nnodes > 0 && heap[queue->nnodes].val != node) {
      node = heap[queue->nnodes].val;
      newgain = heap[queue->nnodes].key;
      oldgain = heap[i].key;

      if (oldgain < newgain) { /* Filter-up */
        while (i > 0) {
          j = (i-1)>>1;
          if (heap[j].key < newgain) {
            heap[i] = heap[j];
            locator[heap[i].val] = i;
            i = j;
          }
          else
            break;
        }
      }
      else { /* Filter down */
        while ((j=2*i+1) < queue->nnodes) {
          if (heap[j].key > newgain) {
            if (j+1 < queue->nnodes && heap[j+1].key > heap[j].key)
              j = j+1;
            heap[i] = heap[j];
            locator[heap[i].val] = i;
            i = j;
          }
          else if (j+1 < queue->nnodes && heap[j+1].key > newgain) {
            j = j+1;
            heap[i] = heap[j];
            locator[heap[i].val] = i;
            i = j;
          }
          else
            break;
        }
      }

      heap[i].key = newgain;
      heap[i].val = node;
      locator[node] = i;
    }

    ASSERT(CheckHeap(queue));
  }

  return 0;
}



/*************************************************************************
* This function deletes a node from a partition and reinserts it with
* an updated gain
**************************************************************************/
int PQueueUpdate(PQueueType *queue, int node, int oldgain, int newgain)
{
  int i, j;
  idxtype *locator;
  KeyValueType *heap;

  if (oldgain == newgain)
    return 0;

  if (queue->type == 1) {
    /* First delete the node and then insert it */
    PQueueDelete(queue, node, oldgain);
    return PQueueInsert(queue, node, newgain);
  }
  else { /* Heap Priority Queue */
    heap = queue->heap;
    locator = queue->locator;

    ASSERT(locator[node] != -1);
    ASSERT(heap[locator[node]].val == node);
    ASSERT(heap[locator[node]].key == oldgain);
    ASSERT(CheckHeap(queue));

    i = locator[node];

    if (oldgain < newgain) { /* Filter-up */
      while (i > 0) {
        j = (i-1)>>1;
        if (heap[j].key < newgain) {
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else
          break;
      }
    }
    else { /* Filter down */
      while ((j=2*i+1) < queue->nnodes) {
        if (heap[j].key > newgain) {
          if (j+1 < queue->nnodes && heap[j+1].key > heap[j].key)
            j = j+1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else if (j+1 < queue->nnodes && heap[j+1].key > newgain) {
          j = j+1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else
          break;
      }
    }

    heap[i].key = newgain;
    heap[i].val = node;
    locator[node] = i;

    ASSERT(CheckHeap(queue));
  }

  return 0;
}



/*************************************************************************
* This function deletes a node from a partition and reinserts it with
* an updated gain
**************************************************************************/
void PQueueUpdateUp(PQueueType *queue, int node, int oldgain, int newgain)
{
  int i, j;
  idxtype *locator;
  ListNodeType *newnode, **buckets;
  KeyValueType *heap;

  if (oldgain == newgain)
    return;

  if (queue->type == 1) {
    ASSERT(oldgain >= -queue->ngainspan && oldgain <= queue->pgainspan);
    ASSERT(newgain >= -queue->ngainspan && newgain <= queue->pgainspan);
    ASSERT(queue->nnodes > 0);

    buckets = queue->buckets;
    newnode = queue->nodes+node;

    /* First delete the node */
    if (newnode->prev != NULL)
      newnode->prev->next = newnode->next;
    else
      buckets[oldgain] = newnode->next;
    if (newnode->next != NULL)
      newnode->next->prev = newnode->prev;

    /* Attach this node in the doubly-linked list */
    newnode->next = buckets[newgain];
    newnode->prev = NULL;
    if (newnode->next != NULL)
      newnode->next->prev = newnode;
    buckets[newgain] = newnode;

    if (queue->maxgain < newgain)
      queue->maxgain = newgain;
  }
  else { /* Heap Priority Queue */
    heap = queue->heap;
    locator = queue->locator;

    ASSERT(locator[node] != -1);
    ASSERT(heap[locator[node]].val == node);
    ASSERT(heap[locator[node]].key == oldgain);
    ASSERT(CheckHeap(queue));


    /* Here we are just filtering up since the newgain is greater than the oldgain */
    i = locator[node];
    while (i > 0) {
      j = (i-1)>>1;
      if (heap[j].key < newgain) {
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
        break;
    }

    heap[i].key = newgain;
    heap[i].val = node;
    locator[node] = i;

    ASSERT(CheckHeap(queue));
  }

}


/*************************************************************************
* This function returns the vertex with the largest gain from a partition
* and removes the node from the bucket list
**************************************************************************/
int PQueueGetMax(PQueueType *queue)
{
  int vtx, i, j, gain, node;
  idxtype *locator;
  ListNodeType *tptr;
  KeyValueType *heap;

  if (queue->nnodes == 0)
    return -1;

  queue->nnodes--;

  if (queue->type == 1) {
    tptr = queue->buckets[queue->maxgain];
    queue->buckets[queue->maxgain] = tptr->next;
    if (tptr->next != NULL) {
      tptr->next->prev = NULL;
    }
    else {
      if (queue->nnodes == 0) {
        queue->maxgain = -queue->ngainspan;
      }
      else
        for (; queue->buckets[queue->maxgain]==NULL; queue->maxgain--);
    }

    return tptr->id;
  }
  else {
    heap = queue->heap;
    locator = queue->locator;

    vtx = heap[0].val;
    locator[vtx] = -1;

    if ((i = queue->nnodes) > 0) {
      gain = heap[i].key;
      node = heap[i].val;
      i = 0;
      while ((j=2*i+1) < queue->nnodes) {
        if (heap[j].key > gain) {
          if (j+1 < queue->nnodes && heap[j+1].key > heap[j].key)
            j = j+1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else if (j+1 < queue->nnodes && heap[j+1].key > gain) {
          j = j+1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else
          break;
      }

      heap[i].key = gain;
      heap[i].val = node;
      locator[node] = i;
    }

    ASSERT(CheckHeap(queue));
    return vtx;
  }
}


/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
int PQueueSeeMax(PQueueType *queue)
{
  int vtx;

  if (queue->nnodes == 0)
    return -1;

  if (queue->type == 1)
    vtx = queue->buckets[queue->maxgain]->id;
  else
    vtx = queue->heap[0].val;

  return vtx;
}


/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
int PQueueGetKey(PQueueType *queue)
{
  int key;

  if (queue->nnodes == 0)
    return -1;

  if (queue->type == 1)
    key = queue->maxgain;
  else
    key = queue->heap[0].key;

  return key;
}




/*************************************************************************
* This functions checks the consistency of the heap
**************************************************************************/
int CheckHeap(PQueueType *queue)
{
  int i, j, nnodes;
  idxtype *locator;
  KeyValueType *heap;

  heap = queue->heap;
  locator = queue->locator;
  nnodes = queue->nnodes;

  if (nnodes == 0)
    return 1;

  ASSERT(locator[heap[0].val] == 0);
  for (i=1; i<nnodes; i++) {
    ASSERTP(locator[heap[i].val] == i, ("%d %d %d %d\n", nnodes, i, heap[i].val, locator[heap[i].val]));
    ASSERTP(heap[i].key <= heap[(i-1)/2].key, ("%d %d %d %d %d\n", i, (i-1)/2, nnodes, heap[i].key, heap[(i-1)/2].key));
  }
  for (i=1; i<nnodes; i++)
    ASSERT(heap[i].key <= heap[0].key);

  for (j=i=0; i<queue->maxnodes; i++) {
    if (locator[i] != -1)
      j++;
  }
  ASSERTP(j == nnodes, ("%d %d\n", j, nnodes));

  return 1;
}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * refine.c
 *
 * This file contains the driving routines for multilevel refinement
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void Refine2Way(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int *tpwgts, float ubfactor)
{

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  Compute2WayPartitionParams(ctrl, graph);

  for (;;) {
    ASSERT(CheckBnd(graph));

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));
    switch (ctrl->RType) {
      case 1:
        Balance2Way(ctrl, graph, tpwgts, ubfactor);
        FM_2WayEdgeRefine(ctrl, graph, tpwgts, 8);
        break;
      default:
        errexit("Unknown refinement type: %d\n", ctrl->RType);
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    Project2WayPartition(ctrl, graph);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}


/*************************************************************************
* This function allocates memory for 2-way edge refinement
**************************************************************************/
void Allocate2WayPartitionMemory(CtrlType *ctrl, GraphType *graph)
{
  int nvtxs;

  nvtxs = graph->nvtxs;

  graph->rdata = idxmalloc(5*nvtxs+2, "Allocate2WayPartitionMemory: rdata");
  graph->pwgts 		= graph->rdata;
  graph->where		= graph->rdata + 2;
  graph->id		= graph->rdata + nvtxs + 2;
  graph->ed		= graph->rdata + 2*nvtxs + 2;
  graph->bndptr		= graph->rdata + 3*nvtxs + 2;
  graph->bndind		= graph->rdata + 4*nvtxs + 2;
}


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void Compute2WayPartitionParams(CtrlType *ctrl, GraphType *graph)
{
  int i, j, nvtxs, nbnd, mincut;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *pwgts;
  idxtype *id, *ed, *where;
  idxtype *bndptr, *bndind;
  int me;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = idxset(2, 0, graph->pwgts);
  id = idxset(nvtxs, 0, graph->id);
  ed = idxset(nvtxs, 0, graph->ed);
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;


  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  nbnd = mincut = 0;
  for (i=0; i<nvtxs; i++) {
    ASSERT(where[i] >= 0 && where[i] <= 1);
    me = where[i];
    pwgts[me] += vwgt[i];

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        id[i] += adjwgt[j];
      else
        ed[i] += adjwgt[j];
    }

    if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
      mincut += ed[i];
      bndptr[i] = nbnd;
      bndind[nbnd++] = i;
    }
  }

  graph->mincut = mincut/2;
  graph->nbnd = nbnd;

  ASSERT(pwgts[0]+pwgts[1] == idxsum(nvtxs, vwgt));
}



/*************************************************************************
* This function projects a partition, and at the same time computes the
* parameters for refinement.
**************************************************************************/
void Project2WayPartition(CtrlType *ctrl, GraphType *graph)
{
  int i, j, k, nvtxs, nbnd, me;
  idxtype *xadj, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *cmap, *where, *id, *ed, *bndptr, *bndind;
  idxtype *cwhere, *cid, *ced, *cbndptr;
  GraphType *cgraph;

  cgraph = graph->coarser;
  cwhere = cgraph->where;
  cid = cgraph->id;
  ced = cgraph->ed;
  cbndptr = cgraph->bndptr;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  Allocate2WayPartitionMemory(ctrl, graph);

  where = graph->where;
  id = idxset(nvtxs, 0, graph->id);
  ed = idxset(nvtxs, 0, graph->ed);
  bndptr = idxset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;


  /* Go through and project partition and compute id/ed for the nodes */
  for (i=0; i<nvtxs; i++) {
    k = cmap[i];
    where[i] = cwhere[k];
    cmap[i] = cbndptr[k];
  }

  for (nbnd=0, i=0; i<nvtxs; i++) {
    me = where[i];

    id[i] = adjwgtsum[i];

    if (xadj[i] == xadj[i+1]) {
      bndptr[i] = nbnd;
      bndind[nbnd++] = i;
    }
    else {
      if (cmap[i] != -1) { /* If it is an interface node. Note that cmap[i] = cbndptr[cmap[i]] */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          if (me != where[adjncy[j]])
            ed[i] += adjwgt[j];
        }
        id[i] -= ed[i];

        if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
          bndptr[i] = nbnd;
          bndind[nbnd++] = i;
        }
      }
    }
  }

  graph->mincut = cgraph->mincut;
  graph->nbnd = nbnd;
  idxcopy(2, cgraph->pwgts, graph->pwgts);

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * separator.c
 *
 * This file contains code for separator extraction
 *
 * Started 8/1/97
 * George
 *
 * $Id$
 *
 */



/*************************************************************************
* This function takes a bisection and constructs a minimum weight vertex
* separator out of it. It uses the node-based separator refinement for it.
**************************************************************************/
void ConstructSeparator(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  int i, j, nvtxs, nbnd;
  idxtype *xadj, *where, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  nbnd = graph->nbnd;
  bndind = graph->bndind;

  where = idxcopy(nvtxs, graph->where, idxwspacemalloc(ctrl, nvtxs));

  /* Put the nodes in the boundary into the separator */
  for (i=0; i<nbnd; i++) {
    j = bndind[i];
    if (xadj[j+1]-xadj[j] > 0)  /* Ignore islands */
      where[j] = 2;
  }

  GKfree((void **) &graph->rdata, LTERM);
  Allocate2WayNodePartitionMemory(ctrl, graph);
  idxcopy(nvtxs, where, graph->where);
  idxwspacefree(ctrl, nvtxs);

  ASSERT(IsSeparable(graph));

  Compute2WayNodePartitionParams(ctrl, graph);

  ASSERT(CheckNodePartitionParams(graph));

  FM_2WayNodeRefine(ctrl, graph, ubfactor, 8);

  ASSERT(IsSeparable(graph));
}



/*************************************************************************
* This function takes a bisection and constructs a minimum weight vertex
* separator out of it. It uses an unweighted minimum-cover algorithm
* followed by node-based separator refinement.
**************************************************************************/
void ConstructMinCoverSeparator0(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  int i, ii, j, jj, k, l, nvtxs, nbnd, bnvtxs[3], bnedges[2], csize;
  idxtype *xadj, *adjncy, *bxadj, *badjncy;
  idxtype *where, *bndind, *bndptr, *vmap, *ivmap, *cover;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  nbnd = graph->nbnd;
  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;

  vmap = idxwspacemalloc(ctrl, nvtxs);
  ivmap = idxwspacemalloc(ctrl, nbnd);
  cover = idxwspacemalloc(ctrl, nbnd);

  if (nbnd > 0) {
    /* Go through the boundary and determine the sizes of the bipartite graph */
    bnvtxs[0] = bnvtxs[1] = bnedges[0] = bnedges[1] = 0;
    for (i=0; i<nbnd; i++) {
      j = bndind[i];
      k = where[j];
      if (xadj[j+1]-xadj[j] > 0) {
        bnvtxs[k]++;
        bnedges[k] += xadj[j+1]-xadj[j];
      }
    }

    bnvtxs[2] = bnvtxs[0]+bnvtxs[1];
    bnvtxs[1] = bnvtxs[0];
    bnvtxs[0] = 0;

    bxadj = idxmalloc(bnvtxs[2]+1, "ConstructMinCoverSeparator: bxadj");
    badjncy = idxmalloc(bnedges[0]+bnedges[1]+1, "ConstructMinCoverSeparator: badjncy");

    /* Construct the ivmap and vmap */
    ASSERT(idxset(nvtxs, -1, vmap) == vmap);
    for (i=0; i<nbnd; i++) {
      j = bndind[i];
      k = where[j];
      if (xadj[j+1]-xadj[j] > 0) {
        vmap[j] = bnvtxs[k];
        ivmap[bnvtxs[k]++] = j;
      }
    }

    /* OK, go through and put the vertices of each part starting from 0 */
    bnvtxs[1] = bnvtxs[0];
    bnvtxs[0] = 0;
    bxadj[0] = l = 0;
    for (k=0; k<2; k++) {
      for (ii=0; ii<nbnd; ii++) {
        i = bndind[ii];
        if (where[i] == k && xadj[i] < xadj[i+1]) {
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            jj = adjncy[j];
            if (where[jj] != k) {
              ASSERT(bndptr[jj] != -1);
              ASSERTP(vmap[jj] != -1, ("%d %d %d\n", jj, vmap[jj], graph->bndptr[jj]));
              badjncy[l++] = vmap[jj];
            }
          }
          bxadj[++bnvtxs[k]] = l;
        }
      }
    }

    ASSERT(l <= bnedges[0]+bnedges[1]);

    MinCover(bxadj, badjncy, bnvtxs[0], bnvtxs[1], cover, &csize);

    IFSET(ctrl->dbglvl, DBG_SEPINFO,
      printf("Nvtxs: %6d, [%5d %5d], Cut: %6d, SS: [%6d %6d], Cover: %6d\n", nvtxs, graph->pwgts[0], graph->pwgts[1], graph->mincut, bnvtxs[0], bnvtxs[1]-bnvtxs[0], csize));

    for (i=0; i<csize; i++) {
      j = ivmap[cover[i]];
      where[j] = 2;
    }

    GKfree((void **) &bxadj, (void **) &badjncy, LTERM);

    for (i=0; i<nbnd; i++)
      bndptr[bndind[i]] = -1;
    for (nbnd=i=0; i<nvtxs; i++) {
      if (where[i] == 2) {
        bndind[nbnd] = i;
        bndptr[i] = nbnd++;
      }
    }
  }
  else {
    IFSET(ctrl->dbglvl, DBG_SEPINFO,
      printf("Nvtxs: %6d, [%5d %5d], Cut: %6d, SS: [%6d %6d], Cover: %6d\n", nvtxs, graph->pwgts[0], graph->pwgts[1], graph->mincut, 0, 0, 0));
  }

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, graph->nbnd);
  idxwspacefree(ctrl, graph->nbnd);
  graph->nbnd = nbnd;


  ASSERT(IsSeparable(graph));
}



/*************************************************************************
* This function takes a bisection and constructs a minimum weight vertex
* separator out of it. It uses an unweighted minimum-cover algorithm
* followed by node-based separator refinement.
**************************************************************************/
void ConstructMinCoverSeparator(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  int i, ii, j, jj, k, l, nvtxs, nbnd, bnvtxs[3], bnedges[2], csize;
  idxtype *xadj, *adjncy, *bxadj, *badjncy;
  idxtype *where, *bndind, *bndptr, *vmap, *ivmap, *cover;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  nbnd = graph->nbnd;
  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;

  vmap = idxwspacemalloc(ctrl, nvtxs);
  ivmap = idxwspacemalloc(ctrl, nbnd);
  cover = idxwspacemalloc(ctrl, nbnd);

  if (nbnd > 0) {
    /* Go through the boundary and determine the sizes of the bipartite graph */
    bnvtxs[0] = bnvtxs[1] = bnedges[0] = bnedges[1] = 0;
    for (i=0; i<nbnd; i++) {
      j = bndind[i];
      k = where[j];
      if (xadj[j+1]-xadj[j] > 0) {
        bnvtxs[k]++;
        bnedges[k] += xadj[j+1]-xadj[j];
      }
    }

    bnvtxs[2] = bnvtxs[0]+bnvtxs[1];
    bnvtxs[1] = bnvtxs[0];
    bnvtxs[0] = 0;

    bxadj = idxmalloc(bnvtxs[2]+1, "ConstructMinCoverSeparator: bxadj");
    badjncy = idxmalloc(bnedges[0]+bnedges[1]+1, "ConstructMinCoverSeparator: badjncy");

    /* Construct the ivmap and vmap */
    ASSERT(idxset(nvtxs, -1, vmap) == vmap);
    for (i=0; i<nbnd; i++) {
      j = bndind[i];
      k = where[j];
      if (xadj[j+1]-xadj[j] > 0) {
        vmap[j] = bnvtxs[k];
        ivmap[bnvtxs[k]++] = j;
      }
    }

    /* OK, go through and put the vertices of each part starting from 0 */
    bnvtxs[1] = bnvtxs[0];
    bnvtxs[0] = 0;
    bxadj[0] = l = 0;
    for (k=0; k<2; k++) {
      for (ii=0; ii<nbnd; ii++) {
        i = bndind[ii];
        if (where[i] == k && xadj[i] < xadj[i+1]) {
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            jj = adjncy[j];
            if (where[jj] != k) {
              ASSERT(bndptr[jj] != -1);
              ASSERTP(vmap[jj] != -1, ("%d %d %d\n", jj, vmap[jj], graph->bndptr[jj]));
              badjncy[l++] = vmap[jj];
            }
          }
          bxadj[++bnvtxs[k]] = l;
        }
      }
    }

    ASSERT(l <= bnedges[0]+bnedges[1]);

    MinCover(bxadj, badjncy, bnvtxs[0], bnvtxs[1], cover, &csize);

    IFSET(ctrl->dbglvl, DBG_SEPINFO,
      printf("Nvtxs: %6d, [%5d %5d], Cut: %6d, SS: [%6d %6d], Cover: %6d\n", nvtxs, graph->pwgts[0], graph->pwgts[1], graph->mincut, bnvtxs[0], bnvtxs[1]-bnvtxs[0], csize));

    for (i=0; i<csize; i++) {
      j = ivmap[cover[i]];
      where[j] = 2;
    }

    GKfree((void **) &bxadj, (void **) &badjncy, LTERM);
  }
  else {
    IFSET(ctrl->dbglvl, DBG_SEPINFO,
      printf("Nvtxs: %6d, [%5d %5d], Cut: %6d, SS: [%6d %6d], Cover: %6d\n", nvtxs, graph->pwgts[0], graph->pwgts[1], graph->mincut, 0, 0, 0));
  }

  /* Prepare to refine the vertex separator */
  idxcopy(nvtxs, graph->where, vmap);
  GKfree((void **) &graph->rdata, LTERM);

  Allocate2WayNodePartitionMemory(ctrl, graph);
  idxcopy(nvtxs, vmap, graph->where);
  idxwspacefree(ctrl, nvtxs+2*graph->nbnd);

  Compute2WayNodePartitionParams(ctrl, graph);

  ASSERT(CheckNodePartitionParams(graph));

  FM_2WayNodeRefine_OneSided(ctrl, graph, ubfactor, 6);

  ASSERT(IsSeparable(graph));
}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * sfm.c
 *
 * This file contains code that implementes an FM-based separator refinement
 *
 * Started 8/1/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function performs a node-based FM refinement
**************************************************************************/
void FM_2WayNodeRefine(CtrlType *ctrl, GraphType *graph, float ubfactor, int npasses)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *mptr, *mind, *moved, *swaps, *perm;
  PQueueType parts[2];
  NRInfoType *rinfo;
  int higain, oldgain, mincut, initcut, mincutorder;
  int pass, to, other, limit;
  int badmaxpwgt, mindiff, newdiff;
  int u[2], g[2];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;
  pwgts = graph->pwgts;
  rinfo = graph->nrinfo;


  i = ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt);
  PQueueInit(ctrl, &parts[0], nvtxs, i);
  PQueueInit(ctrl, &parts[1], nvtxs, i);

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  mptr = idxwspacemalloc(ctrl, nvtxs+1);
  mind = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions: [%6d %6d] Nv-Nb[%6d %6d]. ISep: %6d\n", pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  badmaxpwgt = (int)(ubfactor*(pwgts[0]+pwgts[1]+pwgts[2])/2);

  for (pass=0; pass<npasses; pass++) {
    idxset(nvtxs, -1, moved);
    PQueueReset(&parts[0]);
    PQueueReset(&parts[1]);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(where[i] == 2);
      PQueueInsert(&parts[0], i, vwgt[i]-rinfo[i].edegrees[1]);
      PQueueInsert(&parts[1], i, vwgt[i]-rinfo[i].edegrees[0]);
    }

    ASSERT(CheckNodeBnd(graph, nbnd));
    ASSERT(CheckNodePartitionParams(graph));

    limit = (ctrl->oflags&OFLAG_COMPRESS ? amin(5*nbnd, 400) : amin(2*nbnd, 300));

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = 0;
    mindiff = abs(pwgts[0]-pwgts[1]);
    to = (pwgts[0] < pwgts[1] ? 0 : 1);
    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      u[0] = PQueueSeeMax(&parts[0]);
      u[1] = PQueueSeeMax(&parts[1]);
      if (u[0] != -1 && u[1] != -1) {
        g[0] = vwgt[u[0]]-rinfo[u[0]].edegrees[1];
        g[1] = vwgt[u[1]]-rinfo[u[1]].edegrees[0];

        to = (g[0] > g[1] ? 0 : (g[0] < g[1] ? 1 : pass%2));
        /* to = (g[0] > g[1] ? 0 : (g[0] < g[1] ? 1 : (pwgts[0] < pwgts[1] ? 0 : 1))); */

        if (pwgts[to]+vwgt[u[to]] > badmaxpwgt)
          to = (to+1)%2;
      }
      else if (u[0] == -1 && u[1] == -1) {
        break;
      }
      else if (u[0] != -1 && pwgts[0]+vwgt[u[0]] <= badmaxpwgt) {
        to = 0;
      }
      else if (u[1] != -1 && pwgts[1]+vwgt[u[1]] <= badmaxpwgt) {
        to = 1;
      }
      else
        break;

      other = (to+1)%2;

      higain = PQueueGetMax(&parts[to]);
      if (moved[higain] == -1) /* Delete if it was in the separator originally */
        PQueueDelete(&parts[other], higain, vwgt[higain]-rinfo[higain].edegrees[to]);

      ASSERT(bndptr[higain] != -1);

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

      newdiff = abs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
      if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff)) {
        mincut = pwgts[2];
        mincutorder = nswaps;
        mindiff = newdiff;
      }
      else {
        if (nswaps - mincutorder > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[other]);
          break; /* No further improvement, break out */
        }
      }

      BNDDelete(nbnd, bndind, bndptr, higain);
      pwgts[to] += vwgt[higain];
      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;


      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          oldgain = vwgt[k]-rinfo[k].edegrees[to];
          rinfo[k].edegrees[to] += vwgt[higain];
          if (moved[k] == -1 || moved[k] == -(2+other))
            PQueueUpdate(&parts[other], k, oldgain, oldgain-vwgt[higain]);
        }
        else if (where[k] == other) { /* This vertex is pulled into the separator */
          ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
          BNDInsert(nbnd, bndind, bndptr, k);

          mind[nmind++] = k;  /* Keep track for rollback */
          where[k] = 2;
          pwgts[other] -= vwgt[k];

          edegrees = rinfo[k].edegrees;
          edegrees[0] = edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            if (where[kk] != 2)
              edegrees[where[kk]] += vwgt[kk];
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
              rinfo[kk].edegrees[other] -= vwgt[k];
              if (moved[kk] == -1 || moved[kk] == -(2+to))
                PQueueUpdate(&parts[to], kk, oldgain, oldgain+vwgt[k]);
            }
          }

          /* Insert the new vertex into the priority queue. Only one side! */
          if (moved[k] == -1) {
            PQueueInsert(&parts[to], k, vwgt[k]-edegrees[other]);
            moved[k] = -(2+to);
          }
        }
      }
      mptr[nswaps+1] = nmind;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moved %6d to %3d, Gain: %5d [%5d] [%4d %4d] \t[%5d %5d %5d]\n", higain, to, g[to], g[other], vwgt[u[to]], vwgt[u[other]], pwgts[0], pwgts[1], pwgts[2]));

    }


    /****************************************************************
    * Roll back computation
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      ASSERT(CheckNodePartitionParams(graph));

      to = where[higain];
      other = (to+1)%2;
      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;
      BNDInsert(nbnd, bndind, bndptr, higain);

      edegrees = rinfo[higain].edegrees;
      edegrees[0] = edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2)
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERT(where[k] == 2);
        where[k] = other;
        INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
        BNDDelete(nbnd, bndind, bndptr, k);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2)
            rinfo[kk].edegrees[other] += vwgt[k];
        }
      }
    }

    ASSERT(mincut == pwgts[2]);

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum sep: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, mincutorder, pwgts[0], pwgts[1], nbnd));

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut >= initcut)
      break;
  }

  PQueueFree(ctrl, &parts[0]);
  PQueueFree(ctrl, &parts[1]);

  idxwspacefree(ctrl, nvtxs+1);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function performs a node-based FM refinement
**************************************************************************/
void FM_2WayNodeRefine2(CtrlType *ctrl, GraphType *graph, float ubfactor, int npasses)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *mptr, *mind, *moved, *swaps, *perm;
  PQueueType parts[2];
  NRInfoType *rinfo;
  int higain, oldgain, mincut, initcut, mincutorder;
  int pass, to, other, limit;
  int badmaxpwgt, mindiff, newdiff;
  int u[2], g[2];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;
  pwgts = graph->pwgts;
  rinfo = graph->nrinfo;


  i = ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt);
  PQueueInit(ctrl, &parts[0], nvtxs, i);
  PQueueInit(ctrl, &parts[1], nvtxs, i);

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  mptr = idxwspacemalloc(ctrl, nvtxs+1);
  mind = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions: [%6d %6d] Nv-Nb[%6d %6d]. ISep: %6d\n", pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  badmaxpwgt = (int)(ubfactor*(pwgts[0]+pwgts[1]+pwgts[2])/2);

  for (pass=0; pass<npasses; pass++) {
    idxset(nvtxs, -1, moved);
    PQueueReset(&parts[0]);
    PQueueReset(&parts[1]);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(where[i] == 2);
      PQueueInsert(&parts[0], i, vwgt[i]-rinfo[i].edegrees[1]);
      PQueueInsert(&parts[1], i, vwgt[i]-rinfo[i].edegrees[0]);
    }

    ASSERT(CheckNodeBnd(graph, nbnd));
    ASSERT(CheckNodePartitionParams(graph));

    limit = (ctrl->oflags&OFLAG_COMPRESS ? amin(5*nbnd, 400) : amin(2*nbnd, 300));

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = 0;
    mindiff = abs(pwgts[0]-pwgts[1]);
    to = (pwgts[0] < pwgts[1] ? 0 : 1);
    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      badmaxpwgt = (int)(ubfactor*(pwgts[0]+pwgts[1]+pwgts[2]/2)/2);

      u[0] = PQueueSeeMax(&parts[0]);
      u[1] = PQueueSeeMax(&parts[1]);
      if (u[0] != -1 && u[1] != -1) {
        g[0] = vwgt[u[0]]-rinfo[u[0]].edegrees[1];
        g[1] = vwgt[u[1]]-rinfo[u[1]].edegrees[0];

        to = (g[0] > g[1] ? 0 : (g[0] < g[1] ? 1 : pass%2));
        /* to = (g[0] > g[1] ? 0 : (g[0] < g[1] ? 1 : (pwgts[0] < pwgts[1] ? 0 : 1))); */

        if (pwgts[to]+vwgt[u[to]] > badmaxpwgt)
          to = (to+1)%2;
      }
      else if (u[0] == -1 && u[1] == -1) {
        break;
      }
      else if (u[0] != -1 && pwgts[0]+vwgt[u[0]] <= badmaxpwgt) {
        to = 0;
      }
      else if (u[1] != -1 && pwgts[1]+vwgt[u[1]] <= badmaxpwgt) {
        to = 1;
      }
      else
        break;

      other = (to+1)%2;

      higain = PQueueGetMax(&parts[to]);
      if (moved[higain] == -1) /* Delete if it was in the separator originally */
        PQueueDelete(&parts[other], higain, vwgt[higain]-rinfo[higain].edegrees[to]);

      ASSERT(bndptr[higain] != -1);

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

      newdiff = abs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
      if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff)) {
        mincut = pwgts[2];
        mincutorder = nswaps;
        mindiff = newdiff;
      }
      else {
        if (nswaps - mincutorder > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[other]);
          break; /* No further improvement, break out */
        }
      }

      BNDDelete(nbnd, bndind, bndptr, higain);
      pwgts[to] += vwgt[higain];
      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;


      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          oldgain = vwgt[k]-rinfo[k].edegrees[to];
          rinfo[k].edegrees[to] += vwgt[higain];
          if (moved[k] == -1 || moved[k] == -(2+other))
            PQueueUpdate(&parts[other], k, oldgain, oldgain-vwgt[higain]);
        }
        else if (where[k] == other) { /* This vertex is pulled into the separator */
          ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
          BNDInsert(nbnd, bndind, bndptr, k);

          mind[nmind++] = k;  /* Keep track for rollback */
          where[k] = 2;
          pwgts[other] -= vwgt[k];

          edegrees = rinfo[k].edegrees;
          edegrees[0] = edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            if (where[kk] != 2)
              edegrees[where[kk]] += vwgt[kk];
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
              rinfo[kk].edegrees[other] -= vwgt[k];
              if (moved[kk] == -1 || moved[kk] == -(2+to))
                PQueueUpdate(&parts[to], kk, oldgain, oldgain+vwgt[k]);
            }
          }

          /* Insert the new vertex into the priority queue. Only one side! */
          if (moved[k] == -1) {
            PQueueInsert(&parts[to], k, vwgt[k]-edegrees[other]);
            moved[k] = -(2+to);
          }
        }
      }
      mptr[nswaps+1] = nmind;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moved %6d to %3d, Gain: %5d [%5d] [%4d %4d] \t[%5d %5d %5d]\n", higain, to, g[to], g[other], vwgt[u[to]], vwgt[u[other]], pwgts[0], pwgts[1], pwgts[2]));

    }


    /****************************************************************
    * Roll back computation
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      ASSERT(CheckNodePartitionParams(graph));

      to = where[higain];
      other = (to+1)%2;
      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;
      BNDInsert(nbnd, bndind, bndptr, higain);

      edegrees = rinfo[higain].edegrees;
      edegrees[0] = edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2)
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERT(where[k] == 2);
        where[k] = other;
        INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
        BNDDelete(nbnd, bndind, bndptr, k);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2)
            rinfo[kk].edegrees[other] += vwgt[k];
        }
      }
    }

    ASSERT(mincut == pwgts[2]);

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum sep: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, mincutorder, pwgts[0], pwgts[1], nbnd));

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut >= initcut)
      break;
  }

  PQueueFree(ctrl, &parts[0]);
  PQueueFree(ctrl, &parts[1]);

  idxwspacefree(ctrl, nvtxs+1);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function performs a node-based FM refinement
**************************************************************************/
void FM_2WayNodeRefineEqWgt(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *mptr, *mind, *moved, *swaps, *perm;
  PQueueType parts[2];
  NRInfoType *rinfo;
  int higain, oldgain, mincut, initcut, mincutorder;
  int pass, to, other, limit;
  int mindiff, newdiff;
  int u[2], g[2];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;
  pwgts = graph->pwgts;
  rinfo = graph->nrinfo;


  i = ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt);
  PQueueInit(ctrl, &parts[0], nvtxs, i);
  PQueueInit(ctrl, &parts[1], nvtxs, i);

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  mptr = idxwspacemalloc(ctrl, nvtxs+1);
  mind = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions: [%6d %6d] Nv-Nb[%6d %6d]. ISep: %6d\n", pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  for (pass=0; pass<npasses; pass++) {
    idxset(nvtxs, -1, moved);
    PQueueReset(&parts[0]);
    PQueueReset(&parts[1]);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(where[i] == 2);
      PQueueInsert(&parts[0], i, vwgt[i]-rinfo[i].edegrees[1]);
      PQueueInsert(&parts[1], i, vwgt[i]-rinfo[i].edegrees[0]);
    }

    ASSERT(CheckNodeBnd(graph, nbnd));
    ASSERT(CheckNodePartitionParams(graph));

    limit = (ctrl->oflags&OFLAG_COMPRESS ? amin(5*nbnd, 400) : amin(2*nbnd, 300));

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = 0;
    mindiff = abs(pwgts[0]-pwgts[1]);
    to = (pwgts[0] < pwgts[1] ? 0 : 1);
    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      to = (pwgts[0] < pwgts[1] ? 0 : 1);

      if (pwgts[0] == pwgts[1]) {
        u[0] = PQueueSeeMax(&parts[0]);
        u[1] = PQueueSeeMax(&parts[1]);
        if (u[0] != -1 && u[1] != -1) {
          g[0] = vwgt[u[0]]-rinfo[u[0]].edegrees[1];
          g[1] = vwgt[u[1]]-rinfo[u[1]].edegrees[0];

          to = (g[0] > g[1] ? 0 : (g[0] < g[1] ? 1 : pass%2));
        }
      }
      other = (to+1)%2;

      if ((higain = PQueueGetMax(&parts[to])) == -1)
        break;

      if (moved[higain] == -1) /* Delete if it was in the separator originally */
        PQueueDelete(&parts[other], higain, vwgt[higain]-rinfo[higain].edegrees[to]);

      ASSERT(bndptr[higain] != -1);

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

      newdiff = abs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
      if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff)) {
        mincut = pwgts[2];
        mincutorder = nswaps;
        mindiff = newdiff;
      }
      else {
        if (nswaps - mincutorder > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[other]);
          break; /* No further improvement, break out */
        }
      }

      BNDDelete(nbnd, bndind, bndptr, higain);
      pwgts[to] += vwgt[higain];
      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;


      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          oldgain = vwgt[k]-rinfo[k].edegrees[to];
          rinfo[k].edegrees[to] += vwgt[higain];
          if (moved[k] == -1 || moved[k] == -(2+other))
            PQueueUpdate(&parts[other], k, oldgain, oldgain-vwgt[higain]);
        }
        else if (where[k] == other) { /* This vertex is pulled into the separator */
          ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
          BNDInsert(nbnd, bndind, bndptr, k);

          mind[nmind++] = k;  /* Keep track for rollback */
          where[k] = 2;
          pwgts[other] -= vwgt[k];

          edegrees = rinfo[k].edegrees;
          edegrees[0] = edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            if (where[kk] != 2)
              edegrees[where[kk]] += vwgt[kk];
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
              rinfo[kk].edegrees[other] -= vwgt[k];
              if (moved[kk] == -1 || moved[kk] == -(2+to))
                PQueueUpdate(&parts[to], kk, oldgain, oldgain+vwgt[k]);
            }
          }

          /* Insert the new vertex into the priority queue. Only one side! */
          if (moved[k] == -1) {
            PQueueInsert(&parts[to], k, vwgt[k]-edegrees[other]);
            moved[k] = -(2+to);
          }
        }
      }
      mptr[nswaps+1] = nmind;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moved %6d to %3d, Gain: %5d [%5d] [%4d %4d] \t[%5d %5d %5d]\n", higain, to, g[to], g[other], vwgt[u[to]], vwgt[u[other]], pwgts[0], pwgts[1], pwgts[2]));

    }


    /****************************************************************
    * Roll back computation
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      ASSERT(CheckNodePartitionParams(graph));

      to = where[higain];
      other = (to+1)%2;
      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;
      BNDInsert(nbnd, bndind, bndptr, higain);

      edegrees = rinfo[higain].edegrees;
      edegrees[0] = edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2)
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERT(where[k] == 2);
        where[k] = other;
        INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
        BNDDelete(nbnd, bndind, bndptr, k);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2)
            rinfo[kk].edegrees[other] += vwgt[k];
        }
      }
    }

    ASSERT(mincut == pwgts[2]);

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum sep: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, mincutorder, pwgts[0], pwgts[1], nbnd));

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut >= initcut)
      break;
  }

  PQueueFree(ctrl, &parts[0]);
  PQueueFree(ctrl, &parts[1]);

  idxwspacefree(ctrl, nvtxs+1);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function performs a node-based FM refinement. This is the
* one-way version
**************************************************************************/
void FM_2WayNodeRefine_OneSided(CtrlType *ctrl, GraphType *graph, float ubfactor, int npasses)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *mptr, *mind, *swaps, *perm;
  PQueueType parts;
  NRInfoType *rinfo;
  int higain, oldgain, mincut, initcut, mincutorder;
  int pass, to, other, limit;
  int badmaxpwgt, mindiff, newdiff;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;
  pwgts = graph->pwgts;
  rinfo = graph->nrinfo;

  PQueueInit(ctrl, &parts, nvtxs, ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt));

  perm = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  mptr = idxwspacemalloc(ctrl, nvtxs);
  mind = idxwspacemalloc(ctrl, nvtxs+1);

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions-N1: [%6d %6d] Nv-Nb[%6d %6d]. ISep: %6d\n", pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  badmaxpwgt = (int)(ubfactor*(pwgts[0]+pwgts[1]+pwgts[2])/2);

  to = (pwgts[0] < pwgts[1] ? 1 : 0);
  for (pass=0; pass<npasses; pass++) {
    other = to;
    to = (to+1)%2;

    PQueueReset(&parts);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(where[i] == 2);
      PQueueInsert(&parts, i, vwgt[i]-rinfo[i].edegrees[other]);
    }

    ASSERT(CheckNodeBnd(graph, nbnd));
    ASSERT(CheckNodePartitionParams(graph));

    limit = (ctrl->oflags&OFLAG_COMPRESS ? amin(5*nbnd, 400) : amin(2*nbnd, 300));

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = 0;
    mindiff = abs(pwgts[0]-pwgts[1]);
    for (nswaps=0; nswaps<nvtxs; nswaps++) {

      if ((higain = PQueueGetMax(&parts)) == -1)
        break;

      ASSERT(bndptr[higain] != -1);

      if (pwgts[to]+vwgt[higain] > badmaxpwgt)
        break;  /* No point going any further. Balance will be bad */

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

      newdiff = abs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
      if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff)) {
        mincut = pwgts[2];
        mincutorder = nswaps;
        mindiff = newdiff;
      }
      else {
        if (nswaps - mincutorder > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[other]);
          break; /* No further improvement, break out */
        }
      }

      BNDDelete(nbnd, bndind, bndptr, higain);
      pwgts[to] += vwgt[higain];
      where[higain] = to;
      swaps[nswaps] = higain;


      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          rinfo[k].edegrees[to] += vwgt[higain];
        }
        else if (where[k] == other) { /* This vertex is pulled into the separator */
          ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
          BNDInsert(nbnd, bndind, bndptr, k);

          mind[nmind++] = k;  /* Keep track for rollback */
          where[k] = 2;
          pwgts[other] -= vwgt[k];

          edegrees = rinfo[k].edegrees;
          edegrees[0] = edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            if (where[kk] != 2)
              edegrees[where[kk]] += vwgt[kk];
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
              rinfo[kk].edegrees[other] -= vwgt[k];

              /* Since the moves are one-sided this vertex has not been moved yet */
              PQueueUpdateUp(&parts, kk, oldgain, oldgain+vwgt[k]);
            }
          }

          /* Insert the new vertex into the priority queue. Safe due to one-sided moves */
          PQueueInsert(&parts, k, vwgt[k]-edegrees[other]);
        }
      }
      mptr[nswaps+1] = nmind;


      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moved %6d to %3d, Gain: %5d [%5d] \t[%5d %5d %5d] [%3d %2d]\n",
                       higain, to, (vwgt[higain]-rinfo[higain].edegrees[other]), vwgt[higain], pwgts[0], pwgts[1], pwgts[2], nswaps, limit));

    }


    /****************************************************************
    * Roll back computation
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      ASSERT(CheckNodePartitionParams(graph));
      ASSERT(where[higain] == to);

      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;
      BNDInsert(nbnd, bndind, bndptr, higain);

      edegrees = rinfo[higain].edegrees;
      edegrees[0] = edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2)
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERT(where[k] == 2);
        where[k] = other;
        INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
        BNDDelete(nbnd, bndind, bndptr, k);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2)
            rinfo[kk].edegrees[other] += vwgt[k];
        }
      }
    }

    ASSERT(mincut == pwgts[2]);

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum sep: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, mincutorder, pwgts[0], pwgts[1], nbnd));

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (pass%2 == 1 && (mincutorder == -1 || mincut >= initcut))
      break;
  }

  PQueueFree(ctrl, &parts);

  idxwspacefree(ctrl, nvtxs+1);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function performs a node-based FM refinement
**************************************************************************/
void FM_2WayNodeBalance(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *perm, *moved;
  PQueueType parts;
  NRInfoType *rinfo;
  int higain, oldgain;
  int to, other;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where = graph->where;
  pwgts = graph->pwgts;
  rinfo = graph->nrinfo;

  if (abs(pwgts[0]-pwgts[1]) < (int)((ubfactor-1.0)*(pwgts[0]+pwgts[1])))
    return;
  if (abs(pwgts[0]-pwgts[1]) < 3*idxsum(nvtxs, vwgt)/nvtxs)
    return;

  to = (pwgts[0] < pwgts[1] ? 0 : 1);
  other = (to+1)%2;

  PQueueInit(ctrl, &parts, nvtxs, ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt));

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxset(nvtxs, -1, idxwspacemalloc(ctrl, nvtxs));

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions: [%6d %6d] Nv-Nb[%6d %6d]. ISep: %6d [B]\n", pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  nbnd = graph->nbnd;
  RandomPermute(nbnd, perm, 1);
  for (ii=0; ii<nbnd; ii++) {
    i = bndind[perm[ii]];
    ASSERT(where[i] == 2);
    PQueueInsert(&parts, i, vwgt[i]-rinfo[i].edegrees[other]);
  }

  ASSERT(CheckNodeBnd(graph, nbnd));
  ASSERT(CheckNodePartitionParams(graph));

  /******************************************************
  * Get into the FM loop
  *******************************************************/
  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if ((higain = PQueueGetMax(&parts)) == -1)
      break;

    moved[higain] = 1;

    if (pwgts[other] - rinfo[higain].edegrees[other] < (pwgts[0]+pwgts[1])/2)
      continue;
#ifdef XXX
    if (pwgts[other] - rinfo[higain].edegrees[other] < pwgts[to]+vwgt[higain])
      break;
#endif

    ASSERT(bndptr[higain] != -1);

    pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

    BNDDelete(nbnd, bndind, bndptr, higain);
    pwgts[to] += vwgt[higain];
    where[higain] = to;

    IFSET(ctrl->dbglvl, DBG_MOVEINFO,
          printf("Moved %6d to %3d, Gain: %3d, \t[%5d %5d %5d]\n", higain, to, vwgt[higain]-rinfo[higain].edegrees[other], pwgts[0], pwgts[1], pwgts[2]));


    /**********************************************************
    * Update the degrees of the affected nodes
    ***********************************************************/
    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
        rinfo[k].edegrees[to] += vwgt[higain];
      }
      else if (where[k] == other) { /* This vertex is pulled into the separator */
        ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
        BNDInsert(nbnd, bndind, bndptr, k);

        where[k] = 2;
        pwgts[other] -= vwgt[k];

        edegrees = rinfo[k].edegrees;
        edegrees[0] = edegrees[1] = 0;
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] != 2)
            edegrees[where[kk]] += vwgt[kk];
          else {
            ASSERT(bndptr[kk] != -1);
            oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
            rinfo[kk].edegrees[other] -= vwgt[k];

            if (moved[kk] == -1)
              PQueueUpdateUp(&parts, kk, oldgain, oldgain+vwgt[k]);
          }
        }

        /* Insert the new vertex into the priority queue */
        PQueueInsert(&parts, k, vwgt[k]-edegrees[other]);
      }
    }

    if (pwgts[to] > pwgts[other])
      break;
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("\tBalanced sep: %6d at %4d, PWGTS: [%6d %6d], NBND: %6d\n", pwgts[2], nswaps, pwgts[0], pwgts[1], nbnd));

  graph->mincut = pwgts[2];
  graph->nbnd = nbnd;


  PQueueFree(ctrl, &parts);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function computes the maximum possible gain for a vertex
**************************************************************************/
int ComputeMaxNodeGain(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt)
{
  int i, j, k, max;

  max = 0;
  for (j=xadj[0]; j<xadj[1]; j++)
    max += vwgt[adjncy[j]];

  for (i=1; i<nvtxs; i++) {
    for (k=0, j=xadj[i]; j<xadj[i+1]; j++)
      k += vwgt[adjncy[j]];
    if (max < k)
      max = k;
  }

  return max;
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * srefine.c
 *
 * This file contains code for the separator refinement algortihms
 *
 * Started 8/1/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function is the entry point of the separator refinement
**************************************************************************/
void Refine2WayNode(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, float ubfactor)
{

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  for (;;) {
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RefTmr));
    if (ctrl->RType != 15)
      FM_2WayNodeBalance(ctrl, graph, ubfactor);

    switch (ctrl->RType) {
      case 1:
        FM_2WayNodeRefine(ctrl, graph, ubfactor, 8);
        break;
      case 2:
        FM_2WayNodeRefine_OneSided(ctrl, graph, ubfactor, 8);
        break;
      case 3:
        FM_2WayNodeRefine(ctrl, graph, ubfactor, 8);
        FM_2WayNodeRefine_OneSided(ctrl, graph, ubfactor, 8);
        break;
      case 4:
        FM_2WayNodeRefine_OneSided(ctrl, graph, ubfactor, 8);
        FM_2WayNodeRefine(ctrl, graph, ubfactor, 8);
        break;
      case 5:
        FM_2WayNodeRefineEqWgt(ctrl, graph, 8);
        break;
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    Project2WayNodePartition(ctrl, graph);
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
}


/*************************************************************************
* This function allocates memory for 2-way edge refinement
**************************************************************************/
void Allocate2WayNodePartitionMemory(CtrlType *ctrl, GraphType *graph)
{
  int nvtxs, pad64;

  nvtxs = graph->nvtxs;

  pad64 = (3*nvtxs+3)%2;

  graph->rdata = idxmalloc(3*nvtxs+3+(sizeof(NRInfoType)/sizeof(idxtype))*nvtxs+pad64, "Allocate2WayPartitionMemory: rdata");
  graph->pwgts          = graph->rdata;
  graph->where          = graph->rdata + 3;
  graph->bndptr         = graph->rdata + nvtxs + 3;
  graph->bndind         = graph->rdata + 2*nvtxs + 3;
  graph->nrinfo         = (NRInfoType *)(graph->rdata + 3*nvtxs + 3 + pad64);
}



/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void Compute2WayNodePartitionParams(CtrlType *ctrl, GraphType *graph)
{
  int i, j, nvtxs, nbnd;
  idxtype *xadj, *adjncy, *adjwgt, *vwgt;
  idxtype *where, *pwgts, *bndind, *bndptr, *edegrees;
  NRInfoType *rinfo;
  int me, other;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  rinfo = graph->nrinfo;
  pwgts = idxset(3, 0, graph->pwgts);
  bndind = graph->bndind;
  bndptr = idxset(nvtxs, -1, graph->bndptr);


  /*------------------------------------------------------------
  / Compute now the separator external degrees
  /------------------------------------------------------------*/
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];

    ASSERT(me >=0 && me <= 2);

    if (me == 2) { /* If it is on the separator do some computations */
      BNDInsert(nbnd, bndind, bndptr, i);

      edegrees = rinfo[i].edegrees;
      edegrees[0] = edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (other != 2)
          edegrees[other] += vwgt[adjncy[j]];
      }
    }
  }

  ASSERT(CheckNodeBnd(graph, nbnd));

  graph->mincut = pwgts[2];
  graph->nbnd = nbnd;
}


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void Project2WayNodePartition(CtrlType *ctrl, GraphType *graph)
{
  int i, nvtxs;
  idxtype *cmap, *where, *cwhere;
  GraphType *cgraph;

  cgraph = graph->coarser;
  cwhere = cgraph->where;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;

  Allocate2WayNodePartitionMemory(ctrl, graph);
  where = graph->where;

  /* Project the partition */
  for (i=0; i<nvtxs; i++) {
    where[i] = cwhere[cmap[i]];
    ASSERTP(where[i] >= 0 && where[i] <= 2, ("%d %d %d %d\n", i, cmap[i], where[i], cwhere[cmap[i]]));
  }

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

  Compute2WayNodePartitionParams(ctrl, graph);
}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stat.c
 *
 * This file computes various statistics
 *
 * Started 7/25/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function computes cuts and balance information
**************************************************************************/
void ComputePartitionInfo(GraphType *graph, int nparts, idxtype *where)
{
  int i, j, nvtxs, ncon, mustfree=0;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt, *kpwgts, *tmpptr;
  idxtype *padjncy, *padjwgt, *padjcut;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  adjwgt = graph->adjwgt;

  if (vwgt == NULL) {
    vwgt = graph->vwgt = idxsmalloc(nvtxs, 1, "vwgt");
    mustfree = 1;
  }
  if (adjwgt == NULL) {
    adjwgt = graph->adjwgt = idxsmalloc(xadj[nvtxs], 1, "adjwgt");
    mustfree += 2;
  }

  printf("%d-way Cut: %5d, Vol: %5d, ", nparts, ComputeCut(graph, where), ComputeVolume(graph, where));

  /* Compute balance information */
  kpwgts = idxsmalloc(ncon*nparts, 0, "ComputePartitionInfo: kpwgts");

  for (i=0; i<nvtxs; i++) {
    for (j=0; j<ncon; j++)
      kpwgts[where[i]*ncon+j] += vwgt[i*ncon+j];
  }

  if (ncon == 1) {
    printf("\tBalance: %5.3f out of %5.3f\n",
            1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)),
            1.0*nparts*vwgt[idxamax(nvtxs, vwgt)]/(1.0*idxsum(nparts, kpwgts)));
  }
  else {
    printf("\tBalance:");
    for (j=0; j<ncon; j++)
      printf(" (%5.3f out of %5.3f)",
            1.0*nparts*kpwgts[ncon*idxamax_strd(nparts, kpwgts+j, ncon)+j]/(1.0*idxsum_strd(nparts, kpwgts+j, ncon)),
            1.0*nparts*vwgt[ncon*idxamax_strd(nvtxs, vwgt+j, ncon)+j]/(1.0*idxsum_strd(nparts, kpwgts+j, ncon)));
    printf("\n");
  }


  /* Compute p-adjncy information */
  padjncy = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjncy");
  padjwgt = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");
  padjcut = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");

  idxset(nparts, 0, kpwgts);
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        padjncy[where[i]*nparts+where[adjncy[j]]] = 1;
        padjcut[where[i]*nparts+where[adjncy[j]]] += adjwgt[j];
        if (kpwgts[where[adjncy[j]]] == 0) {
          padjwgt[where[i]*nparts+where[adjncy[j]]]++;
          kpwgts[where[adjncy[j]]] = 1;
        }
      }
    }
    for (j=xadj[i]; j<xadj[i+1]; j++)
      kpwgts[where[adjncy[j]]] = 0;
  }

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjncy+i*nparts);
  printf("Min/Max/Avg/Bal # of adjacent     subdomains: %5d %5d %5.2f %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)],
    1.0*idxsum(nparts, kpwgts)/(1.0*nparts),
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)));

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjcut+i*nparts);
  printf("Min/Max/Avg/Bal # of adjacent subdomain cuts: %5d %5d %5d %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts,
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)));

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjwgt+i*nparts);
  printf("Min/Max/Avg/Bal/Frac # of interface    nodes: %5d %5d %5d %7.3f %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts,
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)), 1.0*idxsum(nparts, kpwgts)/(1.0*nvtxs));

  tmpptr = graph->where;
  graph->where = where;
  for (i=0; i<nparts; i++)
    IsConnectedSubdomain(NULL, graph, i, 1);
  graph->where = tmpptr;

  if (mustfree == 1 || mustfree == 3) {
    free(vwgt);
    graph->vwgt = NULL;
  }
  if (mustfree == 2 || mustfree == 3) {
    free(adjwgt);
    graph->adjwgt = NULL;
  }

  GKfree((void **) &kpwgts,  (void **) &padjncy,
         (void **) &padjwgt, (void **) &padjcut, LTERM);
}


/*************************************************************************
* This function computes cuts and balance information
**************************************************************************/
void ComputePartitionInfoBipartite(GraphType *graph, int nparts, idxtype *where)
{
  int i, j, nvtxs, ncon, mustfree=0;
  idxtype *xadj, *adjncy, *vwgt, *vsize, *adjwgt, *kpwgts;
  idxtype *padjncy, *padjwgt, *padjcut;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  adjwgt = graph->adjwgt;

  if (vwgt == NULL) {
    vwgt = graph->vwgt = idxsmalloc(nvtxs, 1, "vwgt");
    mustfree = 1;
  }
  if (adjwgt == NULL) {
    adjwgt = graph->adjwgt = idxsmalloc(xadj[nvtxs], 1, "adjwgt");
    mustfree += 2;
  }

  printf("%d-way Cut: %5d, Vol: %5d, ", nparts, ComputeCut(graph, where), ComputeVolume(graph, where));

  /* Compute balance information */
  kpwgts = idxsmalloc(ncon*nparts, 0, "ComputePartitionInfo: kpwgts");

  for (i=0; i<nvtxs; i++) {
    for (j=0; j<ncon; j++)
      kpwgts[where[i]*ncon+j] += vwgt[i*ncon+j];
  }

  if (ncon == 1) {
    printf("\tBalance: %5.3f out of %5.3f\n",
            1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)),
            1.0*nparts*vwgt[idxamax(nvtxs, vwgt)]/(1.0*idxsum(nparts, kpwgts)));
  }
  else {
    printf("\tBalance:");
    for (j=0; j<ncon; j++)
      printf(" (%5.3f out of %5.3f)",
            1.0*nparts*kpwgts[ncon*idxamax_strd(nparts, kpwgts+j, ncon)+j]/(1.0*idxsum_strd(nparts, kpwgts+j, ncon)),
            1.0*nparts*vwgt[ncon*idxamax_strd(nvtxs, vwgt+j, ncon)+j]/(1.0*idxsum_strd(nparts, kpwgts+j, ncon)));
    printf("\n");
  }


  /* Compute p-adjncy information */
  padjncy = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjncy");
  padjwgt = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");
  padjcut = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");

  idxset(nparts, 0, kpwgts);
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        padjncy[where[i]*nparts+where[adjncy[j]]] = 1;
        padjcut[where[i]*nparts+where[adjncy[j]]] += adjwgt[j];
        if (kpwgts[where[adjncy[j]]] == 0) {
          padjwgt[where[i]*nparts+where[adjncy[j]]] += vsize[i];
          kpwgts[where[adjncy[j]]] = 1;
        }
      }
    }
    for (j=xadj[i]; j<xadj[i+1]; j++)
      kpwgts[where[adjncy[j]]] = 0;
  }

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjncy+i*nparts);
  printf("Min/Max/Avg/Bal # of adjacent     subdomains: %5d %5d %5d %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts,
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)));

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjcut+i*nparts);
  printf("Min/Max/Avg/Bal # of adjacent subdomain cuts: %5d %5d %5d %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts,
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)));

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjwgt+i*nparts);
  printf("Min/Max/Avg/Bal/Frac # of interface    nodes: %5d %5d %5d %7.3f %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts,
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)), 1.0*idxsum(nparts, kpwgts)/(1.0*nvtxs));


  if (mustfree == 1 || mustfree == 3) {
    free(vwgt);
    graph->vwgt = NULL;
  }
  if (mustfree == 2 || mustfree == 3) {
    free(adjwgt);
    graph->adjwgt = NULL;
  }

  GKfree((void **) &kpwgts,  (void **) &padjncy,
         (void **) &padjwgt, (void **) &padjcut, LTERM);
}



/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void ComputePartitionBalance(GraphType *graph, int nparts, idxtype *where, float *ubvec)
{
  int i, j, nvtxs, ncon;
  idxtype *kpwgts, *vwgt;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  vwgt = graph->vwgt;

  kpwgts = idxsmalloc(nparts, 0, "ComputePartitionInfo: kpwgts");

  if (vwgt == NULL) {
    for (i=0; i<nvtxs; i++)
      kpwgts[where[i]]++;
    ubvec[0] = 1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*nvtxs);
  }
  else {
    for (j=0; j<ncon; j++) {
      idxset(nparts, 0, kpwgts);
      for (i=0; i<graph->nvtxs; i++)
        kpwgts[where[i]] += vwgt[i*ncon+j];

      ubvec[j] = 1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts));
    }
  }

  free(kpwgts);

}


/*************************************************************************
* This function computes the balance of the element partitioning
**************************************************************************/
float ComputeElementBalance(int ne, int nparts, idxtype *where)
{
  int i;
  idxtype *kpwgts;
  float balance;

  kpwgts = idxsmalloc(nparts, 0, "ComputeElementBalance: kpwgts");

  for (i=0; i<ne; i++)
    kpwgts[where[i]]++;

  balance = 1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts));

  free(kpwgts);

  return balance;

}
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * subdomains.c
 *
 * This file contains functions that deal with prunning the number of
 * adjacent subdomains in KMETIS
 *
 * Started 7/15/98
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Random_KWayEdgeRefineMConn(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor, int npasses, int ffactor)
{
  int i, ii, iii, j, k, l, pass, nvtxs, nmoves, nbnd, tvwgt, myndegrees;
  int from, me, to, oldcut, vwgt, gain;
  int maxndoms, nadd;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *itpwgts;
  idxtype *phtable, *pmat, *pmatptr, *ndoms;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndptr = graph->bndptr;
  bndind = graph->bndind;

  where = graph->where;
  pwgts = graph->pwgts;

  pmat = ctrl->wspace.pmat;
  phtable = idxwspacemalloc(ctrl, nparts);
  ndoms = idxwspacemalloc(ctrl, nparts);

  ComputeSubDomainGraph(graph, nparts, pmat, ndoms);

  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d]-[%6d %6d], Balance: %5.3f, Nv-Nb[%6d %6d]. Cut: %6d\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    maxndoms = ndoms[idxamax(nparts, ndoms)];

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (nmoves=iii=0; iii<graph->nbnd; iii++) {
      ii = perm[iii];
      if (ii >= nbnd)
        continue;
      i = bndind[ii];

      myrinfo = graph->rinfo+i;

      if (myrinfo->ed >= myrinfo->id) { /* Total ED is too high */
        from = where[i];
        vwgt = graph->vwgt[i];

        if (myrinfo->id > 0 && pwgts[from]-vwgt < minwgt[from])
          continue;   /* This cannot be moved! */

        myedegrees = myrinfo->edegrees;
        myndegrees = myrinfo->ndegrees;

        /* Determine the valid domains */
        for (j=0; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          phtable[to] = 1;
          pmatptr = pmat + to*nparts;
          for (nadd=0, k=0; k<myndegrees; k++) {
            if (k == j)
              continue;

            l = myedegrees[k].pid;
            if (pmatptr[l] == 0) {
              if (ndoms[l] > maxndoms-1) {
                phtable[to] = 0;
                nadd = maxndoms;
                break;
              }
              nadd++;
            }
          }
          if (ndoms[to]+nadd > maxndoms)
            phtable[to] = 0;
          if (nadd == 0)
            phtable[to] = 2;
        }

        /* Find the first valid move */
        j = myrinfo->id;
        for (k=0; k<myndegrees; k++) {
          to = myedegrees[k].pid;
          if (!phtable[to])
            continue;
          gain = myedegrees[k].ed-j; /* j = myrinfo->id. Allow good nodes to move */
          if (pwgts[to]+vwgt <= maxwgt[to]+ffactor*gain && gain >= 0)
            break;
        }
        if (k == myndegrees)
          continue;  /* break out if you did not find a candidate */

        for (j=k+1; j<myndegrees; j++) {
          to = myedegrees[j].pid;
          if (!phtable[to])
            continue;
          if ((myedegrees[j].ed > myedegrees[k].ed && pwgts[to]+vwgt <= maxwgt[to]) ||
              (myedegrees[j].ed == myedegrees[k].ed &&
               itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid]))
            k = j;
        }

        to = myedegrees[k].pid;

        j = 0;
        if (myedegrees[k].ed-myrinfo->id > 0)
          j = 1;
        else if (myedegrees[k].ed-myrinfo->id == 0) {
          if (/*(iii&7) == 0  ||*/ phtable[myedegrees[k].pid] == 2 || pwgts[from] >= maxwgt[from] || itpwgts[from]*(pwgts[to]+vwgt) < itpwgts[to]*pwgts[from])
            j = 1;
        }
        if (j == 0)
          continue;

        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to'
        *======================================================================*/
        graph->mincut -= myedegrees[k].ed-myrinfo->id;

        IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

        /* Update pmat to reflect the move of 'i' */
        pmat[from*nparts+to] += (myrinfo->id-myedegrees[k].ed);
        pmat[to*nparts+from] += (myrinfo->id-myedegrees[k].ed);
        if (pmat[from*nparts+to] == 0) {
          ndoms[from]--;
          if (ndoms[from]+1 == maxndoms)
            maxndoms = ndoms[idxamax(nparts, ndoms)];
        }
        if (pmat[to*nparts+from] == 0) {
          ndoms[to]--;
          if (ndoms[to]+1 == maxndoms)
            maxndoms = ndoms[idxamax(nparts, ndoms)];
        }

        /* Update where, weight, and ID/ED information of the vertex you moved */
        where[i] = to;
        INC_DEC(pwgts[to], pwgts[from], vwgt);
        myrinfo->ed += myrinfo->id-myedegrees[k].ed;
        SWAP(myrinfo->id, myedegrees[k].ed, j);
        if (myedegrees[k].ed == 0)
          myedegrees[k] = myedegrees[--myrinfo->ndegrees];
        else
          myedegrees[k].pid = from;

        if (myrinfo->ed-myrinfo->id < 0)
          BNDDelete(nbnd, bndind, bndptr, i);

        /* Update the degrees of adjacent vertices */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          ii = adjncy[j];
          me = where[ii];

          myrinfo = graph->rinfo+ii;
          if (myrinfo->edegrees == NULL) {
            myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
            ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
          }
          myedegrees = myrinfo->edegrees;

          ASSERT(CheckRInfo(myrinfo));

          if (me == from) {
            INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

            if (myrinfo->ed-myrinfo->id >= 0 && bndptr[ii] == -1)
              BNDInsert(nbnd, bndind, bndptr, ii);
          }
          else if (me == to) {
            INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

            if (myrinfo->ed-myrinfo->id < 0 && bndptr[ii] != -1)
              BNDDelete(nbnd, bndind, bndptr, ii);
          }

          /* Remove contribution from the .ed of 'from' */
          if (me != from) {
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (myedegrees[k].pid == from) {
                if (myedegrees[k].ed == adjwgt[j])
                  myedegrees[k] = myedegrees[--myrinfo->ndegrees];
                else
                  myedegrees[k].ed -= adjwgt[j];
                break;
              }
            }
          }

          /* Add contribution to the .ed of 'to' */
          if (me != to) {
            for (k=0; k<myrinfo->ndegrees; k++) {
              if (myedegrees[k].pid == to) {
                myedegrees[k].ed += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->ndegrees) {
              myedegrees[myrinfo->ndegrees].pid = to;
              myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
            }
          }

          /* Update pmat to reflect the move of 'i' for domains other than 'from' and 'to' */
          if (me != from && me != to) {
            pmat[me*nparts+from] -= adjwgt[j];
            pmat[from*nparts+me] -= adjwgt[j];
            if (pmat[me*nparts+from] == 0) {
              ndoms[me]--;
              if (ndoms[me]+1 == maxndoms)
                maxndoms = ndoms[idxamax(nparts, ndoms)];
            }
            if (pmat[from*nparts+me] == 0) {
              ndoms[from]--;
              if (ndoms[from]+1 == maxndoms)
                maxndoms = ndoms[idxamax(nparts, ndoms)];
            }

            if (pmat[me*nparts+to] == 0) {
              ndoms[me]++;
              if (ndoms[me] > maxndoms) {
                printf("You just increased the maxndoms: %d %d\n", ndoms[me], maxndoms);
                maxndoms = ndoms[me];
              }
            }
            if (pmat[to*nparts+me] == 0) {
              ndoms[to]++;
              if (ndoms[to] > maxndoms) {
                printf("You just increased the maxndoms: %d %d\n", ndoms[to], maxndoms);
                maxndoms = ndoms[to];
              }
            }
            pmat[me*nparts+to] += adjwgt[j];
            pmat[to*nparts+me] += adjwgt[j];
          }

          ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
          ASSERT(CheckRInfo(myrinfo));

        }
        nmoves++;
      }
    }

    graph->nbnd = nbnd;

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %5d, Vol: %5d, %d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves,
               graph->mincut, ComputeVolume(graph, where), idxsum(nparts, ndoms)));

    if (graph->mincut == oldcut)
      break;
  }

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Greedy_KWayEdgeBalanceMConn(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor, int npasses)
{
  int i, ii, j, k, l, pass, nvtxs, nbnd, tvwgt, myndegrees, oldgain, gain, nmoves;
  int from, me, to, oldcut, vwgt, maxndoms, nadd;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *pwgts, *perm, *bndptr, *bndind, *minwgt, *maxwgt, *moved, *itpwgts;
  idxtype *phtable, *pmat, *pmatptr, *ndoms;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;
  PQueueType queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;

  where = graph->where;
  pwgts = graph->pwgts;

  pmat = ctrl->wspace.pmat;
  phtable = idxwspacemalloc(ctrl, nparts);
  ndoms = idxwspacemalloc(ctrl, nparts);

  ComputeSubDomainGraph(graph, nparts, pmat, ndoms);


  /* Setup the weight intervals of the various subdomains */
  minwgt =  idxwspacemalloc(ctrl, nparts);
  maxwgt = idxwspacemalloc(ctrl, nparts);
  itpwgts = idxwspacemalloc(ctrl, nparts);
  tvwgt = idxsum(nparts, pwgts);
  ASSERT(tvwgt == idxsum(nvtxs, graph->vwgt));

  for (i=0; i<nparts; i++) {
    itpwgts[i] = tpwgts[i]*tvwgt;
    maxwgt[i] = tpwgts[i]*tvwgt*ubfactor;
    minwgt[i] = tpwgts[i]*tvwgt*(1.0/ubfactor);
  }

  perm = idxwspacemalloc(ctrl, nvtxs);
  moved = idxwspacemalloc(ctrl, nvtxs);

  PQueueInit(ctrl, &queue, nvtxs, graph->adjwgtsum[idxamax(nvtxs, graph->adjwgtsum)]);

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%6d %6d]-[%6d %6d], Balance: %5.3f, Nv-Nb[%6d %6d]. Cut: %6d [B]\n",
             pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)], minwgt[0], maxwgt[0],
             1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nvtxs, graph->nbnd,
             graph->mincut));

  for (pass=0; pass<npasses; pass++) {
    ASSERT(ComputeCut(graph, where) == graph->mincut);

    /* Check to see if things are out of balance, given the tolerance */
    for (i=0; i<nparts; i++) {
      if (pwgts[i] > maxwgt[i])
        break;
    }
    if (i == nparts) /* Things are balanced. Return right away */
      break;

    PQueueReset(&queue);
    idxset(nvtxs, -1, moved);

    oldcut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      PQueueInsert(&queue, i, graph->rinfo[i].ed - graph->rinfo[i].id);
      moved[i] = 2;
    }

    maxndoms = ndoms[idxamax(nparts, ndoms)];

    for (nmoves=0;;) {
      if ((i = PQueueGetMax(&queue)) == -1)
        break;
      moved[i] = 1;

      myrinfo = graph->rinfo+i;
      from = where[i];
      vwgt = graph->vwgt[i];

      if (pwgts[from]-vwgt < minwgt[from])
        continue;   /* This cannot be moved! */

      myedegrees = myrinfo->edegrees;
      myndegrees = myrinfo->ndegrees;

      /* Determine the valid domains */
      for (j=0; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        phtable[to] = 1;
        pmatptr = pmat + to*nparts;
        for (nadd=0, k=0; k<myndegrees; k++) {
          if (k == j)
            continue;

          l = myedegrees[k].pid;
          if (pmatptr[l] == 0) {
            if (ndoms[l] > maxndoms-1) {
              phtable[to] = 0;
              nadd = maxndoms;
              break;
            }
            nadd++;
          }
        }
        if (ndoms[to]+nadd > maxndoms)
          phtable[to] = 0;
      }

      for (k=0; k<myndegrees; k++) {
        to = myedegrees[k].pid;
        if (!phtable[to])
          continue;
        if (pwgts[to]+vwgt <= maxwgt[to] || itpwgts[from]*(pwgts[to]+vwgt) <= itpwgts[to]*pwgts[from])
          break;
      }
      if (k == myndegrees)
        continue;  /* break out if you did not find a candidate */

      for (j=k+1; j<myndegrees; j++) {
        to = myedegrees[j].pid;
        if (!phtable[to])
          continue;
        if (itpwgts[myedegrees[k].pid]*pwgts[to] < itpwgts[to]*pwgts[myedegrees[k].pid])
          k = j;
      }

      to = myedegrees[k].pid;

      if (pwgts[from] < maxwgt[from] && pwgts[to] > minwgt[to] && myedegrees[k].ed-myrinfo->id < 0)
        continue;

      /*=====================================================================
      * If we got here, we can now move the vertex from 'from' to 'to'
      *======================================================================*/
      graph->mincut -= myedegrees[k].ed-myrinfo->id;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO, printf("\t\tMoving %6d to %3d. Gain: %4d. Cut: %6d\n", i, to, myedegrees[k].ed-myrinfo->id, graph->mincut));

      /* Update pmat to reflect the move of 'i' */
      pmat[from*nparts+to] += (myrinfo->id-myedegrees[k].ed);
      pmat[to*nparts+from] += (myrinfo->id-myedegrees[k].ed);
      if (pmat[from*nparts+to] == 0) {
        ndoms[from]--;
        if (ndoms[from]+1 == maxndoms)
          maxndoms = ndoms[idxamax(nparts, ndoms)];
      }
      if (pmat[to*nparts+from] == 0) {
        ndoms[to]--;
        if (ndoms[to]+1 == maxndoms)
          maxndoms = ndoms[idxamax(nparts, ndoms)];
      }


      /* Update where, weight, and ID/ED information of the vertex you moved */
      where[i] = to;
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      myrinfo->ed += myrinfo->id-myedegrees[k].ed;
      SWAP(myrinfo->id, myedegrees[k].ed, j);
      if (myedegrees[k].ed == 0)
        myedegrees[k] = myedegrees[--myrinfo->ndegrees];
      else
        myedegrees[k].pid = from;

      if (myrinfo->ed == 0)
        BNDDelete(nbnd, bndind, bndptr, i);

      /* Update the degrees of adjacent vertices */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        me = where[ii];

        myrinfo = graph->rinfo+ii;
        if (myrinfo->edegrees == NULL) {
          myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
          ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
        }
        myedegrees = myrinfo->edegrees;

        ASSERT(CheckRInfo(myrinfo));

        oldgain = (myrinfo->ed-myrinfo->id);

        if (me == from) {
          INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

          if (myrinfo->ed > 0 && bndptr[ii] == -1)
            BNDInsert(nbnd, bndind, bndptr, ii);
        }
        else if (me == to) {
          INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

          if (myrinfo->ed == 0 && bndptr[ii] != -1)
            BNDDelete(nbnd, bndind, bndptr, ii);
        }

        /* Remove contribution from the .ed of 'from' */
        if (me != from) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == from) {
              if (myedegrees[k].ed == adjwgt[j])
                myedegrees[k] = myedegrees[--myrinfo->ndegrees];
              else
                myedegrees[k].ed -= adjwgt[j];
              break;
            }
          }
        }

        /* Add contribution to the .ed of 'to' */
        if (me != to) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (myedegrees[k].pid == to) {
              myedegrees[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            myedegrees[myrinfo->ndegrees].pid = to;
            myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
          }
        }

        /* Update pmat to reflect the move of 'i' for domains other than 'from' and 'to' */
        if (me != from && me != to) {
          pmat[me*nparts+from] -= adjwgt[j];
          pmat[from*nparts+me] -= adjwgt[j];
          if (pmat[me*nparts+from] == 0) {
            ndoms[me]--;
            if (ndoms[me]+1 == maxndoms)
              maxndoms = ndoms[idxamax(nparts, ndoms)];
          }
          if (pmat[from*nparts+me] == 0) {
            ndoms[from]--;
            if (ndoms[from]+1 == maxndoms)
              maxndoms = ndoms[idxamax(nparts, ndoms)];
          }

          if (pmat[me*nparts+to] == 0) {
            ndoms[me]++;
            if (ndoms[me] > maxndoms) {
              printf("You just increased the maxndoms: %d %d\n", ndoms[me], maxndoms);
              maxndoms = ndoms[me];
            }
          }
          if (pmat[to*nparts+me] == 0) {
            ndoms[to]++;
            if (ndoms[to] > maxndoms) {
              printf("You just increased the maxndoms: %d %d\n", ndoms[to], maxndoms);
              maxndoms = ndoms[to];
            }
          }
          pmat[me*nparts+to] += adjwgt[j];
          pmat[to*nparts+me] += adjwgt[j];
        }

        /* Update the queue */
        if (me == to || me == from) {
          gain = myrinfo->ed-myrinfo->id;
          if (moved[ii] == 2) {
            if (myrinfo->ed > 0)
              PQueueUpdate(&queue, ii, oldgain, gain);
            else {
              PQueueDelete(&queue, ii, oldgain);
              moved[ii] = -1;
            }
          }
          else if (moved[ii] == -1 && myrinfo->ed > 0) {
            PQueueInsert(&queue, ii, gain);
            moved[ii] = 2;
          }
        }

        ASSERT(myrinfo->ndegrees <= xadj[ii+1]-xadj[ii]);
        ASSERT(CheckRInfo(myrinfo));
      }
      nmoves++;
    }

    graph->nbnd = nbnd;

    IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Balance: %5.3f, Nb: %6d. Nmoves: %5d, Cut: %6d, %d\n",
               pwgts[idxamin(nparts, pwgts)], pwgts[idxamax(nparts, pwgts)],
               1.0*nparts*pwgts[idxamax(nparts, pwgts)]/tvwgt, graph->nbnd, nmoves, graph->mincut,idxsum(nparts, ndoms)));
  }

  PQueueFree(ctrl, &queue);

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}




/*************************************************************************
* This function computes the subdomain graph
**************************************************************************/
void PrintSubDomainGraph(GraphType *graph, int nparts, idxtype *where)
{
  int i, j, k, me, nvtxs, total, max;
  idxtype *xadj, *adjncy, *adjwgt, *pmat;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  pmat = idxsmalloc(nparts*nparts, 0, "ComputeSubDomainGraph: pmat");

  for (i=0; i<nvtxs; i++) {
    me = where[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] != me)
        pmat[me*nparts+where[k]] += adjwgt[j];
    }
  }

  /* printf("Subdomain Info\n"); */
  total = max = 0;
  for (i=0; i<nparts; i++) {
    for (k=0, j=0; j<nparts; j++) {
      if (pmat[i*nparts+j] > 0)
        k++;
    }
    total += k;

    if (k > max)
      max = k;
/*
    printf("%2d -> %2d  ", i, k);
    for (j=0; j<nparts; j++) {
      if (pmat[i*nparts+j] > 0)
        printf("[%2d %4d] ", j, pmat[i*nparts+j]);
    }
    printf("\n");
*/
  }
  printf("Total adjacent subdomains: %d, Max: %d\n", total, max);

  free(pmat);
}



/*************************************************************************
* This function computes the subdomain graph
**************************************************************************/
void ComputeSubDomainGraph(GraphType *graph, int nparts, idxtype *pmat, idxtype *ndoms)
{
  int i, j, k, me, nvtxs, ndegrees;
  idxtype *xadj, *adjncy, *adjwgt, *where;
  RInfoType *rinfo;
  EDegreeType *edegrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  rinfo = graph->rinfo;

  idxset(nparts*nparts, 0, pmat);

  for (i=0; i<nvtxs; i++) {
    if (rinfo[i].ed > 0) {
      me = where[i];
      ndegrees = rinfo[i].ndegrees;
      edegrees = rinfo[i].edegrees;

      k = me*nparts;
      for (j=0; j<ndegrees; j++)
        pmat[k+edegrees[j].pid] += edegrees[j].ed;
    }
  }

  for (i=0; i<nparts; i++) {
    ndoms[i] = 0;
    for (j=0; j<nparts; j++) {
      if (pmat[i*nparts+j] > 0)
        ndoms[i]++;
    }
  }

}





/*************************************************************************
* This function computes the subdomain graph
**************************************************************************/
void EliminateSubDomainEdges(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts)
{
  int i, ii, j, k, me, other, nvtxs, total, max, avg, totalout, nind, ncand, ncand2, target, target2, nadd;
  int min, move, cpwgt, tvwgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt, *pwgts, *where, *maxpwgt, *pmat, *ndoms, *mypmat, *otherpmat, *ind;
  KeyValueType *cand, *cand2;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = graph->pwgts;  /* We assume that this is properly initialized */

  maxpwgt = idxwspacemalloc(ctrl, nparts);
  ndoms = idxwspacemalloc(ctrl, nparts);
  otherpmat = idxwspacemalloc(ctrl, nparts);
  ind = idxwspacemalloc(ctrl, nvtxs);
  pmat = ctrl->wspace.pmat;

  cand = (KeyValueType *)GKmalloc(nparts*sizeof(KeyValueType), "EliminateSubDomainEdges: cand");
  cand2 = (KeyValueType *)GKmalloc(nparts*sizeof(KeyValueType), "EliminateSubDomainEdges: cand");

  /* Compute the pmat matrix and ndoms */
  ComputeSubDomainGraph(graph, nparts, pmat, ndoms);


  /* Compute the maximum allowed weight for each domain */
  tvwgt = idxsum(nparts, pwgts);
  for (i=0; i<nparts; i++)
    maxpwgt[i] = 1.25*tpwgts[i]*tvwgt;


  /* Get into the loop eliminating subdomain connections */
  for (;;) {
    total = idxsum(nparts, ndoms);
    avg = total/nparts;
    max = ndoms[idxamax(nparts, ndoms)];

    /* printf("Adjacent Subdomain Stats: Total: %3d, Max: %3d, Avg: %3d [%5d]\n", total, max, avg, idxsum(nparts*nparts, pmat)); */

    if (max < 1.4*avg)
      break;

    me = idxamax(nparts, ndoms);
    mypmat = pmat + me*nparts;
    totalout = idxsum(nparts, mypmat);

    /*printf("Me: %d, TotalOut: %d,\n", me, totalout);*/

    /* Sort the connections according to their cut */
    for (ncand2=0, i=0; i<nparts; i++) {
      if (mypmat[i] > 0) {
        cand2[ncand2].key = mypmat[i];
        cand2[ncand2++].val = i;
      }
    }
    ikeysort(ncand2, cand2);

    move = 0;
    for (min=0; min<ncand2; min++) {
      if (cand2[min].key > totalout/(2*ndoms[me]))
        break;

      other = cand2[min].val;

      /*printf("\tMinOut: %d to %d\n", mypmat[other], other);*/

      idxset(nparts, 0, otherpmat);

      /* Go and find the vertices in 'other' that are connected in 'me' */
      for (nind=0, i=0; i<nvtxs; i++) {
        if (where[i] == other) {
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            if (where[adjncy[j]] == me) {
              ind[nind++] = i;
              break;
            }
          }
        }
      }

      /* Go and construct the otherpmat to see where these nind vertices are connected to */
      for (cpwgt=0, ii=0; ii<nind; ii++) {
        i = ind[ii];
        cpwgt += vwgt[i];

        for (j=xadj[i]; j<xadj[i+1]; j++)
          otherpmat[where[adjncy[j]]] += adjwgt[j];
      }
      otherpmat[other] = 0;

      for (ncand=0, i=0; i<nparts; i++) {
        if (otherpmat[i] > 0) {
          cand[ncand].key = -otherpmat[i];
          cand[ncand++].val = i;
        }
      }
      ikeysort(ncand, cand);

      /*
       * Go through and the select the first domain that is common with 'me', and
       * does not increase the ndoms[target] higher than my ndoms, subject to the
       * maxpwgt constraint. Traversal is done from the mostly connected to the least.
       */
      target = target2 = -1;
      for (i=0; i<ncand; i++) {
        k = cand[i].val;

        if (mypmat[k] > 0) {
          if (pwgts[k] + cpwgt > maxpwgt[k])  /* Check if balance will go off */
            continue;

          for (j=0; j<nparts; j++) {
            if (otherpmat[j] > 0 && ndoms[j] >= ndoms[me]-1 && pmat[nparts*j+k] == 0)
              break;
          }
          if (j == nparts) { /* No bad second level effects */
            for (nadd=0, j=0; j<nparts; j++) {
              if (otherpmat[j] > 0 && pmat[nparts*k+j] == 0)
                nadd++;
            }

            /*printf("\t\tto=%d, nadd=%d, %d\n", k, nadd, ndoms[k]);*/
            if (target2 == -1 && ndoms[k]+nadd < ndoms[me]) {
              target2 = k;
            }
            if (nadd == 0) {
              target = k;
              break;
            }
          }
        }
      }
      if (target == -1 && target2 != -1)
        target = target2;

      if (target == -1) {
        /* printf("\t\tCould not make the move\n");*/
        continue;
      }

      /*printf("\t\tMoving to %d\n", target);*/

      /* Update the partition weights */
      INC_DEC(pwgts[target], pwgts[other], cpwgt);

      MoveGroupMConn(ctrl, graph, ndoms, pmat, nparts, target, nind, ind);

      move = 1;
      break;
    }

    if (move == 0)
      break;
  }

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);

  GKfree((void **) &cand, (void **) &cand2, LTERM);
}


/*************************************************************************
* This function moves a collection of vertices and updates their rinfo
**************************************************************************/
void MoveGroupMConn(CtrlType *ctrl, GraphType *graph, idxtype *ndoms, idxtype *pmat,
                    int nparts, int to, int nind, idxtype *ind)
{
  int i, ii, iii, j, k, nvtxs, nbnd;
  int from, me;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *bndptr, *bndind;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  nbnd = graph->nbnd;

  for (iii=0; iii<nind; iii++) {
    i = ind[iii];
    from = where[i];

    myrinfo = graph->rinfo+i;
    if (myrinfo->edegrees == NULL) {
      myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += xadj[i+1]-xadj[i];
      myrinfo->ndegrees = 0;
    }
    myedegrees = myrinfo->edegrees;

    /* find the location of 'to' in myrinfo or create it if it is not there */
    for (k=0; k<myrinfo->ndegrees; k++) {
      if (myedegrees[k].pid == to)
        break;
    }
    if (k == myrinfo->ndegrees) {
      myedegrees[k].pid = to;
      myedegrees[k].ed = 0;
      myrinfo->ndegrees++;
    }

    graph->mincut -= myedegrees[k].ed-myrinfo->id;

    /* Update pmat to reflect the move of 'i' */
    pmat[from*nparts+to] += (myrinfo->id-myedegrees[k].ed);
    pmat[to*nparts+from] += (myrinfo->id-myedegrees[k].ed);
    if (pmat[from*nparts+to] == 0)
      ndoms[from]--;
    if (pmat[to*nparts+from] == 0)
      ndoms[to]--;

    /* Update where, weight, and ID/ED information of the vertex you moved */
    where[i] = to;
    myrinfo->ed += myrinfo->id-myedegrees[k].ed;
    SWAP(myrinfo->id, myedegrees[k].ed, j);
    if (myedegrees[k].ed == 0)
      myedegrees[k] = myedegrees[--myrinfo->ndegrees];
    else
      myedegrees[k].pid = from;

    if (myrinfo->ed-myrinfo->id < 0 && bndptr[i] != -1)
      BNDDelete(nbnd, bndind, bndptr, i);

    /* Update the degrees of adjacent vertices */
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ii = adjncy[j];
      me = where[ii];

      myrinfo = graph->rinfo+ii;
      if (myrinfo->edegrees == NULL) {
        myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
        ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
      }
      myedegrees = myrinfo->edegrees;

      ASSERT(CheckRInfo(myrinfo));

      if (me == from) {
        INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

        if (myrinfo->ed-myrinfo->id >= 0 && bndptr[ii] == -1)
          BNDInsert(nbnd, bndind, bndptr, ii);
      }
      else if (me == to) {
        INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

        if (myrinfo->ed-myrinfo->id < 0 && bndptr[ii] != -1)
          BNDDelete(nbnd, bndind, bndptr, ii);
      }

      /* Remove contribution from the .ed of 'from' */
      if (me != from) {
        for (k=0; k<myrinfo->ndegrees; k++) {
          if (myedegrees[k].pid == from) {
            if (myedegrees[k].ed == adjwgt[j])
              myedegrees[k] = myedegrees[--myrinfo->ndegrees];
            else
              myedegrees[k].ed -= adjwgt[j];
            break;
          }
        }
      }

      /* Add contribution to the .ed of 'to' */
      if (me != to) {
        for (k=0; k<myrinfo->ndegrees; k++) {
          if (myedegrees[k].pid == to) {
            myedegrees[k].ed += adjwgt[j];
            break;
          }
        }
        if (k == myrinfo->ndegrees) {
          myedegrees[myrinfo->ndegrees].pid = to;
          myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
        }
      }

      /* Update pmat to reflect the move of 'i' for domains other than 'from' and 'to' */
      if (me != from && me != to) {
        pmat[me*nparts+from] -= adjwgt[j];
        pmat[from*nparts+me] -= adjwgt[j];
        if (pmat[me*nparts+from] == 0)
          ndoms[me]--;
        if (pmat[from*nparts+me] == 0)
          ndoms[from]--;

        if (pmat[me*nparts+to] == 0)
          ndoms[me]++;
        if (pmat[to*nparts+me] == 0)
          ndoms[to]++;

        pmat[me*nparts+to] += adjwgt[j];
        pmat[to*nparts+me] += adjwgt[j];
      }

      ASSERT(CheckRInfo(myrinfo));
    }

    ASSERT(CheckRInfo(graph->rinfo+i));
  }

  graph->nbnd = nbnd;

}




/*************************************************************************
* This function finds all the connected components induced by the
* partitioning vector in wgraph->where and tries to push them around to
* remove some of them
**************************************************************************/
void EliminateComponents(CtrlType *ctrl, GraphType *graph, int nparts, float *tpwgts, float ubfactor)
{
  int i, ii, j, jj, k, me, nvtxs, tvwgt, first, last, nleft, ncmps, cwgt, target, deltawgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt, *where, *pwgts, *maxpwgt;
  idxtype *cpvec, *touched, *perm, *todo, *cind, *cptr, *npcmps;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  adjwgt = graph->adjwgt;

  where = graph->where;
  pwgts = graph->pwgts;

  touched = idxset(nvtxs, 0, idxwspacemalloc(ctrl, nvtxs));
  cptr = idxwspacemalloc(ctrl, nvtxs);
  cind = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  todo = idxwspacemalloc(ctrl, nvtxs);
  maxpwgt = idxwspacemalloc(ctrl, nparts);
  cpvec = idxwspacemalloc(ctrl, nparts);
  npcmps = idxset(nparts, 0, idxwspacemalloc(ctrl, nparts));

  for (i=0; i<nvtxs; i++)
    perm[i] = todo[i] = i;

  /* Find the connected componends induced by the partition */
  ncmps = -1;
  first = last = 0;
  nleft = nvtxs;
  while (nleft > 0) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      ASSERT(touched[todo[0]] == 0);
      i = todo[0];
      cind[last++] = i;
      touched[i] = 1;
      me = where[i];
      npcmps[me]++;
    }

    i = cind[first++];
    k = perm[i];
    j = todo[k] = todo[--nleft];
    perm[j] = k;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] == me && !touched[k]) {
        cind[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  /* printf("I found %d components, for this %d-way partition\n", ncmps, nparts); */

  if (ncmps > nparts) { /* There are more components than processors */
    /* First determine the max allowed load imbalance */
    tvwgt = idxsum(nparts, pwgts);
    for (i=0; i<nparts; i++)
      maxpwgt[i] = ubfactor*tpwgts[i]*tvwgt;

    deltawgt = 5;

    for (i=0; i<ncmps; i++) {
      me = where[cind[cptr[i]]];  /* Get the domain of this component */
      if (npcmps[me] == 1)
        continue;  /* Skip it because it is contigous */

      /*printf("Trying to move %d from %d\n", i, me); */

      /* Determine the weight of the block to be moved and abort if too high */
      for (cwgt=0, j=cptr[i]; j<cptr[i+1]; j++)
        cwgt += vwgt[cind[j]];

      if (cwgt > .30*pwgts[me])
        continue;  /* Skip the component if it is over 30% of the weight */

      /* Determine the connectivity */
      idxset(nparts, 0, cpvec);
      for (j=cptr[i]; j<cptr[i+1]; j++) {
        ii = cind[j];
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++)
          cpvec[where[adjncy[jj]]] += adjwgt[jj];
      }
      cpvec[me] = 0;

      target = -1;
      for (j=0; j<nparts; j++) {
        if (cpvec[j] > 0 && (cwgt < deltawgt || pwgts[j] + cwgt < maxpwgt[j])) {
          if (target == -1 || cpvec[target] < cpvec[j])
            target = j;
        }
      }

      /* printf("\tMoving it to %d [%d]\n", target, cpvec[target]);*/

      if (target != -1) {
        /* Assign all the vertices of 'me' to 'target' and update data structures */
        INC_DEC(pwgts[target], pwgts[me], cwgt);
        npcmps[me]--;

        MoveGroup(ctrl, graph, nparts, target, i, cptr, cind);
      }
    }

  }

  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nparts);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}


/*************************************************************************
* This function moves a collection of vertices and updates their rinfo
**************************************************************************/
void MoveGroup(CtrlType *ctrl, GraphType *graph, int nparts, int to, int gid, idxtype *ptr, idxtype *ind)
{
  int i, ii, iii, j, k, nvtxs, nbnd;
  int from, me;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *bndptr, *bndind;
  EDegreeType *myedegrees;
  RInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  nbnd = graph->nbnd;

  for (iii=ptr[gid]; iii<ptr[gid+1]; iii++) {
    i = ind[iii];
    from = where[i];

    myrinfo = graph->rinfo+i;
    if (myrinfo->edegrees == NULL) {
      myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
      ctrl->wspace.cdegree += xadj[i+1]-xadj[i];
      myrinfo->ndegrees = 0;
    }
    myedegrees = myrinfo->edegrees;

    /* find the location of 'to' in myrinfo or create it if it is not there */
    for (k=0; k<myrinfo->ndegrees; k++) {
      if (myedegrees[k].pid == to)
        break;
    }
    if (k == myrinfo->ndegrees) {
      myedegrees[k].pid = to;
      myedegrees[k].ed = 0;
      myrinfo->ndegrees++;
    }

    graph->mincut -= myedegrees[k].ed-myrinfo->id;


    /* Update where, weight, and ID/ED information of the vertex you moved */
    where[i] = to;
    myrinfo->ed += myrinfo->id-myedegrees[k].ed;
    SWAP(myrinfo->id, myedegrees[k].ed, j);
    if (myedegrees[k].ed == 0)
      myedegrees[k] = myedegrees[--myrinfo->ndegrees];
    else
      myedegrees[k].pid = from;

    if (myrinfo->ed-myrinfo->id < 0 && bndptr[i] != -1)
      BNDDelete(nbnd, bndind, bndptr, i);

    /* Update the degrees of adjacent vertices */
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ii = adjncy[j];
      me = where[ii];

      myrinfo = graph->rinfo+ii;
      if (myrinfo->edegrees == NULL) {
        myrinfo->edegrees = ctrl->wspace.edegrees+ctrl->wspace.cdegree;
        ctrl->wspace.cdegree += xadj[ii+1]-xadj[ii];
      }
      myedegrees = myrinfo->edegrees;

      ASSERT(CheckRInfo(myrinfo));

      if (me == from) {
        INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);

        if (myrinfo->ed-myrinfo->id >= 0 && bndptr[ii] == -1)
          BNDInsert(nbnd, bndind, bndptr, ii);
      }
      else if (me == to) {
        INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);

        if (myrinfo->ed-myrinfo->id < 0 && bndptr[ii] != -1)
          BNDDelete(nbnd, bndind, bndptr, ii);
      }

      /* Remove contribution from the .ed of 'from' */
      if (me != from) {
        for (k=0; k<myrinfo->ndegrees; k++) {
          if (myedegrees[k].pid == from) {
            if (myedegrees[k].ed == adjwgt[j])
              myedegrees[k] = myedegrees[--myrinfo->ndegrees];
            else
              myedegrees[k].ed -= adjwgt[j];
            break;
          }
        }
      }

      /* Add contribution to the .ed of 'to' */
      if (me != to) {
        for (k=0; k<myrinfo->ndegrees; k++) {
          if (myedegrees[k].pid == to) {
            myedegrees[k].ed += adjwgt[j];
            break;
          }
        }
        if (k == myrinfo->ndegrees) {
          myedegrees[myrinfo->ndegrees].pid = to;
          myedegrees[myrinfo->ndegrees++].ed = adjwgt[j];
        }
      }

      ASSERT(CheckRInfo(myrinfo));
    }

    ASSERT(CheckRInfo(graph->rinfo+i));
  }

  graph->nbnd = nbnd;

}

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * timing.c
 *
 * This file contains routines that deal with timing Metis
 *
 * Started 7/24/97
 * George
 *
 * $Id$
 *
 */




/*************************************************************************
* This function clears the timers
**************************************************************************/
void InitTimers(CtrlType *ctrl)
{
  cleartimer(ctrl->TotalTmr);
  cleartimer(ctrl->InitPartTmr);
  cleartimer(ctrl->MatchTmr);
  cleartimer(ctrl->ContractTmr);
  cleartimer(ctrl->CoarsenTmr);
  cleartimer(ctrl->UncoarsenTmr);
  cleartimer(ctrl->RefTmr);
  cleartimer(ctrl->ProjectTmr);
  cleartimer(ctrl->SplitTmr);
  cleartimer(ctrl->SepTmr);
  cleartimer(ctrl->AuxTmr1);
  cleartimer(ctrl->AuxTmr2);
  cleartimer(ctrl->AuxTmr3);
  cleartimer(ctrl->AuxTmr4);
  cleartimer(ctrl->AuxTmr5);
  cleartimer(ctrl->AuxTmr6);
}



/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(CtrlType *ctrl)
{
  printf("\nTiming Information -------------------------------------------------");
  printf("\n Multilevel: \t\t %7.3f", gettimer(ctrl->TotalTmr));
  printf("\n     Coarsening: \t\t %7.3f", gettimer(ctrl->CoarsenTmr));
  printf("\n            Matching: \t\t\t %7.3f", gettimer(ctrl->MatchTmr));
  printf("\n            Contract: \t\t\t %7.3f", gettimer(ctrl->ContractTmr));
  printf("\n     Initial Partition: \t %7.3f", gettimer(ctrl->InitPartTmr));
  printf("\n   Construct Separator: \t %7.3f", gettimer(ctrl->SepTmr));
  printf("\n     Uncoarsening: \t\t %7.3f", gettimer(ctrl->UncoarsenTmr));
  printf("\n          Refinement: \t\t\t %7.3f", gettimer(ctrl->RefTmr));
  printf("\n          Projection: \t\t\t %7.3f", gettimer(ctrl->ProjectTmr));
  printf("\n     Splitting: \t\t %7.3f", gettimer(ctrl->SplitTmr));
  printf("\n          AUX1: \t\t %7.3f", gettimer(ctrl->AuxTmr1));
  printf("\n          AUX2: \t\t %7.3f", gettimer(ctrl->AuxTmr2));
  printf("\n          AUX3: \t\t %7.3f", gettimer(ctrl->AuxTmr3));
  printf("\n********************************************************************\n");
}


/*************************************************************************
* This function returns the seconds
**************************************************************************/
double seconds(void)
{
  return((double) clock()/CLOCKS_PER_SEC);
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id$
 */




/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s", out1);

  fprintf(stdout, "%s", out2);
  fflush(stdout);

  abort();
}



#ifndef DMALLOC
/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
int *imalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (int *)GKmalloc(sizeof(int)*n, msg);
}


/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *idxmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)GKmalloc(sizeof(idxtype)*n, msg);
}


/*************************************************************************
* The following function allocates an array of float
**************************************************************************/
float *fmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (float *)GKmalloc(sizeof(float)*n, msg);
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *ismalloc(int n, int ival, char *msg)
{
  if (n == 0)
    return NULL;

  return iset(n, ival, (int *)GKmalloc(sizeof(int)*n, msg));
}



/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *idxsmalloc(int n, idxtype ival, char *msg)
{
  if (n == 0)
    return NULL;

  return idxset(n, ival, (idxtype *)GKmalloc(sizeof(idxtype)*n, msg));
}


/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(int nbytes, char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL)
    errexit("***Memory allocation failed for %s. Requested size: %d bytes", msg, nbytes);

  return ptr;
}
#endif

/*************************************************************************
* This function is my wrapper around free, allows multiple pointers
**************************************************************************/
void GKfree(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    free(*ptr1);
  *ptr1 = NULL;

  va_start(plist, ptr1);

  /* while ((int)(ptr = va_arg(plist, void **)) != -1) { */
  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
      free(*ptr);
    *ptr = NULL;
  }

  va_end(plist);
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
int *iset(int n, int val, int *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *idxset(int n, idxtype val, idxtype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
float *sset(int n, float val, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int iamax(int n, int *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax(int n, idxtype *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax_strd(int n, idxtype *x, int incx)
{
  int i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int samax(int n, float *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
int samax2(int n, float *x)
{
  int i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int idxamin(int n, idxtype *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int samin(int n, float *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum(int n, idxtype *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum_strd(int n, idxtype *x, int incx)
{
  int i, sum = 0;

  for (i=0; i<n; i++, x+=incx) {
    sum += *x;
  }

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void idxadd(int n, idxtype *x, idxtype *y)
{
  for (n--; n>=0; n--)
    y[n] += x[n];
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int charsum(int n, char *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int isum(int n, int *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum(int n, float *x)
{
  int i;
  float sum = 0.0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum_strd(int n, float *x, int incx)
{
  int i;
  float sum = 0.0;

  for (i=0; i<n; i++, x+=incx)
    sum += *x;

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void sscale(int n, float alpha, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] *= alpha;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float snorm2(int n, float *v)
{
  int i;
  float partial = 0;

  for (i = 0; i<n; i++)
    partial += v[i] * v[i];

  return sqrt(partial);
}



/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float sdot(int n, float *x, float *y)
{
  int i;
  float partial = 0;

  for (i = 0; i<n; i++)
    partial += x[i] * y[i];

  return partial;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
void saxpy(int n, float alpha, float *x, int incx, float *y, int incy)
{
  int i;

  for (i=0; i<n; i++, x+=incx, y+=incy)
    *y += alpha*(*x);
}




/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i
**************************************************************************/
void RandomPermute(int n, idxtype *p, int flag)
{
  int i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n <= 4)
    return;

  for (i=0; i<n; i+=16) {
    u = RandomInRangeFast(n-4);
    v = RandomInRangeFast(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
}



/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
int ispow2(int a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(int seed)
{
  if (seed == -1) {
#ifndef __VC__
    srand48(7654321L);
#endif
    srand(4321);
  }
  else {
#ifndef __VC__
    srand48(seed);
#endif
    srand(seed);
  }
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int log2(int a)
{
  int i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}

