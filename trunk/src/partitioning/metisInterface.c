/*
       ******************************************************************
       *                                                                *
       * File:          metisInterface.c                                *
       * Author:        Edwin van der Weide                             *
       * Starting date: 02-12-2003                                      *
       * Last modified: 07-07-2005                                      *
       *                                                                *
       ******************************************************************
*/

#include "metis.h"

/*
       ******************************************************************
       *                                                                *
       * metisInterface is an interface between the application         *
       * program, either written in Fortran or C, and the               *
       * multi-constraint graph partitioning routines of metis.         *
       *                                                                *
       ******************************************************************
*/

void metisInterface(int *n, int *ncon, idxtype *xadj, idxtype *adjncy,
                    idxtype *vwgt, idxtype *adjwgt, int *wgtflag,
                    int *numflag, int *nparts, float *ubvec,
                    int *options, int *edgecut, idxtype *part)
{
  /* According to the Metis manual, the graph partitioner should be  */
  /* used when the number of partitions is greater than 8. Otherwise */
  /* the recursive bisection is to be preferred.                     */
  /* Take care of the special case nparts == 1 here.                 */

  if(*nparts == 1)
  {
    int i;
    for(i=0; i<(*n); i++) part[i] = 0;
  }
  else if(*nparts > 8)
    METIS_mCPartGraphKway(n, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag,
                          numflag, nparts, ubvec, options, edgecut,
                          part);
  else
    METIS_mCPartGraphRecursive(n, ncon, xadj, adjncy, vwgt, adjwgt,
                               wgtflag, numflag, nparts, options,
                               edgecut, part);
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void METISINTERFACE(int *n, int *ncon, idxtype *xadj, idxtype *adjncy,
                    idxtype *vwgt, idxtype *adjwgt, int *wgtflag,
                    int *numflag, int *nparts, float *ubvec,
                    int *options, int *edgecut, idxtype *part)
{
  metisInterface(n, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag,
                 numflag, nparts, ubvec, options, edgecut, part);
}

void metisinterface_(int *n, int *ncon, idxtype *xadj, idxtype *adjncy,
                     idxtype *vwgt, idxtype *adjwgt, int *wgtflag,
                     int *numflag, int *nparts, float *ubvec,
                     int *options, int *edgecut, idxtype *part)
{
  metisInterface(n, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag,
                 numflag, nparts, ubvec, options, edgecut, part);
}

void metisinterface(int *n, int *ncon, idxtype *xadj, idxtype *adjncy,
                    idxtype *vwgt, idxtype *adjwgt, int *wgtflag,
                    int *numflag, int *nparts, float *ubvec,
                    int *options, int *edgecut, idxtype *part)
{
  metisInterface(n, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag,
                 numflag, nparts, ubvec, options, edgecut, part);
}

void metisinterface__(int *n, int *ncon, idxtype *xadj, idxtype *adjncy,
                      idxtype *vwgt, idxtype *adjwgt, int *wgtflag,
                      int *numflag, int *nparts, float *ubvec,
                      int *options, int *edgecut, idxtype *part)
{
  metisInterface(n, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag,
                 numflag, nparts, ubvec, options, edgecut, part);
}
