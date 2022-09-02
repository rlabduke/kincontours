/*			    readpoints.h                        */
/* Author: J. Michael Word              Date Written: 1/28/99   */

/* Note: usage requires linking to the math library -lm         */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999-2002 J. Michael Word                       */
/*****************************************************************/

#ifndef READPOINTS_H
#define READPOINTS_H 1

#include "nrmat.h"

typedef struct val3Dnode_t {
   real_val            x, y, z;
   real_val            v;       /* value assigned to the point */
   struct val3Dnode_t* next;    /* linked list pointer */
} val3Dnode;

#define MAPINFO_MAX_TITLES (20)
#define XPLOR_MAPTYPE (1)

typedef struct mapInfo_t {
   double   axis[3];          /* unit cell axes   */
   double   angle[3];         /* unit cell angles */
   int      grid[3];          /* num gaps in unit cell        */
   int      origin[3];        /* data starts at this position */
   int      extent[3];        /* data ends here               */
   int      fast, med, slow;  /* which dimension is which     */
   int      ntitles;          /* number of title lines */
   char*    title[MAPINFO_MAX_TITLES];
   double   gapsz[3];         /* grid spacing */
   double   cf[6];            /* grid to xyz conversion factors */
   int      orthog;           /* are axes orthogonal? */
} mapInfo;

typedef struct coordlist_t {
   real_val   xmin, xmax; /* data ranges */
   real_val   ymin, ymax;
   real_val   zmin, zmax;
   real_val   vmin, vmax;
   real_val   vmean, vsdev; /* mean and sqrt(variance) of values */

   long       ncoords;  /* number of values    */
   int        ndim;     /* either 2 or 3       */
   val3Dnode* head;     /* linked list pointer */
   mapInfo*   mapinfo;  /* density map header  */
} val3Dlist;

val3Dlist* load_coords(FILE *inf, int ndim, int vfirst, char *delim);
void destroy3Dlist(val3Dlist* vl);
void print3Dlist(FILE *outf, val3Dlist* vl, int showVals);

val3Dnode* read_pointval(FILE *inf, int ndim, int vfirst, char *delim);

val3Dlist* load_xplor_map(FILE *inf);
int read_xplor_header(FILE *inf, char* buffer, int bufsz, val3Dlist* vl);
int read_xplor_slices(FILE *inf, char* buffer, int bufsz, val3Dlist* vl,
                                                            double *sum);
val3Dnode* xyzv_to_node(double x, double y, double z, double v);
int update_mapinfo_factors(mapInfo* minfo);
int mapindex_to_xyz(val3Dlist* vl, int ia, int ib, int ic,
                     double *x, double *y, double *z);
void destroyMapinfo(mapInfo* minfo);

#endif
