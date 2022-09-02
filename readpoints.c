/*			    readpoints.c                        */
/* Author: J. Michael Word              Date Written: 1/27/99   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999-2002 J. Michael Word                       */
/*****************************************************************/

/* Modifications:                                               */
/*  3/ 1/99 - jmw - made it read 1D data as well as 2D and 3D   */
/*  8/ 1/02 - jmw - reads xplor maps                            */
/*  2/13/03 - jmw - fixed (- vs *) bug reading  xplor files     */
/*  8/30/22 - mgp - changed format ‘%d’ to %ld in printf w/ long int  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "readpoints.h"

#define BUFMAX 1000
char InputBuf[BUFMAX];

val3Dnode* read_pointval(FILE *inf, int ndim, int vfirst, char *delim) {
   char *word = NULL, *bufend = NULL;
   real_val col1 = 0.0, col2 = 0.0, col3 = 0.0, col4 = 0.0;
   int inputcols = 0, extraChars = 0;
   val3Dnode* v3D = NULL;
   
   if (ndim > 3) { ndim = 3; } /* force dimensions to 1, 2 or 3 */
   if (ndim < 1) { ndim = 1; }

   if (delim == NULL) { delim = " \t\n\r,:"; }

   do { /* loop-back point to try next line */

      bufend = &(InputBuf[BUFMAX-2]);
      *bufend = '\0';
      if (fgets(InputBuf, BUFMAX, inf) == NULL) {
	 return NULL;
      }
      extraChars = (*bufend != '\0' && *bufend != '\n');

      inputcols = 0;
      if ((word = strtok(InputBuf, delim)) != NULL) {
	 if (word[0] != '#') { /* pound sign marks a comment */
	    col1 = atof(word); inputcols = 1;
	    if ((word = strtok(NULL, delim)) != NULL) {
	       col2 = atof(word); inputcols = 2;
	       if ( (ndim >= 2)
	        && ((word = strtok(NULL, delim)) != NULL)) {
		  col3 = atof(word); inputcols = 3;
		  if (ndim >= 3) {
		     if ((word = strtok(NULL, delim)) != NULL) {
			col4 = atof(word); inputcols = 4;
		     }
		  }
	       }
	    }
	 }
      }
      if (inputcols >= (ndim+1)) {
	 v3D = (val3Dnode*)malloc(sizeof(val3Dnode));

	 if (v3D == NULL) {
	    fprintf(stderr, "ERROR: out of memory in read_pointval()\n");
	 }
	 else {
	    if (ndim == 1) {
	       if (vfirst) { /* two data orders for the 1D case */
		  v3D->v = col1;
		  v3D->x = col2;
	       }
	       else {
		  v3D->v = col2;
		  v3D->x = col1;
	       }
	       v3D->y = v3D->z = 0.0;
	    }
	    else if (vfirst) { /* two different data orders for 2D & 3D */
	       v3D->v = col1;
	       v3D->x = col2;
	       v3D->y = col3;
	       v3D->z = (ndim == 3) ? col4 : 0.0;
	    }
	    else {
	       v3D->x = col1;
	       v3D->y = col2;
	       v3D->z = (ndim == 3) ? col3 : 0.0;
	       v3D->v = (ndim == 3) ? col4 : col3;
	    }
	    v3D->next = NULL;
	 }
      }
      else if (inputcols > 0) { /* partial data */
	 fprintf(stderr, "WARNING: skipping partial input: %s\n",
	       InputBuf);
      }

      if (extraChars && v3D) { /* consume characters up-to newline */
	 int ch;
	 while(((ch = fgetc(inf)) != EOF) && (ch != '\n')) {
	    /* do nothing */
	 }
      }

      /* try again if we did not get a full set of data */
   } while (inputcols < (ndim+1));

   return v3D;
}

val3Dlist* load_coords(FILE *inf, int ndim, int vfirst, char *delim) {
   val3Dlist* vl  = NULL;
   val3Dnode* v3D = NULL;
   double vdev = 0.0, vsum = 0.0, vssq = 0.0;

   if (ndim > 3) { ndim = 3; } /* force dimensions to 1, 2 or 3 */
   if (ndim < 1) { ndim = 1; }

   vl = (val3Dlist*)malloc(sizeof(val3Dlist));
   if (vl == NULL) {
      fprintf(stderr, "ERROR: out of memory in load_coords()\n");
   }
   else {
      vl->xmin  = vl->xmax  = 0.0;
      vl->ymin  = vl->ymax  = 0.0;
      vl->zmin  = vl->zmax  = 0.0;
      vl->vmin  = vl->vmax  = 0.0;
      vl->vmean = vl->vsdev = 0.0;

      vl->head    = NULL;
      vl->mapinfo = NULL;
      vl->ncoords = 0;
      vl->ndim    = ndim;
      
      v3D = read_pointval(inf, ndim, vfirst, delim);

      vsum = 0.0;
      if (v3D != NULL) { /* initialize min & max to first record */
	 vl->xmin = v3D->x;
	 vl->xmax = v3D->x;
	 vl->ymin = v3D->y;
	 vl->ymax = v3D->y;
	 vl->zmin = v3D->z;
	 vl->zmax = v3D->z;
	 vl->vmin = v3D->v;
	 vl->vmax = v3D->v;
	 vl->vmean = v3D->v;
	 vl->vsdev = 0.0;
      }
      while(v3D != NULL) {
	    /* gather input range information */
	 if (vl->xmin > v3D->x) { vl->xmin = v3D->x; }
	 if (vl->xmax < v3D->x) { vl->xmax = v3D->x; }
	 if (vl->ymin > v3D->y) { vl->ymin = v3D->y; }
	 if (vl->ymax < v3D->y) { vl->ymax = v3D->y; }
	 if (vl->zmin > v3D->z) { vl->zmin = v3D->z; }
	 if (vl->zmax < v3D->z) { vl->zmax = v3D->z; }
	 if (vl->vmin > v3D->v) { vl->vmin = v3D->v; }
	 if (vl->vmax < v3D->v) { vl->vmax = v3D->v; }

	 vsum += v3D->v;

	 v3D->next = vl->head; /* link into list */
	 vl->head  = v3D;
	 vl->ncoords++;

	 v3D = read_pointval(inf, ndim, vfirst, delim); /* next rec */
      }
      
      vl->vmean = vsum / vl->ncoords;
      if (vl->ncoords > 1) {
	 vssq = 0.0;
	 for(v3D = vl->head; v3D; v3D = v3D->next) {
	    vdev = v3D->v - vl->vmean;
	    vssq += vdev*vdev;
	 }
	 vl->vsdev = sqrt((vssq)/(vl->ncoords - 1));
      }
   }
   
   return vl;
}

int read_xplor_header(FILE *inf, char* inbuf, int bufsz, val3Dlist* vl) {
   char *rc = NULL;
   int   i = 0, n = 0, na, amin, amax, nb, bmin, bmax, nc, cmin, cmax;
   double a, b, c, alpha, beta, gamma;

   /* first line blank so we ignore */
   if (! fgets(inbuf, bufsz, inf)) {
      fprintf(stderr, "ERROR: end of data reading header\n");
      return 0;
   }

   /* number of title lines */
   if ((! fgets(inbuf, bufsz, inf))
   ||  (! sscanf(inbuf, "%d",  &n)) ) {
      fprintf(stderr, "ERROR: end of data reading num title lines\n");
      return 0;
   }

   for(i = 0; i < n; i++) {
      if (! fgets(inbuf, bufsz, inf)) {
         fprintf(stderr, "ERROR: end of data reading title lines\n");
         return 0;
      }
      if (i < MAPINFO_MAX_TITLES) { /* save (trimmed) title lines */
         char* endp = NULL;
	 for (endp = inbuf + strlen(inbuf); endp != inbuf; endp--) {
	    if (strchr(" \t\r\n", *endp) == NULL) {
	       endp[1] = '\0';
	       break;
	    }
	 }
         vl->mapinfo->ntitles  = i+1;
         vl->mapinfo->title[i] = strdup(inbuf);
      }
   }

   /* read bounds */
   if ((! fgets(inbuf, bufsz, inf))
   ||  (sscanf(inbuf, "%d %d %d %d %d %d %d %d %d",
              &na, &amin, &amax, &nb, &bmin, &bmax, &nc, &cmin, &cmax) != 9) ) {
      fprintf(stderr, "ERROR: end of data reading bounds\n");
      return 0;
   }
   vl->mapinfo->grid[0]   = na;
   vl->mapinfo->grid[1]   = nb;
   vl->mapinfo->grid[2]   = nc;
   vl->mapinfo->origin[0] = amin;
   vl->mapinfo->origin[1] = bmin;
   vl->mapinfo->origin[2] = cmin;
   vl->mapinfo->extent[0] = amax;
   vl->mapinfo->extent[1] = bmax;
   vl->mapinfo->extent[2] = cmax;

   /* read unit cell */
   if ((! fgets(inbuf, bufsz, inf))
   ||  (sscanf(inbuf, "%12lf%12lf%12lf%12lf%12lf%12lf",
              &a, &b, &c, &alpha, &beta, &gamma) != 6) ) {
      fprintf(stderr, "ERROR: end of data reading unit cell\n");
      return 0;
   }
   vl->mapinfo->axis[0] = a;
   vl->mapinfo->axis[1] = b;
   vl->mapinfo->axis[2] = c;
   vl->mapinfo->angle[0]= alpha;
   vl->mapinfo->angle[1]= beta;
   vl->mapinfo->angle[2]= gamma;

   /* mode */
   if (! fgets(inbuf, bufsz, inf)) {
      fprintf(stderr, "ERROR: end of data reading mode\n");
      return 0;
   }
   if(((inbuf[0] != 'Z') && (inbuf[0] != 'z'))
   || ((inbuf[1] != 'Y') && (inbuf[1] != 'y'))
   || ((inbuf[2] != 'X') && (inbuf[2] != 'x'))) {
      fprintf(stderr, "ERROR: not an XPLOR \"ZYX\" file\n");
      return 0;
   }
   vl->mapinfo->fast = 0;
   vl->mapinfo->med  = 1;
   vl->mapinfo->slow = 2;

   update_mapinfo_factors(vl->mapinfo);

   return 1;
}

#define PER_XPLOR_LINE (6)

int read_xplor_slices(FILE *inf, char* inbuf, int bufsz, val3Dlist* vl,
                                                            double *sumptr) {
   val3Dnode* v3D = NULL;
   double vsum = 0.0, value[PER_XPLOR_LINE];
   double tx = 0.0, ty = 0.0, tz = 0.0;
   double cp = 0.0, cq = 0.0, cr = 0.0;
   double cs = 0.0, ct = 0.0, cu = 0.0;
   int slow    = 0, med    = 0, fast    = 0;
   int oslow   = 0, omed   = 0, ofast   = 0;
   int eslow   = 0, emed   = 0, efast   = 0;
   int islow   = 0, imed   = 0, ifast   = 0;
   int numslow = 0, nummed = 0, numfast = 0;
   int *slowptr = NULL, *medptr = NULL, *fastptr = NULL;

   int ksect = 0, per = 0, numfullrows = 0, numinremainder = 0, orthog = 0;
   int i = 0, j = 0, k = 0;
   int abc[3]; /* coordinate in a, b, c order */
   long nslice = 0L;

   cp = vl->mapinfo->cf[0];
   cq = vl->mapinfo->cf[1];
   cr = vl->mapinfo->cf[2];
   cs = vl->mapinfo->cf[3];
   ct = vl->mapinfo->cf[4];
   cu = vl->mapinfo->cf[5];
   orthog = vl->mapinfo->orthog;

   slow = vl->mapinfo->slow;
   med  = vl->mapinfo->med;
   fast = vl->mapinfo->fast;

   ofast = vl->mapinfo->origin[fast];
   omed  = vl->mapinfo->origin[med];
   oslow = vl->mapinfo->origin[slow];

   efast = vl->mapinfo->extent[fast];
   emed  = vl->mapinfo->extent[med];
   eslow = vl->mapinfo->extent[slow];

   numfast = efast - ofast + 1;
   nummed  = emed  - omed  + 1;
   numslow = eslow - oslow + 1;

   ksect = -1;

   nslice = numfast * nummed;

   numfullrows    = nslice / PER_XPLOR_LINE;
   numinremainder = nslice % PER_XPLOR_LINE;

   slowptr = abc + vl->mapinfo->slow;
   medptr  = abc + vl->mapinfo->med;
   fastptr = abc + vl->mapinfo->fast;
   
   vsum = 0.0;
   for(islow = oslow; islow <= eslow; islow++) {
      *slowptr = islow;

      if ((! fgets(inbuf, bufsz, inf))
      ||  (! sscanf(inbuf, "%d", &ksect)) ) {
         fprintf(stderr, "ERROR: end of data reading ksect number\n");
         return 0;
      }

      *medptr  = imed  = omed;
      *fastptr = ifast = ofast;

      for(j = 0; j < numfullrows; j++) {
         if ((! fgets(inbuf, bufsz, inf))
         ||  (sscanf(inbuf, "%12lf%12lf%12lf%12lf%12lf%12lf",
              &value[0], &value[1], &value[2], &value[3], &value[4], &value[5])
	         != PER_XPLOR_LINE) ) {
            fprintf(stderr, "ERROR: end of data reading ksect %d\n", ksect);
            return 0;
         }
         for(k = 0; k < PER_XPLOR_LINE; k++) {
	    if (orthog) {
	       tx = abc[0] * cp;
	       ty = abc[1] * cs;
	       tz = abc[2] * cu;
	    }
	    else {
	       tx = abc[0] * cp + abc[1] * cq + abc[2] * cr;
	       ty = abc[1] * cs + abc[2] * ct;
	       tz = abc[2] * cu;
	    }
	    vsum += value[k];
	    if (v3D = xyzv_to_node(tx, ty, tz, value[k])) {
	       v3D->next = vl->head; /* link into list */
	       vl->head  = v3D;
	       vl->ncoords++;
	    }
            *fastptr = ++ifast;
	    if(ifast > efast) {
               *fastptr =ifast = ofast;
	       *medptr = ++imed;
	    }
	 }
      }
      if (numinremainder > 0) {
         if ((! fgets(inbuf, bufsz, inf))
         ||  (sscanf(inbuf, "%12lf%12lf%12lf%12lf%12lf%12lf",
              &value[0], &value[1], &value[2], &value[3], &value[4], &value[5])
	         < numinremainder) ) {
            fprintf(stderr, "ERROR: end of data reading ksect %d remainder\n", ksect);
            return 0;
         }
         for(k = 0; k < numinremainder; k++) {
	    if (orthog) {
	       tx = abc[0] * cp;
	       ty = abc[1] * cs;
	       tz = abc[2] * cu;
	    }
	    else {
	       tx = abc[0] * cp + abc[1] * cq + abc[2] * cr;
	       ty = abc[1] * cs + abc[2] * ct;
	       tz = abc[2] * cu;
	    }
	    vsum += value[k];
	    if (v3D = xyzv_to_node(tx, ty, tz, value[k])) {
	       v3D->next = vl->head; /* link into list */
	       vl->head  = v3D;
	       vl->ncoords++;
	    }
            *fastptr = ++ifast;
	    if(ifast > efast) {
               *fastptr =ifast = ofast;
	       *medptr = ++imed;
	    }
	 }
      }
   }

   *sumptr = vsum;
   return 1;
}

val3Dnode* xyzv_to_node(double x, double y, double z, double v) {
   val3Dnode* v3D = NULL;

   v3D = (val3Dnode*)malloc(sizeof(val3Dnode));

   if (v3D == NULL) {
      fprintf(stderr, "ERROR: out of memory in xyzv_to_node()\n");
   }
   else {
      v3D->x = x;
      v3D->y = y;
      v3D->z = z;
      v3D->v = v;
      v3D->next = NULL;
   }
   return v3D;
}

#define DEG_TO_RAD   (0.0174532925)
#define CHOP_CUTOFF  (0.000000001)
#define REALABS(x)  (((x)<0.0)?-(x):(x))

int update_mapinfo_factors(mapInfo* minfo) {
   double cosalpha = 0.0, cosbeta = 0.0, cosgamma = 0.0, singamma = 0.0;
   double scsq = 0.0, vol = 0.0;
   double a_gap = 0.0, b_gap = 0.0, c_gap = 0.0;
   if (minfo) {
      cosalpha = cos(minfo->angle[0] * DEG_TO_RAD);
      cosbeta  = cos(minfo->angle[1] * DEG_TO_RAD);
      cosgamma = cos(minfo->angle[2] * DEG_TO_RAD);
      singamma = sin(minfo->angle[2] * DEG_TO_RAD);

      scsq = 1.0 - cosalpha * cosalpha
                 - cosbeta  * cosbeta
                 - cosgamma * cosgamma
                 + 2.0 * cosalpha * cosbeta * cosgamma;
      if (REALABS(scsq) < CHOP_CUTOFF) { scsq = 0.0; }

      vol = minfo->axis[0] * minfo->axis[1] * minfo->axis[2] * sqrt(scsq);
      a_gap = minfo->axis[0] / minfo->grid[0];
      b_gap = minfo->axis[1] / minfo->grid[1];
      c_gap = minfo->axis[2] / minfo->grid[2];
      
      /* output */
      minfo->gapsz[0] = a_gap;
      minfo->gapsz[1] = b_gap;
      minfo->gapsz[2] = c_gap;
      minfo->cf[0]  = a_gap;
      minfo->cf[1]  = b_gap * cosgamma;
      minfo->cf[2]  = c_gap * cosbeta;
      minfo->cf[3]  = b_gap * singamma;
      minfo->cf[4]  = c_gap * (cosalpha - cosbeta * cosgamma) / singamma;
      minfo->cf[5]  = c_gap * vol / (minfo->axis[0] *
                                     minfo->axis[1] *
			             minfo->axis[2] * singamma);
      minfo->orthog = ((REALABS(minfo->cf[1]) +
                        REALABS(minfo->cf[2]) +
			REALABS(minfo->cf[4])) < CHOP_CUTOFF);
      return 1;
   }
   return 0;
}

int mapindex_to_xyz(val3Dlist* vl, int ia, int ib, int ic,
                     double *x, double *y, double *z) {
   double cp = 0.0, cq = 0.0, cr = 0.0;
   double cs = 0.0, ct = 0.0, cu = 0.0;

   if (vl && vl->mapinfo) {
      cp = vl->mapinfo->cf[0];
      cq = vl->mapinfo->cf[1];
      cr = vl->mapinfo->cf[2];
      cs = vl->mapinfo->cf[3];
      ct = vl->mapinfo->cf[4];
      cu = vl->mapinfo->cf[5];

      *x = ia * cp + ib * cq + ic * cr;
      *y = ib * cs + ic * ct;
      *z = ic * cu;

      return 1;
   }
   return 0;
}


#define BUFFERSZ (256)

val3Dlist* load_xplor_map(FILE *inf) {
   val3Dlist* vl  = NULL;
   val3Dnode* v3D = NULL;
   double vdev = 0.0, vsum = 0.0, vssq = 0.0;

   vl = (val3Dlist*)malloc(sizeof(val3Dlist));
   if (vl == NULL) {
      fprintf(stderr, "ERROR: out of memory in load_xplor_map()\n");
   }
   else {
      vl->xmin  = vl->xmax  = 0.0;
      vl->ymin  = vl->ymax  = 0.0;
      vl->zmin  = vl->zmax  = 0.0;
      vl->vmin  = vl->vmax  = 0.0;
      vl->vmean = vl->vsdev = 0.0;

      vl->head    = NULL;
      vl->mapinfo = NULL;
      vl->ncoords = 0;
      vl->ndim    = 3;

      vl->mapinfo = (mapInfo*)malloc(sizeof(mapInfo));
      if (vl->mapinfo == NULL) {
         fprintf(stderr, "ERROR: out of mapinfo memory in load_xplor_map()\n");
	 destroy3Dlist(vl);
	 return NULL;
      }

      /* initialize mapinfo */
      vl->mapinfo->axis[0] = 0.0;
      vl->mapinfo->axis[1] = 0.0;
      vl->mapinfo->axis[2] = 0.0;
      vl->mapinfo->angle[0]= 0.0;
      vl->mapinfo->angle[1]= 0.0;
      vl->mapinfo->angle[2]= 0.0;
      vl->mapinfo->grid[0]   = 0;
      vl->mapinfo->grid[1]   = 0;
      vl->mapinfo->grid[2]   = 0;
      vl->mapinfo->origin[0] = 0;
      vl->mapinfo->origin[1] = 0;
      vl->mapinfo->origin[2] = 0;
      vl->mapinfo->extent[0] = 0;
      vl->mapinfo->extent[1] = 0;
      vl->mapinfo->extent[2] = 0;
      vl->mapinfo->fast = 0;
      vl->mapinfo->med  = 1;
      vl->mapinfo->slow = 2;
      vl->mapinfo->ntitles = 0;
      vl->mapinfo->title[0] = NULL;
      vl->mapinfo->gapsz[0] = 0.0;
      vl->mapinfo->gapsz[1] = 0.0;
      vl->mapinfo->gapsz[2] = 0.0;
      vl->mapinfo->cf[0]   = 0.0;
      vl->mapinfo->cf[1]   = 0.0;
      vl->mapinfo->cf[2]   = 0.0;
      vl->mapinfo->cf[3]   = 0.0;
      vl->mapinfo->cf[4]   = 0.0;
      vl->mapinfo->cf[5]   = 0.0;

      if ((! read_xplor_header(inf, InputBuf, BUFMAX, vl))
      ||  (! read_xplor_slices(inf, InputBuf, BUFMAX, vl, &vsum))) {
         fprintf(stderr, "ERROR: stopped loading XPLOR map\n");
	 destroy3Dlist(vl);
	 return NULL;
      }

      vl->vmean = vsum / vl->ncoords;
      if (vl->ncoords > 1) {
	 vssq = 0.0;
         if (vl->head != NULL) { /* initialize min & max to first record */
	    vl->xmin = vl->head->x;
	    vl->xmax = vl->head->x;
	    vl->ymin = vl->head->y;
	    vl->ymax = vl->head->y;
	    vl->zmin = vl->head->z;
	    vl->zmax = vl->head->z;
	    vl->vmin = vl->head->v;
	    vl->vmax = vl->head->v;
	    vl->vmean = vl->head->v;
	    vl->vsdev = 0.0;
         }
	 for(v3D = vl->head; v3D; v3D = v3D->next) {
	    if (vl->xmin > v3D->x) { vl->xmin = v3D->x; }
	    if (vl->xmax < v3D->x) { vl->xmax = v3D->x; }
	    if (vl->ymin > v3D->y) { vl->ymin = v3D->y; }
	    if (vl->ymax < v3D->y) { vl->ymax = v3D->y; }
	    if (vl->zmin > v3D->z) { vl->zmin = v3D->z; }
	    if (vl->zmax < v3D->z) { vl->zmax = v3D->z; }
	    if (vl->vmin > v3D->v) { vl->vmin = v3D->v; }
	    if (vl->vmax < v3D->v) { vl->vmax = v3D->v; }

	    vdev = v3D->v - vl->vmean;
	    vssq += vdev*vdev;
	 }
	 vl->vsdev = sqrt((vssq)/(vl->ncoords - 1));
      }
   }

   return vl;
}

void destroy3Dlist(val3Dlist* vl) {
   val3Dnode *curr = NULL;
   val3Dnode *nxt = NULL;

   if (vl != NULL) {

      if (vl->head != NULL) {

	 curr = vl->head;
	 while (curr) {
	    nxt = curr->next;
	    curr->next = NULL;
	    free(curr);
	    curr = nxt;
	 }
	 vl->head = NULL;
	 destroyMapinfo(vl->mapinfo);
	 vl->mapinfo = NULL;
	 free(vl);
      }
   }
}

void destroyMapinfo(mapInfo* minfo) {
   int i = 0;

   if (minfo != NULL) {
      for (i = 0; i < minfo->ntitles; i++) {
	 if (minfo->title[i] != NULL) {
	    free(minfo->title[i]);
	    minfo->title[i] = NULL;
	 }
      }
      free(minfo);
   }
}

void print3Dlist(FILE *outf, val3Dlist* vl, int showVals) {
   val3Dnode* v3D = NULL;
   int i;

   if (vl == NULL) {
      fprintf(outf, "*val3Dlist: NULL\n");
   }
   else if (vl->head == NULL) {
      fprintf(outf, "*val3Dlist: EMPTY\n");
   }
   else {
      fprintf(outf, "*val3Dlist: %ld values, %dD\n",
	 vl->ncoords, vl->ndim);
      fprintf(outf, "   x range: %g .. %g\n", vl->xmin, vl->xmax);
      if (vl->ndim >= 2) {
	 fprintf(outf, "   y range: %g .. %g\n", vl->ymin, vl->ymax);
      }
      if (vl->ndim >= 3) {
	 fprintf(outf, "   z range: %g .. %g\n", vl->zmin, vl->zmax);
      }
      fprintf(outf, "   v range: %g .. %g\n", vl->vmin, vl->vmax);
      fprintf(outf, "   v  mean: %g, sdev: %g\n", vl->vmean, vl->vsdev);
      if (showVals) {
	 i = 0;
	 for(v3D = vl->head; v3D; v3D = v3D->next) {
	    if (vl->ndim == 3) {
	       fprintf(outf, "   %d) [%g, %g, %g] = %g\n", ++i,
		  v3D->x, v3D->y, v3D->z, v3D->v);
	    }
	    else {
	       fprintf(outf, "   %d) [%g, %g] = %g\n", ++i,
		  v3D->x, v3D->y, v3D->v);
	    }
	 }
      }
   }
}
