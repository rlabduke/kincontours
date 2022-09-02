/*			    kin3Dcont.c                         */
/* Author: J. Michael Word              Date Written: 1/27/99   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999-2003 J. Michael Word                       */
/*****************************************************************/

/* Modifications:                                               */
/*  2/ 3/99 - jmw - made input robust to blanks, added # as     */
/*                  comment and supported processing one pt/val */
/*  2/18/99 - jmw - added -wrap flag to support periodic bounds */
/*  2/22/99 - jmw - added -sampled flag for sampled datasets    */
/*                  also considers average for default contour  */
/*  2/24/99 - jmw - default contours +/- 1 sigma smoothed data  */
/*  2/26/99 - jmw - modified -dump function                     */
/*  3/ 6/99 - jmw - fixed minor formatting problem in caption   */
/*  4/30/02 - jmw - added -nosmooth flag                        */
/*  7/24/02 - jmw - added -alignminxyz flag                     */
/*  7/25/02 - jmw - added -SNAP, -SL, -SM, and fixed -SP        */
/*  7/29/02 - jmw - added -XPLOR and -NOPERSP                   */
/*  8/ 1/02 - jmw - got xplor non-orthogonal unit cells to work */
/*                  and added map info to caption               */
/*  2/12/03 - jmw - fixed bug with non-orthogonal maps          */
/*  2/13/03 - jmw - added clipping!                             */
/*  8/30/22 - mgp - changed format ‘%d’ to %ld in printf w/ long int  */
/*  8/30/22 - mgp - added prototype for remap_slices                  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "connect.h"
#include "readpoints.h"

#define MAXLEVELS 300

#define DEFAULT_DENS  1.0
#define DEFAULT_FSDEV 2.0
#define CHOP_EPSILON  0.00000000001

/* span is filter cuttoff in standard deviations */
#define DEFAULT_SPAN 2.0

#define PROGRAM_NAME "kin3Dcont"
//#define VERSION_STR "kin3Dcont: version 1.12, Feb 13, 2003, Copyright 1999-2003, J. Michael Word"
#define VERSION_STR "kin3Dcont: version 1.13, Aug 30, 2022, Copyright 1999-2022, J. Michael Word and Duke University"

typedef struct progparms_t {
   int       verbose;        /* print out descriptive info? */
   int       usesmoothing;   /* do the smoothing thing? */
   char*     infname;        /* input file name */

   int       nlevels;        /* contour level & style */
   ctree_lvl levels[MAXLEVELS]; /* what we use (may have been transformed) */
   ctree_lvl inlvls[MAXLEVELS]; /* what was input on the command line      */
   char*     colors[MAXLEVELS];
   int       sigmalevels;    /* true if levels input as sigma */

   int       wrap;           /* do the edges wrap arround? */

   coord_val dx,  dy,  dz;   /* lattice spacing  */
   coord_val sx,  sy,  sz;   /* stddev of filter */
   coord_val filtspan;       /* filter cuttoff   */

   int       sampled;        /* sampled data (not normalized) */
   int       alignGridToMin; /* grid aligned with min value */
   int       snap;           /* if data very close to grid make it on grid */
   
   int        vfirst;        /* value first (else last) */
   char*      delim;         /* special delimiter       */
   val3Dlist* vl;            /* list of input points    */

   int       writekin;       /* display @kinemage */
   int       writegrp;       /* display @group    */
   int       writesub;       /* display @subgroup */
   int       dominant;       /* make the grp or subgrp dominant? */
   int       perspective;    /* use the perspective flag? */
   int       lens;           /* add lens keyword? */
   int       axis;           /* draw box edges? */
   char*     grpname;        /* group or subgroup name */

   int       dumplattice;    /* write lattice values for debugging */

   coord_val xmin, ymin, zmin; /* bounds for data when */
   coord_val xmax, ymax, zmax; /* points wrap around   */

   int       clipout;                /* clip output around a point */
   coord_val xout, yout, zout, rout; /* output center and radius */
   coord_val rout_sq;                /* and radius (squared) */

   /* parameters calculated from input */

   int       ifx,  ify, ifz; /* integer filter radii */
   int       nfx,  nfy, nfz; /* number of filter points   */
   int       nx,   ny,  nz;  /* number of lattice points  */
   coord_val x0,   y0,  z0;  /* coord of lattice[0][0][0] */

   real_val ***lat;          /* three dimensional data lattice */
   real_val ***filt;         /* three dimensional filter */
} progparms;

typedef struct clip_box_t {
   coord_val xout, yout, zout, rout; /* output center and radius */
   coord_val rout_sq;                /* and radius (squared) */

   coord_val xb, yb, zb; /* bounding box beginning (real) */
   coord_val xe, ye, ze; /* bounding box end       (real) */

   int xo,   yo,   zo;   /* bounding box origin (int) */
   int xf,   yf,   zf;   /* bounding box finish (int) */
   int xrgn, yrgn, zrgn; /* bounding box region (int) */
} clip_box;

val3Dlist* processCommandline(int argc, char *argv[], progparms *parms);
void draw_axis_box(progparms *parms);
void kin_draw_box(coord_val xmin, coord_val ymin, coord_val zmin,
                  coord_val xmax, coord_val ymax, coord_val zmax);
void kin_draw_new_line(coord_val x1,coord_val y1,coord_val z1,
                       coord_val x2,coord_val y2,coord_val z2);
void kin_continue_line(coord_val x1,coord_val y1,coord_val z1,
                       coord_val x2,coord_val y2,coord_val z2);
void kin_rombus_draw_cell(val3Dlist* vl,
                          int a_lo, int b_lo, int c_lo,
                          int a_hi, int b_hi, int c_hi);
void kin_rombus_draw_new_line(val3Dlist* vl,
                          int a_beg, int b_beg, int c_beg,
                          int a_end, int b_end, int c_end);
void kin_rombus_continue_line(val3Dlist* vl,
                          int a_beg, int b_beg, int c_beg,
                          int a_end, int b_end, int c_end);
void analyze_data(progparms *parms);
int insert_data_into_lattice(progparms *parms);
void splat(progparms *parms, real_val wt, coord_val x, coord_val y, coord_val z);
void placemask(progparms *parms, real_val wt, int lx, int ly, int lz);
void placemaskwrap(progparms *parms, real_val wt, int lx, int ly, int lz);
void placesample(progparms *parms, real_val wt, int lx, int ly, int lz);
void input_sigma_to_contour_levels(progparms *parms);
void calc_default_contour_levels(progparms *parms);
void convert_sigma_contour_levels(progparms *parms);
void contour_lattice(progparms *parms);
real_val*** build_lattice(progparms *parms);
real_val*** build_filter(progparms *parms);
void destroy_lattice(progparms *parms);
void destroy_filter(progparms *parms);
void describeInput(progparms *parms);
void helpInfo(int fullhelp);
int compArgStr(char *str, char *arg, int min);
void dup_periodic_edges(progparms *parms);
void build_clip_box(progparms *parms, clip_box *box);
void xy_contour(progparms *parms, cont_3d_info *cp);
void xz_contour(progparms *parms, cont_3d_info *cp);
void yz_contour(progparms *parms, cont_3d_info *cp);
void dump_lattice(progparms *parms);

int remap_slices(cont_3d_info *cp, slice_info *sp);

int
main(int argc, char *argv[]) {
   progparms parms;

   if (processCommandline(argc, argv, &parms)) {
      analyze_data(&parms);
   }
   return 0;
}

void
analyze_data(progparms *parms) {

   if (parms) {
      if (parms->dumplattice) {
	 /* override -dominant flag to suppress -dump header*/
	 if (! (parms->dominant)) {
	    fprintf(stdout, "#dump from %s\n", VERSION_STR);
	    fprintf(stdout, "#input: %s, %ld coords\n", parms->infname, parms->vl->ncoords);
	    fprintf(stdout, "#lattice gap: [%g, %g, %g]\n", parms->dx, parms->dy, parms->dz);
            if (parms->usesmoothing) {
	       fprintf(stdout, "#filter stddev: [%g, %g, %g]\n", parms->sx, parms->sy, parms->sz);
            }
            else {
               fprintf(stdout, "#no smoothing\n");
            }
	    if (parms->wrap) {
	       fprintf(stdout, "#bounding box: x[%g .. %g] y[%g .. %g] z[%g .. %g] *wrap*\n",
		  parms->xmin,  parms->xmax - parms->dx,
		  parms->ymin,  parms->ymax - parms->dy,
		  parms->zmin,  parms->zmax - parms->dz);
	    }
	    else {
	       fprintf(stdout, "#bounding box: x[%g .. %g] y[%g .. %g] z[%g .. %g]\n",
		  parms->xmin,  parms->xmax,
		  parms->ymin,  parms->ymax,
		  parms->zmin,  parms->zmax);
	    }
	 }
      }
      else {
	 if (parms->writekin) {
	    fprintf(stdout, "@kinemage 1\n");
	 }
	 fprintf(stdout, "@caption %s\n", VERSION_STR);
	 fprintf(stdout, " input: %s, %ld coords\n", parms->infname, parms->vl->ncoords);
	 fprintf(stdout, " lattice gap: [%g, %g, %g]", parms->dx, parms->dy, parms->dz);
	 if (parms->usesmoothing) {
	    fprintf(stdout, ", filter stddev: [%g, %g, %g]", parms->sx, parms->sy, parms->sz);
         }
	 else {
	    fprintf(stdout, ", no smoothing");
	 }
	 if (parms->alignGridToMin) {
	    fprintf(stdout, ", alignminxyz");
         }
	 if (parms->snap) {
	    fprintf(stdout, ", snaptogrid");
         }
	 if (parms->sampled) {
	    fprintf(stdout, ", sampled");
         }
	 fprintf(stdout, "\n");

	 if (parms->wrap) {
	    fprintf(stdout, " wrap box: x[%g .. %g] y[%g .. %g] z[%g .. %g]\n",
	            parms->xmin,  parms->xmax,
	            parms->ymin,  parms->ymax,
	            parms->zmin,  parms->zmax);
	 }

	 if (parms->sigmalevels) {
	    int n = 0;
	    fprintf(stdout, " levels input as sigma values:");
            for (n = 0; n < parms->nlevels; n++) {
	       fprintf(stdout, " %g %s", parms->inlvls[n], parms->colors[n]);
	    }
	    fprintf(stdout, "\n");
	 }

         if (parms->vl && parms->vl->mapinfo) {
	    int thetitle = 0;
	    fprintf(stdout, " from map file:\n");
            for (thetitle = 0;
	         thetitle < parms->vl->mapinfo->ntitles;
	         thetitle++ ) {
               if (parms->vl->mapinfo->title[thetitle]) {
	          fprintf(stdout, " --%s\n",
		          parms->vl->mapinfo->title[thetitle]);
               }
            }
	    fprintf(stdout, "    grid[%d, %d, %d]\n",
                    parms->vl->mapinfo->grid[0],
                    parms->vl->mapinfo->grid[1],
                    parms->vl->mapinfo->grid[2]);
	    fprintf(stdout, "    origin[%d, %d, %d]\n",
                    parms->vl->mapinfo->origin[0],
		    parms->vl->mapinfo->origin[1],
		    parms->vl->mapinfo->origin[2]);
	    fprintf(stdout, "    extent[%d, %d, %d]\n",
		    parms->vl->mapinfo->extent[0],
		    parms->vl->mapinfo->extent[1],
		    parms->vl->mapinfo->extent[2]);
	    fprintf(stdout, "    cell_axes[%.6g, %.6g, %.6g]\n",
		    parms->vl->mapinfo->axis[0],
		    parms->vl->mapinfo->axis[1],
		    parms->vl->mapinfo->axis[2]);
	    fprintf(stdout, "    cell_angles[%.6g, %.6g, %.6g]\n",
		    parms->vl->mapinfo->angle[0],
		    parms->vl->mapinfo->angle[1],
		    parms->vl->mapinfo->angle[2]);
	    fprintf(stdout, "    fast,med,slow dim [%d, %d, %d]\n",
		    parms->vl->mapinfo->fast + 1,
		    parms->vl->mapinfo->med + 1,
		    parms->vl->mapinfo->slow + 1);
         }
	 if (parms->perspective) { fprintf(stdout, "@perspective\n"); }
	 if (parms->writegrp) {
	    fprintf(stdout, "@group {%s}%s%s\n",
	       parms->grpname, (parms->dominant?" dominant":""),
	       (parms->lens?" lens":""));
	 }
	 else if (parms->writesub) {
	    fprintf(stdout, "@subgroup {%s}%s%s\n",
	       parms->grpname, (parms->dominant?" dominant":""),
	       (parms->lens?" lens":""));
	 }
	 if (parms->axis) { draw_axis_box(parms); }
      }

      /* create and initialize lattice */
      if (build_lattice(parms)) {

	 if (build_filter(parms)) {

	    insert_data_into_lattice(parms);

	    destroy3Dlist(parms->vl);  parms->vl   = NULL;
	    destroy_filter(parms);     parms->filt = NULL;

	    if (parms->dumplattice) {
	       dump_lattice(parms);
	    }
	    else {
	       contour_lattice(parms);
	    }
	 }

	 destroy_lattice(parms); parms->lat = NULL;
      }
   }
}

void helpInfo(int fullhelp) {
   fprintf(stderr, "\nSyntax: %s [-flags] [inputfilename | -] >> out.kin\n\n",
                    PROGRAM_NAME);
if (fullhelp) {

   fprintf(stderr, "Read x, y, z, value or value, x, y, z sample data and\n");
   fprintf(stderr, "place on a 3D grid using a gaussian spot function and\n");
   fprintf(stderr, "generate one or more 3D contour grids in kinemage format.\n\n");
   fprintf(stderr, "Flags:\n");
}
else {
   fprintf(stderr, "Key Flags: (use -help to see more)\n");
}
   fprintf(stderr, "\n(contour level specification)\n");
   fprintf(stderr, "  -Level #.# color [...]          contour at one or more levels\n");
   fprintf(stderr, "  -Multi #.# #.# #.# color [...]  contour From To By Color\n");
   fprintf(stderr, "     where [...] represents optional additional level specifications\n");
   fprintf(stderr, "  -SLevel #.# color [...]         same as -Level except levels are given\n");
   fprintf(stderr, "                                  in standard deviations from the mean\n");
   fprintf(stderr, "  -SMulti #.# #.# #.# color [...] same as -Level except levels are given\n");
   fprintf(stderr, "                                  in standard deviations from the mean\n");
   fprintf(stderr, "\n(discrete lattice and spot filter)\n");
   fprintf(stderr, "  -Grain#.#         spacing between lattice points (default %g)\n", DEFAULT_DENS);
if (fullhelp) {
   fprintf(stderr, "  -GXYZ #.# #.# #.# spacing between lattice points in each dim\n");
}
   fprintf(stderr, "  -Stddev#.#        spot filter stddev (default %g)\n", DEFAULT_FSDEV);
if (fullhelp) {
   fprintf(stderr, "  -SXYZ #.# #.# #.# spot filter stddev in each dim\n");
}
   fprintf(stderr, "  -WRAP #.# #.# #.# #.# #.# #.#  bounds of box where edges wrap\n");
   fprintf(stderr, "                     (6 numbers: xmin xmax ymin ymax zmin zmax)\n");
   fprintf(stderr, "  -SNAPtogrid       force to grid data points which are very close to grid\n");
   fprintf(stderr, "  -SAMPled          input data is a uniform grid of samples\n");
   fprintf(stderr, "  -ALIGNminxyz      sampled and snapped with grid aligned to min xyz\n");
   fprintf(stderr, "\n(input source)\n");
   fprintf(stderr, "  -XPLORascii       XPLOR format ascii map file\n");
   fprintf(stderr, "  -FIRST            value in first column then x,y,z (default)\n");
   fprintf(stderr, "  -LAST             value in last column after x,y,z\n");
   fprintf(stderr, "  -CLIP #.# #.# #.# #.#  clip output around a point\n");
   fprintf(stderr, "                     (4 numbers: x,y,z and halfwidth, in that order)\n");
if (fullhelp) {
   fprintf(stderr, "  -                 a dash represents the standard input stream\n");
   fprintf(stderr, "  -INput filename   alternate way to identify the input file\n");
   fprintf(stderr, "  by default, input is read from standard input\n");
   fprintf(stderr, "\n(other)\n");
   fprintf(stderr, "  -KINemage         writes @kinemage in header\n");
   fprintf(stderr, "  -GROUP            writes @group in header\n");
   fprintf(stderr, "  -SUBgroup         writes @subgroup in header\n");
   fprintf(stderr, " (-group and -subgroup are mutually exculsive.)\n");
   fprintf(stderr, "  -NAMe string      identifier for @group or @subgroup\n");
   fprintf(stderr, "  -NOPERSPective    do not include @perspective statement\n");
   fprintf(stderr, "  -DOMinant         add dominant keyword to @group or @subgroup\n");
   fprintf(stderr, "                    (supresses header when used with -dump)\n");
   fprintf(stderr, "  -LENS             add lens keyword to @group or @subgroup\n");
   fprintf(stderr, "  -NOAXIS           do not draw axis lines marking the box edges\n");
   fprintf(stderr, "  -T_               Use _ as the field separator character(s).\n");
   fprintf(stderr, "                    For example -t@ would separate on at-signs.\n");
   fprintf(stderr, "                    In the input, multiple consecutive separators count as one.\n");
   fprintf(stderr, "(The default delimiters are: blank, tab, newline, return, comma and colon)\n");
   fprintf(stderr, "  -NOSMooth         do not attempt to smooth the input data\n");
   fprintf(stderr, "  -SPan#.#          spot filter span in stddev (default %g)\n", DEFAULT_SPAN);
   fprintf(stderr, "  -DUMPlattice      output lattice data for debugging\n\n");
   fprintf(stderr, "  -Verbose          displays extra info about the calculation\n");
   fprintf(stderr, "  -Quiet            suppresses extra output\n");
   fprintf(stderr, "\n  -Help             prints this help message\n");
   fprintf(stderr, "\nExample: %s -g1.5 -s3 -l 125 sea mydata >>my.kin\n", PROGRAM_NAME);
   fprintf(stderr, "generates a single 3D contour with a spacing of 1.5 and a filter stddev of 3\n");
}
else {
   fprintf(stderr, "\n  -Help             prints full help message (additional info!)\n");
}
   fprintf(stderr, "\n%s\n", VERSION_STR);
   exit(1);
}

void describeInput(progparms *parms) {

   fprintf(stderr, "input file: \"%s\"\n", parms->infname);
   fprintf(stderr, "    levels: n=%d\n", parms->nlevels);
   fprintf(stderr, "     delta: %g, %g, %g\n", parms->dx, parms->dy, parms->dz);
   if (parms->usesmoothing) {
      fprintf(stderr, "filt. sdev: %g, %g, %g\n", parms->sx, parms->sy, parms->sz);
      fprintf(stderr, "filt. span: %g stddev\n",  parms->filtspan);

      fprintf(stderr, "filter dim: %d .. %d, %d .. %d, %d .. %d (%d total)\n",
         -parms->ifx, parms->ifx,
         -parms->ify, parms->ify,
         -parms->ifz, parms->ifz, parms->nfx*parms->nfy*parms->nfz);
   }
   else {
      fprintf(stderr, " no smoothing\n");
   }

   fprintf(stderr, "  data dim: %d, %d, %d (%d total)\n",
      parms->nx,  parms->ny,  parms->nz, parms->nx*parms->ny*parms->nz);
   fprintf(stderr, "    origin: %g, %g, %g\n", parms->x0,  parms->y0,  parms->z0);
   if (parms->wrap) {
      fprintf(stderr, "   *wrap x: %g .. %g\n", parms->xmin,  parms->xmax);
      fprintf(stderr, "   *wrap y: %g .. %g\n", parms->ymin,  parms->ymax);
      fprintf(stderr, "   *wrap z: %g .. %g\n", parms->zmin,  parms->zmax);
   }

   print3Dlist(stderr, parms->vl, 0);
}

int compArgStr(char *str, char *arg, int min) {
   int i, max;
   char s, a;

   if (!str || !arg) return 0;

   max = strlen(arg);

   for(i=0; i<max; i++) {
      s = toupper(str[i]);
      a = toupper(arg[i]);

      if (i >= min && (s == '\0' || s == '.'
         || s == '+' || s == '-' || isdigit(s))) {
	 break; /* good ending point */
      }
      else if (s != a) {
	 i = 0; /* failed to match */
	 break;
      }
   }

   return i;
}

val3Dlist*
processCommandline(int argc, char *argv[], progparms *parms) {
   FILE *inf = stdin;
   char *p = NULL;
   int i = 0, n = 0, nfiles = 0;
   coord_val rval = 0.0;
   int mapinput = 0;

   if (parms == NULL) {
      fprintf(stderr, "ERROR: processCommandline(..., parms=NULL)\n");
      return NULL;
   }

   /* initialize parameters */

   parms->verbose = 0;
   parms->usesmoothing = 1;
   parms->infname = "-- none --";
   parms->nlevels = 0;
   parms->levels[0] = 0.0;
   parms->inlvls[0] = 0.0;
   parms->colors[0] = NULL;
   parms->sigmalevels = 0;
   parms->wrap = 0;
   parms->dx  = parms->dy  = parms->dz  = DEFAULT_DENS;
   parms->sx  = parms->sy  = parms->sz  = DEFAULT_FSDEV;
   parms->filtspan = DEFAULT_SPAN;
   parms->sampled = 0;
   parms->alignGridToMin = 0;
   parms->snap = 0;

   parms->vfirst = 1;
   parms->delim  = " \t\n\r,:"; /* make sure help screen matches this */
   parms->vl     = NULL;

   parms->writekin = 0;
   parms->writegrp = 0;
   parms->writesub = 0;
   parms->dominant = 0;
   parms->perspective = 1;
   parms->lens     = 0;
   parms->axis     = 1;
   parms->grpname  = "cont";

   parms->dumplattice = 0;

   parms->xmin = parms->ymin = parms->zmin = 0.0;
   parms->xmax = parms->ymax = parms->zmax = 0.0;

   parms->clipout = 0;
   parms->xout = parms->yout = parms->zout = 0.0;
   parms->rout    = FLT_MAX;
   parms->rout_sq = FLT_MAX;
   
   /* zero out stuff to be calculated from data */

   parms->ifx = parms->ify = parms->ifz = 0;
   parms->nfx = parms->nfy = parms->nfz = 0;
   parms->nx  = parms->ny  = parms->nz  = 0;
   parms->x0  = parms->y0  = parms->z0  = 0.0;

   parms->lat  = NULL; /* lattice goes here */
   parms->filt = NULL; /*  filter goes here */

   nfiles = 0;
   for (i = 1; i < argc; i++) {
      p = argv[i];
      if (p[0] == '-') {
	 if (p[1] == '\0') {
	    if (nfiles == 0) {
	       nfiles++;
	       parms->infname = "-";
	       inf = stdin;
	    }
	    else {
	       fprintf(stderr, "ERROR: stray '-' in command line\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "HELP", 1)){
	    helpInfo(1);
	 }
	 else if(compArgStr(p+1, "INput", 2)){
	    if (nfiles == 0) {
	       nfiles++;
	       if (++i < argc) {
		  parms->infname = argv[i];
		  inf = fopen(parms->infname, "r");
		  if (inf == NULL) {
		     fprintf(stderr, "ERROR: could not open file, \"%s\"\n",
		                     parms->infname);
		     return NULL;
		  }
	       }
	       else {
		  fprintf(stderr, "ERROR: -INput name missing\n");
		  helpInfo(0);
	       }
	    }
	    else {
	       fprintf(stderr, "ERROR: extra -INput flag in command line\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "Verbose", 1)){
	    parms->verbose = 1;
	 }
	 else if(compArgStr(p+1, "Quiet", 1)){
	    parms->verbose = 0;
	 }
	 else if(compArgStr(p+1, "FIRST", 5)){
	    parms->vfirst = 1;
	 }
	 else if(compArgStr(p+1, "LAST", 4)){
	    parms->vfirst = 0;
	 }
	 else if(compArgStr(p+1, "XPLORascii", 5)){
            mapinput = XPLOR_MAPTYPE;
	 }
	 else if(compArgStr(p+1, "SAMPled", 4)){
	    parms->sampled = 1;
	 }
	 else if(compArgStr(p+1, "ALIGNminxyz", 5)){
	    parms->sampled = 1;
            parms->alignGridToMin = 1;
            parms->wrap = 0;
	    parms->snap = 1;
	 }
	 else if(compArgStr(p+1, "SNAPtogrid", 4)){
	    parms->snap = 1;
	 }
	 else if(compArgStr(p+1, "KINemage", 3)){
	    parms->writekin = 1;
	 }
	 else if(compArgStr(p+1, "GROUP", 5)){
	    parms->writegrp = 1;
	    parms->writesub = 0;
	 }
	 else if(compArgStr(p+1, "SUBgroup", 3)){
	    parms->writegrp = 0;
	    parms->writesub = 1;
	 }
	 else if(compArgStr(p+1, "NAMe", 3)){
	    if (++i < argc) {
	       parms->grpname = argv[i];
	    }
	    else {
	       fprintf(stderr, "ERROR: -NAMe string missing\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "DOMinant", 3)){
	    parms->dominant = 1;
	 }
	 else if(compArgStr(p+1, "NOPERSPective", 7)){
            parms->perspective = 0;
	 }
	 else if(compArgStr(p+1, "LENS", 4)){
	    parms->lens = 1;
	 }
	 else if(compArgStr(p+1, "NOAXIS", 6)){
	    parms->axis = 0;
	 }
	 else if(compArgStr(p+1, "T", 1)){
	    parms->delim = p+2;
	 }
	 else if(n=compArgStr(p+1, "Grain", 1)){
            rval = atof(p+n+1);
	    if (rval < 0.0) { rval = -rval; }
	    parms->dx  = parms->dy  = parms->dz = rval;
	 }
	 else if(compArgStr(p+1, "GXYZ", 4)){
	    if (i+3 < argc) {
	       parms->dx = atof(argv[i+1]);
	       parms->dy = atof(argv[i+2]);
	       parms->dz = atof(argv[i+3]);
	       if (parms->dx < 0.0) { parms->dx = -parms->dx; }
	       if (parms->dy < 0.0) { parms->dy = -parms->dy; }
	       if (parms->dz < 0.0) { parms->dz = -parms->dz; }
	       i += 3;
	    }
	    else {
	       fprintf(stderr, "ERROR: -GXYZ flag requires 3 numbers\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "NOSMooth", 4)){
	    parms->usesmoothing = 0;
	 }
	 else if(n=compArgStr(p+1, "Stddev", 1)){
            rval = atof(p+n+1);
	    if (rval < 0.0) { rval = -rval; }
	    parms->sx  = parms->sy  = parms->sz = rval;
	 }
	 else if(compArgStr(p+1, "SXYZ", 4)){
	    if (i+3 < argc) {
	       parms->sx = atof(argv[i+1]);
	       parms->sy = atof(argv[i+2]);
	       parms->sz = atof(argv[i+3]);
	       if (parms->sx < 0.0) { parms->sx = -parms->sx; }
	       if (parms->sy < 0.0) { parms->sy = -parms->sy; }
	       if (parms->sz < 0.0) { parms->sz = -parms->sz; }
	       i += 3;
	    }
	    else {
	       fprintf(stderr, "ERROR: -RXYZ  flag requires 3 numbers\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "Levels",  1) ||
	         compArgStr(p+1, "SLevels", 2)){
	    if(compArgStr(p+1, "SLevels",     2)) {
	       parms->sigmalevels = 1;
	    }
	    while (i+2 < argc) {
	       n = (parms->nlevels)++;
	       if (n > MAXLEVELS) {
		  fprintf(stderr, "ERROR: levels exceed limit(%d)\n", MAXLEVELS);
		  return NULL;
	       }
	       parms->inlvls[n] = parms->levels[n] = atof(argv[i+1]);
	       parms->colors[n] = argv[i+2];
	       i += 2;

	       if ((i+2 >= argc)
	        || ((argv[i+1][0] != '-') &&
		    (argv[i+1][0] != '+') &&
		    (argv[i+1][0] != '.') && !isdigit(argv[i+1][0]))
	        || ((argv[i+1][0] == '-') && !isdigit(argv[i+1][1]))) {
		  /* If not enough room or not the way a number begins */
		  /* or starts with a minus but has no digits... */
 		  /* we are done with list of level info. */
		  break;
	       }
	    }
	 }
	 else if(compArgStr(p+1, "Multilevel",  1) ||
	         compArgStr(p+1, "SMultilevel", 2)){
	    double lcv = 0.0, lto = 0.0, lby = 0.0;

	    if(compArgStr(p+1, "SMultilevel", 2)) {
	       parms->sigmalevels = 1;
	    }
	    while (i+4 < argc) {
	       lcv = atof(argv[i+1]);
	       lto = atof(argv[i+2]);
	       lby = atof(argv[i+3]);

	       while (lcv <= (lto + CHOP_EPSILON)) {
		  n = (parms->nlevels)++;
		  if (n > MAXLEVELS) {
		     fprintf(stderr, "ERROR: levels exceed limit(%d)\n", MAXLEVELS);
		     return NULL;
		  }
		  
		  if (lcv > -CHOP_EPSILON && lcv < CHOP_EPSILON) {
		     lcv = 0.0;   /* chop levels very close to zero */
		  }
		  parms->inlvls[n] = parms->levels[n] = lcv;
		  parms->colors[n] = argv[i+4];
		  lcv += lby;
	       }

	       i += 4;

	       if ((i+4 >= argc)
	        || ((argv[i+1][0] != '-') &&
		    (argv[i+1][0] != '+') &&
		    (argv[i+1][0] != '.') && !isdigit(argv[i+1][0]))
	        || ((argv[i+1][0] == '-') && !isdigit(argv[i+1][1]))) {
		  /* If not enough room or not the way a number begins */
		  /* or starts with a minus but has no digits... */
 		  /* we are done with list of level info. */
		  break;
	       }
	    }
	 }
	 else if(n=compArgStr(p+1, "SPan", 2)){
            rval = atof(p+n+1);
	    parms->filtspan = rval;
	 }
	 else if(compArgStr(p+1, "DUMPlattice", 4)){
	    parms->dumplattice = 1;
	 }
	 else if(compArgStr(p+1, "WRAP", 4)){
	    if (i+6 < argc) {
	       parms->wrap = 1;
	       parms->xmin = atof(argv[i+1]);
	       parms->xmax = atof(argv[i+2]);
	       parms->ymin = atof(argv[i+3]);
	       parms->ymax = atof(argv[i+4]);
	       parms->zmin = atof(argv[i+5]);
	       parms->zmax = atof(argv[i+6]);
	       /* make min no greater than max */
	       if (parms->xmin > parms->xmax) {
		  rval = parms->xmin;
		  parms->xmin = parms->xmax;
		  parms->xmax = rval;
	       }
	       if (parms->ymin > parms->ymax) {
		  rval = parms->ymin;
		  parms->ymin = parms->ymax;
		  parms->ymax = rval;
	       }
	       if (parms->xmin > parms->zmax) {
		  rval = parms->zmin;
		  parms->zmin = parms->zmax;
		  parms->zmax = rval;
	       }
	       i += 6;
	    }
	    else {
	       fprintf(stderr, "ERROR: -WRAP flag requires 6 numbers (xmin,xmax,ymin,ymax,zmin,zmax)\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "CLIP", 4)){
	    if (i+4 < argc) {
	       parms->xout = atof(argv[i+1]);
	       parms->yout = atof(argv[i+2]);
	       parms->zout = atof(argv[i+3]);
	       parms->rout = atof(argv[i+4]);

	       parms->clipout = 1;
	       parms->rout_sq   = parms->rout*parms->rout;
	       i += 4;
	    }
	    else {
	       fprintf(stderr, "ERROR: -CLIP flag requires 4 numbers (x,y,z,halfwidth)\n");
	       helpInfo(0);
	    }
	 }
         else {
	    fprintf(stderr,"ERROR: unrecognized flag, \"%s\"\n", p);
	    helpInfo(0);
	 }
      }
      else {
	 if (nfiles == 0) {
	    nfiles++;
	    parms->infname = p;
	    inf = fopen(parms->infname, "r");
	    if (inf == NULL) {
	       fprintf(stderr, "ERROR: could not open file, \"%s\"\n",
	                       parms->infname);
	       return NULL;
	    }
	 }
	 else {
	    fprintf(stderr,"ERROR: extra parameter, \"%s\"\n", p);
	    helpInfo(0);
	 }
      }
   }
   if (nfiles < 1) {
      fprintf(stderr,"ERROR: no input filename or -\n");
      helpInfo(0);
   }

   if ((parms->dx < CHOP_EPSILON)||(parms->dy < CHOP_EPSILON)||(parms->dz < CHOP_EPSILON)) {
      fprintf(stderr, "ERROR: grid spacing must be greater than zero: %g, %g, %g\n",
                       parms->dx, parms->dy, parms->dz);
   }
   else if (mapinput == XPLOR_MAPTYPE) {
      parms->vl = load_xplor_map(inf);
      if (parms->vl && parms->vl->mapinfo) {
         parms->dx = parms->vl->mapinfo->gapsz[0];
         parms->dy = parms->vl->mapinfo->gapsz[1];
         parms->dz = parms->vl->mapinfo->gapsz[2];
      }

   }
   else {
      parms->vl = load_coords(inf, 3, parms->vfirst, parms->delim);
   }

   if (inf != stdin) { fclose(inf); }

   if (! parms->usesmoothing) {
      parms->sx  = parms->sy  = parms->sz = 0.0;
   }
   
   if (parms->vl) {
      if (parms->vl->ncoords > 0) { /* work out the calculated parms */
	 double dval = 0.0;

	 /* size the spot filter (integer radii) */
	 parms->ifx = ceil(parms->filtspan * parms->sx/parms->dx);
	 parms->ify = ceil(parms->filtspan * parms->sy/parms->dy);
	 parms->ifz = ceil(parms->filtspan * parms->sz/parms->dz);

	 parms->nfx = 2*parms->ifx + 1; /* total filter dim */
	 parms->nfy = 2*parms->ify + 1;
	 parms->nfz = 2*parms->ifz + 1;

/* these deal with where the range is whole number of grain val */
/* either go over */
#define WIGGLE_ROOM ( 0.0001)
/* or round with a bias towards upward (needed in a world with roundoff) */
#define SHRINKAGE   (-0.25 )

	 if (parms->wrap) {
	    /* period of the lattice + edge */
	    parms->nx = ceil((parms->xmax - parms->xmin)/parms->dx) + 1;
	    parms->ny = ceil((parms->ymax - parms->ymin)/parms->dy) + 1;
	    parms->nz = ceil((parms->zmax - parms->zmin)/parms->dz) + 1;

	    /* fine tune the lattice spacing (for integer period) */
	    parms->dx = (parms->xmax - parms->xmin)/(parms->nx - 1);
	    parms->dy = (parms->ymax - parms->ymin)/(parms->ny - 1);
	    parms->dz = (parms->zmax - parms->zmin)/(parms->nz - 1);

	    /* lattice origin */
	    parms->x0 = parms->xmin;
	    parms->y0 = parms->ymin;
	    parms->z0 = parms->zmin;
	 }
	 else if (parms->alignGridToMin) {
	    /* period of the lattice + filter + border */
	    n = ceil(SHRINKAGE + (parms->vl->xmax - parms->vl->xmin)/parms->dx) + 1;
	    parms->nx = n + (parms->nfx - 1) + 2;
	    n = ceil(SHRINKAGE + (parms->vl->ymax - parms->vl->ymin)/parms->dy) + 1;
	    parms->ny = n + (parms->nfy - 1) + 2;
	    n = ceil(SHRINKAGE + (parms->vl->zmax - parms->vl->zmin)/parms->dz) + 1;
	    parms->nz = n + (parms->nfz - 1) + 2;

	    /* lattice origin */
	    parms->xmin = parms->vl->xmin - (0.5*(parms->nfx - 1) + 1.0)*parms->dx;
	    parms->ymin = parms->vl->ymin - (0.5*(parms->nfy - 1) + 1.0)*parms->dy;
	    parms->zmin = parms->vl->zmin - (0.5*(parms->nfz - 1) + 1.0)*parms->dz;

	    parms->xmax = parms->vl->xmax + (0.5*(parms->nfx - 1) + 1.0)*parms->dx;
	    parms->ymax = parms->vl->ymax + (0.5*(parms->nfy - 1) + 1.0)*parms->dy;
	    parms->zmax = parms->vl->zmax + (0.5*(parms->nfz - 1) + 1.0)*parms->dz;

	    parms->x0 = parms->xmin;
	    parms->y0 = parms->ymin;
	    parms->z0 = parms->zmin;

	    /* fine tune the lattice spacing (for integer period) */
	    parms->dx = (parms->xmax - parms->xmin)/(parms->nx - 1);
	    parms->dy = (parms->ymax - parms->ymin)/(parms->ny - 1);
	    parms->dz = (parms->zmax - parms->zmin)/(parms->nz - 1);
	 }
	 else { /* no wrap - center box on data */
	    /* size the lattice + filter */
	    n = ceil(WIGGLE_ROOM +
	             (parms->vl->xmax - parms->vl->xmin)/parms->dx) + 1;
	    parms->nx = n + (parms->nfx - 1);
	    n = ceil(WIGGLE_ROOM +
	             (parms->vl->ymax - parms->vl->ymin)/parms->dy) + 1;
	    parms->ny = n + (parms->nfy - 1);
	    n = ceil(WIGGLE_ROOM +
	             (parms->vl->zmax - parms->vl->zmin)/parms->dz) + 1;
	    parms->nz = n + (parms->nfz - 1);

	    /* position the origin */
	    dval = parms->vl->xmin + 0.5*(parms->vl->xmax - parms->vl->xmin);
	    parms->x0 = dval - 0.5*(parms->nx - 1)*parms->dx;

	    dval = parms->vl->ymin + 0.5*(parms->vl->ymax - parms->vl->ymin);
	    parms->y0 = dval - 0.5*(parms->ny - 1)*parms->dy;

	    dval = parms->vl->zmin + 0.5*(parms->vl->zmax - parms->vl->zmin);
	    parms->z0 = dval - 0.5*(parms->nz - 1)*parms->dz;

	    parms->xmin = parms->x0; /* bounds of lattice */
	    parms->ymin = parms->y0;
	    parms->zmin = parms->z0;
	    parms->xmax = parms->x0 + parms->dx*(parms->nx - 1);
	    parms->ymax = parms->y0 + parms->dy*(parms->ny - 1);
	    parms->zmax = parms->z0 + parms->dz*(parms->nz - 1);
	 }
      }
      else {
	 fprintf(stderr, "WARNING: no data read from file, \"%s\"\n",
	    parms->infname);
      }
   }
   else {
      fprintf(stderr, "ERROR: problems reading file, \"%s\"\n",
	 parms->infname);
   }

   if (parms->verbose) { describeInput(parms); }
   
   return parms->vl;
}

real_val*** build_lattice(progparms *parms) {
   int ix=0, iy=0, iz=0, i = 0;
   coord_val rval = 0.0;

   if (parms == NULL) { return NULL; }

   parms->lat  = matrix3D(0, parms->nx -1,
                          0, parms->ny -1,
                          0, parms->nz -1);
   if (parms->lat) {
      for (ix = 0; ix < parms->nx; ix++) {
	 for (iy = 0; iy < parms->ny; iy++) {
	    for (iz = 0; iz < parms->nz; iz++) {
	       parms->lat[ix][iy][iz] = 0.0; /* initialize to zero */
	    }
	 }
      }
   }

#ifdef DEBUG_LATTICE_GRID
   rval = parms->x0;
   fprintf(stderr, "X grid points: ", rval);
   for (i = 0; i < parms->nx; i++) {
      if (i%8 == 0) { fprintf(stderr, "\n"); }
      fprintf(stderr, "%g ", rval);
      rval += parms->dx;
   }
   fprintf(stderr, "\n");

   rval = parms->y0;
   fprintf(stderr, "Y grid points: ", rval);
   for (i = 0; i < parms->ny; i++) {
      if (i%8 == 0) { fprintf(stderr, "\n"); }
      fprintf(stderr, "%g ", rval);
      rval += parms->dy;
   }
   fprintf(stderr, "\n");

   rval = parms->z0;
   fprintf(stderr, "Z grid points: ", rval);
   for (i = 0; i < parms->nz; i++) {
      if (i%8 == 0) { fprintf(stderr, "\n"); }
      fprintf(stderr, "%g ", rval);
      rval += parms->dz;
   }
   fprintf(stderr, "\n");
#endif

   return parms->lat;
}

void destroy_lattice(progparms *parms) {
   free_matrix3D(parms->lat, 0, 0, 0);
}

real_val*** build_filter(progparms *parms) {
   int ix=0, iy=0, iz=0;
   real_val sumval=0.0, scalef=0.0;
   real_val qx=0.0, qy=0.0, qz=0.0;
   
   if (parms == NULL) { return NULL; }

   parms->filt = matrix3D(-parms->ifx, parms->ifx,
                          -parms->ify, parms->ify,
                          -parms->ifz, parms->ifz);

   if (parms->filt == NULL) { return NULL; }

   if (parms->sx > -CHOP_EPSILON && parms->sx < CHOP_EPSILON) {
      qx = 1.0;
   }
   else {
      qx = (parms->dx * parms->dx) / (2.0 * parms->sx * parms->sx);
   }
   if (parms->sy > -CHOP_EPSILON && parms->sy < CHOP_EPSILON) {
      qy = 1.0;
   }
   else {
      qy = (parms->dy * parms->dy) / (2.0 * parms->sy * parms->sy);
   }
   if (parms->sz > -CHOP_EPSILON && parms->sz < CHOP_EPSILON) {
      qz = 1.0;
   }
   else {
      qz = (parms->dz * parms->dz) / (2.0 * parms->sz * parms->sz);
   }

   sumval = 0.0;
   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 for (iz = -parms->ifz; iz <= parms->ifz; iz++) {
	    parms->filt[ix][iy][iz]
		  = exp(-(  ix*ix*qx + iy*iy*qy + iz*iz*qz ));
	    sumval += parms->filt[ix][iy][iz];
	 }
      }
   }

   scalef = 1.0/sumval;
   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 for (iz = -parms->ifz; iz <= parms->ifz; iz++) {
	    parms->filt[ix][iy][iz] *= scalef;
	 }
      }
   }

   if (parms->verbose) { /* display a section of the filter */
      fprintf(stderr, "100 X Y>");
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 fprintf(stderr, "%7d ", iy);
      }
      fprintf(stderr, "\n");
      for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
	 fprintf(stderr, "%7d ", ix);
	 for (iy = -parms->ify; iy <= parms->ify; iy++) {
	    fprintf(stderr, "%7.3f ", parms->filt[ix][iy][0]*100.0);
	 }
	 fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
   }

   return parms->filt;
}

void destroy_filter(progparms *parms) {
   if (parms && parms->filt) {
      free_matrix3D(parms->filt, -parms->ifx, -parms->ify, -parms->ifz);
   }
}

int insert_data_into_lattice(progparms *parms) {
   val3Dnode* v3D = NULL;
   int n = 0;
   float pctfact = 0.0, pervol = 0.0;

   pctfact = 100.0/parms->vl->ncoords;

   /* normalization factor */
   if (parms->sampled) { pervol = 1.0; }
   else {
      pervol = 1.0/(parms->dx * parms->dy * parms->dz);
   }

   if (parms && parms->vl && parms->vl->head) {
      if (parms->verbose) {
	 fprintf(stderr, "inserting data into lattice\n");
      }
      for(v3D = parms->vl->head; v3D; v3D = v3D->next) {
	 n++;
	 splat(parms, pervol * v3D->v, v3D->x, v3D->y, v3D->z);

	 if (parms->verbose) {
	    fprintf(stderr, "   %3.0f%% \r", n * pctfact);
	 }
      }
      if (parms->verbose) {
	 fprintf(stderr, "complete.     \n");
      }
   }
   return n;
}

void splat(progparms *parms, real_val wt, coord_val x, coord_val y, coord_val z) {
   int lx = 0, ly = 0, lz = 0;
   coord_val deltax = 0.0, deltay = 0.0, deltaz = 0.0;
   coord_val gridx  = 0.0, gridy  = 0.0, gridz  = 0.0;
   coord_val fractx = 0.0, fracty = 0.0, fractz = 0.0;

   deltax = (x - parms->x0)/parms->dx;
   deltay = (y - parms->y0)/parms->dy;
   deltaz = (z - parms->z0)/parms->dz;

   gridx = floor(deltax);
   gridy = floor(deltay);
   gridz = floor(deltaz);

   lx = (int)gridx;
   ly = (int)gridy;
   lz = (int)gridz;

   fractx = deltax - gridx;
   fracty = deltay - gridy;
   fractz = deltaz - gridz;

   if (parms->snap) {
      if (fractx < 0.01) { fractx = 0.0; }
      if (fracty < 0.01) { fracty = 0.0; }
      if (fractz < 0.01) { fractz = 0.0; }
      if (fractx > 0.99) { fractx = 1.0; }
      if (fracty > 0.99) { fracty = 1.0; }
      if (fractz > 0.99) { fractz = 1.0; }
   }

   if (parms->usesmoothing) {
      if (parms->wrap) {
         /* skip duplicate points past wrap edge */
         if ((x < parms->xmax) && (y < parms->ymax) && (z < parms->zmax)) {
	    placemaskwrap(parms, wt*(1.0-fractx)*(1.0-fracty)*(1.0-fractz), lx,   ly,   lz);
	    placemaskwrap(parms, wt*(    fractx)*(1.0-fracty)*(1.0-fractz), lx+1, ly,   lz);
	    placemaskwrap(parms, wt*(1.0-fractx)*(    fracty)*(1.0-fractz), lx,   ly+1, lz);
	    placemaskwrap(parms, wt*(    fractx)*(    fracty)*(1.0-fractz), lx+1, ly+1, lz);
	    placemaskwrap(parms, wt*(1.0-fractx)*(1.0-fracty)*(    fractz), lx,   ly,   lz+1);
	    placemaskwrap(parms, wt*(    fractx)*(1.0-fracty)*(    fractz), lx+1, ly,   lz+1);
	    placemaskwrap(parms, wt*(1.0-fractx)*(    fracty)*(    fractz), lx,   ly+1, lz+1);
	    placemaskwrap(parms, wt*(    fractx)*(    fracty)*(    fractz), lx+1, ly+1, lz+1);
         }
      }
      else {
         placemask(parms, wt*(1.0-fractx)*(1.0-fracty)*(1.0-fractz), lx,   ly,   lz);
         placemask(parms, wt*(    fractx)*(1.0-fracty)*(1.0-fractz), lx+1, ly,   lz);
         placemask(parms, wt*(1.0-fractx)*(    fracty)*(1.0-fractz), lx,   ly+1, lz);
         placemask(parms, wt*(    fractx)*(    fracty)*(1.0-fractz), lx+1, ly+1, lz);
         placemask(parms, wt*(1.0-fractx)*(1.0-fracty)*(    fractz), lx,   ly,   lz+1);
         placemask(parms, wt*(    fractx)*(1.0-fracty)*(    fractz), lx+1, ly,   lz+1);
         placemask(parms, wt*(1.0-fractx)*(    fracty)*(    fractz), lx,   ly+1, lz+1);
         placemask(parms, wt*(    fractx)*(    fracty)*(    fractz), lx+1, ly+1, lz+1);
      }
   }
   else { /* no filter required */

      placesample(parms, wt*(1.0-fractx)*(1.0-fracty)*(1.0-fractz), lx,   ly,   lz);
      placesample(parms, wt*(    fractx)*(1.0-fracty)*(1.0-fractz), lx+1, ly,   lz);
      placesample(parms, wt*(1.0-fractx)*(    fracty)*(1.0-fractz), lx,   ly+1, lz);
      placesample(parms, wt*(    fractx)*(    fracty)*(1.0-fractz), lx+1, ly+1, lz);
      placesample(parms, wt*(1.0-fractx)*(1.0-fracty)*(    fractz), lx,   ly,   lz+1);
      placesample(parms, wt*(    fractx)*(1.0-fracty)*(    fractz), lx+1, ly,   lz+1);
      placesample(parms, wt*(1.0-fractx)*(    fracty)*(    fractz), lx,   ly+1, lz+1);
      placesample(parms, wt*(    fractx)*(    fracty)*(    fractz), lx+1, ly+1, lz+1);
   }
}

void placemask(progparms *parms, real_val wt, int lx, int ly, int lz) {
   int ix = 0, iy = 0, iz = 0;

   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 for (iz = -parms->ifx; iz <= parms->ifz; iz++) {
	    parms->lat[lx+ix][ly+iy][lz+iz] += wt * parms->filt[ix][iy][iz];
	 }
      }
   }
}

void placemaskwrap(progparms *parms, real_val wt, int lx, int ly, int lz) {
   int ix = 0, iy = 0, iz = 0;
   int mx = 0, my = 0, mz = 0;

   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      mx = lx + ix;
      while (mx < 0) { mx += (parms->nx - 1); }
      while (mx >= (parms->nx - 1)) { mx -= (parms->nx - 1); }
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 my = ly + iy;
	 while (my < 0) { my += (parms->ny - 1); }
	 while (my >= (parms->ny - 1)) { my -= (parms->ny - 1); }
	 for (iz = -parms->ifx; iz <= parms->ifz; iz++) {
	    mz = lz + iz;
	    while (mz < 0) { mz += (parms->nz - 1); }
	    while (mz >= (parms->nz - 1)) { mz -= (parms->nz - 1); }
	    parms->lat[mx][my][mz] += wt * parms->filt[ix][iy][iz];
	 }
      }
   }
}

void placesample(progparms *parms, real_val wt, int lx, int ly, int lz) {
   parms->lat[lx][ly][lz] += wt;
}

void contour_lattice(progparms *parms) {
   cont_3d_info *cp=NULL;

   if (parms->nlevels < 1) {
      calc_default_contour_levels(parms);
   }
   else if (parms->sigmalevels) { /* convert levels from sigma to absolute */
      input_sigma_to_contour_levels(parms);
   }

   cp = build_cont_info(parms->levels, parms->nlevels,
                        parms->colors, 3);

   if (cp) {
      if (parms->wrap) {
	 if (parms->verbose) { fprintf(stderr, "duplicating edges\n"); }
	 dup_periodic_edges(parms);
      }

      if (parms->verbose) { fprintf(stderr, "X-Y contours\n"); }

      xy_contour(parms, cp);

      if (parms->verbose) { fprintf(stderr, "X-Z contours\n"); }

      xz_contour(parms, cp);

      if (parms->verbose) { fprintf(stderr, "Y-Z contours\n"); }

      yz_contour(parms, cp);

      if (parms->verbose) { fprintf(stderr, "writing to kin file.\n"); }

      simple_kin_contour_output(stdout, cp, 0, 1);

      destroy_cont_info(cp);
   }
}

void dup_periodic_edges(progparms *parms) {
   int ix = 0, iy = 0, iz = 0;
   int mx = 0, my = 0, mz = 0;

   mx = parms->nx - 1;
   my = parms->ny - 1;
   mz = parms->nz - 1;

   for (iz = 0; iz < mz; iz++) {
      /* copy x-y edges on each slice */
      for (ix = 0; ix < mx; ix++) {
	 parms->lat[ix][my][iz] = parms->lat[ix][0][iz];
      }
      for (iy = 0; iy < my; iy++) {
	 parms->lat[mx][iy][iz] = parms->lat[0][iy][iz];
      }
      parms->lat[mx][my][iz] = parms->lat[0][0][iz];
   }

   /* the max z plane is the same as the 0th one */
   for (ix = 0; ix <= mx; ix++) {
      for (iy = 0; iy <= my; iy++) {
	 parms->lat[ix][iy][mz] = parms->lat[ix][iy][0];
      }
   }
}

void build_clip_box(progparms *parms, clip_box *box) {
   int      lb      = 0,  le      = 0,   t = 0;
   coord_val b      = 0.0, e      = 0.0;
   coord_val deltab = 0.0, deltae = 0.0;
   coord_val gridb  = 0.0, gride  = 0.0;

   box->xout    = parms->xout;
   box->yout    = parms->yout;
   box->zout    = parms->zout;
   box->rout    = parms->rout;
   box->rout_sq = parms->rout_sq;

   if (! parms->clipout) {
     box->xb = parms->x0;
     box->yb = parms->y0;
     box->zb = parms->z0;

     box->xe = parms->x0 + (parms->nx - 1) * parms->dx;
     box->ye = parms->y0 + (parms->ny - 1) * parms->dy;
     box->ze = parms->z0 + (parms->nz - 1) * parms->dz;

     box->xo = box->xrgn = 0; /* all of it */
     box->yo = box->yrgn = 0;
     box->zo = box->zrgn = 0;
     
     box->xf = parms->nx - 1;
     box->yf = parms->ny - 1;
     box->zf = parms->nz - 1;
     return; /* quit while you're ahead */
   }
/* -------------- x ----------------------*/
   b = parms->xout - parms->rout;
   e = parms->xout + parms->rout;

   deltab = (b - parms->x0)/parms->dx;
   deltae = (e - parms->x0)/parms->dx;

   gridb = floor(deltab);
   gride =  ceil(deltae);

   lb = (int)gridb;
   le = (int)gride;
   
   if (lb < 0)               { lb = 0; }
   if (le > (parms->nx - 1)) { le = parms->nx - 1; }
   if (lb > le) { t = le; le = lb; lb = t; }

   box->xb = parms->x0 + lb * parms->dx;
   box->xe = parms->x0 + le * parms->dx;

   box->xo   = lb;
   box->xrgn = le - lb;
   box->xf   = le;

/* -------------- y ----------------------*/
   b = parms->yout - parms->rout;
   e = parms->yout + parms->rout;

   deltab = (b - parms->y0)/parms->dy;
   deltae = (e - parms->y0)/parms->dy;

   gridb = floor(deltab);
   gride =  ceil(deltae);

   lb = (int)gridb;
   le = (int)gride;
   
   if (lb < 0)               { lb = 0; }
   if (le > (parms->ny - 1)) { le = parms->ny - 1; }
   if (lb > le) { t = le; le = lb; lb = t; }

   box->yb = parms->y0 + lb * parms->dy;
   box->ye = parms->y0 + le * parms->dy;

   box->yo   = lb;
   box->yrgn = le - lb;
   box->yf   = le;

/* -------------- z ----------------------*/
   b = parms->zout - parms->rout;
   e = parms->zout + parms->rout;

   deltab = (b - parms->z0)/parms->dz;
   deltae = (e - parms->z0)/parms->dz;

   gridb = floor(deltab);
   gride =  ceil(deltae);

   lb = (int)gridb;
   le = (int)gride;
   
   if (lb < 0)               { lb = 0; }
   if (le > (parms->nz - 1)) { le = parms->nz - 1; }
   if (lb > le) { t = le; le = lb; lb = t; }

   box->zb = parms->z0 + lb * parms->dz;
   box->ze = parms->z0 + le * parms->dz;

   box->zo   = lb;
   box->zrgn = le - lb;
   box->zf   = le;
}

void xy_contour(progparms *parms, cont_3d_info *cp) {
   int i = 0, rc = 0;
   int r = 0, c = 0, co = 0, ro = 0, crgn = 0, rrgn = 0;
   coord_val xspan = 0.0, yspan = 0.0;
   ctree_data data;
   ctree_val **rawdata = NULL;
   ctree_val *work1 = NULL, *work2 = NULL;
   coord_val slice_pos = 0.0;
   slice_info *sp=NULL;
   int dimmap[3];
   clip_box clipbox;
   
   dimmap[0] = X_DIM_NUM;
   dimmap[1] = Y_DIM_NUM;
   dimmap[2] = Z_DIM_NUM;

   build_clip_box(parms, &clipbox);

   rawdata = matrix(0, parms->ny - 1, 0, parms->nx - 1);
   work1   = vector(0, parms->nx - 1);
   work2   = vector(0, parms->nx - 1);

   if (rawdata == NULL || work1 == NULL || work2 == NULL) {
      fprintf(stderr, "ERROR: out of memory xy_contour(%d x %d)\n",
	 parms->nx, parms->ny);
      return; /* should we just exit? */
   }

   yspan = (parms->ny - 1) * parms->dy;
   xspan = (parms->nx - 1) * parms->dx;
   r = parms->ny;
   c = parms->nx;
   co = clipbox.xo; ro = clipbox.yo;
   crgn = clipbox.xrgn; rrgn = clipbox.yrgn;

#ifdef DEBUG_CLIP
fprintf(stderr, "DEBUG: xy(%g, %g) cr(%d, %d)\n       cro(%d, %d) crrgn(%d, %d)\n",
               parms->x0, parms->y0, c, r, co, ro, crgn, rrgn);
#endif
   rc = ctree_compose_rawdata(&data, &(rawdata[0][0]),
	    parms->x0, parms->y0, xspan, yspan,
	    c, r, co, ro, crgn, rrgn, &(work1[0]), &(work2[0]), NULL);

   if (rc) {
      slice_pos = clipbox.zb;

      sp = build_slice_info(parms->levels, parms->nlevels, dimmap,
	 parms->dx, parms->dy, 1, slice_pos, 1);

      if (sp) {
	 for (i = clipbox.zo; i <= clipbox.zf; i++) { /* each slice */

	    for (c = 0; c < parms->nx; c++) { /* load array */
	       for (r = 0; r < parms->ny; r++) {
		  rawdata[r][c] = parms->lat[c][r][i];
	       }
	    }

	    ctree_plot_region(&data, parms->levels, parms->nlevels,
		  basic_data_func, connecting_plot_func, sp);
	    slice_all_done(sp);

	    remap_slices(cp, sp);

	    slice_pos += parms->dz;
	    update_slice_info(sp, dimmap, parms->dx, parms->dy,
	       1, slice_pos);
	 }
      }

      destroy_slice_info(sp);
   }
   else {
      fprintf(stderr, "ERROR: bad rawdata parameters in xy_contour()\n");
   }
   free_vector(work2, 0);
   free_vector(work1, 0);
   free_matrix(rawdata, 0, 0);
}

void xz_contour(progparms *parms, cont_3d_info *cp) {
   int i = 0, rc = 0;
   int r = 0, c = 0, co = 0, ro = 0, crgn = 0, rrgn = 0;
   coord_val xspan = 0.0, zspan = 0.0;
   ctree_data data;
   ctree_val **rawdata = NULL;
   ctree_val *work1 = NULL, *work2 = NULL;
   coord_val slice_pos = 0.0;
   slice_info *sp=NULL;
   int dimmap[3];
   clip_box clipbox;

   dimmap[0] = X_DIM_NUM;
   dimmap[1] = Z_DIM_NUM;
   dimmap[2] = Y_DIM_NUM;

   build_clip_box(parms, &clipbox);

   rawdata = matrix(0, parms->nz - 1, 0, parms->nx - 1);
   work1   = vector(0, parms->nx - 1);
   work2   = vector(0, parms->nx - 1);

   if (rawdata == NULL || work1 == NULL || work2 == NULL) {
      fprintf(stderr, "ERROR: out of memory xz_contour(%d x %d)\n",
	 parms->nx, parms->nz);
      return; /* should we just exit? */
   }

   zspan = (parms->nz - 1) * parms->dz;
   xspan = (parms->nx - 1) * parms->dx;
   r = parms->nz;
   c = parms->nx;
   co = clipbox.xo; ro = clipbox.zo;
   crgn = clipbox.xrgn; rrgn = clipbox.zrgn;

#ifdef DEBUG_CLIP
fprintf(stderr, "DEBUG: xz(%g, %g) cr(%d, %d)\n       cro(%d, %d) crrgn(%d, %d)\n",
               parms->x0, parms->z0, c, r, co, ro, crgn, rrgn);
#endif
   rc = ctree_compose_rawdata(&data, &(rawdata[0][0]),
	    parms->x0, parms->z0, xspan, zspan,
	    c, r, co, ro, crgn, rrgn, &(work1[0]), &(work2[0]), NULL);

   if (rc) {
      slice_pos = clipbox.yb;

      sp = build_slice_info(parms->levels, parms->nlevels, dimmap,
	 parms->dx, parms->dz, 1, slice_pos, 1);

      if (sp) {
	 for (i = clipbox.yo; i <= clipbox.yf; i++) { /* each slice */

	    for (c = 0; c < parms->nx; c++) { /* load array */
	       for (r = 0; r < parms->nz; r++) {
		  rawdata[r][c] = parms->lat[c][i][r];
	       }
	    }

	    ctree_plot_region(&data, parms->levels, parms->nlevels,
		  basic_data_func, connecting_plot_func, sp);
	    slice_all_done(sp);

	    remap_slices(cp, sp);

	    slice_pos += parms->dy;
	    update_slice_info(sp, dimmap, parms->dx, parms->dz,
	       1, slice_pos);
	 }
      }

      destroy_slice_info(sp);
   }
   else {
      fprintf(stderr, "ERROR: bad rawdata parameters in xz_contour()\n");
   }
   free_vector(work2, 0);
   free_vector(work1, 0);
   free_matrix(rawdata, 0, 0);
}

void yz_contour(progparms *parms, cont_3d_info *cp) {
   int i = 0, rc = 0;
   int r = 0, c = 0, co = 0, ro = 0, crgn = 0, rrgn = 0;
   coord_val yspan = 0.0, zspan = 0.0;
   ctree_data data;
   ctree_val **rawdata = NULL;
   ctree_val *work1 = NULL, *work2 = NULL;
   coord_val slice_pos = 0.0;
   slice_info *sp=NULL;
   int dimmap[3];
   clip_box clipbox;

   dimmap[0] = Y_DIM_NUM;
   dimmap[1] = Z_DIM_NUM;
   dimmap[2] = X_DIM_NUM;

   build_clip_box(parms, &clipbox);

   rawdata = matrix(0, parms->nz - 1, 0, parms->ny - 1);
   work1   = vector(0, parms->ny - 1);
   work2   = vector(0, parms->ny - 1);

   if (rawdata == NULL || work1 == NULL || work2 == NULL) {
      fprintf(stderr, "ERROR: out of memory yz_contour(%d x %d)\n",
	 parms->ny, parms->nz);
      return; /* should we just exit? */
   }

   zspan = (parms->nz - 1) * parms->dz;
   yspan = (parms->ny - 1) * parms->dy;
   r = parms->nz;
   c = parms->ny;
   co = clipbox.yo; ro = clipbox.zo;
   crgn = clipbox.yrgn; rrgn = clipbox.zrgn;

#ifdef DEBUG_CLIP
fprintf(stderr, "DEBUG: yz(%g, %g) cr(%d, %d)\n       cro(%d, %d) crrgn(%d, %d)\n",
               parms->y0, parms->z0, c, r, co, ro, crgn, rrgn);
#endif
   rc = ctree_compose_rawdata(&data, &(rawdata[0][0]),
	    parms->y0, parms->z0, yspan, zspan,
	    c, r, co, ro, crgn, rrgn, &(work1[0]), &(work2[0]), NULL);

   if (rc) {
      slice_pos = clipbox.xb;

      sp = build_slice_info(parms->levels, parms->nlevels, dimmap,
	 parms->dy, parms->dz, 1, slice_pos, 1);

      if (sp) {
	 for (i = clipbox.xo; i <= clipbox.xf; i++) { /* each slice */

	    for (c = 0; c < parms->ny; c++) { /* load array */
	       for (r = 0; r < parms->nz; r++) {
		  rawdata[r][c] = parms->lat[i][c][r];
	       }
	    }

	    ctree_plot_region(&data, parms->levels, parms->nlevels,
		  basic_data_func, connecting_plot_func, sp);
	    slice_all_done(sp);

	    remap_slices(cp, sp);

	    slice_pos += parms->dx;
	    update_slice_info(sp, dimmap, parms->dy, parms->dz,
	       1, slice_pos);
	 }
      }

      destroy_slice_info(sp);
   }
   else {
      fprintf(stderr, "ERROR: bad rawdata parameters in yz_contour()\n");
   }
   free_vector(work2, 0);
   free_vector(work1, 0);
   free_matrix(rawdata, 0, 0);
}

void draw_axis_box(progparms *parms) {

   if (parms->vl && parms->vl->mapinfo) {

      fprintf(stdout, "@vectorlist {unit_cell} color= gray width= 1 off\n");
      kin_rombus_draw_cell(parms->vl,
                      0, 0, 0,
                      parms->vl->mapinfo->grid[0],
		      parms->vl->mapinfo->grid[1],
		      parms->vl->mapinfo->grid[2]);

      fprintf(stdout, "@vectorlist {data_range} color= yellow width= 1 off\n");
      kin_rombus_draw_cell(parms->vl,
                      parms->vl->mapinfo->origin[0],
		      parms->vl->mapinfo->origin[1],
		      parms->vl->mapinfo->origin[2],
		      parms->vl->mapinfo->extent[0],
		      parms->vl->mapinfo->extent[1],
		      parms->vl->mapinfo->extent[2]);
   }
   else {
      fprintf(stdout, "@vectorlist {axis} color= white width= 1 off\n");

      kin_draw_box(parms->xmin, parms->ymin, parms->zmin,
                   parms->xmax, parms->ymax, parms->zmax);
   }
}

void kin_rombus_draw_cell(val3Dlist* vl,
                          int a_lo, int b_lo, int c_lo,
                          int a_hi, int b_hi, int c_hi) {

   kin_rombus_draw_new_line(vl, a_lo, b_lo, c_lo, a_hi, b_lo, c_lo); /* c_lo square */
   kin_rombus_continue_line(vl, a_hi, b_lo, c_lo, a_hi, b_hi, c_lo);
   kin_rombus_continue_line(vl, a_hi, b_hi, c_lo, a_lo, b_hi, c_lo);
   kin_rombus_continue_line(vl, a_lo, b_hi, c_lo, a_lo, b_lo, c_lo);

   kin_rombus_continue_line(vl, a_lo, b_lo, c_lo, a_lo, b_lo, c_hi); /* connector */

   kin_rombus_continue_line(vl, a_lo, b_lo, c_hi, a_hi, b_lo, c_hi); /* c_hi square */
   kin_rombus_continue_line(vl, a_hi, b_lo, c_hi, a_hi, b_hi, c_hi);
   kin_rombus_continue_line(vl, a_hi, b_hi, c_hi, a_lo, b_hi, c_hi);
   kin_rombus_continue_line(vl, a_lo, b_hi, c_hi, a_lo, b_lo, c_hi);

   kin_rombus_draw_new_line(vl, a_hi, b_lo, c_lo, a_hi, b_lo, c_hi); /* missing pieces */
   kin_rombus_draw_new_line(vl, a_hi, b_hi, c_lo, a_hi, b_hi, c_hi);
   kin_rombus_draw_new_line(vl, a_lo, b_hi, c_lo, a_lo, b_hi, c_hi);
}

void kin_rombus_draw_new_line(val3Dlist* vl,
                          int a_beg, int b_beg, int c_beg,
                          int a_end, int b_end, int c_end) {
   double x1, y1, z1, x2, y2, z2;

   mapindex_to_xyz(vl, a_beg, b_beg, c_beg, &x1, &y1, &z1);
   mapindex_to_xyz(vl, a_end, b_end, c_end, &x2, &y2, &z2);

   kin_draw_new_line(x1, y1, z1, x2, y2, z2);
}

void kin_rombus_continue_line(val3Dlist* vl,
                          int a_beg, int b_beg, int c_beg,
                          int a_end, int b_end, int c_end) {
   double x1, y1, z1, x2, y2, z2;

   mapindex_to_xyz(vl, a_beg, b_beg, c_beg, &x1, &y1, &z1);
   mapindex_to_xyz(vl, a_end, b_end, c_end, &x2, &y2, &z2);

   kin_continue_line(x1, y1, z1, x2, y2, z2);
}

void kin_draw_box(coord_val xmin, coord_val ymin, coord_val zmin,
                  coord_val xmax, coord_val ymax, coord_val zmax) {

   kin_draw_new_line(xmin, ymin, zmin, xmax, ymin, zmin); /* zmin square */
   kin_continue_line(xmax, ymin, zmin, xmax, ymax, zmin);
   kin_continue_line(xmax, ymax, zmin, xmin, ymax, zmin);
   kin_continue_line(xmin, ymax, zmin, xmin, ymin, zmin);

   kin_continue_line(xmin, ymin, zmin, xmin, ymin, zmax); /* connector */

   kin_continue_line(xmin, ymin, zmax, xmax, ymin, zmax); /* zmax square */
   kin_continue_line(xmax, ymin, zmax, xmax, ymax, zmax);
   kin_continue_line(xmax, ymax, zmax, xmin, ymax, zmax);
   kin_continue_line(xmin, ymax, zmax, xmin, ymin, zmax);

   kin_draw_new_line(xmax, ymin, zmin, xmax, ymin, zmax); /* missing pieces */
   kin_draw_new_line(xmax, ymax, zmin, xmax, ymax, zmax);
   kin_draw_new_line(xmin, ymax, zmin, xmin, ymax, zmax);
}

void kin_draw_new_line(coord_val x1, coord_val y1, coord_val z1,
                       coord_val x2, coord_val y2, coord_val z2) {

   fprintf(stdout, "{%g, %g, %g} P %.3f %.3f %.3f\n",
	       x1, y1, z1, x1, y1, z1);

   kin_continue_line(x1, y1, z1, x2, y2, z2);
}

void kin_continue_line(coord_val x1, coord_val y1, coord_val z1,
                       coord_val x2, coord_val y2, coord_val z2) {
   coord_val xav, yav, zav;
   xav = 0.5*(x1+x2);
   yav = 0.5*(y1+y2);
   zav = 0.5*(z1+z2);

   fprintf(stdout, "{%g, %g, %g} L %.3f %.3f %.3f\n",
	       xav, yav, zav, xav, yav, zav);
   fprintf(stdout, "{%g, %g, %g} L %.3f %.3f %.3f\n",
	       x2, y2, z2, x2, y2, z2);
}

void input_sigma_to_contour_levels(progparms *parms) {
   if (parms->sigmalevels) { /* sanity check */
      convert_sigma_contour_levels(parms);
   }
}

void calc_default_contour_levels(progparms *parms) {
   int nlev = 0;

   parms->inlvls[nlev] = parms->levels[nlev] = -1.0;
   parms->colors[nlev] = "grey";
   nlev++;
   if (parms->verbose) {
      fprintf(stderr, "auto level %d s.d.: %g %s\n", nlev,
                       parms->levels[nlev-1], parms->colors[nlev-1]);
   }

   parms->inlvls[nlev] = parms->levels[nlev] = 0.0;
   parms->colors[nlev] = "brown";
   nlev++;
   if (parms->verbose) {
      fprintf(stderr, "auto level %d s.d.: %g %s\n", nlev,
                       parms->levels[nlev-1], parms->colors[nlev-1]);
   }

   parms->inlvls[nlev] = parms->levels[nlev] = +1.0;
   parms->colors[nlev] = "orange";
   nlev++;
   if (parms->verbose) {
      fprintf(stderr, "auto level %d s.d.: %g %s\n", nlev,
                       parms->levels[nlev-1], parms->colors[nlev-1]);
   }
   parms->nlevels = nlev;

   parms->sigmalevels = 1;

   convert_sigma_contour_levels(parms);
}

void convert_sigma_contour_levels(progparms *parms) {
   int ix = 0, iy = 0, iz = 0;
   int bx = 0, by = 0, bz = 0, ex = 0, ey = 0, ez = 0;
   int nlev = 0, numsamp = 0;
   double sum = 0.0, dev = 0.0, ssq = 0.0;
   double mean = 0.0, sdev = 0.0, lvl = 0.0, slvl = 0.0;

   if (parms->wrap) {
      bx = by = bz = 0;
      ex = parms->nx - 1;
      ey = parms->ny - 1;
      ez = parms->nz - 1;
   }
   else {
      bx = 1 + ((parms->nfx - 1)/2); /* ignore the ragged edges */
      by = 1 + ((parms->nfy - 1)/2);
      bz = 1 + ((parms->nfz - 1)/2);
      ex = parms->nx - bx;
      ey = parms->ny - bx;
      ez = parms->nz - bx;
   }

   numsamp = (ex - bx)*(ey - by)*(ez - bz);
   if (numsamp < 1) {
      fprintf(stderr, "ERROR: %d points when edges ignored\n", numsamp);
      return; /* should we just exit? */
   }

   sum = 0.0;
   for (ix = bx; ix < ex; ix++) {
      for (iy = by; iy < ey; iy++) {
	 for (iz = bz; iz < ez; iz++) {
	    sum += parms->lat[ix][iy][iz];
	 }
      }
   }
   mean = sum/numsamp;

   ssq = 0.0;
   for (ix = bx; ix < ex; ix++) {
      for (iy = by; iy < ey; iy++) {
	 for (iz = bz; iz < ez; iz++) {
	    dev = parms->lat[ix][iy][iz] - mean;
	    ssq += dev*dev;
	 }
      }
   }
   sdev = sqrt(ssq/(numsamp - 1));

   if (parms->verbose) {
      fprintf(stderr, "data stats on %d samples: mean %g, s.d. %g\n", numsamp, mean, sdev);
   }

   /* --------------- */

   for (nlev = 0; nlev < parms->nlevels; nlev++) {
      slvl = parms->levels[nlev];
      lvl = mean + slvl*sdev;
      if (lvl > -CHOP_EPSILON && lvl < CHOP_EPSILON) {
         lvl = 0.0;   /* chop levels very close to zero */
      }
      parms->levels[nlev] = lvl;
      if (parms->verbose) {
         fprintf(stderr, "converting level %d: %g s.d. => %g\n", nlev+1, slvl, lvl);
      }
   }
}

void dump_lattice(progparms *parms) {
   int ix = 0, iy = 0, iz = 0, nx = 0, ny = 0, nz = 0;
   double xpos= 0.0, ypos= 0.0, zpos= 0.0;

   if (parms->wrap) {
      nx = parms->nx - 1;
      ny = parms->ny - 1;
      nz = parms->nz - 1;
   }
   else {
      nx = parms->nx;
      ny = parms->ny;
      nz = parms->nz;
   }

   xpos = parms->x0;
   for (ix = 0; ix < nx; ix++) {
      ypos = parms->y0;
      for (iy = 0; iy < ny; iy++) {
	 zpos = parms->z0;
	 for (iz = 0; iz < nz; iz++) {

	    fprintf(stdout, "%g %g %g %g\n",
	       parms->lat[ix][iy][iz], xpos, ypos, zpos);

	    zpos += parms->dz;
	 }
	 ypos += parms->dy;
      }
      xpos += parms->dx;
   }

}
