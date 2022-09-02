/*			    kin2Dcont.c                         */
/* Author: J. Michael Word              Date Written: 2/12/99   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

/* Modifications:                                               */
/*  2/12/99 - jmw - created as a modified kin3Dcont             */
/*  2/16/99 - jmw - added -scale and -ps flags                  */
/*  2/18/99 - jmw - added -wrap flag to support periodic bounds */
/*  2/22/99 - jmw - added -sampled flag for sampled datasets    */
/*                  also considers average for default contour  */
/*  2/24/99 - jmw - default contours +/- 1 sigma smoothed data  */
/*  2/26/99 - jmw - modified -dump function                     */
/*  3/ 6/99 - jmw - fixed minor formatting problem in caption   */
/*  4/26/02 - jmw - added -nosmooth flag                        */
/*  4/30/02 - jmw - modified @caption output for -nosmooth flag */
/*  7/24/02 - jmw - added -alignminxyz flag                     */
/*  7/25/02 - jmw - added -SNAP, -SL, -SM, and fixed -SP        */
/*  8/30/22 - mgp - changed format ‘%d’ to %ld in printf w/ long int  */
/*  8/30/22 - mgp - added prototype for remap_slices                  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "connect.h"
#include "readpoints.h"

#define MAXLEVELS 300

#define DEFAULT_DENS  1.0
#define DEFAULT_FSDEV 2.0
#define CHOP_EPSILON  0.00000000001

/* span is filter cuttoff in standard deviations */
#define DEFAULT_SPAN   2.0

#define DEFAULT_ZSCALE 1.0

#define PROGRAM_NAME "kin2Dcont"
//#define VERSION_STR "kin2Dcont: version 1.8 7/25/02, Copyright 1999-2002, J. Michael Word"
#define VERSION_STR "kin2Dcont: version 1.81 8/30/02, Copyright 1999-2022, J. Michael Word and Duke University"

typedef struct progparms_t {
   int       verbose; /* print out descriptive info? */
   int       usesmoothing; /* do the smoothing thing? */
   char*     infname; /* input file name */

   int       nlevels; /* contour level & style */
   ctree_lvl levels[MAXLEVELS];
   char*     colors[MAXLEVELS];
   int       sigmalevels; /* true if levels input as sigma */

   int       wrap; /* do the edges wrap arround? */

   coord_val dx,  dy;   /* lattice spacing  */
   coord_val sx,  sy;   /* stddev of filter */
   coord_val filtspan;  /* filter cuttoff   */

   int       sampled; /* sampled data (not normalized) */
   int       alignGridToMin; /* grid aligned with min value */
   int       snap; /* if data very close to grid make it on grid */
   
   int        vfirst; /* value first (else last) */
   char*      delim;  /* special delimiter       */
   val3Dlist* vl;     /* list of input points    */

   int       writekin; /* display @kinemage */
   int       writegrp; /* display @group    */
   int       writesub; /* display @subgroup */
   int       dominant; /* make the grp or subgrp dominant? */
   int       lens;     /* add lens keyword? */
   int       axis;     /* draw box edges? */
   char*     grpname;  /* group or subgroup name */
   int       sliceZ;   /* do we show levels in Z? */
   int       dumplattice; /* write lattice values for debugging */
   float     zscale;   /* controls output of layers in kin file */
   int       psoutput; /* write PostScript output (vs .kin) */

   coord_val xmin, ymin; /* bounds for data when */
   coord_val xmax, ymax; /* points wrap around   */

   /* parameters calculated from input */
   int       ifx,  ify;  /* integer filter radii      */
   int       nfx,  nfy;  /* number of filter points   */
   int       nx,   ny;   /* number of lattice points  */
   coord_val x0,   y0;   /* coord of lattice[0][0]    */

   real_val **lat;  /* two dimensional data lattice */
   real_val **filt; /* two dimensional filter */
} progparms;

val3Dlist* processCommandline(int argc, char *argv[], progparms *parms);
void kin_draw_axis_box(progparms *parms);
void kin_draw_line(coord_val x1,coord_val y1,coord_val x2,coord_val y2);
void ps_draw_axis_box(progparms *parms);
void analyze_data(progparms *parms);
int insert_data_into_lattice(progparms *parms);
void splat(progparms *parms, real_val wt, coord_val x, coord_val y);
void placemask(progparms *parms, real_val wt, int lx, int ly);
void placemaskwrap(progparms *parms, real_val wt, int lx, int ly);
void placesample(progparms *parms, real_val wt, int lx, int ly);
void input_sigma_to_contour_levels(progparms *parms);
void calc_default_contour_levels(progparms *parms);
void convert_sigma_contour_levels(progparms *parms);
void contour_lattice(progparms *parms);
real_val** build_lattice(progparms *parms);
real_val** build_filter(progparms *parms);
void destroy_lattice(progparms *parms);
void destroy_filter(progparms *parms);
void describeInput(progparms *parms);
void helpInfo(int fullhelp);
int compArgStr(char *str, char *arg, int min);
void dup_periodic_edges(progparms *parms);
void xy_contour(progparms *parms, cont_3d_info *cp);
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
	    fprintf(stdout, "#lattice gap: [%g, %g]\n", parms->dx, parms->dy);
            if (parms->usesmoothing) {
	       fprintf(stdout, "#filter stddev: [%g, %g]\n", parms->sx, parms->sy);
            }
            else {
               fprintf(stdout, "#no smoothing\n");
            }
	    if (parms->wrap) {
	       fprintf(stdout, "#bounding rect: x[%g .. %g] y[%g .. %g] *wrap*\n",
		  parms->xmin,  parms->xmax - parms->dx,
		  parms->ymin,  parms->ymax - parms->dy);
	    }
	    else {
	       fprintf(stdout, "#bounding rect: x[%g .. %g] y[%g .. %g]\n",
		  parms->xmin,  parms->xmax,
		  parms->ymin,  parms->ymax);
	    }
	 }
      }
      else if (parms->psoutput) {
	 simple_2Dps_procset(stdout, 1);
	 if (parms->axis) { ps_draw_axis_box(parms); }
      }
      else {
	 if (parms->writekin) {
	    fprintf(stdout, "@kinemage 1\n");
	 }
	 fprintf(stdout, "@caption %s\n", VERSION_STR);
	 fprintf(stdout, " input: %s, %ld coords\n", parms->infname, parms->vl->ncoords);
	 fprintf(stdout, " lattice gap: [%g, %g],", parms->dx, parms->dy);
	 if (parms->usesmoothing) {
	    fprintf(stdout, " filter stddev: [%g, %g]\n", parms->sx, parms->sy);
         }
	 else {
	    fprintf(stdout, " no smoothing\n");
	 }
	 if (parms->wrap) {
	    fprintf(stdout, " wrap rect: x[%g .. %g] y[%g .. %g]\n",
	       parms->xmin,  parms->xmax,
	       parms->ymin,  parms->ymax);
	 }
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
	 if (parms->axis) { kin_draw_axis_box(parms); }
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

   fprintf(stderr, "Read x, y, value or value, x, y sample data and\n");
   fprintf(stderr, "place on a 2D grid using a gaussian spot function and\n");
   fprintf(stderr, "generate one or more 2D contour levels in kinemage format.\n\n");
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
   fprintf(stderr, "  -GXY #.# #.#      spacing between lattice points in each dim\n");
}
   fprintf(stderr, "  -Stddev#.#        spot filter stddev (default %g)\n", DEFAULT_FSDEV);
if (fullhelp) {
   fprintf(stderr, "  -SXY #.# #.#      spot filter stddev in each dim\n");
}
   fprintf(stderr, "  -WRAP #.# #.# #.# #.#  bounds of rectangle where edges wrap\n");
   fprintf(stderr, "                         (4 numbers: xmin xmax ymin ymax)\n");
   fprintf(stderr, "  -SNAPtogrid       force to grid data points which are very close to grid\n");
   fprintf(stderr, "  -SAMPled          input data is a uniform grid of samples\n");
   fprintf(stderr, "  -ALIGNminxyz      sampled and snapped with grid aligned to min xyz\n");
   fprintf(stderr, "\n(input source)\n");
   fprintf(stderr, "  -FIRST            value in first column then x,y (default)\n");
   fprintf(stderr, "  -LAST             value in last column after x,y\n");
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
   fprintf(stderr, "  -DOMinant         add dominant keyword to @group or @subgroup\n");
   fprintf(stderr, "                    (supresses header when used with -dump)\n");
   fprintf(stderr, "  -LENS             add lens keyword to @group or @subgroup\n");
   fprintf(stderr, "  -NOAXIS           do not draw axis lines marking the box edges\n");
   fprintf(stderr, "  -PS               output in PostScript rather than kinemage format\n");
   fprintf(stderr, "  -FLAT             make the contour levels flat on the 2D plane\n");
   fprintf(stderr, "  -T_               Use _ as the field separator character(s).\n");
   fprintf(stderr, "                    For example -t@ would separate on at-signs.\n");
   fprintf(stderr, "                    In the input, multiple consecutive separators count as one.\n");
   fprintf(stderr, "(The default delimiters are: blank, tab, newline, return, comma and colon)\n");
   fprintf(stderr, "  -NOSMooth         do not attempt to smooth the input data\n");
   fprintf(stderr, "  -SPan#.#          spot filter span in stddev (default %g)\n", DEFAULT_SPAN);
   fprintf(stderr, "  -SCALE#.#         Z dimension scale factor (default %g)\n", DEFAULT_ZSCALE);
   fprintf(stderr, "  -DUMPlattice      output lattice data for debugging\n\n");
   fprintf(stderr, "  -Verbose          displays extra info about the calculation\n");
   fprintf(stderr, "  -Quiet            suppresses extra output\n");
   fprintf(stderr, "\n  -Help             prints this help message\n");
   fprintf(stderr, "\nExample: %s -g5 -s10 -l 125 sea mydata >>my.kin\n", PROGRAM_NAME);
   fprintf(stderr, "generates a single 2D contour with a spacing of 5 and a filter stddev of 10\n");
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
   fprintf(stderr, "     delta: %g, %g\n", parms->dx, parms->dy);
   if (parms->usesmoothing) {
      fprintf(stderr, "filt. sdev: %g, %g\n", parms->sx, parms->sy);
      fprintf(stderr, "filt. span: %g stddev\n",  parms->filtspan);

      fprintf(stderr, "filter dim: %d .. %d, %d .. %d (%d total)\n",
         -parms->ifx, parms->ifx,
         -parms->ify, parms->ify, parms->nfx*parms->nfy);
   }
   else {
      fprintf(stderr, " no smoothing\n");
   }

   fprintf(stderr, "  data dim: %d, %d (%d total)\n",
      parms->nx,  parms->ny, parms->nx*parms->ny);
   fprintf(stderr, "    origin: %g, %g\n", parms->x0,  parms->y0);
   if (parms->wrap) {
      fprintf(stderr, "   *wrap x: %g .. %g\n", parms->xmin,  parms->xmax);
      fprintf(stderr, "   *wrap y: %g .. %g\n", parms->ymin,  parms->ymax);
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
   parms->colors[0] = NULL;
   parms->sigmalevels = 0;
   parms->wrap = 0;
   parms->dx  = parms->dy  = DEFAULT_DENS;
   parms->sx  = parms->sy  = DEFAULT_FSDEV;
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
   parms->lens     = 0;
   parms->axis     = 1;
   parms->sliceZ   = 1;
   parms->grpname  = "cont";

   parms->dumplattice = 0;
   parms->zscale = DEFAULT_ZSCALE;
   parms->psoutput = 0;

   parms->xmin = parms->ymin = 0.0;
   parms->xmax = parms->ymax = 0.0;

   /* zero out stuff to be calculated from data */

   parms->ifx = parms->ify = 0;
   parms->nfx = parms->nfy = 0;
   parms->nx  = parms->ny  = 0;
   parms->x0  = parms->y0  = 0.0;

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
	 else if(compArgStr(p+1, "PS", 2)){
	    parms->psoutput = 1;
	    parms->writekin = 0;
	    parms->writegrp = 0;
	    parms->writesub = 0;
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
	 else if(compArgStr(p+1, "LENS", 4)){
	    parms->lens = 1;
	 }
	 else if(compArgStr(p+1, "NOAXIS", 6)){
	    parms->axis = 0;
	 }
	 else if(compArgStr(p+1, "FLAT", 4)){
	    parms->sliceZ = 0;
	 }
	 else if(compArgStr(p+1, "T", 1)){
	    parms->delim = p+2;
	 }
	 else if(n=compArgStr(p+1, "Grain", 1)){
            rval = atof(p+n+1);
	    if (rval < 0.0) { rval = -rval; }
	    parms->dx = parms->dy = rval;
	 }
	 else if(compArgStr(p+1, "GXY", 3)){
	    if (i+2 < argc) {
	       parms->dx = atof(argv[i+1]);
	       parms->dy = atof(argv[i+2]);
	       if (parms->dx < 0.0) { parms->dx = -parms->dx; }
	       if (parms->dy < 0.0) { parms->dy = -parms->dy; }
	       i += 2;
	    }
	    else {
	       fprintf(stderr, "ERROR: -GXY flag requires 2 numbers\n");
	       helpInfo(0);
	    }
	 }
	 else if(compArgStr(p+1, "NOSMooth", 4)){
	    parms->usesmoothing = 0;
	 }
	 else if(n=compArgStr(p+1, "Stddev", 1)){
            rval = atof(p+n+1);
	    if (rval < 0.0) { rval = -rval; }
	    parms->sx  = parms->sy  = rval;
	 }
	 else if(compArgStr(p+1, "SXY", 3)){
	    if (i+2 < argc) {
	       parms->sx = atof(argv[i+1]);
	       parms->sy = atof(argv[i+2]);
	       if (parms->sx < 0.0) { parms->sx = -parms->sx; }
	       if (parms->sy < 0.0) { parms->sy = -parms->sy; }
	       i += 2;
	    }
	    else {
	       fprintf(stderr, "ERROR: -RXY  flag requires 2 numbers\n");
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
	       parms->levels[n] = atof(argv[i+1]);
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
		  parms->levels[n] = lcv;
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
	 else if(n=compArgStr(p+1, "SCALE", 5)){
            rval = atof(p+n+1);
	    parms->zscale = rval;
	 }
	 else if(compArgStr(p+1, "DUMPlattice", 4)){
	    parms->dumplattice = 1;
	 }
	 else if(compArgStr(p+1, "WRAP", 4)){
	    if (i+4 < argc) {
	       parms->wrap = 1;
	       parms->xmin = atof(argv[i+1]);
	       parms->xmax = atof(argv[i+2]);
	       parms->ymin = atof(argv[i+3]);
	       parms->ymax = atof(argv[i+4]);
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
	       i += 4;
	    }
	    else {
	       fprintf(stderr, "ERROR: -WRAP flag requires 4 numbers\n");
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

   if ((parms->dx < CHOP_EPSILON) || (parms->dy < CHOP_EPSILON)) {
      fprintf(stderr, "ERROR: grid spacing must be greater than zero: %g, %g\n",
                       parms->dx, parms->dy);
   }
   else {
      parms->vl = load_coords(inf, 2, parms->vfirst, parms->delim);
   }
   if (inf != stdin) { fclose(inf); }
   

   if (! parms->usesmoothing) {
      parms->sx  = parms->sy  = 0.0;
   }

   if (parms->vl) {
      if (parms->vl->ncoords > 0) { /* work out the calculated parms */
	 double dval = 0.0;

	 /* size the spot filter (integer radii) */
	 parms->ifx = ceil(parms->filtspan * parms->sx/parms->dx);
	 parms->ify = ceil(parms->filtspan * parms->sy/parms->dy);

	 parms->nfx = 2*parms->ifx + 1; /* total filter dim */
	 parms->nfy = 2*parms->ify + 1;

/* these deal with where the range is whole number of grain val */
/* either go over */
#define WIGGLE_ROOM ( 0.0001)
/* or round with a bias towards upward (needed in a world with roundoff) */
#define SHRINKAGE   (-0.25 )

	 if (parms->wrap) {
	    /* period of the lattice + edge */
	    parms->nx = ceil((parms->xmax - parms->xmin)/parms->dx) + 1;
	    parms->ny = ceil((parms->ymax - parms->ymin)/parms->dy) + 1;

	    /* fine tune the lattice spacing (for integer period) */
	    parms->dx = (parms->xmax - parms->xmin)/(parms->nx - 1);
	    parms->dy = (parms->ymax - parms->ymin)/(parms->ny - 1);

	    /* lattice origin */
	    parms->x0 = parms->xmin;
	    parms->y0 = parms->ymin;
	 }
	 else if (parms->alignGridToMin) {
	    /* period of the lattice + filter + border */
	    n = ceil(SHRINKAGE + (parms->vl->xmax - parms->vl->xmin)/parms->dx) + 1;
	    parms->nx = n + (parms->nfx - 1) + 2;
	    n = ceil(SHRINKAGE + (parms->vl->ymax - parms->vl->ymin)/parms->dy) + 1;
	    parms->ny = n + (parms->nfy - 1) + 2;

	    /* lattice origin */
	    /* lattice origin */
	    parms->xmin = parms->vl->xmin - (0.5*(parms->nfx - 1) + 1.0)*parms->dx;
	    parms->ymin = parms->vl->ymin - (0.5*(parms->nfy - 1) + 1.0)*parms->dy;

	    parms->xmax = parms->vl->xmax + (0.5*(parms->nfx - 1) + 1.0)*parms->dx;
	    parms->ymax = parms->vl->ymax + (0.5*(parms->nfy - 1) + 1.0)*parms->dy;

	    parms->x0 = parms->xmin;
	    parms->y0 = parms->ymin;

	    /* fine tune the lattice spacing (for integer period) */
	    parms->dx = (parms->xmax - parms->xmin)/(parms->nx - 1);
	    parms->dy = (parms->ymax - parms->ymin)/(parms->ny - 1);
	 }
	 else { /* no wrap - center box on data */
	    /* size the lattice */
	    n = ceil(WIGGLE_ROOM +
	             (parms->vl->xmax - parms->vl->xmin)/parms->dx) + 1;
	    parms->nx = n + (parms->nfx - 1);
	    n = ceil(WIGGLE_ROOM +
	             (parms->vl->ymax - parms->vl->ymin)/parms->dy) + 1;
	    parms->ny = n + (parms->nfy - 1);

	    /* position the origin */
	    dval = parms->vl->xmin + 0.5*(parms->vl->xmax - parms->vl->xmin);
	    parms->x0 = dval - 0.5*(parms->nx - 1)*parms->dx;

	    dval = parms->vl->ymin + 0.5*(parms->vl->ymax - parms->vl->ymin);
	    parms->y0 = dval - 0.5*(parms->ny - 1)*parms->dy;

	    parms->xmin = parms->x0; /* bounds of lattice */
	    parms->ymin = parms->y0;
	    parms->xmax = parms->x0 + parms->dx*(parms->nx - 1);
	    parms->ymax = parms->y0 + parms->dy*(parms->ny - 1);
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

real_val** build_lattice(progparms *parms) {
   int ix=0, iy=0, i = 0;
   coord_val rval = 0.0;

   if (parms == NULL) { return NULL; }

   parms->lat  = matrix(0, parms->nx -1,
                        0, parms->ny -1);
   if (parms->lat) {
      for (ix = 0; ix < parms->nx; ix++) {
	 for (iy = 0; iy < parms->ny; iy++) {
	    parms->lat[ix][iy] = 0.0; /* initialize to zero */
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
#endif

   return parms->lat;
}

void destroy_lattice(progparms *parms) {
   free_matrix(parms->lat, 0, 0);
}

real_val** build_filter(progparms *parms) {
   int ix=0, iy=0;
   real_val sumval=0.0, scalef=0.0;
   real_val qx=0.0, qy=0.0;
   
   if (parms == NULL) { return NULL; }

   parms->filt = matrix(-parms->ifx, parms->ifx,
                        -parms->ify, parms->ify);

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

   sumval = 0.0;
   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	    parms->filt[ix][iy]
		  = exp(-(  ix*ix*qx + iy*iy*qy ));
	    sumval += parms->filt[ix][iy];
      }
   }

   scalef = 1.0/sumval;
   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	    parms->filt[ix][iy] *= scalef;
      }
   }

   if (parms->verbose) { /* display the filter */
      fprintf(stderr, "100 X Y>");
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 fprintf(stderr, "%7d ", iy);
      }
      fprintf(stderr, "\n");
      for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
	 fprintf(stderr, "%7d ", ix);
	 for (iy = -parms->ify; iy <= parms->ify; iy++) {
	    fprintf(stderr, "%7.3f ", parms->filt[ix][iy]*100.0);
	 }
	 fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
   }

   return parms->filt;
}

void destroy_filter(progparms *parms) {
   if (parms && parms->filt) {
      free_matrix(parms->filt, -parms->ifx, -parms->ify);
   }
}

int insert_data_into_lattice(progparms *parms) {
   val3Dnode* v3D = NULL;
   int n = 0;
   float pctfact = 0.0, perarea = 0.0;

   pctfact = 100.0/parms->vl->ncoords;

   /* normalization factor */
   if (parms->sampled) { perarea = 1.0; }
   else {
      perarea = 1.0/(parms->dx * parms->dy);
   }

   if (parms && parms->vl && parms->vl->head) {
      if (parms->verbose) {
	 fprintf(stderr, "inserting data into lattice\n");
      }
      for(v3D = parms->vl->head; v3D; v3D = v3D->next) {
	 n++;
	 splat(parms, perarea * v3D->v, v3D->x, v3D->y);

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

void splat(progparms *parms, real_val wt, coord_val x, coord_val y) {
   int lx = 0, ly = 0;
   coord_val deltax = 0.0, deltay = 0.0;
   coord_val gridx  = 0.0, gridy  = 0.0;
   coord_val fractx = 0.0, fracty = 0.0;

   deltax = (x - parms->x0)/parms->dx;
   deltay = (y - parms->y0)/parms->dy;

   gridx = floor(deltax);
   gridy = floor(deltay);

   lx = (int)gridx;
   ly = (int)gridy;

   fractx = deltax - gridx;
   fracty = deltay - gridy;

   
   if (parms->snap) {
      if (fractx < 0.01) { fractx = 0.0; }
      if (fracty < 0.01) { fracty = 0.0; }
      if (fractx > 0.99) { fractx = 1.0; }
      if (fracty > 0.99) { fracty = 1.0; }
   }

   if (parms->usesmoothing) {
      if (parms->wrap) {
         /* skip duplicate points past wrap edge */
         if ((x < parms->xmax) && (y < parms->ymax)) {
	    placemaskwrap(parms, wt*(1.0-fractx)*(1.0-fracty), lx,   ly);
	    placemaskwrap(parms, wt*(    fractx)*(1.0-fracty), lx+1, ly);
	    placemaskwrap(parms, wt*(1.0-fractx)*(    fracty), lx,   ly+1);
	    placemaskwrap(parms, wt*(    fractx)*(    fracty), lx+1, ly+1);
         }
      }
      else {
         placemask(parms, wt*(1.0-fractx)*(1.0-fracty), lx,   ly);
         placemask(parms, wt*(    fractx)*(1.0-fracty), lx+1, ly);
         placemask(parms, wt*(1.0-fractx)*(    fracty), lx,   ly+1);
         placemask(parms, wt*(    fractx)*(    fracty), lx+1, ly+1);
      }
   }
   else { /* no filter required */
      placesample(parms, wt*(1.0-fractx)*(1.0-fracty), lx,   ly);
      placesample(parms, wt*(    fractx)*(1.0-fracty), lx+1, ly);
      placesample(parms, wt*(1.0-fractx)*(    fracty), lx,   ly+1);
      placesample(parms, wt*(    fractx)*(    fracty), lx+1, ly+1);
   }
}

void placemask(progparms *parms, real_val wt, int lx, int ly) {
   int ix = 0, iy = 0;

   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	    parms->lat[lx+ix][ly+iy] += wt * parms->filt[ix][iy];
      }
   }
}

void placemaskwrap(progparms *parms, real_val wt, int lx, int ly) {
   int ix = 0, iy = 0;
   int mx = 0, my = 0;

   for (ix = -parms->ifx; ix <= parms->ifx; ix++) {
      mx = lx + ix;
      while (mx < 0) { mx += (parms->nx - 1); }
      while (mx >= (parms->nx - 1)) { mx -= (parms->nx - 1); }
      for (iy = -parms->ify; iy <= parms->ify; iy++) {
	 my = ly + iy;
	 while (my < 0) { my += (parms->ny - 1); }
	 while (my >= (parms->ny - 1)) { my -=  (parms->ny - 1); }
	 parms->lat[mx][my] += wt * parms->filt[ix][iy];
      }
   }
}

void placesample(progparms *parms, real_val wt, int lx, int ly) {
   parms->lat[lx][ly] += wt;
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
                        parms->colors, 2);

   if (cp) {
      if (parms->wrap) {
	 if (parms->verbose) { fprintf(stderr, "duplicating edges\n"); }
	 dup_periodic_edges(parms);
      }

      if (parms->verbose) { fprintf(stderr, "calculating contours\n"); }

      xy_contour(parms, cp);

      if (parms->psoutput) {
	 if (parms->verbose) { fprintf(stderr, "writing to ps file.\n"); }

	 simple_2Dps_contour_output(stdout, cp, 0);
	 simple_2Dps_showpage(stdout);
      }
      else {
	 if (parms->verbose) { fprintf(stderr, "writing to kin file.\n"); }

	 scaled_kin_contour_output(stdout, cp, 0, 1, parms->zscale);
      }

      destroy_cont_info(cp);
   }
}

void dup_periodic_edges(progparms *parms) {
   int ix = 0, iy = 0;
   int mx = 0, my = 0;

   mx = parms->nx - 1;
   my = parms->ny - 1;

   for (ix = 0; ix < mx; ix++) {
      parms->lat[ix][my] = parms->lat[ix][0];
   }
   for (iy = 0; iy < my; iy++) {
      parms->lat[mx][iy] = parms->lat[0][iy];
   }
   parms->lat[mx][my] = parms->lat[0][0];
}

void xy_contour(progparms *parms, cont_3d_info *cp) {
   int rc = 0;
   int r = 0, c = 0, co = 0, ro = 0, crgn = 0, rrgn = 0;
   coord_val xspan = 0.0, yspan = 0.0;
   ctree_data data;
   ctree_val **rawdata = NULL;
   ctree_val *work1 = NULL, *work2 = NULL;
   slice_info *sp=NULL;
   int dimmap[3];

   dimmap[0] = X_DIM_NUM;
   dimmap[1] = Y_DIM_NUM;
   dimmap[2] = Z_DIM_NUM;

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
   co   = 0;   ro = 0;
   crgn = 0; rrgn = 0; /* i.e. calculate region */

   rc = ctree_compose_rawdata(&data, &(rawdata[0][0]),
	    parms->x0, parms->y0, xspan, yspan,
	    c, r, co, ro, crgn, rrgn, &(work1[0]), &(work2[0]), NULL);

   if (rc) {
      sp = build_slice_info(parms->levels, parms->nlevels, dimmap,
	 parms->dx, parms->dy, (! parms->sliceZ), 0.0, 1);

      if (sp) {
	 for (c = 0; c < parms->nx; c++) { /* load array */
	    for (r = 0; r < parms->ny; r++) {
	       rawdata[r][c] = parms->lat[c][r];
	    }
	 }

	 ctree_plot_region(&data, parms->levels, parms->nlevels,
		  basic_data_func, connecting_plot_func, sp);
	 slice_all_done(sp);

	 remap_slices(cp, sp);

	 update_slice_info(sp, dimmap, parms->dx, parms->dy,
	       (! parms->sliceZ), 0.0);
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

void kin_draw_axis_box(progparms *parms) {
   fprintf(stdout, "@vectorlist {axis} color= white width= 1\n");

   kin_draw_line(parms->xmin, parms->ymin, parms->xmax, parms->ymin);
   kin_draw_line(parms->xmin, parms->ymax, parms->xmax, parms->ymax);
   kin_draw_line(parms->xmin, parms->ymin, parms->xmin, parms->ymax);
   kin_draw_line(parms->xmax, parms->ymin, parms->xmax, parms->ymax);
}

void kin_draw_line(coord_val x1,coord_val y1,coord_val x2,coord_val y2) {
   coord_val xav, yav, zav;
   xav = 0.5*(x1+x2);
   yav = 0.5*(y1+y2);

   fprintf(stdout, "{%g, %g} P %.3f %.3f 0\n",
	       x1, y1, x1, y1);
   fprintf(stdout, "{%g, %g} L %.3f %.3f 0\n",
	       xav, yav, xav, yav);
   fprintf(stdout, "{%g, %g} L %.3f %.3f 0\n",
	       x2, y2, x2, y2);
}

void ps_draw_axis_box(progparms *parms) {
   fprintf(stdout, "newpath %.3f %.3f M ",
	       parms->xmin, parms->ymin);
   fprintf(stdout, "%.3f %.3f L closepath stroke\n",
	       parms->xmax, parms->ymin);
   fprintf(stdout, "newpath %.3f %.3f M ",
	       parms->xmin, parms->ymax);
   fprintf(stdout, "%.3f %.3f L closepath stroke\n",
	       parms->xmax, parms->ymax);

   fprintf(stdout, "newpath %.3f %.3f M ",
	       parms->xmin, parms->ymin);
   fprintf(stdout, "%.3f %.3f L closepath stroke\n",
	       parms->xmin, parms->ymax);
   fprintf(stdout, "newpath %.3f %.3f M ",
	       parms->xmax, parms->ymin);
   fprintf(stdout, "%.3f %.3f L closepath stroke\n",
	       parms->xmax, parms->ymax);
}

void input_sigma_to_contour_levels(progparms *parms) {
   if (parms->sigmalevels) { /* sanity check */
      convert_sigma_contour_levels(parms);
   }
}

void calc_default_contour_levels(progparms *parms) {
   int nlev = 0;

   parms->levels[nlev] = -1.0;
   parms->colors[nlev] = "grey";
   nlev++;
   if (parms->verbose) {
      fprintf(stderr, "auto level %d s.d.: %g %s\n", nlev,
                       parms->levels[nlev-1], parms->colors[nlev-1]);
   }

   parms->levels[nlev] = 0.0;
   parms->colors[nlev] = "brown";
   nlev++;
   if (parms->verbose) {
      fprintf(stderr, "auto level %d s.d.: %g %s\n", nlev,
                       parms->levels[nlev-1], parms->colors[nlev-1]);
   }

   parms->levels[nlev] = +1.0;
   parms->colors[nlev] = "orange";
   nlev++;
   if (parms->verbose) {
      fprintf(stderr, "auto level %d s.d.: %g %s\n", nlev,
                       parms->levels[nlev-1], parms->colors[nlev-1]);
   }
   parms->nlevels = nlev;

   convert_sigma_contour_levels(parms);
}


void convert_sigma_contour_levels(progparms *parms) {
   int ix = 0, iy = 0;
   int bx = 0, by = 0, ex = 0, ey = 0;
   int nlev = 0, numsamp = 0;
   double sum = 0.0, dev = 0.0, ssq = 0.0;
   double mean = 0.0, sdev = 0.0, lvl = 0.0, slvl = 0.0;

   if (parms->wrap) {
      bx = by = 0;
      ex = parms->nx - 1;
      ey = parms->ny - 1;
   }
   else {
      bx = 1 + ((parms->nfx - 1)/2); /* ignore the ragged edges */
      by = 1 + ((parms->nfy - 1)/2);
      ex = parms->nx - bx;
      ey = parms->ny - by;
   }

   numsamp = (ex - bx)*(ey - by);
   if (numsamp < 1) {
      fprintf(stderr, "ERROR: %d points when edges ignored\n", numsamp);
      return; /* should we just exit? */
   }

   sum = 0.0;
   for (ix = bx; ix < ex; ix++) {
      for (iy = by; iy < ey; iy++) {
	    sum += parms->lat[ix][iy];
      }
   }
   mean = sum/numsamp;

   ssq = 0.0;
   for (ix = bx; ix < ex; ix++) {
      for (iy = by; iy < ey; iy++) {
	    dev = parms->lat[ix][iy] - mean;
	    ssq += dev*dev;
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
   int ix = 0, iy = 0, nx = 0, ny = 0;
   double xpos= 0.0, ypos= 0.0;

   if (parms->wrap) {
      nx = parms->nx - 1;
      ny = parms->ny - 1;
   }
   else {
      nx = parms->nx;
      ny = parms->ny;
   }

   xpos = parms->x0;
   for (ix = 0; ix < nx; ix++) {
      ypos = parms->y0;
      for (iy = 0; iy < ny; iy++) {

	 fprintf(stdout, "%g %g %g\n",
	       parms->lat[ix][iy], xpos, ypos);

	 ypos += parms->dy;
      }
      xpos += parms->dx;
   }

}
