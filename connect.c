/*			    connect.c				*/
/* Author: J. Michael Word		Date Written: 1/19/99	*/

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include "connect.h"

void print_cont_info(FILE *fp, cont_3d_info *cp);
void print_slice_info(FILE *fp, slice_info *sp);
int remap_slices(cont_3d_info *cp, slice_info *sp);
void remap_line_segment(line_segment* lsp, int map[]);
void destroy_line_seg(line_segment* lsp);
void reverse_line_segment(line_segment* lsp);
line_segment*
join_segments(line_segment* ls1, line_segment* ls2,
	       int selection1, int selection2);
int line_seg_ends_match(line_segment* ls1, line_segment* ls2,
		  int *selection1, int *selection2,
		  coord_val epsilon0, coord_val epsilon1);
int check_for_loop(line_segment* lsp, coord_val epsilon0, coord_val epsilon1);
int line_seg_too_far_back(line_segment* lsp,
	 coord_val r1, coord_val r2, coord_val cuttoff);
int insert_one_line_segment(int level, slice_info* sp, line_segment* lsp);
void print_line_segment(FILE *fp, line_segment* lsp);

cont_3d_info*
build_cont_info(ctree_lvl level[], int nlevels, char* props[], int ndim) {
   cont_3d_info *cp = NULL;
   level_3d_segs *level_segs = NULL;
   int i = 0;
   char *prop_str = NULL;

   if (nlevels < 1) {
      fprintf(stderr, "ERROR: build_cont_info() nlevels(%d) < 1\n", nlevels);
      return NULL;
   }

   cp = (cont_3d_info *) malloc(sizeof(cont_3d_info));
   if (cp == NULL) {
      fprintf(stderr, "ERROR: out of memory in build_cont_info() at 1\n");
      return NULL;
   }
   cp->nlevels = nlevels;

   if (ndim < 2) {
      fprintf(stderr, "ERROR: build_cont_info() ndim(%d) resetting to 2\n", ndim);
      ndim = 2;
   }
   if (ndim > 3) {
      fprintf(stderr, "ERROR: build_cont_info() ndim(%d) resetting to 3\n", ndim);
      ndim = 3;
   }
   cp->ndim = ndim;     /* number of dimensions to be output */

   level_segs = (level_3d_segs *) malloc(nlevels * sizeof(level_3d_segs));
   if (level_segs == NULL) {
      fprintf(stderr, "ERROR: out of memory in build_cont_info() at 2\n");
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) cont_3d_info\n", cp);
#endif
      free(cp);
      return NULL;
   }
   cp->level = level_segs;

   for(i=0; i < nlevels; i++) {
      level_segs[i].level_num = i+1;
      level_segs[i].iso_val   = level[i];
      prop_str = (props ? props[i] : NULL);
      level_segs[i].property  = strdup(prop_str?prop_str:"");
      level_segs[i].lines     = NULL;
   }

   return cp;
}

/* print all the data in the cont info */
void simple_contour_output(FILE *fp, cont_3d_info *cp) {
   int i = 0, j = 0, ndim = 0;
   line_segment *curr_seg = NULL;
   coordinate_node *curr = NULL;
   char *prop_str = NULL;
   char *linetype = NULL;

   if (cp == NULL || cp->level == NULL) { return; }
   ndim = cp->ndim;

   for(i=0; i < cp->nlevels; i++) { /* loop over each level */
      prop_str = cp->level[i].property;
      fprintf(fp, "#level%d {%g} %s\n",
	    i+1, cp->level[i].iso_val, (prop_str?prop_str:""));

      for(curr_seg = cp->level[i].lines; curr_seg; curr_seg = curr_seg->next) {

	 if (curr_seg->loop) {
	    fprintf(fp, "loop\n");
	 }

	 j = 0;
	 for (curr = curr_seg->head; curr; curr = curr->next) {
	    if (j == 0) {
	       if (curr->next == NULL) { /* point */
		  fprintf(fp, (ndim == 2)?"point %g %g\n":"point %g %g %g\n",
			curr->d[0], curr->d[1], curr->d[2]);
	       }
	       else { /* move */
		  fprintf(fp, (ndim == 2)?"move %g %g\n":"move %g %g %g\n",
			curr->d[0], curr->d[1], curr->d[2]);
	       }
	    }
	    else { /* draw */
		  fprintf(fp, (ndim == 2)?"draw %g %g\n":"draw %g %g %g\n",
			curr->d[0], curr->d[1], curr->d[2]);
	    }
	    j++;
	 }
	 fprintf(fp, "\n");
      }
   }
}

/* PostScript procedure definitions for contour routine */
void simple_2Dps_procset(FILE *fp, int dobang) {
   if (dobang) {
      fprintf(fp, "%%!PS-Adobe-2.0\n");
   }

   fprintf(fp, "/setcmykcolor where { pop } {\n");
   fprintf(fp, "      /setcmykcolor {\n");
   fprintf(fp, "         1 sub 4 1 roll\n");
   fprintf(fp, "         3 {\n");
   fprintf(fp, "           3 index add neg dup 0 lt { pop 0 } if\n");
   fprintf(fp, "           3 1 roll\n");
   fprintf(fp, "         } repeat\n");
   fprintf(fp, "         setrgbcolor pop\n");
   fprintf(fp, "      } bdef\n");
   fprintf(fp, "   } ifelse\n");

   fprintf(fp, "/M /moveto load def\n");
   fprintf(fp, "/L /lineto load def\n");
   fprintf(fp, "/point {2 copy moveto lineto} bind def\n");

   fprintf(fp, "/color_           { 0.0  0.0  0.0  1.0 } def\n");
   fprintf(fp, "/color_black      { 0.0  0.0  0.0  1.0 } def\n");
   fprintf(fp, "/color_red        { 0.15 0.80 0.80 0.0 } def\n");
   fprintf(fp, "/color_green      { 0.80 0.20 0.80 0.0 } def\n");
   fprintf(fp, "/color_blue       { 0.80 0.60 0.00 0.0 } def\n");
   fprintf(fp, "/color_cyan       { 0.80 0.20 0.00 0.0 } def\n");
   fprintf(fp, "/color_yellow     { 0.00 0.10 0.80 0.0 } def\n");
   fprintf(fp, "/color_magenta    { 0.30 0.90 0.00 0.0 } def\n");
   fprintf(fp, "/color_white      { 0.00 0.00 0.00 0.0 } def\n");
   fprintf(fp, "/color_pink       { 0.00 0.50 0.15 0.0 } def\n");
   fprintf(fp, "/color_orange     { 0.00 0.50 0.60 0.0 } def\n");
   fprintf(fp, "/color_purple     { 0.50 0.90 0.00 0.0 } def\n");
   fprintf(fp, "/color_sky        { 0.80 0.33 0.00 0.0 } def\n");
   fprintf(fp, "/color_brown      { 0.15 0.30 0.30 0.0 } def\n");
   fprintf(fp, "/color_gray       { 0.0  0.0  0.0 0.35 } def\n");
   fprintf(fp, "/color_grey       { 0.0  0.0  0.0 0.35 } def\n");
   fprintf(fp, "/color_gold       { 0.00 0.30 0.80 0.0 } def\n");
   fprintf(fp, "/color_yellowtint { 0.00 0.05 0.50 0.0 } def\n");
   fprintf(fp, "/color_sea        { 0.70 0.00 0.50 0.0 } def\n");
   fprintf(fp, "/color_pinktint   { 0.00 0.30 0.10 0.0 } def\n");
   fprintf(fp, "/color_bluetint   { 0.35 0.10 0.00 0.0 } def\n");
   fprintf(fp, "/color_greentint  { 0.50 0.00 0.50 0.0 } def\n");
   fprintf(fp, "/color_hotpink    { 0.10 0.80 0.10 0.0 } def\n");
}

/* print all the data in the cont info in PostScript format */
void simple_2Dps_contour_output(FILE *fp, cont_3d_info *cp, int doevery) {
   int i = 0, j = 0, ndim = 0;
   line_segment *curr_seg = NULL;
   coordinate_node *curr = NULL;
   char *prop_str = NULL;
   char *linetype = NULL;

   if (cp == NULL || cp->level == NULL) { return; }
   ndim = cp->ndim;

   for(i=0; i < cp->nlevels; i++) { /* loop over each level */
      if ((cp->level[i].lines != NULL) || doevery) {
	 prop_str = cp->level[i].property;
	 fprintf(fp, "%%level%d iso_val: %g\ncolor_%s setcmykcolor\n",
	    i+1, cp->level[i].iso_val, (prop_str?prop_str:"grey"));

	 for(curr_seg = cp->level[i].lines; curr_seg;
	     curr_seg = curr_seg->next) {

	    fprintf(fp, "newpath\n");
	    
	    j = 0;
	    curr = curr_seg->head;
	    if (curr_seg->loop && curr && curr->next) {
	       curr = curr->next; /* skip the duplicate endpoints of a loop */
	    }
	    while (curr) {
	       if (j == 0) {
		  if (curr->next == NULL) { /* point */
		     fprintf(fp, "%g %g point\n", curr->d[0], curr->d[1]);
		  }
		  else { /* move */
		     fprintf(fp, "%g %g M\n", curr->d[0], curr->d[1]);
		  }
	       }
	       else { /* draw */
		  fprintf(fp, "%g %g L\n", curr->d[0], curr->d[1]);
	       }
	       j++;
	       curr = curr->next;
	    }

	    if (curr_seg->loop) {
	       fprintf(fp, "closepath stroke\n");
	    }
	    else {
	       fprintf(fp, "stroke\n");
	    }
	 }
      }
   }
}

/* PostScript to end and display the page */
void simple_2Dps_showpage(FILE *fp) {
   fprintf(fp, "showpage\n");
}

/* print all the data in the cont info in kinemage format */
void simple_kin_contour_output(FILE *fp, cont_3d_info *cp,
                               int doevery, int thinline) {
   int i = 0, j = 0, ndim = 0;
   line_segment *curr_seg = NULL;
   coordinate_node *curr = NULL;
   char *prop_str = NULL;
   char *linetype = NULL;
   ctree_val iso_val;

   if (cp == NULL || cp->level == NULL) { return; }
   ndim = cp->ndim;

   for(i=0; i < cp->nlevels; i++) { /* loop over each level */
      if ((cp->level[i].lines != NULL) || doevery) {
	 iso_val = cp->level[i].iso_val;
	 prop_str = cp->level[i].property;
	 fprintf(fp, "@vectorlist {%g} color= %s%s\n",
	    iso_val, (prop_str?prop_str:"grey"),
	    (thinline?" width= 1":""));

	 for(curr_seg = cp->level[i].lines; curr_seg;
	     curr_seg = curr_seg->next) {

	    j = 0;
	    for (curr = curr_seg->head; curr; curr = curr->next) {
	       if (j == 0) {
		  if (curr->next == NULL) { /* point */
		     fprintf(fp, "{level %g}P %0.3f %0.3f %0.3f {\"}L %0.3f %0.3f %0.3f\n",iso_val,
			curr->d[0], curr->d[1], curr->d[2],
			curr->d[0], curr->d[1], curr->d[2]);
		  }
		  else { /* move */
		     fprintf(fp, "{level %g}P %0.3f %0.3f %0.3f\n", iso_val,
			curr->d[0], curr->d[1], curr->d[2]);
		  }
	       }
	       else { /* draw */
		  fprintf(fp, "{\"}L %0.3f %0.3f %0.3f\n",
			curr->d[0], curr->d[1], curr->d[2]);
	       }
	       j++;
	    }
	 }
      }
   }
}

/* print all the data in the cont info in kinemage format with a zscale */
void scaled_kin_contour_output(FILE *fp, cont_3d_info *cp,
	int doevery, int thinline, float zscale) {
   int i = 0, j = 0, ndim = 0;
   line_segment *curr_seg = NULL;
   coordinate_node *curr = NULL;
   char *prop_str = NULL;
   char *linetype = NULL;
   ctree_val iso_val;

   if (cp == NULL || cp->level == NULL) { return; }
   ndim = cp->ndim;

   for(i=0; i < cp->nlevels; i++) { /* loop over each level */
      if ((cp->level[i].lines != NULL) || doevery) {
	 iso_val = cp->level[i].iso_val;
	 prop_str = cp->level[i].property;
	 fprintf(fp, "@vectorlist {%g} color= %s%s\n",
	    iso_val, (prop_str?prop_str:"grey"),
	    (thinline?" width= 1":""));

	 for(curr_seg = cp->level[i].lines; curr_seg;
	     curr_seg = curr_seg->next) {

	    j = 0;
	    for (curr = curr_seg->head; curr; curr = curr->next) {
	       if (j == 0) {
		  if (curr->next == NULL) { /* point */
		     fprintf(fp, "{level %g}P %0.3f %0.3f %0.3f {\"}L %0.3f %0.3f %0.3f\n",iso_val,
			curr->d[0], curr->d[1], curr->d[2]*zscale,
			curr->d[0], curr->d[1], curr->d[2]*zscale);
		  }
		  else { /* move */
		     fprintf(fp, "{level %g}P %0.3f %0.3f %0.3f\n", iso_val,
			curr->d[0], curr->d[1], curr->d[2]*zscale);
		  }
	       }
	       else { /* draw */
		  fprintf(fp, "{\"}L %0.3f %0.3f %0.3f\n",
			curr->d[0], curr->d[1], curr->d[2]*zscale);
	       }
	       j++;
	    }
	 }
      }
   }
}

/* print all the data in the cont info */
void print_cont_info(FILE *fp, cont_3d_info *cp) {
   int i = 0, j = 0;
   line_segment *curr;
   char *prop_str = NULL;

   if (cp == NULL) {
      fprintf(fp, "cont_info NULL\n");
      return;
   }
   fprintf(fp, " cont_info:0x%p\n", cp);
   fprintf(fp, "      ndim: %d\n", cp->ndim);
   fprintf(fp, "   nlevels: %d\n", cp->nlevels);

   if (cp->level == NULL) {
      fprintf(fp, "     level: pointer NULL\n");
      return;
   }

   for(i=0; i < cp->nlevels; i++) {
      fprintf(fp, "     level[%d]\n", i);
      fprintf(fp, "     level_num: %d\n", cp->level[i].level_num);
      fprintf(fp, "       iso_val: %g\n", cp->level[i].iso_val);
      prop_str = cp->level[i].property;
      fprintf(fp, "      property: \"%s\"\n", (prop_str?prop_str:""));
      fprintf(fp, "         lines: 0x%p\n", cp->level[i].lines);
      curr = cp->level[i].lines;
      j = 0;
      while (curr) {
	 fprintf(fp, " ** line segment %d\n", ++j);
	 print_line_segment(fp, curr);
	 curr = curr->next;
      }
   }
}

/* move all the line segments from sp to cp and remap */
int remap_slices(cont_3d_info *cp, slice_info *sp) {
   line_segment  *curr_seg   = NULL;
   line_segment  *next_seg   = NULL;
   int i = 0, nlevels = 0;

   if (cp == NULL || sp == NULL) {
      fprintf(stderr, "ERROR: remap_slices(0x%p, 0x%p) NULL pointer\n", cp, sp);
      return 0;
   }

   slice_all_done(sp); /* consolidate line segments */

   if (cp->nlevels != sp->nlevels) {
      fprintf(stderr, "ERROR: remap_slices() number of levels differ(%d vs %d)\n",
	 cp->nlevels, sp->nlevels);
      return 0;
   }

   nlevels = sp->nlevels;

   for(i=0; i < nlevels; i++) {

      /* unhook open segs for this level */
      curr_seg = sp->level[i].done;
      sp->level[i].done = NULL;

      while (curr_seg) {
	 next_seg = curr_seg->next; /* remember next */

	 /* add curr to lines */
	 curr_seg->next = cp->level[i].lines;
	 cp->level[i].lines = curr_seg;

	 remap_line_segment(curr_seg, sp->level[i].map);

	 curr_seg = next_seg; /* continue */
      }
   }
   return 1;
}

/* alter the order of the dimensions for a line segment */
void remap_line_segment(line_segment* lsp, int map[]) {
   coordinate_node *curr = NULL;
   int i = 0;
   coord_val x, y, z;

   if (lsp != NULL && lsp->head != NULL) {

      for (curr = lsp->head; curr; curr = curr->next) {
	 x = curr->d[0];
	 y = curr->d[1];
	 z = curr->d[2];
	 curr->d[map[0]] = x;
	 curr->d[map[1]] = y;
	 curr->d[map[2]] = z;
      }
   }
}

/* function to kill all the data in the slice info */
void destroy_cont_info(cont_3d_info *cp) {
   line_segment  *curr_seg   = NULL;
   line_segment  *next_seg   = NULL;
   int i = 0, nlevels = 0;

   if (cp != NULL) {
      nlevels = cp->nlevels;

      for(i=0; i < nlevels; i++) {

	 if (cp->level[i].property != NULL) {
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(\"%s\") level property\n", cp->level[i].property);
#endif
	       free(cp->level[i].property);
	 }

	 curr_seg = cp->level[i].lines;
	 while (curr_seg) {
	    next_seg = curr_seg->next;
	    destroy_line_seg(curr_seg);
	    curr_seg = next_seg;
	 }
	 cp->level[i].lines = NULL;
      }
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) cont level[]\n", cp->level);
#endif
      free(cp->level);
      cp->level = NULL;
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) cont_3d_info\n", cp);
#endif
      free(cp);
   }
}

slice_info*
build_slice_info(ctree_lvl level[], int nlevels, int map[],
		  coord_val dx, coord_val dy,
		  int set_last_val, ctree_lvl last_val, int drop1pts) {
   slice_info *sp = NULL;
   level_2d_segs *level_segs = NULL;
   int i = 0;

   if (nlevels < 1) {
      fprintf(stderr, "ERROR: build_slice_info() nlevels(%d) < 1\n", nlevels);
      return NULL;
   }

   sp = (slice_info *) malloc(sizeof(slice_info));
   if (sp == NULL) {
      fprintf(stderr, "ERROR: out of memory in build_slice_info() at 1\n");
      return NULL;
   }
   sp->nlevels = nlevels;
      /* conservatively, we slightly expand the row delta cuttoff value */
   sp->row_delta = 1.1           * ((dy < 0.0) ? -dy : dy);
   sp->epsilon0  = SCALE_EPSILON * ((dx < 0.0) ? -dx : dx);
   sp->epsilon1  = SCALE_EPSILON * ((dy < 0.0) ? -dy : dy);
   sp->drop1pts  = drop1pts;

   level_segs = (level_2d_segs *) malloc(nlevels * sizeof(level_2d_segs));
   if (level_segs == NULL) {
      fprintf(stderr, "ERROR: out of memory in build_slice_info() at 2\n");
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) slice_info\n", sp);
#endif
      free(sp);
      return NULL;
   }
   sp->level = level_segs;

   for(i=0; i < nlevels; i++) {
      level_segs[i].level_num = i+1;
      level_segs[i].iso_val   = level[i];
      if (set_last_val) {
	 level_segs[i].coord_last_val = last_val; /* fixed third coordinate */
      }
      else {
	 level_segs[i].coord_last_val = level[i]; /* third coordinate == level */
      }
      level_segs[i].map[0]    = map[0];
      level_segs[i].map[1]    = map[1];
      level_segs[i].map[2]    = map[2];
      level_segs[i].open      = NULL;
      level_segs[i].done      = NULL;
   }

   return sp;
}

/* modify the coordinate system data for the slice information */
void
update_slice_info(slice_info *sp, int map[],
		  coord_val dx, coord_val dy,
		  int set_last_val, ctree_lvl last_val) {
   level_2d_segs *level_segs = NULL;
   int i = 0, nlevels = 0;

   if (sp == NULL) {
      fprintf(stderr, "ERROR: update_slice_info(0x%p) NULL slice\n", sp);
      return;
   }
   nlevels = sp->nlevels;

      /* conservatively, we slightly expand the row delta cuttoff value */
   sp->row_delta = 1.1           * ((dy < 0.0) ? -dy : dy);
   sp->epsilon0  = SCALE_EPSILON * ((dx < 0.0) ? -dx : dx);
   sp->epsilon1  = SCALE_EPSILON * ((dy < 0.0) ? -dy : dy);

   level_segs = sp->level;

   for(i=0; i < nlevels; i++) {
      if (set_last_val) {      /* fixed third coordinate */
	 level_segs[i].coord_last_val = last_val;
      }
      else {                   /* third coordinate == level */
	 level_segs[i].coord_last_val = level_segs[i].iso_val;
      }
      level_segs[i].map[0]    = map[0];
      level_segs[i].map[1]    = map[1];
      level_segs[i].map[2]    = map[2];
#ifdef WARN_PROBLEMS
      if (level_segs[i].open != NULL) {
	 fprintf(stderr, "WARNING: update_slice_info(0x%p) [%d].open(0x%p) list not empty\n",
	    sp, i, level_segs[i].open);
      }
      if (level_segs[i].done != NULL) {
	 fprintf(stderr, "WARNING: update_slice_info(0x%p) [%d].done(0x%p) list not empty\n",
	    sp, i, level_segs[i].done);
      }
#endif
   }
}

/* print all the data in the slice info */
void print_slice_info(FILE *fp, slice_info *sp) {
   int i = 0, j = 0;
   line_segment *curr;

   if (sp == NULL) {
      fprintf(fp, "slice_info NULL\n");
      return;
   }
   fprintf(fp, "slice_info:0x%p\n", sp);
   fprintf(fp, "   nlevels: %d\n", sp->nlevels);
   fprintf(fp, " row_delta: %g\n", sp->row_delta);
   fprintf(fp, "  epsilon0: %g\n", sp->epsilon0);
   fprintf(fp, "  epsilon1: %g\n", sp->epsilon1);

   if (sp->level == NULL) {
      fprintf(fp, "     level: pointer NULL\n");
      return;
   }

   for(i=0; i < sp->nlevels; i++) {
      fprintf(fp, "     level[%d]\n", i);
      fprintf(fp, "     level_num: %d\n", sp->level[i].level_num);
      fprintf(fp, "       iso_val: %g\n", sp->level[i].iso_val);
      fprintf(fp, "coord_last_val: %g\n", sp->level[i].coord_last_val);
      fprintf(fp, "           map: [%d, %d, %d]\n",
	 sp->level[i].map[0], sp->level[i].map[1], sp->level[i].map[2]);
      fprintf(fp, "          open: 0x%p\n", sp->level[i].open);
      curr = sp->level[i].open;
      j = 0;
      while (curr) {
	 fprintf(fp, " ** line segment %d\n", ++j);
	 print_line_segment(fp, curr);
	 curr = curr->next;
      }
      fprintf(fp, "          done: 0x%p\n", sp->level[i].done);
      curr = sp->level[i].done;
      j = 0;
      while (curr) {
	 fprintf(fp, " ** line segment %d\n", ++j);
	 print_line_segment(fp, curr);
	 curr = curr->next;
      }
   }
}

/* move all the open segments to the done list */
void slice_all_done(slice_info *sp) {
   line_segment  *curr_seg   = NULL;
   line_segment  *next_seg   = NULL;
   int i = 0, nlevels = 0;

   if (sp != NULL) {
      nlevels = sp->nlevels;

      for(i=0; i < nlevels; i++) {

	 /* unhook open segs for this level */
	 curr_seg = sp->level[i].open;
	 sp->level[i].open = NULL;

	 while (curr_seg) {
	    next_seg = curr_seg->next; /* remember next */

	     /* add curr to done */
	    curr_seg->next = sp->level[i].done;
	    sp->level[i].done = curr_seg;

	    curr_seg = next_seg; /* continue */
	 }
      }
   }
}

/* function to kill all the data in the slice info */
void destroy_slice_info(slice_info *sp) {
   line_segment  *curr_seg   = NULL;
   line_segment  *next_seg   = NULL;
   int i = 0, nlevels = 0;

   if (sp != NULL) {
      nlevels = sp->nlevels;

      for(i=0; i < nlevels; i++) {

	 curr_seg = sp->level[i].open;
	 while (curr_seg) {
	    next_seg = curr_seg->next;
	    destroy_line_seg(curr_seg);
	    curr_seg = next_seg;
	 }
	 sp->level[i].open = NULL;

	 curr_seg = sp->level[i].done;
	 while (curr_seg) {
	    next_seg = curr_seg->next;
	    destroy_line_seg(curr_seg);
	    curr_seg = next_seg;
	 }
	 sp->level[i].done = NULL;
      }
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) slice level[]\n", sp->level);
#endif
      free(sp->level);
      sp->level = NULL;
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) slice_info\n", sp);
#endif
      free(sp);
   }
}

/* delete all the components of a line segment */
void destroy_line_seg(line_segment* lsp) {

   if (lsp != NULL) {
      if (lsp->head != NULL) {
	 coordinate_node *curr_coord = NULL;
	 coordinate_node *next_coord = NULL;

	 curr_coord = lsp->head;

	 while (curr_coord) {
	    next_coord = curr_coord->next;

	    curr_coord->next = NULL;
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) coord => (%g, %g, %g)\n",
	    curr_coord, curr_coord->d[0], curr_coord->d[1], curr_coord->d[2]);
#endif
	    free(curr_coord); /* delete the coordinate */
	    curr_coord = next_coord;
	 }
      }
      lsp->head = NULL;
      lsp->tail = NULL;
      lsp->next = NULL;
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) line_segment\n", lsp);
#endif
      free(lsp); /* delete the line segment header */
   }
}

/* invert the direction of a line segment */
void reverse_line_segment(line_segment* lsp) {
   coordinate_node *cdr     = NULL;
   coordinate_node *oldlist = NULL;

   if (lsp != NULL && lsp->head != NULL) {

      oldlist = lsp->head->next;

      lsp->head->next = NULL;
      lsp->tail = lsp->head;

      while (oldlist) {
	 cdr           = oldlist->next;
	 oldlist->next = lsp->head;
	 lsp->head     = oldlist;
	 oldlist       = cdr;
      }
   }
}

/* combine two line segments in a specified orientation, returning only one segment */
line_segment*
join_segments(line_segment* ls1, line_segment* ls2, int selection1, int selection2) {

   if (ls1 == NULL) {
      fprintf(stderr, "ERROR: join_segments() list 1 NULL, list 2 = 0x%p\n", ls2);
      return ls2;
   }
   if (ls2 == NULL) {
      fprintf(stderr, "ERROR: join_segments() list 1 = 0x%p, list 2 NULL\n", ls1);
      return ls1;
   }
   if (ls1 == ls2) {
      fprintf(stderr, "ERROR: join_segments() duplicate lists = 0x%p\n", ls2);
      return ls1;
   }

   if ((selection1 == SELECT_TAIL) && (selection2 == SELECT_HEAD)) {
      /* the order is the way we prefer: tail of 1 -> head of 2 */
      /* *** we do not do anything here *** */
   }
   else if ((selection1 == SELECT_TAIL) && (selection2 == SELECT_TAIL)) {
      /* make the tail of 2 into the head of 2 */
      reverse_line_segment(ls2);
   }
   else if ((selection1 == SELECT_HEAD) && (selection2 == SELECT_HEAD)) {
      /* make the head of 1 into the tail of 1 */
      reverse_line_segment(ls1);
   }
   else if ((selection1 == SELECT_HEAD) && (selection2 == SELECT_TAIL)) {
      coordinate_node *hd1 = NULL, *tl1 = NULL;

      /* swap the two lists */
      hd1 = ls1->head;
      tl1 = ls1->tail;
      ls1->head = ls2->head;
      ls1->tail = ls2->tail;
      ls2->head = hd1;
      ls2->tail = tl1;
   }
   else {
      fprintf(stderr, "ERROR: join_segments(..., %d, %d) unknown selection type\n",
			   selection1, selection2);
      fprintf(stderr, "       NOTE: trying to join TAIL to HEAD...\n");
      
      /* fall through */
   }

   if (ls1->tail == NULL) {
      fprintf(stderr, "ERROR: join_segments() list 1 coords missing => 0x%p\n", ls1->head);
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) line_segment\n", ls1);
#endif
      free(ls1); /* kill header of list 1 */
      return ls2;
   }
   else if (ls2->head == NULL) {
      fprintf(stderr, "ERROR: join_segments() list 2 coords missing => 0x%p\n", ls2->tail);
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) line_segment\n", ls2);
#endif
      free(ls2); /* kill header of list 2 */
      return ls1;
   }
   else { /* tail and head are ok */

      ls1->tail->next = ls2->head->next; /* join the lists */

      if (ls1->tail->next == NULL) { /* if list 2 is only one coordinate long */
	 ls2->tail = ls1->tail;
      }
      ls2->head->next = NULL;
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) coordinate\n", ls2->head);
#endif
      free(ls2->head); /* skip the head of list 2 (it should be a duplicate) */

      ls1->tail = ls2->tail; /* set the new tail */

      ls2->head = NULL;
      ls2->tail = NULL;
#ifdef DEBUG_FREE
      fprintf(stderr, "DEBUG: free(0x%p) line_segment\n", ls2);
#endif
      free(ls2); /* kill the old list 2 */

      ls1->next = NULL; /* make extra sure this segment is unconnected at this point */
   }
   return ls1;
}

#define WITHIN_EPSILON(av, bv, eps) (((av)>(bv)) ? \
                        (((av)-(bv))<=(eps)) : (((bv)-(av))<=(eps)))

/* do the ends of the two line segments match? */
/* if so, which ends (selections are updated)  */
int line_seg_ends_match(line_segment* ls1, line_segment* ls2,
		  int *selection1, int *selection2,
		  coord_val epsilon0, coord_val epsilon1) {
   int rc = 0;

   if (ls1 == NULL || ls2 == NULL) {
      fprintf(stderr, "ERROR: line_seg_ends_match(0x%p, 0x%p) NULL segments\n", ls1, ls2);
      rc=0;
   }
   else if (ls1 == ls2) {
      fprintf(stderr, "ERROR: line_seg_ends_match() duplicate lists = 0x%p\n", ls2);
      rc=0;
   }
   else if (ls1->head == NULL || ls1->tail == NULL
    || ls2->head == NULL || ls2->tail == NULL) {
      fprintf(stderr, "ERROR: line_seg_ends_match([0x%p, 0x%p], [0x%p, 0x%p]) NULL head/tail\n",
      ls1->head, ls1->tail, ls2->head, ls2->tail);
      rc=0;
   }
   else if (WITHIN_EPSILON(ls1->head->d[0], ls2->head->d[0], epsilon0)
        &&  WITHIN_EPSILON(ls1->head->d[1], ls2->head->d[1], epsilon1)) {
      *selection1 = SELECT_HEAD;
      *selection2 = SELECT_HEAD;
      rc=1;
   }
   else if (WITHIN_EPSILON(ls1->head->d[0], ls2->tail->d[0], epsilon0)
        &&  WITHIN_EPSILON(ls1->head->d[1], ls2->tail->d[1], epsilon1)) {
      *selection1 = SELECT_HEAD;
      *selection2 = SELECT_TAIL;
      rc=1;
   }
   else if (WITHIN_EPSILON(ls1->tail->d[0], ls2->head->d[0], epsilon0)
        &&  WITHIN_EPSILON(ls1->tail->d[1], ls2->head->d[1], epsilon1)) {
      *selection1 = SELECT_TAIL;
      *selection2 = SELECT_HEAD;
      rc=1;
   }
   else if (WITHIN_EPSILON(ls1->tail->d[0], ls2->tail->d[0], epsilon0)
        &&  WITHIN_EPSILON(ls1->tail->d[1], ls2->tail->d[1], epsilon1)) {
      *selection1 = SELECT_TAIL;
      *selection2 = SELECT_TAIL;
      rc=1;
   }
   else { rc = 0; }  /* no match */

#ifdef DEBUG_LINK
if (rc) {
fprintf(stderr, "DEBUG: connect dim 0 (%g - %g) = %g, e = %g\n",
   ls1->head->d[0], ls2->tail->d[0],
   ls1->head->d[0]- ls2->tail->d[0], epsilon0);
fprintf(stderr, "DEBUG: connect dim 1 (%g - %g) = %g, e = %g\n\n",
   ls1->head->d[1], ls2->tail->d[1],
   ls1->head->d[1]- ls2->tail->d[1], epsilon1);
}
#endif

   return rc;
}

/* Is the list circular? If so set the loop flag */
int check_for_loop(line_segment* lsp, coord_val epsilon0, coord_val epsilon1) {
   int rc = 0;

   if (lsp == NULL) {
      fprintf(stderr, "ERROR: check_for_loop() NULL pointer\n");
      rc=0;
   }
   else if (lsp->head == NULL || lsp->tail == NULL) {
      fprintf(stderr, "ERROR: check_for_loop([0x%p, 0x%p]) NULL head/tail\n",
	 lsp->head, lsp->tail);
      rc=0;
   }
   else if (lsp->head->next == NULL) {
      rc=0; /* single points do not loop */
   }
   else if (WITHIN_EPSILON(lsp->head->d[0], lsp->tail->d[0], epsilon0)
        &&  WITHIN_EPSILON(lsp->head->d[1], lsp->tail->d[1], epsilon1)) {

#ifdef DEBUG_LINK
fprintf(stderr, "DEBUG: loop dim 0 (%g - %g) = %g, e = %g\n",
   lsp->head->d[0], lsp->tail->d[0],
   lsp->head->d[0]- lsp->tail->d[0], epsilon0);
fprintf(stderr, "DEBUG: loop dim 1 (%g - %g) = %g, e = %g\n\n",
   lsp->head->d[1], lsp->tail->d[1],
   lsp->head->d[1]- lsp->tail->d[1], epsilon1);
#endif

      rc=1;
      lsp->loop = 1; /* line segment makes a loop */
   }
   else {
      rc = 0;
      lsp->loop = 0; /* no match, must be linear */
   }

   return rc;
}

/* We can optimize the search by putting line segments which can not     */
/* be reached into the done list. Relies on the new segs end row values. */
int line_seg_too_far_back(line_segment* lsp,
	       coord_val r1, coord_val r2, coord_val cuttoff) {
   coord_val headrow, tailrow;

   if (lsp != NULL && lsp->head != NULL) {

      headrow = lsp->head->d[ROW_DIM_NUM];
      if (lsp->tail != NULL) {
	 tailrow = lsp->tail->d[ROW_DIM_NUM];
      }
      else {  /* Should be impossible. Must be corrupted line_segment.*/
	 fprintf(stderr, "ERROR: line_seg_too_far_back(0x%p) NULL tail\n", lsp);
	 return 0;
      }

      if (   (! WITHIN_EPSILON(r1, headrow, cuttoff))
          && (! WITHIN_EPSILON(r1, tailrow, cuttoff))
          && (! WITHIN_EPSILON(r2, headrow, cuttoff))
          && (! WITHIN_EPSILON(r2, tailrow, cuttoff))) {
			   /* out of range -- too many rows back */
	 return 1;
      }

   }
   return 0;
}

/* store a line segment into the slice structure at a given level */

int insert_one_line_segment(int level, slice_info* sp, line_segment* lsp) {
   level_2d_segs *level_segs = NULL;

   if (lsp == NULL) {
      fprintf(stderr, "ERROR: insert_one_line_segment() line segment NULL\n");
      return 0;
   }
   else if (lsp->head == NULL) {
      fprintf(stderr, "ERROR: insert_one_line_segment(0x%p) NULL head\n", lsp);
      return 0;
   }

   lsp->next = NULL; /* initialize link */

   level_segs = &(sp->level[level]);
   
   if (level_segs->open == NULL) {
      level_segs->open = lsp; /* add as the head of the empty open list */
   }
   else {
      line_segment *prev_ls = NULL, *curr_ls = NULL;
      int isamatch, isaloop, toofarback, sel1, sel2;
      coord_val headrow, tailrow;

      /* Remember the row positions of the new segment to identify stale line segs. */
      /* Assumes that line segments are inserted in a sweep across rows             */
      /* so we can determine if an old line seg is too far back to connect.         */
      headrow = lsp->head->d[ROW_DIM_NUM];
      if (lsp->tail != NULL) {
	 tailrow = lsp->tail->d[ROW_DIM_NUM];
      }
      else {  /* Should be impossible. Must be corrupted line_segment.*/
	 fprintf(stderr, "ERROR: insert_one_line_segment(0x%p) NULL tail\n", lsp);
	 return 0;
      }
      
      /* must loop over the list of open segments */

      curr_ls = level_segs->open;
      while (curr_ls != NULL) {

	 /* look for connecting segments */

	 isamatch = line_seg_ends_match(lsp, curr_ls,
			   &sel1, &sel2, sp->epsilon0, sp->epsilon1);
	 if (isamatch) {

	    if (prev_ls == NULL) { /* unlink matching segment */
	       level_segs->open = curr_ls->next;
	    }
	    else {
	       prev_ls->next = curr_ls->next;
	    }
	    curr_ls->next = NULL;
	    
	    lsp = join_segments(lsp, curr_ls, sel1, sel2);
	    if (lsp == NULL) { /* problem while storing segment */
	       return 0;       /* signal we have a problem      */
	    }

	    isaloop = check_for_loop(lsp, sp->epsilon0, sp->epsilon1);
	    if (isaloop) { /* loops are placed on the done list */
	       lsp->next = level_segs->done;
	       level_segs->done = lsp;

	       curr_ls = NULL; /* stop combing through the list  */
	       lsp     = NULL; /* signal we have stored the list */
	    }
	    else { /* linear */
	       if (prev_ls == NULL) { /* continue scanning */
		  curr_ls = level_segs->open;
	       }
	       else {
		  curr_ls = prev_ls->next;
	       }
	    }

	 }
	 else {	  /* line segment does NOT connect with current seg */

	    /* is the current segment too many rows back? */

	    toofarback = line_seg_too_far_back(curr_ls, headrow, tailrow, sp->row_delta);
	    if (toofarback) {
	       if (prev_ls == NULL) { /* unlink stale segment */
		  level_segs->open = curr_ls->next;
	       }
	       else {
		  prev_ls->next = curr_ls->next;
	       }
	       curr_ls->next = NULL;

	       curr_ls->next = level_segs->done; /* add to the done list */
	       level_segs->done = curr_ls;

	       if (prev_ls == NULL) { /* continue combing through the list */
		  curr_ls = level_segs->open;
	       }
	       else {
		  curr_ls = prev_ls->next;
	       }
	    }
	    else { /* line end is sufficiently nearby */
		  prev_ls = curr_ls;
		  curr_ls = curr_ls->next;
	    }
	 } /* end -- no match */

      } /* end while curr_ls */

      if (lsp != NULL) { /* line seg has not been put on done list */
	 lsp->next = level_segs->open;/* must add to the open list */
	 level_segs->open = lsp;
      }    
   }

   return 1; /* signal success */
}

/* write out diagnostic information about a line segment */
void print_line_segment(FILE *fp, line_segment* lsp) {
   coordinate_node *curr = NULL;
   int i = 0;

   if (lsp == NULL) {
      fprintf(fp, "NULL line segment pointer\n");
      return;
   }
   else {
      fprintf(fp, "segment:0x%p %s next(0x%p)\n",
	 lsp, ((lsp->loop)?"loop":"linear"), lsp->next);
   }
   if (lsp->head == NULL) {
      fprintf(fp, "  NULL line segment HEAD\n");
      return;
   }

   for (curr = lsp->head; curr; curr = curr->next) {
      i++;
      fprintf(fp, "   %3d: 0x%p (%g, %g, %g)\n", i,
	    curr, curr->d[0], curr->d[1], curr->d[2]);
   }

   if (lsp->tail == NULL) {
      fprintf(fp, "  NULL line segment TAIL\n");
   }
   else {
      fprintf(fp, "  TAIL: 0x%p (%g, %g, %g)\n", lsp->tail,
	    lsp->tail->d[0], lsp->tail->d[1], lsp->tail->d[2]);
   }
}

/* simple mapping from raw data to a row of data */
int
basic_data_func(ctree_data *data, int r, int cf, int cl,
                                                ctree_val *outvec) {
   float *fp = (float *)(data->rawdata);
   int i = 0, cnt = 0;

   if (fp == NULL) { return 0; }

   for (i = cf; i <= cl; i++) {
      outvec[cnt++] = (ctree_val)(*(fp + r*(data->cols) + i));
   }

   return cnt;
}

/* currently, we are not displaying single points */
void
connecting_plot_func(ctree_disp_list dl[],
		     int ndl, int style, void *user_ptr) {
   int i = 0, ncoords = 0;
   slice_info      *sp = (slice_info *)(user_ptr);
   line_segment    *lsp = NULL;
   coordinate_node *cnp = NULL;
   coord_val xprev = -9999.137, yprev = -9999.314; /* dummy values */

   if (ndl < 1)   { return; } /* bail out if list is empty */

   if (ndl == 1 && sp->drop1pts) { return; } /* single point not wanted */

   if (sp == NULL) {
      fprintf(stderr, "ERROR: connecting_plot_func() slice info missing\n");
      return;
   }
   
   if (sp->nlevels <= style) {
      fprintf(stderr, "ERROR: connecting_plot_func() style(%d) >= nlevels(%d)\n",
	       style, sp->nlevels);
      return;
   }

   for (i = 0; i < ndl; i++) {
      switch (dl[i].c) {
      case CT_POINT: /* single points are like move-tos */
	 if (sp->drop1pts) { break; }
	 else {
	       /* fall through (no break) */
	 }
      case CT_MOVE:
	 if (lsp != NULL) {
	    if (ncoords == 1 && sp->drop1pts) {
	       destroy_line_seg(lsp);
	    }
	    else {
	       if (! insert_one_line_segment(style, sp, lsp)) {
		  /* problem trying to insert */
		  return;     /* do not go on */
	       }
	    }
	    lsp = NULL;
	    ncoords = 0;
	 }

	 lsp = (line_segment *) malloc(sizeof(line_segment));
	 if (lsp == NULL) {
	    fprintf(stderr, "ERROR: out of memory in connecting_plot_func() at 1\n");
	    return;
	 }
	 lsp->next = NULL; /* initialize */
	 lsp->head = NULL;
	 lsp->tail = NULL;
	 lsp->loop = 0;

	 cnp = (coordinate_node *) malloc(sizeof(coordinate_node));
	 if (cnp == NULL) {
	    fprintf(stderr, "ERROR: out of memory in connecting_plot_func() at 2\n");
	    return;
	 }
	 cnp->d[0] = dl[i].x; /* store coordinates */
	 cnp->d[1] = dl[i].y;
	 cnp->d[2] = sp->level[style].coord_last_val;
	 
	 cnp->next = NULL;
	 lsp->head = cnp;  /* setup this coord as new list first item */
	 lsp->tail = cnp;

	 ncoords = 1;
	 xprev = dl[i].x; /* remember we were here */
	 yprev = dl[i].y;
	 break;
      case CT_DRAW:
	 if (lsp == NULL) {
	    fprintf(stderr,
	       "ERROR: DRAW with no line seg struct in connecting_plot_func()\n");
	    return;
	 }

	 /* make sure it is not a duplicate point */
	 if ((! WITHIN_EPSILON(xprev, dl[i].x, sp->epsilon0))
	  || (! WITHIN_EPSILON(yprev, dl[i].y, sp->epsilon1)) ) {

	    cnp = (coordinate_node *) malloc(sizeof(coordinate_node));
	    if (cnp == NULL) {
	       fprintf(stderr, "ERROR: out of memory in connecting_plot_func() at 2\n");
	       return;
	    }
	    cnp->d[0] = dl[i].x; /* store coordinates */
	    cnp->d[1] = dl[i].y;
	    cnp->d[2] = sp->level[style].coord_last_val;

	    cnp->next = lsp->head;
	    lsp->head = cnp;    /* add to the head of the list */

	    ncoords++;
	    xprev = dl[i].x; /* remember we were here */
	    yprev = dl[i].y;
	 }
	 break;
      default:	break;
      }
   }
   if (lsp != NULL) { /* insert any outstanding segment */
      if (ncoords == 1 && sp->drop1pts) {
	 destroy_line_seg(lsp);
      }
      else {
	 insert_one_line_segment(style, sp, lsp);
      }
      lsp = NULL;
      ncoords = 0;
   }
}
