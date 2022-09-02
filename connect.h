/*			   connect.h				*/
/* Author: J. Michael Word              Date Written: 1/19/99	*/
/*                                                              */
/* Purpose: definitions for connecting contour segments		*/
/*                                                              */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef CONNECT_H
#define CONNECT_H 1

#include "contree.h"
#include "connect.h"

#define SCALE_EPSILON 0.00001

/* absolute array position of each axis */
#define   X_DIM_NUM 0
#define   Y_DIM_NUM 1
#define   Z_DIM_NUM 2
#define ROW_DIM_NUM 1

/* end selected for linking segment */
#define SELECT_HEAD 0
#define SELECT_TAIL 1

/* node in a linked list of coordinates */
/* used in both 2d and 3d point lists   */
typedef struct coordinate_node_t {
   coord_val d[3];                  /* coordinates         */
   struct coordinate_node_t *next;  /* linked list pointer */
} coordinate_node;

/* connected line segment info can be  */
/* linked into a list of line segments */
typedef struct line_seg_t {
   coordinate_node   *head;  /* head of line seg */
   coordinate_node   *tail;  /*  end of line seg */
   struct line_seg_t *next;  /* linked list ptr  */
   int                loop;  /* is this line seg circular? */
} line_segment;

/* line segment information for a given level for once 2d slice */
typedef struct {
   int        level_num;      /* contour level number */
   ctree_val  iso_val;        /* value for contours in this level   */
   coord_val  coord_last_val; /* coordinate d[2] gets this value    */
   /* (usually either the iso-value or the 3rd dimension value) */
   int        map[3];         /* dimension index of each coordinate */
   line_segment *open;        /* incomplete segments   */
   line_segment *done;        /* segs we are done with */
} level_2d_segs;

/* line segment information for an entire 2d slice */
typedef struct {
   coord_val row_delta;          /* helps work out when line seg is done */
   coord_val epsilon0, epsilon1; /* how close is close enough to be "=="?*/
   int       drop1pts;           /* do we ignore single points?          */
   int       nlevels;            /* how many levels are being maintained */
   level_2d_segs *level;         /* array of line segment lists          */
} slice_info;

/* The 3d info structures are used to store the transformed,   */
/* re-mapped coordinates after all the segments are connected. */

/* line segment information for a given level for once 2d slice */
typedef struct {
   int        level_num;    /* contour level number */
   ctree_val  iso_val;      /* value for contours in this level */
   char      *property;     /* this usually holds the color name */
   line_segment *lines; /* segs we are done with */
} level_3d_segs;

/* 3d contour line segment package    */
/* The segment data in slice_info     */
/* gets transformed into this format  */
/* and the dimensions get re-mapped.  */
/* This is used even if the output is */
/* just a 2d map.                     */
typedef struct {
   int ndim;              /* number of dimensions being output    */
   int nlevels;           /* how many levels are being maintained */
   level_3d_segs *level;  /* array of line segment lists          */
} cont_3d_info;

/* ------------------------- prototypes ------------------------- */

void simple_contour_output(FILE *fp, cont_3d_info *cp);
void simple_2Dps_procset(FILE *fp, int dobang);
void simple_2Dps_contour_output(FILE *fp, cont_3d_info *cp, int doevery);
void simple_2Dps_showpage(FILE *fp);
void simple_kin_contour_output(FILE *fp, cont_3d_info *cp,
                               int doevery, int thinline);
void scaled_kin_contour_output(FILE *fp, cont_3d_info *cp,
                               int doevery, int thinline, float zscale);

cont_3d_info*
build_cont_info(ctree_lvl level[], int nlevels, char* props[], int ndim);

void destroy_cont_info(cont_3d_info *cp);

slice_info*
build_slice_info(ctree_lvl level[], int nlevels, int map[],
		  coord_val dx, coord_val dy,
		  int set_last_val, ctree_lvl last_val, int drop1pts);

void
update_slice_info(slice_info *sp, int map[],
		  coord_val dx, coord_val dy,
		  int set_last_val, ctree_lvl last_val);

void slice_all_done(slice_info *sp);

void destroy_slice_info(slice_info *sp);

/* ------------------------- prototypes ------------------------- */

int basic_data_func(ctree_data *data, int r, int cf, int cl,
                                                ctree_val *outvec);

void connecting_plot_func(ctree_disp_list dl[],
		     int ndl, int style, void *user_ptr);
#endif
