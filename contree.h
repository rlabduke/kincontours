/*			   contree.h				*/
/* Author: J. Michael Word              Date Written: 12/3/91	*/
/*                                                              */
/* Purpose: definitions for contour tree routines		*/
/* Ref: A Decomposable Algorithm for Contour Surface Display	*/
/*      Generation, Michael J. Zyda, ACM Transactions on	*/
/*	Graphics, Vol. 7, No. 2, April 1988, Pages 129-148.	*/

#ifndef CONTREE_H
#define CONTREE_H 1

/* three types of drawing primatives generated */
#define CT_MOVE  'm'
#define CT_DRAW  'd'
#define CT_POINT 'p'
typedef char ctree_command;

/* data values and levels are real numbers */
typedef float ctree_val;
typedef ctree_val ctree_lvl;

/* coordinates are real numbers */
typedef float coord_val;

/* display list is an array of coordinates and drawing primatives */
typedef struct {
        coord_val x, y;         /* coordinates        */
        ctree_command c;        /* plotting command   */
} ctree_disp_list;
         /* can't generate more than this many for a single cell */
#define DLIST_MAX 10

/* node in a contour tree data structure */
typedef struct {
        short n;                        /* node number     */
        ctree_command comm;             /* drawing command */
        short d[3];                     /* descendents     */
} ctree_node;

/* data value and coordinate for corner or center of a cell */
typedef struct {
        ctree_val v;            	/* grid value   */
        coord_val x, y;                 /* coordinates  */
} ctree_locval;

/* a single 2x2 cell is what we focus our attention on to contour */
typedef struct {                 /* One contour cell              */
        ctree_locval node[6];    /* Ignore zeroth element.        */
        ctree_val minval, maxval;/* to check against contour lvls */
        unsigned char index;     /* tells us which contour type   */
}ctree_cell;

/* cell nodes numbered:

	4 ------------- 3

	|       5 (avg) |

	1 ------------- 2
*/

/* a cell can have this many trees */
#define CT_DIMAX 2

/* data used to organize the process of contouring a rect region */
typedef struct ctrd {
        void *rawdata;          /* Raw data pointer             */
        int cols, rows;         /* size of entire plane         */
        int c_offset, r_offset; /* region to be plotted         */
        int rgn_c_last;
	int rgn_r_last;         /* region column and row end pt */
        coord_val x, y;         /* coordinates of region tl     */
        coord_val dx, dy;       /* dist between grid points     */
	ctree_val *work1,*work2;/* work vectors, 'col rgn' long */
	void *user_data;        /* supplimental data for plots	*/
        struct ctrd *next;      /* link to other datasets       */
} ctree_data;

#define VALID_REGION(dat) ((((dat)->rgn_c_last) < (dat)->cols)\
                        && (((dat)->rgn_r_last) < (dat)->rows))

                                 /* user supplied function types */
typedef void (*ctree_plot_func)(ctree_disp_list *, int, int,void *);
typedef int (*ctree_data_func)(ctree_data *, int, int, int,
                                                     ctree_val *);

                                          /* other functions */
int ctree_compose_rawdata(ctree_data *, void *,
      coord_val, coord_val, coord_val, coord_val,
      int, int, int, int, int, int,
      ctree_val *, ctree_val *, void *);

int ctree_plot_region(ctree_data *, ctree_lvl *level, int,
      ctree_data_func, ctree_plot_func, void *);

void ctree_make(void);
void ctree_printree(FILE *outf, char *bases_name,
      char *leaves_name);

#endif
