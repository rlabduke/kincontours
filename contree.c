/*			    contree.c				*/
/* Author: J. Michael Word		Date Written: 11/19/91	*/
/*								*/
/* Purpose: Generate countors from a 2D array of data values	*/
/*								*/
/* Ref: A Decomposable Algorithm for Contour Surface Display	*/
/*      Generation, Michael J. Zyda, ACM Transactions on	*/
/*	Graphics, Vol. 7, No. 2, April 1988, Pages 129-148.	*/

#include <stdio.h>
#include "contree.h"
#include "ctdata.h"

unsigned char ctree_compose_cell(ctree_cell *,
   coord_val, coord_val, coord_val, coord_val,
   ctree_lvl, ctree_lvl);
void ctree_visit(int, int, ctree_cell *, ctree_lvl,
   ctree_disp_list*, int*);
void ctree_visit_subtree(int, int, ctree_cell *, ctree_lvl,
   ctree_disp_list*, int*);

#define EPSILON (0.000000001)

/* ctree_compose_rawdata() - Setup rawdata pointers and parms    */
/*                           used in plotting a rect region.     */
/* returns (true/false) whether the region described makes sense */
/*   data         - output data structure being built            */
/*   raw          - pointer to raw data (input to data fetch)    */
/*   x0, y0       - user unit beginning position (lower left)    */
/*   xspan, yspan - size in user units                           */
/*   c, r         - number of cols and rows in the data          */
/*   co, ro       - starting col and row                         */
/*   crgn, rrgn   - size of region to be contoured (num pt spans)*/
/*                  (if zero the default is the rest of data)    */
/*   work1, work2 - work vectors (at least crgn long)            */
/*   user         - extra data to data fetch routine             */
/* NOTE: the data will be plotted in rows but what is a row      */
/*       and what is a col is somewhat a matter of opinion and   */
/*       is specified here.                                      */
int
ctree_compose_rawdata(ctree_data *data, void *raw,
      coord_val x0, coord_val y0, coord_val xspan, coord_val yspan,
      int c, int r, int co, int ro, int crgn, int rrgn,
      ctree_val *work1, ctree_val *work2,
      void *user) {

   data->rawdata = raw;
   data->cols = c;
   data->rows = r;
   data->c_offset = co;
   data->r_offset = ro;
   data->rgn_c_last = (crgn > 1) ? co + crgn : c - 1; /* default */
   data->rgn_r_last = (rrgn > 1) ? ro + rrgn : r - 1; /* to rest */
   data->user_data = user;
   data->next = NULL;

   if (r <= 1 || c <= 1 || co < 0 || ro < 0)
      return 0;

   data->dx = xspan/(c - 1);
   data->dy = yspan/(r - 1);
   data->x = x0 + (data->dx)*co;
   data->y = y0 + (data->dy)*ro;

   data->work1 = work1;
   data->work2 = work2;

   return VALID_REGION(data);
}

/* ctree_compose_cell() - set-up cell data parameters	*/
unsigned char
ctree_compose_cell(ctree_cell *cell,
      coord_val x0, coord_val y0, coord_val dx, coord_val dy,
      ctree_lvl llevel, ctree_lvl ulevel) {

   ctree_locval *n = cell->node;
   int low_a=1, high_a=2, low_b=3, high_b=4; /* sort temp */

   cell->index = 0;

   /* values for corners already in cell */
                         /* sort to determine low and high values */
   if (n[1].v > n[2].v) { low_a = 2; high_a = 1; }
   if (n[3].v > n[4].v) { low_b = 4; high_b = 3; }
   cell->minval  = n[(n[low_a].v < n[low_b].v)? low_a : low_b].v;
   cell->maxval = n[(n[high_a].v > n[high_b].v)? high_a : high_b].v;

                               /* check if we need go any further */
   if (ulevel < cell->minval || llevel > cell->maxval)
      return 0;         /* NOTE: we know zero has no subtrees */

                           /* center point is just an average */
   n[5].v = (n[1].v + n[2].v + n[3].v + n[4].v)/4.0;

                           /* compute index from order relations */
   cell->index  = ((n[4].v > n[5].v) << 7)/* bit 7: on if node 4 > center*/
                | ((n[3].v > n[5].v) << 6)/* bit 6: on if node 3 > center*/
                | ((n[2].v > n[5].v) << 5)/* etc.  */
                | ((n[1].v > n[5].v) << 4)
                | ((n[4].v > n[1].v) << 3)
                | ((n[3].v > n[4].v) << 2)
                | ((n[2].v > n[3].v) << 1)
                |  (n[1].v > n[2].v);

   n[1].x = n[4].x = x0;     /* x and y coordinates for 4 corners */
   n[2].x = n[3].x = x0 + dx;

   n[1].y = n[2].y = y0;
   n[3].y = n[4].y = y0 + dy;

   n[5].x = x0 + dx*0.5;     /* center x and y coordinates */
   n[5].y = y0 + dy*0.5;

   return cell->index;
}

/* ctree_plot_region() - Output contours for a subset of the data.*/
/*                       Stops and returns zero if data func      */
/*                       returns zero.                            */
/* data - structure describing the data to be contoured           */
/* level - countour levels array                                  */
/* nlevels - number of contour levels to be plotted               */
/* data_func - function to fetch data (returns num items)         */
/*      takes parms (*data, row, begincol, endcol, outputarray[]) */
/* plot_func - function to plot contours                          */
/*      takes parms (*displaylist, dlcount, levelnum, *user_data) */
/* user_plot_data - supplemental data passed to plot routine      */
/* NOTE: the data func returns rows of data[bgncol..endcol]       */
/*       If your data works best the other way, just modify data  */
/*       since which is the row and which is the col is a matter  */
/*       of opinion.                                              */
int
ctree_plot_region(ctree_data *data, ctree_lvl level[], int nlevels,
      ctree_data_func data_func, ctree_plot_func plot_func,
      void *user_plot_data) {
   ctree_cell cell;
   ctree_disp_list dlist[DLIST_MAX];
   unsigned char idx;
   ctree_lvl llevel, ulevel;
   int ir, ic, root, i, ioff, rc, ndl=0;
   coord_val x, y, dx, dy;
   ctree_val *prev_row, *next_row, *temp_row;

                               /* zero out unused array elements */
   cell.node[0].x=cell.node[0].y=cell.node[0].v=0.0;

   llevel = ulevel = level[0];	/* determine min and max levels */
   for (i=0; i<nlevels; i++) {
      if (ulevel < level[i]) ulevel = level[i];
      if (llevel > level[i]) llevel = level[i];
   }
   dx = data->dx;
   dy = data->dy;

   prev_row = data->work1;
   next_row = data->work2;

   rc = (*data_func)(data, data->r_offset,
                     data->c_offset, data->rgn_c_last, next_row);

   if (rc == 0) return 0; /* row function wants us to stop */

   y = data->y;
   for (ir = data->r_offset; ir < data->rgn_r_last; ir++, y+=dy) {
      temp_row  = prev_row;
      prev_row = next_row;
      next_row = temp_row;

      rc = (*data_func)(data, ir+1, data->c_offset,
                              data->rgn_c_last, next_row);
      if (rc == 0) return 0; /* request to stop */

      x = data->x;
      for (ioff = 0, ic = data->c_offset;
            ic < data->rgn_c_last; ic++, ioff++, x+=dx) {

					/* extract 2x2 data */
         cell.node[1].v = prev_row[ioff];
         cell.node[2].v = prev_row[ioff + 1];
         cell.node[3].v = next_row[ioff + 1];
         cell.node[4].v = next_row[ioff];

         idx = ctree_compose_cell(&cell, x, y, dx, dy,
                                                   llevel, ulevel);

                                            /* do contours if any */
         if ((root = Contree[idx][0]) >= 0) {
            for (i = 0; i<nlevels; i++)
               if (level[i] >= cell.minval
                && level[i] <= cell.maxval) {
                  ndl = 0; /* reset display list */
                  ctree_visit(root, root, &cell, level[i], dlist,&ndl);
                  (*plot_func)(dlist, ndl, i, user_plot_data);
               }
                                     /* some require two subtrees */
               if ((root = Contree[idx][1]) >= 0)
                  for (i = 0; i<nlevels; i++)
                     if (level[i] >= cell.minval
                      && level[i] <= cell.maxval) {
                        ndl = 0; /* reset display list */
                        ctree_visit(root,root,&cell,level[i],dlist,&ndl);
                        (*plot_func)(dlist, ndl, i, user_plot_data);
               }
         } /* end first if root */
      } /* next across */
   } /* next down */

   return 1;
}

/* ctree_visit() - generate contour for a single cell */
void
ctree_visit(int node, int ancestor, ctree_cell *cell,
      ctree_lvl level, ctree_disp_list dl[], int *ndl) {
   ctree_locval *n, *a;
   coord_val delta, fract;

   if (node < 0)
      return;

   n = &(cell->node[CT_leaf[node].n]);
   a = &(cell->node[CT_leaf[ancestor].n]);

   if ((n->v <= level && level < a->v)
	 || (n->v == level && node == ancestor)) {
                                      /* edge contains the level */

                                      /* interpolate edge */
      delta = (a->v - n->v);
      fract = (-EPSILON < delta && delta < EPSILON)
                          ? 0.5 : ((level - n->v) / delta);

                              /* store plot info in display list */
      dl[*ndl].x     = n->x + fract*(a->x - n->x);
      dl[*ndl].y     = n->y + fract*(a->y - n->y);
      dl[(*ndl)++].c = CT_leaf[node].comm;

                               /* check for equal valued edges */
      ctree_visit_subtree(CT_leaf[node].d[0],node,cell,level,dl,ndl);
      ctree_visit_subtree(CT_leaf[node].d[1],node,cell,level,dl,ndl);
      ctree_visit_subtree(CT_leaf[node].d[2],node,cell,level,dl,ndl);
   }
   else {                     /* visit children */
      ctree_visit(CT_leaf[node].d[0], node, cell, level, dl, ndl);
      ctree_visit(CT_leaf[node].d[1], node, cell, level, dl, ndl);
      ctree_visit(CT_leaf[node].d[2], node, cell, level, dl, ndl);
   }
}

/* ctree_visit_subtree() - check for equal valued edges */
void
ctree_visit_subtree(int node,int ancestor, ctree_cell *cell, ctree_lvl level,
      ctree_disp_list dl[], int *ndl) {
   ctree_locval *n, *a;

   if (node < 0)
      return;

   n = &(cell->node[CT_leaf[node].n]);

   if (n->v == level) {	/* found adjacent edge with same value  */
      a = &(cell->node[CT_leaf[ancestor].n]);

                         /* store move-draw info in display list */
      dl[*ndl].x     = a->x;
      dl[*ndl].y     = a->y;
      dl[(*ndl)++].c = CT_MOVE;

      dl[*ndl].x     = n->x;
      dl[*ndl].y     = n->y;
      dl[(*ndl)++].c = CT_DRAW;
   }
   if (n->v >= level) {
      ctree_visit_subtree(CT_leaf[node].d[0],node,cell,level,dl,ndl);
      ctree_visit_subtree(CT_leaf[node].d[1],node,cell,level,dl,ndl);
      ctree_visit_subtree(CT_leaf[node].d[2],node,cell,level,dl,ndl);
   }
}
