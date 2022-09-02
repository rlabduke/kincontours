/*			    nrmat.c                               */
/* Matrix allocation routines from Numerical Recipes in C, 2nd ed */
/* By: Press, Teukolsky, Vetterling, and Flannery                 */

/* Modifications:                                                     */
/*  8/30/22 - mgp - changed format ‘%d’ to %ld in printf w/ long int  */


#include <stdio.h>
#include <stdlib.h>
#include "nrmat.h"

/* vector() transcribed from Numerical Recipes */
/* allocates a real vector with range[nl..nh] */

real_val* vector(long nl, long nh) {
   real_val *v;
   
   v = (real_val *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(real_val)));
   if (!v) {
      fprintf(stderr,"ERROR: out of memory in vector() [%ld]\n", nh-nl+1);
      return NULL;
   }
   
   return v - nl + NR_END;
}

/* ivector() transcribed from Numerical Recipes */
/* allocates an integer vector with range[nl..nh] */

int* ivector(long nl, long nh) {
   int *v;
   
   v = (int *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) {
      fprintf(stderr,"ERROR: out of memory in ivector() [%ld]\n", nh-nl+1);
      return NULL;
   }
   
   return v - nl + NR_END;
}

/* matrix() transcribed from Numerical Recipes */
/* allocates a real matrix with range[nrl..nrh][ncl..nch] */

real_val** matrix(long nrl, long nrh, long ncl, long nch) {
   long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
   real_val **m;
   
   /* allocate pointers to rows */
   m = (real_val **)malloc((size_t)((nrow+NR_END)*sizeof(real_val *)));
   if (!m) {
      fprintf(stderr,"ERROR: out of memory in matrix() [%ld]\n", nrow);
      return NULL;
   }
   m += NR_END;
   m -= nrl;
   
   /* allocate rows and set pointers to them */
   m[nrl] = (real_val *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(real_val)));
   if (!m[nrl]) {
      fprintf(stderr,"ERROR: out of memory in matrix() [%ld][%ld]\n", nrow, ncol);
      free((void*) (m + nrl - NR_END));
      return NULL;
   }
   m[nrl] += NR_END;
   m[nrl] -= ncl;
   
   for(i=nrl+1;i<=nrh;i++) { m[i] = m[i-1] + ncol; }

   /* return pointer to array of pointers to rows */
   return m;
}

/* matrix3D() transcribed from f3tensor() in Numerical Recipes */
/* allocates a real matrix with range[nrl..nrh][ncl..nch][ndl..ndh] */

real_val*** matrix3D(long nrl,long nrh,long ncl,long nch,long ndl,long ndh){
   long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
   real_val ***t;

   /* allocate pointers to pointers to rows */
   t = (real_val ***)malloc((size_t)((nrow+NR_END)*sizeof(real_val **)));
   if (!t) {
      fprintf(stderr,"ERROR: out of memory in matrix3D() [%ld]\n", nrow);
      return NULL;
   }
   t += NR_END;
   t -= nrl;
   
   /* allocate pointers to rows and set pointers to them */
   t[nrl] = (real_val **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(real_val *)));
   if (!t[nrl]) {
      fprintf(stderr,"ERROR: out of memory in matrix3D() [%ld][%ld]\n", nrow,ncol);
      free((void*) (t + nrl - NR_END));
      return NULL;
   }
   t[nrl] += NR_END;
   t[nrl] -= ncl;

   /* allocate rows and set pointers to them */
   t[nrl][ncl] = (real_val *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(real_val)));
   if (!t[nrl][ncl]) {
      fprintf(stderr,"ERROR: out of memory in matrix3D() [%ld][%ld][%ld]\n", nrow,ncol,ndep);
      free((void*) (t[nrl] + ncl - NR_END));
      free((void*) (t + nrl - NR_END));
      return NULL;
   }
   t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;
   
   for(j=ncl+1;j<=nch;j++) { t[nrl][j] = t[nrl][j-1] + ndep; }
   for(i=nrl+1;i<=nrh;i++) {
      t[i] = t[i-1] + ncol;
      t[i][ncl] = t[i-1][ncl] + ncol*ndep;
      for(j=ncl+1;j<=nch;j++) { t[i][j] = t[i][j-1] + ndep; }
   }

   /* return pointer to array of pointers to rows */
   return t;
}

void free_vector(real_val *v, long nl) {
   free((void*) (v + nl - NR_END));
}

void free_ivector(int *v, long nl) {
   free((void*) (v + nl - NR_END));
}

void free_matrix(real_val **m, long nrl, long ncl) {
   free((void*) (m[nrl] + ncl - NR_END));
   free((void*) (m + nrl - NR_END));
}

void free_matrix3D(real_val ***t, long nrl, long ncl, long ndl) {
   free((void*) (t[nrl][ncl] + ndl - NR_END));
   free((void*) (t[nrl] + ncl - NR_END));
   free((void*) (t + nrl - NR_END));
}
