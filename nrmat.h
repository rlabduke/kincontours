/*			        nrmat.h                           */
/* Matrix allocation routines from Numerical Recipes in C, 2nd ed */
/* By: Press, Teukolsky, Vetterling, and Flannery                 */

#ifndef NRMAT_H
#define NRMAT_H 1

typedef float real_val;

/* we can set NR_END to 1 if base pointer b does not equal (b-1)+1 */
#define NR_END 0

real_val* vector(long nl, long nh);
int* ivector(long nl, long nh);
real_val** matrix(long nrl, long nrh, long ncl, long nch);
real_val*** matrix3D(long nrl,long nrh,long ncl,long nch,long ndl,long ndh);
void free_vector(real_val *v, long nl);
void free_ivector(int *v, long nl);
void free_matrix(real_val **m, long nrl, long ncl);
void free_matrix3D(real_val ***t, long nrl, long ncl, long ndl);

#endif
