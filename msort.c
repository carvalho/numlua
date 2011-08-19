#include <lua.h>
#include <complex.h>
#include "numlua.h"

/* Adapted from ssort.f (translated by f2c) in SLATEC: */
/* ***PURPOSE  Sort an array and optionally make the same interchanges in */
/*            an auxiliary array.  The array may be sorted in increasing */
/*            or decreasing order.  A slightly modified QUICKSORT */
/*            algorithm is used. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  N6A2B */
/* ***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I) */
/* ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */

#define SORTDIM (31)

/* assume mx is not a section and complex */
int sort1c (nl_Matrix *mx) {
  int i, j, k, l, m, ij;
  int il[SORTDIM], iu[SORTDIM];
  int n = mx->size;
  int s = mx->stride;
  lua_Number r;
  nl_Complex t, tt;
  nl_Complex *x = (nl_Complex *) mx->data;

  if (n < 1) return 0;
  m = 1;
  i = 0;
  j = (n - 1) * s;
  r = .375;

L20:
  if (i == j) goto L60;
  if (r <= .5898437)
    r += .0390625;
  else
    r += -.21875;

L30:
  k = i;
  /* Select a central element of the array and save it in location T */
  ij = i + ((int) ((j - i) / s * r)) * s;
  t = x[ij];
  /* If first element of array is greater than T, interchange with T */
  if (cgt(x[i], t)) {
    x[ij] = x[i]; x[i] = t; t = x[ij];
  }
  l = j;
  /* If last element of array is less than than T, interchange with T */
  if (clt(x[j], t)) {
    x[ij] = x[j]; x[j] = t; t = x[ij];
    /* If first element of array is greater than T, interchange with T */
    if (cgt(x[i], t)) {
      x[ij] = x[i]; x[i] = t; t = x[ij];
    }
  }

  /* Find an element in the second half of the array which is smaller than T */
L40:
  l -= s;
  if (cgt(x[l], t)) goto L40;

  /* Find an element in the first half of the array which is greater than T */
L50:
  k += s;
  if (clt(x[k], t)) goto L50;
  /* Interchange these elements */
  if (k <= l) {
    tt = x[l]; x[l] = x[k]; x[k] = tt;
    goto L40;
  }
  /* Save upper and lower subscripts of the array yet to be sorted */
  if (l - i > j - k) {
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    ++m;
  }
  else {
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
  }
  goto L70;

  /* Begin again on another portion of the unsorted array */
L60:
  --m;
  if (m == 0) return 0;
  i = il[m - 1];
  j = iu[m - 1];

L70:
  if (j - i > 0) goto L30;
  if (i == 0) goto L20;
  i -= s;

L80:
  i += s;
  if (i == j) goto L60;
  t = x[i + s];
  if (!cgt(x[i], t)) goto L80;
  k = i;

L90:
  x[k + s] = x[k];
  k -= s;
  if (clt(t, x[k])) goto L90;
  x[k + s] = t;
  goto L80;
}

/* assume mx is not a section and not complex */
int sort1d (nl_Matrix *mx) {
  int i, j, k, l, m, ij;
  int il[SORTDIM], iu[SORTDIM];
  int n = mx->size;
  int s = mx->stride;
  lua_Number r, t, tt;
  lua_Number *x = mx->data;

  if (n < 1) return 0;
  m = 1;
  i = 0;
  j = (n - 1) * s;
  r = .375;

L20:
  if (i == j) goto L60;
  if (r <= .5898437)
    r += .0390625;
  else
    r += -.21875;

L30:
  k = i;
  /* Select a central element of the array and save it in location T */
  ij = i + ((int) ((j - i) / s * r)) * s;
  t = x[ij];
  /* If first element of array is greater than T, interchange with T */
  if (x[i] > t) {
    x[ij] = x[i]; x[i] = t; t = x[ij];
  }
  l = j;
  /* If last element of array is less than than T, interchange with T */
  if (x[j] < t) {
    x[ij] = x[j]; x[j] = t; t = x[ij];
    /* If first element of array is greater than T, interchange with T */
    if (x[i] > t) {
      x[ij] = x[i]; x[i] = t; t = x[ij];
    }
  }

  /* Find an element in the second half of the array which is smaller than T */
L40:
  l -= s;
  if (x[l] > t) goto L40;

  /* Find an element in the first half of the array which is greater than T */
L50:
  k += s;
  if (x[k] < t) goto L50;
  /* Interchange these elements */
  if (k <= l) {
    tt = x[l]; x[l] = x[k]; x[k] = tt;
    goto L40;
  }
  /* Save upper and lower subscripts of the array yet to be sorted */
  if (l - i > j - k) {
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    ++m;
  }
  else {
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
  }
  goto L70;

  /* Begin again on another portion of the unsorted array */
L60:
  --m;
  if (m == 0) return 0;
  i = il[m - 1];
  j = iu[m - 1];

L70:
  if (j - i > 0) goto L30;
  if (i == 0) goto L20;
  i -= s;

L80:
  i += s;
  if (i == j) goto L60;
  t = x[i + s];
  if (x[i] <= t) goto L80;
  k = i;

L90:
  x[k + s] = x[k];
  k -= s;
  if (t < x[k]) goto L90;
  x[k + s] = t;
  goto L80;
}



/* x and y are size-consistent, y is simple vector */
int sort2c (nl_Matrix *mx, nl_Matrix *my) {
  int i, j, k, l, m, ij;
  int ix, jx, kx, lx, ijx;
  int il[SORTDIM], iu[SORTDIM];
  int n = mx->size;
  int s = mx->stride;
  lua_Number r, ty, tty;
  nl_Complex t, tt;
  nl_Complex *x = (nl_Complex *) mx->data;
  lua_Number *y = my->data;

  if (n < 1) return 0;
  m = 1;
  i = 0; ix = 0;
  j = n - 1; jx = j * s;
  r = .375;

L110:
  if (i == j) goto L150;
  if (r <= .5898437)
    r += .0390625;
  else
    r += -.21875;

L120:
  k = i; kx = ix;
  /* Select a central element of the array and save it in location T */
  ij = i + (int) ((j - i) * r); ijx = ij * s;
  t = x[ijx];
  ty = y[ij];
  /* If first element of array is greater than T, interchange with T */
  if (cgt(x[ix], t)) {
    x[ijx] = x[ix]; x[ix] = t; t = x[ijx];
    y[ij] = y[i]; y[i] = ty; ty = y[ij];
  }
  l = j; lx = jx;
  /* If last element of array is less than T, interchange with T */
  if (clt(x[jx], t)) {
    x[ijx] = x[jx]; x[jx] = t; t = x[ijx];
    y[ij] = y[j]; y[j] = ty; ty = y[ij];
    /* If first element of array is greater than T, interchange with T */
    if (cgt(x[ix], t)) {
      x[ijx] = x[ix]; x[ix] = t; t = x[ijx];
      y[ij] = y[i]; y[i] = ty; ty = y[ij];
    }
  }

  /* Find an element in the second half of the array which is smaller than T */
L130:
  --l; lx -= s;
  if (cgt(x[lx], t)) goto L130;

  /* Find an element in the first half of the array which is greater than T */
L140:
  ++k; kx += s;
  if (clt(x[kx], t)) goto L140;
  /* Interchange these elements */
  if (k <= l) {
    tt = x[lx]; x[lx] = x[kx]; x[kx] = tt;
    tty = y[l]; y[l] = y[k]; y[k] = tty;
    goto L130;
  }
  /* Save upper and lower subscripts of the array yet to be sorted */
  if (l - i > j - k) {
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k; ix = kx;
    ++m;
  }
  else {
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l; jx = lx;
    ++m;
  }
  goto L160;

  /* Begin again on another portion of the unsorted array */
L150:
  --m;
  if (m == 0) return 0;
  i = il[m - 1]; ix = i * s;
  j = iu[m - 1]; jx = j * s;

L160:
  if (j - i > 0) goto L120;
  if (i == 0) goto L110;
  --i; ix -= s;

L170:
  ++i; ix += s;
  if (i == j) goto L150;
  t = x[ix + s];
  ty = y[i + 1];
  if (!cgt(x[ix], t)) goto L170;
  k = i; kx = ix;

L180:
  x[kx + s] = x[kx];
  y[k + 1] = y[k];
  --k; kx -= s;
  if (clt(t, x[kx])) goto L180;
  x[kx + s] = t;
  y[k + 1] = ty;
  goto L170;
}


/* x and y are size-consistent, y is simple vector */
int sort2d (nl_Matrix *mx, nl_Matrix *my) {
  int i, j, k, l, m, ij;
  int ix, jx, kx, lx, ijx;
  int il[SORTDIM], iu[SORTDIM];
  int n = mx->size;
  int s = mx->stride;
  lua_Number r, t, tt, ty, tty;
  lua_Number *x = mx->data;
  lua_Number *y = my->data;

  if (n < 1) return 0;
  m = 1;
  i = 0; ix = 0;
  j = n - 1; jx = j * s;
  r = .375;

L110:
  if (i == j) goto L150;
  if (r <= .5898437)
    r += .0390625;
  else
    r += -.21875;

L120:
  k = i; kx = ix;
  /* Select a central element of the array and save it in location T */
  ij = i + (int) ((j - i) * r); ijx = ij * s;
  t = x[ijx];
  ty = y[ij];
  /* If first element of array is greater than T, interchange with T */
  if (x[ix] > t) {
    x[ijx] = x[ix]; x[ix] = t; t = x[ijx];
    y[ij] = y[i]; y[i] = ty; ty = y[ij];
  }
  l = j; lx = jx;
  /* If last element of array is less than T, interchange with T */
  if (x[jx] < t) {
    x[ijx] = x[jx]; x[jx] = t; t = x[ijx];
    y[ij] = y[j]; y[j] = ty; ty = y[ij];
    /* If first element of array is greater than T, interchange with T */
    if (x[ix] > t) {
      x[ijx] = x[ix]; x[ix] = t; t = x[ijx];
      y[ij] = y[i]; y[i] = ty; ty = y[ij];
    }
  }

  /* Find an element in the second half of the array which is smaller than T */
L130:
  --l; lx -= s;
  if (x[lx] > t) goto L130;

  /* Find an element in the first half of the array which is greater than T */
L140:
  ++k; kx += s;
  if (x[kx] < t) goto L140;
  /* Interchange these elements */
  if (k <= l) {
    tt = x[lx]; x[lx] = x[kx]; x[kx] = tt;
    tty = y[l]; y[l] = y[k]; y[k] = tty;
    goto L130;
  }
  /* Save upper and lower subscripts of the array yet to be sorted */
  if (l - i > j - k) {
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k; ix = kx;
    ++m;
  }
  else {
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l; jx = lx;
    ++m;
  }
  goto L160;

  /* Begin again on another portion of the unsorted array */
L150:
  --m;
  if (m == 0) return 0;
  i = il[m - 1]; ix = i * s;
  j = iu[m - 1]; jx = j * s;

L160:
  if (j - i > 0) goto L120;
  if (i == 0) goto L110;
  --i; ix -= s;

L170:
  ++i; ix += s;
  if (i == j) goto L150;
  t = x[ix + s];
  ty = y[i + 1];
  if (x[ix] <= t) goto L170;
  k = i; kx = ix;

L180:
  x[kx + s] = x[kx];
  y[k + 1] = y[k];
  --k; kx -= s;
  if (t < x[kx]) goto L180;
  x[kx + s] = t;
  y[k + 1] = ty;
  goto L170;
}

