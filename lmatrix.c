/* {=================================================================
 *
 * lmatrix.c
 * Multidimensional matrix library for NumericLua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lua.h>
#include <lauxlib.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include "numlua.h"

/* Constants */
static int one = 1;
static int two = 2;
static lua_Number minusone = -1;
static lua_Number oned = 1.;
static nl_Complex onec = 1.;


/* TODO [wishlist]:
 *
 * _cshift_(m, shift [, what]): if `m` is 2D, `shift` can be vector depending
 * on `what`
 *
 * _eoshift_(m, shift [, boundary] [, what]): if `m` is 2D, `shift` and
 * `boundary` can be vector depending on `what`
 *
 * `concat` is a _builder_ from lower equal dimensions to a higher dimension
 * Options: m1, m2, ..., mn with dims (d1, d2, ..., dp)
 *   (i) Increase dim on dim k: (d1, d2, d{k-1}, n, dk, ..., dp) ["join"]
 *  (ii) Across dim k: (d1, d2, ..., sum(dk), ..., dp) ["concat"]
 *   - particular cases: rowcat (dim=2, k=1), colcat (dim=2, k=2), c (dim=1)
 *
 * advanced linear algebra: `schur` (_gees), `hess` (_gehrd), `qz`
 */


static int matrix_mt_ = 0;
#define MATRIX_MT ((void *)&matrix_mt_)
#define SECTION_STUB ((nl_Section *)&matrix_mt_)
#define DATA_STUB ((lua_Number *)&matrix_mt_)

/* error messages and checks */
#define CONF_ERROR "dimensions are not conformable"
#define DOMC_ERROR "domains are not consistent"

#define checkvector(L,m,a) \
  luaL_argcheck(L, (m)->ndims==1, (a), "vector expected")
#define checkrealvector(L,m,a) \
  luaL_argcheck(L, (m)->ndims==1 && !(m)->iscomplex, \
      (a), "real vector expected")
#define checkarray(L,m,a) \
  luaL_argcheck(L, (m)->ndims==2, (a), "array expected")
#define check2dmatrix(L,m,a) \
  luaL_argcheck(L, (m)->ndims <= 2, (a), "two-dimensional matrix expected")
#define checksquare(L,m,a) \
  luaL_argcheck(L, (m)->ndims==2 && (m)->dim[0]==(m)->dim[1], \
      (a), "square matrix expected")

#define checksection(L,m,a) \
  luaL_argcheck(L, !(m)->section, (a), "sections are not allowed")
#define check2dsection(L,m,a) \
  luaL_argcheck(L, (m)->stride==1 \
      && (!(m)->section \
        || ((m)->section[0].step==1 && (m)->section[1].step==1)), \
      (a), "only simple array sections are allowed")
/* note: simple array sections: unitary steps */


/* {=====================================================================
 *    Auxiliary
 * ======================================================================} */

/* i must be != 0 */
#define CIRC(i,n) ((i)>0 ? (((i)-1)%(n)+1) : ((i)+1)%(n)+(n))

#define conformable(a,b) \
  ((a)->size==(b)->size && (a)->iscomplex==(b)->iscomplex)

#define LD(m,i) (((m)->section) ? (m)->section[i].ld : (m)->dim[i])
#define STEP(m,i) (((m)->section) ? (m)->section[i].step : 1)

#define DSHIFT(m,i) ((m)->data + nl_msshift((m),(i)))
#define CSHIFT(m,i) (CPX((m)->data) + nl_msshift((m),(i)))

/* overall stride in section `m` from element-order index `eo` (zero-based)
 * note that nl_msshift(m, 0) == 0 */
NUMLUA_API int nl_msshift (nl_Matrix *m, int eo) {
  int i, d, shift = 0;
  int stride = m->stride;
  for (i = 0; i < m->ndims; i++) {
    d = eo % m->dim[i];
    shift += (d * m->section[i].step) * stride;
    stride *= m->section[i].ld;
    eo = (eo - d) / m->dim[i];
  }
  return shift;
}

/* [[ Internal ]] */

/* meant for use on functions that have MT as env */
static nl_Matrix *tomatrix (lua_State *L, int narg) {
  nl_Matrix *p = NULL;
  if (lua_type(L, narg) == LUA_TUSERDATA /* userdata? */
      && lua_getmetatable(L, narg)) {  /* has metatable? */
    if (lua_rawequal(L, -1, lua_upvalueindex(1))) /* MT == upvalue? */
      p = (nl_Matrix *) lua_touserdata(L, narg);
    lua_pop(L, 1); /* remove metatable */
  }
  return p;
}

static nl_Matrix *checkmatrix (lua_State *L, int narg) {
  nl_Matrix *p = tomatrix(L, narg);
  if (p == NULL) nl_typeerror(L, narg, MATRIX_LIBNAME);
  return p;
}

/* assumes data block is at top of stack if data != NULL or creates a new data
 * block otherwise */
#define blocksize(c,n) ((n) * ((c) ? sizeof(nl_Complex) : sizeof(lua_Number)))
static nl_Matrix *pushmatrix (lua_State *L, int iscomplex, int ndims,
    int *dim, int stride, int size, nl_Section *section, lua_Number *data) {
  lua_Number *mdata = (data != NULL) ? data
    : lua_newuserdata(L, blocksize(iscomplex, size));
  nl_Matrix *m = lua_newuserdata(L, sizeof(nl_Matrix)
      + (ndims - 1) * sizeof(int)
      + ((section) ? ndims * sizeof(nl_Section) : 0));
  int i;
  lua_pushvalue(L, lua_upvalueindex(1)); /* metatable */
  lua_pushvalue(L, -2);
  lua_pushvalue(L, -4); /* data block */
  lua_rawset(L, -3); /* mt[matrix] = data */
  m->ndims = ndims;
  m->stride = stride;
  m->size = size;
  m->iscomplex = iscomplex;
  m->section = NULL;
  m->data = mdata;
  if (dim) {
    for (i = 0; i < ndims; i++)
      m->dim[i] = dim[i];
  }
  if (section) {
    m->section = (nl_Section *)((int *)(m + 1) + ndims - 1);
    if (section != SECTION_STUB) { /* fill? */
      for (i = 0; i < ndims; i++) {
        m->section[i].ld = section[i].ld;
        m->section[i].step = section[i].step;
      }
    }
  }
  lua_setmetatable(L, -2);
  if (mdata != DATA_STUB)
    lua_replace(L, -2); /* pop data block */
  return m;
}

#define pushframe(L,m) \
  pushmatrix(L, (m)->iscomplex, (m)->ndims, (m)->dim, 1, (m)->size, \
      NULL, NULL)

/* NOTE: setdata operates column-wise */
static void setdatatoscalar (int iscomplex, int size, nl_Complex s,
    int stride, int shift, lua_Number *data) {
  int i;
  if (iscomplex) {
    nl_Complex *e = CPX(data) + shift;
    for (i = 0; i < size; i++, e += stride)
      *e = s;
  }
  else {
    lua_Number *e = data + shift;
    lua_Number rs = creal(s);
    for (i = 0; i < size; i++, e += stride)
      *e = rs;
  }
}

static void setdatatovector (nl_Matrix *v,
    int stride, int shift, lua_Number *data) {
  int i;
  if (v->section) {
    if (v->iscomplex) {
      nl_Complex *e = CPX(data) + shift;
      for (i = 0; i < v->size; i++, e += stride)
        *e = *(CSHIFT(v, i));
    }
    else {
      lua_Number *e = data + shift;
      for (i = 0; i < v->size; i++, e += stride)
        *e = *(DSHIFT(v, i));
    }
  }
  else { /* regular copy */
    if (v->stride >= 0) {
      if (v->iscomplex)
        ZCOPY(&v->size, CPX(v->data), &v->stride, CPX(data) + shift, &stride);
      else
        DCOPY(&v->size, v->data, &v->stride, data + shift, &stride);
    }
    else {
      if (v->iscomplex) {
        nl_Complex *e = CPX(data) + shift;
        for (i = 0; i < v->size; i++, e += stride)
          *e = CPX(v->data)[i * v->stride];
      }
      else {
        lua_Number *e = data + shift;
        for (i = 0; i < v->size; i++, e += stride)
          *e = v->data[i * v->stride];
      }
    }
  }
}

/* diagonal entries are not set; assumes matrix is 2d section */
static void settriangtoscalar (nl_Complex vc, char what, nl_Matrix *m) {
  int i, j, s;
  int n = (m->dim[0] < m->dim[1]) ? m->dim[0] : m->dim[1]; /* min */
  int ldm = LD(m, 0);
  if (m->iscomplex) {
    if (what == 'l' || what == 'L') { /* lower? */
      for (i = 0; i < n; i++) {
        s = i * (ldm + 1) + 1;
        for (j = 0; j < m->dim[0] - i - 1; j++)
          CPX(m->data)[(s + j) * m->stride] = vc;
      }
    }
    else { /* upper */
      for (i = 0; i < n; i++) {
        s = i * (ldm + 1) + ldm;
        for (j = 0; j < m->dim[1] - i - 1; j++)
          CPX(m->data)[(s + j * ldm) * m->stride] = vc;
      }
    }
  }
  else {
    lua_Number rvc = creal(vc);
    if (what == 'l' || what == 'L') { /* lower? */
      for (i = 0; i < n; i++) {
        s = i * (ldm + 1) + 1;
        for (j = 0; j < m->dim[0] - i - 1; j++)
          m->data[(s + j) * m->stride] = rvc;
      }
    }
    else { /* upper */
      for (i = 0; i < n; i++) {
        s = i * (ldm + 1) + ldm;
        for (j = 0; j < m->dim[1] - i - 1; j++)
          m->data[(s + j * ldm) * m->stride] = rvc;
      }
    }
  }
}

/* destination (`m`) should have no larger dimensions than source (`v`) */
static void settriangtovector (nl_Matrix *v, char what, nl_Matrix *m) {
  int i, l, sv, sm;
  int n = (m->dim[0] < m->dim[1]) ? m->dim[0] : m->dim[1]; /* min */
  int ldv = LD(v, 0), ldm = LD(m, 0);
  if (what == 'l' || what == 'L') { /* lower? */
    l = m->dim[0] - 1; /* except diagonal */
    if (m->iscomplex)
      for (i = 0; i < n; i++, l--) {
        sv = (i * (ldv + 1) + 1) * v->stride;
        sm = (i * (ldm + 1) + 1) * m->stride;
        ZCOPY(&l, CPX(v->data) + sv, &v->stride,
            CPX(m->data) + sm, &m->stride);
      }
    else
      for (i = 0; i < n; i++, l--) {
        sv = (i * (ldv + 1) + 1) * v->stride;
        sm = (i * (ldm + 1) + 1) * m->stride;
        DCOPY(&l, v->data + sv, &v->stride, m->data + sm, &m->stride);
      }
  }
  else { /* upper */
    int stv = ldv * v->stride;
    int stm = ldm * m->stride;
    l = m->dim[1] - 1; /* except diagonal */
    if (m->iscomplex)
      for (i = 0; i < n; i++, l--) {
        sv = (i * (ldv + 1) + ldv) * v->stride;
        sm = (i * (ldm + 1) + ldm) * m->stride;
        ZCOPY(&l, CPX(v->data) + sv, &stv, CPX(m->data) + sm, &stm);
      }
    else
      for (i = 0; i < n; i++, l--) {
        sv = (i * (ldv + 1) + ldv) * v->stride;
        sm = (i * (ldm + 1) + ldm) * m->stride;
        DCOPY(&l, v->data + sv, &stv, m->data + sm, &stm);
      }
  }
}

/* sets submatrix of m given by (n, stride, size, shift) to value v at position
 * narg in stack, where v can be a number/complex or a conformable matrix to
 * the submatrix */
static void settoarg (lua_State *L, nl_Matrix *m, int n, int stride, int size,
    int shift, int narg) {
  int i, iscomplex;
  nl_Complex vc = nl_tocomplex(L, narg, &iscomplex);
  nl_Matrix *s = NULL;
  if (m->section) { /* setup subsection frame? */
    if (n > 0) {
      int ndims = m->ndims - n;
      s = lua_newuserdata(L, sizeof(nl_Matrix) + (ndims - 1) * sizeof(int));
      s->ndims = ndims;
      s->stride = stride;
      for (i = 0; i < ndims; i++)
        s->dim[i] = m->dim[n + i];
      s->section = m->section + n; /* remaining sections */
    }
    else s = m;
  }
  /* number/complex? */
  if (iscomplex) {
    if (m->section) {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data) + shift;
        for (i = 0; i < size; i++)
          e[nl_msshift(s, i)] = vc;
      }
      else {
        lua_Number *e = m->data + shift;
        lua_Number rvc = creal(vc);
        for (i = 0; i < size; i++)
          e[nl_msshift(s, i)] = rvc;
      }
    }
    else /* regular copy */
      setdatatoscalar(m->iscomplex, size, vc, stride, shift, m->data);
  }
  /* conformable matrix? */
  else if (lua_type(L, narg) == LUA_TUSERDATA) {
    nl_Matrix *v = checkmatrix(L, narg);
    luaL_argcheck(L, size == v->size && m->iscomplex == v->iscomplex, narg,
        CONF_ERROR);
    if (m->section) {
      if (v->section) {
        if (m->iscomplex) {
          nl_Complex *e = CPX(m->data) + shift;
          for (i = 0; i < size; i++)
            e[nl_msshift(s, i)] = *(CSHIFT(v, i));
        }
        else {
          lua_Number *e = m->data + shift;
          for (i = 0; i < size; i++)
            e[nl_msshift(s, i)] = *(DSHIFT(v, i));
        }
      }
      else {
        if (m->iscomplex) {
          nl_Complex *e = CPX(m->data) + shift;
          nl_Complex *f = CPX(v->data);
          for (i = 0; i < size; i++, f += v->stride)
            e[nl_msshift(s, i)] = *f;
        }
        else {
          lua_Number *e = m->data + shift;
          lua_Number *f = v->data;
          for (i = 0; i < size; i++, f += v->stride)
            e[nl_msshift(s, i)] = *f;
        }
      }
    }
    else
      setdatatovector(v, stride, shift, m->data);
  }
  if (m->section) lua_pop(L, 1); /* remove frame */
}


/* {=====================================================================
 *    API
 * ======================================================================} */

NUMLUA_API nl_Matrix *nl_tomatrix (lua_State *L, int narg) {
  nl_Matrix *p = NULL;
  if (lua_type(L, narg) == LUA_TUSERDATA /* userdata? */
      && lua_getmetatable(L, narg)) { /* has metatable? */
    lua_pushlightuserdata(L, MATRIX_MT);
    lua_rawget(L, LUA_REGISTRYINDEX);
    if (lua_rawequal(L, -1, -2)) /* right MT? */
      p = lua_touserdata(L, narg);
    lua_pop(L, 2); /* MTs */
  }
  return p;
}

NUMLUA_API nl_Matrix *nl_checkmatrix (lua_State *L, int narg) {
  nl_Matrix *p = nl_tomatrix(L, narg);
  if (p == NULL) nl_typeerror(L, narg, MATRIX_LIBNAME);
  return p;
}

/* assumes data block is at top of stack, cannot build sections */
NUMLUA_API nl_Matrix *nl_pushmatrix (lua_State *L, int iscomplex, int ndims,
    int *dim, int stride, int size, lua_Number *data) {
  lua_Number *mdata = (data != NULL) ? data
    : lua_newuserdata(L, blocksize(iscomplex, size));
  nl_Matrix *m = lua_newuserdata(L, sizeof(nl_Matrix)
      + (ndims - 1) * sizeof(int));
  int i;
  nl_getmetatable(L, MATRIX_MT);
  lua_pushvalue(L, -2);
  lua_pushvalue(L, -4); /* data block */
  lua_rawset(L, -3); /* mt[matrix] = data */
  m->ndims = ndims;
  m->stride = stride;
  m->size = size;
  m->iscomplex = iscomplex;
  m->section = NULL;
  m->data = mdata;
  if (dim) {
    for (i = 0; i < ndims; i++)
      m->dim[i] = dim[i];
  }
  lua_setmetatable(L, -2);
  if (mdata != DATA_STUB)
    lua_replace(L, -2); /* pop data block */
  return m;
}


/* {=======   Linear Algebra (API)   =======} */

/* returns the transpose of `m` in the stack */
static nl_Matrix *nl_transpose (lua_State *L, nl_Matrix *m, int hermitian) {
  nl_Matrix *t;
  if (m->ndims == 1) {
    t = pushframe(L, m);
    setdatatovector(m, 1, 0, t->data);
    if (hermitian && t->iscomplex) /* take conjugate? */
      DSCAL(&t->size, &minusone, t->data + 1, &two);
  }
  else {
    int i, j;
    int stride = m->stride * STEP(m, 0);
    int ld = m->stride * LD(m, 0) * STEP(m, 1);
    t = pushmatrix(L, m->iscomplex, m->ndims, NULL, 1, m->size, NULL, NULL);
    t->dim[0] = m->dim[1];
    t->dim[1] = m->dim[0];
    if (m->iscomplex) {
      nl_Complex *e, *f;
      for (j = 0; j < m->dim[1]; j++) { /* each column of m */
        e = CPX(m->data) + j * ld;
        f = CPX(t->data) + j;
        for (i = 0; i < m->dim[0]; i++) {
          *f = (hermitian) ? conj(*e) : *e;
          e += stride;
          f += t->dim[0];
        }
      }
    }
    else {
      lua_Number *e, *f;
      for (j = 0; j < m->dim[1]; j++) { /* each column of m */
        e = m->data + j * ld;
        f = t->data + j;
        for (i = 0; i < m->dim[0]; i++) {
          *f = *e;
          e += stride;
          f += t->dim[0];
        }
      }
    }
  }
  return t;
}


static void nl_trmul (nl_Matrix *x, nl_Matrix *a,
    char uplo, int invert, char trans, char side) {
  char diag = 'N'; /* not unit triangular */
  int n;
  int lda = LD(a, 0);
  if (x->ndims == 1) { /* vector? */
    int incx = x->stride;
    n = x->size;
    if (x->iscomplex) {
      if (invert)
        ZTRSV(&uplo, &trans, &diag, &n, CPX(a->data), &lda,
            CPX(x->data), &incx, 1, 1, 1);
      else
        ZTRMV(&uplo, &trans, &diag, &n, CPX(a->data), &lda,
            CPX(x->data), &incx, 1, 1, 1);
    }
    else {
      if (invert)
        DTRSV(&uplo, &trans, &diag, &n, a->data, &lda,
            x->data, &incx, 1, 1, 1);
      else
        DTRMV(&uplo, &trans, &diag, &n, a->data, &lda,
            x->data, &incx, 1, 1, 1);
    }
  }
  else { /* matrix */
    int m = x->dim[0];
    int ldb = LD(x, 0);
    n = x->dim[1];
    if (x->iscomplex) {
      if (invert)
        ZTRSM(&side, &uplo, &trans, &diag, &m, &n, &onec,
            CPX(a->data), &lda, CPX(x->data), &ldb, 1, 1, 1, 1);
      else
        ZTRMM(&side, &uplo, &trans, &diag, &m, &n, &onec,
            CPX(a->data), &lda, CPX(x->data), &ldb, 1, 1, 1, 1);
    }
    else {
      if (invert)
        DTRSM(&side, &uplo, &trans, &diag, &m, &n, (lua_Number *) &onec,
            a->data, &lda, x->data, &ldb, 1, 1, 1, 1);
      else
        DTRMM(&side, &uplo, &trans, &diag, &m, &n, (lua_Number *) &onec,
            a->data, &lda, x->data, &ldb, 1, 1, 1, 1);
    }
  }
}


static void nl_hemul (nl_Matrix *x, nl_Matrix *a, int inner, char what,
    lua_Number alpha) {
  int n = x->dim[0]; /* order */
  int ldx = LD(x, 0);
  int i, l, p;
  char uplo = 'L';
  if (what == 'u' || what == 'U') uplo = 'U';
  if (a->ndims == 1) { /* vector? */
    int inca = a->stride;
    if (x->iscomplex)
      ZHER(&uplo, &n, &alpha, CPX(a->data), &inca, CPX(x->data), &ldx, 1);
    else
      DSYR(&uplo, &n, &alpha, a->data, &inca, x->data, &ldx, 1);
  }
  else {
    char trans = inner ? 'C' : 'N';
    int k = inner ? a->dim[0] : a->dim[1];
    int lda = LD(a, 0);
    if (x->iscomplex)
      ZHERK(&uplo, &trans, &n, &k, &alpha, CPX(a->data), &lda, &oned,
          CPX(x->data), &ldx, 1, 1);
    else
      DSYRK(&uplo, &trans, &n, &k, &alpha, a->data, &lda, &oned,
          x->data, &ldx, 1, 1);
  }
  if (what == 'f' || what == 'F') { /* reflect lower to upper triangle? */
    if (x->iscomplex) {
      for (i = 1; i < n; i++) {
        l = n - i; p = (i - 1) * (ldx + 1);
        ZCOPY(&l, CPX(x->data) + p + 1, &one, CPX(x->data) + p + ldx, &ldx);
      }
    }
    else {
      for (i = 1; i < n; i++) {
        l = n - i; p = (i - 1) * (ldx + 1);
        DCOPY(&l, x->data + p + 1, &one, x->data + p + ldx, &ldx);
      }
    }
  }
}


static void nl_mmul (nl_Matrix *c, nl_Matrix *a, nl_Matrix *b,
    char transa, char transb, nl_Complex alpha) {
  int m, n;
  if (b->ndims == 1) { /* b vector? */
    if (a->ndims == 1) { /* rank 1 operation? */
      int ldc = LD(c, 0);
      int inca = a->stride;
      int incb = b->stride;
      m = c->dim[0]; n = c->dim[1];
      if (c->iscomplex)
        ZGERC(&m, &n, &alpha, CPX(a->data), &inca, CPX(b->data), &incb,
            CPX(c->data), &ldc);
      else
        DGER(&m, &n, (lua_Number *) &alpha, a->data, &inca, b->data, &incb,
            c->data, &ldc);
    }
    else {
      int lda = LD(a, 0);
      int incb = b->stride;
      int incc = c->stride;
      m = a->dim[0]; n = a->dim[1];
      if (c->iscomplex)
        ZGEMV(&transa, &m, &n, &alpha, CPX(a->data), &lda,
            CPX(b->data), &incb, &onec, CPX(c->data), &incc, 1);
      else
        DGEMV(&transa, &m, &n, (lua_Number *) &alpha, a->data, &lda,
            b->data, &incb, &oned, c->data, &incc, 1);
    }
  }
  else { /* a and b are matrices */
    int lda = LD(a, 0);
    int ldb = LD(b, 0);
    int ldc = LD(c, 0);
    int k = (transa == 'n' || transa == 'N') ? a->dim[1] : a->dim[0];
    m = c->dim[0]; n = c->dim[1];
    if (c->iscomplex)
      ZGEMM(&transa, &transb, &m, &n, &k, &alpha, CPX(a->data), &lda,
          CPX(b->data), &ldb, &onec, CPX(c->data), &ldc, 1, 1);
    else
      DGEMM(&transa, &transb, &m, &n, &k, (lua_Number *) &alpha,
          a->data, &lda, b->data, &ldb, &oned, c->data, &ldc, 1, 1);
  }
}


#define SPPOW(x,p) \
  do { \
    if ((x) != 0) { \
      a = fabs(x); \
      if (scale < a) { \
        s = 1 + s * pow(scale / a, (p)); \
        scale = a; \
      } \
      else \
        s += pow(a / scale, (p)); \
    } \
  } while (0)

#define SNPOW(x,p) \
  do { \
    if ((x) != 0) { \
      a = fabs(x); \
      if (scale > a) { \
        s = 1 + s * pow(scale / a, (p)); \
        scale = a; \
      } \
      else \
        s += pow(a / scale, (p)); \
    } \
  } while (0)

/* slight generalization of LAPACK's dlange */
/* `argm` is only referenced if what == 'm'/'M', and returns the argmax */
static lua_Number nl_norm (nl_Matrix *m, char what, lua_Number p, int *argm) {
  lua_Number a, norm = 0;
  int i;
  if (what == 0 && p == 0) { /* return which(m, nil, "#")? */
    int n = 0;
    if (m->section) {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        for (i = 0; i < m->size; i++, e = CSHIFT(m, i))
          if (cabs(*e) != 0) n++;
      }
      else {
        lua_Number *e = m->data;
        for (i = 0; i < m->size; i++, e = DSHIFT(m, i))
          if (*e != 0) n++;
      }
    }
    else {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        for (i = 0; i < m->size; i++, e += m->stride)
          if (cabs(*e) != 0) n++;
      }
      else {
        lua_Number *e = m->data;
        for (i = 0; i < m->size; i++, e += m->stride)
          if (*e != 0) n++;
      }
    }
    return n;
  }
  switch (what) {
    case 'm': case 'M': { /* max(abs(A(i,j))) */
      if (m->section) {
        if (m->iscomplex) {
          nl_Complex *e = CPX(m->data);
          for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
            a = cabs(*e);
            if (norm < a) {
              *argm = i + 1;
              norm = a;
            }
          }
        }
        else {
          lua_Number *e = m->data;
          for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
            a = fabs(*e);
            if (norm < a) {
              *argm = i + 1;
              norm = a;
            }
          }
        }
      }
      else {
        if (m->iscomplex) {
          i = IZAMAX(&m->size, CPX(m->data), &m->stride) - 1;
          *argm = i + 1;
          norm = cabs(CPX(m->data)[i * m->stride]);
        }
        else {
          i = IDAMAX(&m->size, m->data, &m->stride) - 1;
          *argm = i + 1;
          norm = fabs(m->data[i * m->stride]);
        }
      }
      break;
    }
    case 'o': case 'O': { /* one norm: max column sum */
      if (m->ndims == 1) { /* vector? */
        if (m->iscomplex)
          norm = DZASUM(&m->size, CPX(m->data), &m->stride);
        else
          norm = DASUM(&m->size, m->data, &m->stride);
      }
      else {
        int stride = m->stride * STEP(m, 0);
        int ld = m->stride * LD(m, 0) * STEP(m, 1);
        if (m->iscomplex) {
          nl_Complex *e = CPX(m->data);
          for (i = 0; i < m->dim[1]; i++) { /* for each column */
            a = DZASUM(&m->dim[0], e, &stride);
            if (norm < a) norm = a;
            e += ld;
          }
        }
        else {
          lua_Number *e = m->data;
          for (i = 0; i < m->dim[1]; i++) { /* for each column */
            a = DASUM(&m->dim[0], e, &stride);
            if (norm < a) norm = a;
            e += ld;
          }
        }
      }
      break;
    }
    case 'i': case 'I': { /* sup norm: max row sum */
      if (m->ndims == 1) { /* vector? */
        if (m->iscomplex) {
          i = IZAMAX(&m->size, CPX(m->data), &m->stride) - 1;
          *argm = i + 1;
          norm = cabs(CPX(m->data)[i * m->stride]);
        }
        else {
          i = IDAMAX(&m->size, m->data, &m->stride) - 1;
          *argm = i + 1;
          norm = fabs(m->data[i * m->stride]);
        }
      }
      else {
        int stride = m->stride * LD(m, 0) * STEP(m, 1);
        int ld = m->stride * STEP(m, 0);
        if (m->iscomplex) {
          nl_Complex *e = CPX(m->data);
          for (i = 0; i < m->dim[1]; i++) { /* for each row */
            a = DZASUM(&m->dim[0], e, &stride);
            if (norm < a) {
              *argm = i + 1;
              norm = a;
            }
            e += ld;
          }
        }
        else {
          lua_Number *e = m->data;
          for (i = 0; i < m->dim[1]; i++) { /* for each row */
            a = DASUM(&m->dim[0], e, &stride);
            if (norm < a) {
              *argm = i + 1;
              norm = a;
            }
            e += ld;
          }
        }
      }
      break;
    }
    default: { /* p norm */
      lua_Number scale = 1, s = 0;
      if (m->section) {
        if (m->iscomplex) {
          nl_Complex *e = CPX(m->data);
          if (p > 0) {
            for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
              SPPOW(creal(*e), p);
              SPPOW(cimag(*e), p);
            }
          }
          else {
            for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
              SNPOW(creal(*e), p);
              SNPOW(cimag(*e), p);
            }
          }
        }
        else {
          lua_Number *e = m->data;
          if (p > 0) {
            for (i = 0; i < m->size; i++, e = DSHIFT(m, i))
              SPPOW(*e, p);
          }
          else {
            for (i = 0; i < m->size; i++, e = DSHIFT(m, i))
              SNPOW(*e, p);
          }
        }
        norm = scale * pow(s, 1 / p);
      }
      else {
        if (p == 2) {
          if (m->iscomplex)
            norm = DZNRM2(&m->size, CPX(m->data), &m->stride);
          else
            norm = DNRM2(&m->size, m->data, &m->stride);
        }
        else {
          if (m->iscomplex) {
            nl_Complex *e = CPX(m->data);
            if (p > 0) {
              for (i = 0; i < m->size; i++, e += m->stride) {
                SPPOW(creal(*e), p);
                SPPOW(cimag(*e), p);
              }
            }
            else {
              for (i = 0; i < m->size; i++, e += m->stride) {
                SNPOW(creal(*e), p);
                SNPOW(cimag(*e), p);
              }
            }
          }
          else {
            lua_Number *e = m->data;
            if (p > 0) {
              for (i = 0; i < m->size; i++, e += m->stride)
                SPPOW(*e, p);
            }
            else {
              for (i = 0; i < m->size; i++, e += m->stride)
                SNPOW(*e, p);
            }
          }
          norm = scale * pow(s, 1 / p);
        }
      }
    }
  }
  return norm;
}


static int nl_chol (nl_Matrix *a, char uplo) {
  int n = a->dim[0], lda = LD(a, 0);
  int info;
  if (a->iscomplex)
    ZPOTRF(&uplo, &n, CPX(a->data), &lda, &info, 1);
  else
    DPOTRF(&uplo, &n, a->data, &lda, &info, 1);
  return info;
}


static int nl_lu (nl_Matrix *a, nl_Buffer *ipiv) {
  int m = a->dim[0], n = a->dim[1], lda = LD(a, 0);
  int info;
  if (a->iscomplex)
    ZGETRF(&m, &n, CPX(a->data), &lda, ipiv->data.bint, &info);
  else
    DGETRF(&m, &n, a->data, &lda, ipiv->data.bint, &info);
  return info;
}


/* return the reciprocal of the condition number of `m` using 1-norm
 * what = 'd' for diagonal, 'u'/'l' for triangular, 'p' for symmetric posdef
 * or 'g' for general (default). For 'p' and 'g', assumes that the argument is
 * actually a factorized (Cholesky or LU, resp.) matrix.
 * `ipiv` is only referenced when what=='g'/'G', and it should have size at
 * least `m->dim[0]` */
static lua_Number nl_rcond (lua_State *L, nl_Matrix *m, char what,
    int *ipiv, int *info) {
  lua_Number norm, rcond;
  nl_Buffer *work, *xwork;
  int n = m->dim[0], lda = LD(m, 0);
  char cnorm = '1';
  int rfinfo;
  if (what == 'd' || what == 'D') { /* diagonal? */
    /* note that 1 / ||M^-1||_1 = min_{i=1,n} { M[i][i] } */
    int i, index;
    int stride = (lda + 1) * m->stride; /* stride on diagonal */
    if (m->iscomplex) {
      nl_Complex *e = CPX(m->data);
      index = (IZAMAX(&n, CPX(m->data), &stride) - 1) * stride;
      rcond = norm = cabs(CPX(m->data)[index]);
      for (i = 0; i < n; i++, e += stride) {
        lua_Number a = cabs(*e);
        if (rcond > a) rcond = a;
        if (rcond == 0) break;
      }
    }
    else {
      lua_Number *e = m->data;
      index = (IDAMAX(&n, m->data, &stride) - 1) * stride;
      rcond = norm = fabs(m->data[index]);
      for (i = 0; i < n; i++, e += stride) {
        lua_Number a = fabs(*e);
        if (rcond > a) rcond = a;
        if (rcond == 0) break;
      }
    }
    *info = 0;
    return rcond / norm;
  }
  switch (what) {
    case 'l': case 'L': case 'u': case 'U': { /* triangular? */
      char diag = 'N';
      xwork = nl_getbuffer(L, n);
      if (what == 'l') what = 'L';
      if (what == 'u') what = 'U';
      if (m->iscomplex) {
        work = nl_getbuffer(L, 4 * n);
        ZTRCON(&cnorm, &what, &diag, &m->dim[0], CPX(m->data), &lda,
            &rcond, CPX(work->data.bnum), xwork->data.bnum, info, 1, 1, 1);
      }
      else {
        work = nl_getbuffer(L, 3 * n);
        DTRCON(&cnorm, &what, &diag, &m->dim[0], m->data, &lda,
            &rcond, work->data.bnum, xwork->data.bint, info, 1, 1, 1);
      }
      rfinfo = 0;
      break;
    }
    case 'p': case 'P': { /* symmetric posdef? */
      char uplo = 'L';
      if (m->iscomplex) {
        ZPOTRF(&uplo, &n, CPX(m->data), &lda, &rfinfo, 1);
        if (rfinfo == 0) {
          norm = nl_norm(m, 'o', 0, NULL);
          work = nl_getbuffer(L, 4 * n);
          xwork = nl_getbuffer(L, n);
          ZPOCON(&uplo, &n, CPX(m->data), &lda, &norm, &rcond,
              CPX(work->data.bnum), xwork->data.bnum, info, 1);
        }
      }
      else {
        DPOTRF(&uplo, &n, m->data, &lda, &rfinfo, 1);
        if (rfinfo == 0) {
          norm = nl_norm(m, 'o', 0, NULL);
          work = nl_getbuffer(L, 3 * n);
          xwork = nl_getbuffer(L, n);
          DPOCON(&uplo, &n, m->data, &lda, &norm, &rcond,
              work->data.bnum, xwork->data.bint, info, 1);
        }
      }
      break;
    }
    default: { /* general */
      if (m->iscomplex) {
        ZGETRF(&n, &n, CPX(m->data), &lda, ipiv, &rfinfo);
        if (rfinfo == 0) {
          norm = nl_norm(m, 'o', 0, NULL);
          work = nl_getbuffer(L, 4 * n);
          xwork = nl_getbuffer(L, 2 * n);
          ZGECON(&cnorm, &n, CPX(m->data), &lda, &norm, &rcond,
              CPX(work->data.bnum), xwork->data.bnum, info, 1);
        }
      }
      else {
        DGETRF(&n, &n, m->data, &lda, ipiv, &rfinfo);
        if (rfinfo == 0) {
          norm = nl_norm(m, 'o', 0, NULL);
          work = nl_getbuffer(L, 4 * n);
          xwork = nl_getbuffer(L, n);
          DGECON(&cnorm, &n, m->data, &lda, &norm, &rcond,
              work->data.bnum, xwork->data.bint, info, 1);
        }
      }
    }
  }
  if (rfinfo == 0) {
    nl_freebuffer(work);
    nl_freebuffer(xwork);
  }
  return rcond;
}


/* Hints:
 *  - when inverting, rcond is computed by default;
 *  - rcond and inv decompose the input automatically; if that is not desired,
 *    you have to decompose the input and use 'L'/'U' for rcond and inv;
 */

/* if rcond is not null, receives rcond of m; returns info */
static int nl_inv (lua_State *L, nl_Matrix *m, char what, lua_Number *rcond) {
  int n = m->dim[0], lda = LD(m, 0), info = 0;
  switch (what) {
    case 'd': case 'D': { /* diagonal? */
      int i, stride = (lda + 1) * m->stride;
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        for (i = 0; i < n; i++, e += stride) {
          if (cabs(*e) > 0)
            *e = 1. / *e;
          else {
            info = i + 1; break;
          }
        }
      }
      else {
        lua_Number *e = m->data;
        for (i = 0; i < n; i++, e += stride) {
          if (fabs(*e) > 0)
            *e = 1. / *e;
          else {
            info = i + 1; break;
          }
        }
      }
      if (info == 0 && rcond)
        *rcond = nl_rcond(L, m, what, NULL, &info);
      return info;
    }
    case 'l': case 'L': case 'u': case 'U': { /* triangular? */
      char diag = 'N';
      if (what == 'l') what = 'L';
      if (what == 'u') what = 'U';
      if (m->iscomplex)
        ZTRTRI(&what, &diag, &n, CPX(m->data), &lda, &info, 1, 1);
      else
        DTRTRI(&what, &diag, &n, m->data, &lda, &info, 1, 1);
      if (info == 0 && rcond)
        *rcond = nl_rcond(L, m, what, NULL, &info);
      return info;
    }
    case 'p': case 'P': { /* symmetric posdef? */
      char uplo = 'L';
      if (m->iscomplex) {
        if (rcond)
          *rcond = nl_rcond(L, m, what, NULL, &info);
        else
          ZPOTRF(&uplo, &n, CPX(m->data), &lda, &info, 1);
        if (info == 0)
          ZPOTRI(&uplo, &n, CPX(m->data), &lda, &info, 1);
      }
      else {
        if (rcond)
          *rcond = nl_rcond(L, m, what, NULL, &info);
        else
          DPOTRF(&uplo, &n, m->data, &lda, &info, 1);
        if (info == 0)
          DPOTRI(&uplo, &n, m->data, &lda, &info, 1);
      }
      return info;
    }
    default: { /* general: note that a LU decomposition is computed */
      nl_Buffer *ipiv = nl_getbuffer(L, n);
      nl_Buffer *work;
      nl_Complex qwork;
      int lwork = -1; /* query */
      if (m->iscomplex) {
        if (rcond)
          *rcond = nl_rcond(L, m, what, ipiv->data.bint, &info);
        else
          ZGETRF(&n, &n, CPX(m->data), &lda, ipiv->data.bint, &info);
        if (info == 0) {
          ZGETRI(&n, CPX(m->data), &lda, ipiv->data.bint,
              &qwork, &lwork, &info);
          lua_number2int(lwork, creal(qwork));
          work = nl_getbuffer(L, 2 * lwork);
          ZGETRI(&n, CPX(m->data), &lda, ipiv->data.bint,
              CPX(work->data.bnum), &lwork, &info);
          nl_freebuffer(work);
        }
      }
      else {
        if (rcond)
          *rcond = nl_rcond(L, m, what, ipiv->data.bint, &info);
        else
          DGETRF(&n, &n, m->data, &lda, ipiv->data.bint, &info);
        if (info == 0) {
          DGETRI(&n, m->data, &lda, ipiv->data.bint,
              (lua_Number *) &qwork, &lwork, &info);
          lua_number2int(lwork, creal(qwork));
          work = nl_getbuffer(L, lwork);
          DGETRI(&n, m->data, &lda, ipiv->data.bint,
              work->data.bnum, &lwork, &info);
          nl_freebuffer(work);
        }
      }
      nl_freebuffer(ipiv);
      return info;
    }
  }
}


/* pushes results in stack depending on `what`:
 *  [only S is returned]
 *    what = 'l': the first min(m, n) columns of U are stored in A;
 *    what = 'r': the first min(m, n) rows of V^H are stored in A;
 *    what = 'n': A is copied, only singular values are returned;
 *  [U, S, V]
 *    what = 'a': A is copied, U, S, V^H are returned
 */
static int nl_svd (lua_State *L, nl_Matrix *a, char what) {
  int m = a->dim[0], n = a->dim[1];
  int mn = (m < n) ? m : n; /* min(m, n) */
  int info, lwork = -1;
  char jobu, jobvt;
  nl_Buffer *work, *buf = NULL;
  nl_Complex qwork;
  nl_Matrix *u = NULL, *v = NULL;
  nl_Matrix *s = pushmatrix(L, 0, 1, &mn, 1, mn, NULL, NULL);
  switch (what) {
    case 'l': case 'L': jobu = 'O'; jobvt = 'N'; break;
    case 'r': case 'R': jobu = 'N'; jobvt = 'O'; break;
    case 'n': case 'N': jobu = 'N'; jobvt = 'N'; break;
    default: /* what = 'A' */
      jobu = jobvt = 'A';
      u = pushmatrix(L, a->iscomplex, 2, NULL, 1, m * m, NULL, NULL);
      u->dim[0] = u->dim[1] = m;
      lua_insert(L, -2); /* u before s */
      v = pushmatrix(L, a->iscomplex, 2, NULL, 1, n * n, NULL, NULL);
      v->dim[0] = v->dim[1] = n;
  }
  if (jobu != 'O' && jobvt != 'O') { /* copy a? */
    buf = nl_getbuffer(L, a->iscomplex ? 2 * a->size : a->size);
    setdatatovector(a, 1, 0, buf->data.bnum); /* buf := a */
  }
  if (a->iscomplex) {
    nl_Buffer *rwork = nl_getbuffer(L, 5 * mn);
    nl_Complex *data = (jobu != 'O' && jobvt != 'O')
      ? CPX(buf->data.bnum) : CPX(a->data);
    /* query lwork */
    ZGESVD(&jobu, &jobvt, &m, &n, data, &m, s->data, NULL, &m, NULL, &n,
        &qwork, &lwork, rwork->data.bnum, &info, 1, 1);
    lua_number2int(lwork, creal(qwork));
    work = nl_getbuffer(L, 2 * lwork);
    /* compute svd */
    ZGESVD(&jobu, &jobvt, &m, &n, data, &m, s->data,
        u ? CPX(u->data) : NULL, &m, v ? CPX(v->data) : NULL, &n,
        CPX(work->data.bnum), &lwork, rwork->data.bnum, &info, 1, 1);
    nl_freebuffer(rwork);
  }
  else {
    lua_Number *data = (jobu != 'O' && jobvt != 'O')
      ? buf->data.bnum : a->data;
    /* query lwork */
    DGESVD(&jobu, &jobvt, &m, &n, data, &m, s->data, NULL, &m, NULL, &n,
        (lua_Number *) &qwork, &lwork, &info, 1, 1);
    lua_number2int(lwork, creal(qwork));
    work = nl_getbuffer(L, lwork);
    /* compute svd */
    DGESVD(&jobu, &jobvt, &m, &n, data, &m, s->data,
        u ? u->data : NULL, &m, v ? v->data : NULL, &n,
        work->data.bnum, &lwork, &info, 1, 1);
  }
  nl_freebuffer(work);
  if (jobu != 'O' && jobvt != 'O') nl_freebuffer(buf);
  return info;
}


/* nl_qr: takes matrix `r`, not a section, and a pivot buffer `pvt` and
 * computes a QR decomposition in-place, returning Q and storing R. If `pvt`
 * is NULL, does not permute `r` */
static int nl_qr (lua_State *L, nl_Matrix *r, nl_Buffer *pvt) {
  int m = r->dim[0], n = r->dim[1];
  int mn = (m < n) ? m : n; /* min(m, n) */
  int info, lwork = -1;
  nl_Buffer *work, *tau;
  nl_Complex qwork;
  nl_Matrix *q = pushmatrix(L, r->iscomplex, 2, NULL, 1, m * m, NULL, NULL);
  q->dim[0] = q->dim[1] = m;
  if (r->iscomplex) {
    tau = nl_getbuffer(L, 2 * mn);
    if (pvt != NULL) { /* permute? */
      nl_Buffer *rwork = nl_getbuffer(L, 2 * n);
      /* query lwork */
      ZGEQP3(&m, &n, CPX(r->data), &m, pvt->data.bint,
          CPX(tau->data.bnum), &qwork, &lwork, rwork->data.bnum, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, 2 * lwork);
      /* compute qr */
      ZGEQP3(&m, &n, CPX(r->data), &m, pvt->data.bint,
          CPX(tau->data.bnum), CPX(work->data.bnum), &lwork,
          rwork->data.bnum, &info);
      nl_freebuffer(rwork);
    }
    else {
      /* query lwork */
      ZGEQRF(&m, &n, CPX(r->data), &m, CPX(tau->data.bnum),
          &qwork, &lwork, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, 2 * lwork);
      /* compute qr */
      ZGEQRF(&m, &n, CPX(r->data), &m, CPX(tau->data.bnum),
          CPX(work->data.bnum), &lwork, &info);
    }
    if (info == 0) {
      /* setup q and r */
      settriangtovector(r, 'L', q); /* q.L = r.L */
      settriangtoscalar(0, 'L', r); /* r.L = 0 */
      /* build q */
      if (n > m) n = m;
      ZUNGQR(&m, &n, &mn, CPX(q->data), &m, CPX(tau->data.bnum),
          CPX(work->data.bnum), &lwork, &info);
    }
  }
  else {
    tau = nl_getbuffer(L, mn);
    if (pvt != NULL) { /* permute? */
      /* query lwork */
      DGEQP3(&m, &n, r->data, &m, pvt->data.bint, tau->data.bnum,
          (lua_Number *) &qwork, &lwork, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, lwork);
      /* compute qr */
      DGEQP3(&m, &n, r->data, &m, pvt->data.bint, tau->data.bnum,
          work->data.bnum, &lwork, &info);
    }
    else {
      /* query lwork */
      DGEQRF(&m, &n, r->data, &m, tau->data.bnum,
          (lua_Number *) &qwork, &lwork, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, lwork);
      /* compute qr */
      DGEQRF(&m, &n, r->data, &m, tau->data.bnum,
          work->data.bnum, &lwork, &info);
    }
    if (info == 0) {
      /* setup q and r */
      settriangtovector(r, 'L', q); /* q.L = r.L */
      settriangtoscalar(0, 'L', r); /* r.L = 0 */
      /* build q */
      if (n > m) n = m;
      DORGQR(&m, &n, &mn, q->data, &m, tau->data.bnum,
          work->data.bnum, &lwork, &info);
    }
  }
  nl_freebuffer(work);
  nl_freebuffer(tau);
  return info;
}


/* Computes the eigenvalues and optionally the eigenvectors of matrix `a`, not
 * a section, Hermitian if `hermitian` is true.
 * If `what` == 'l'|'L', left eigenvectors are computed,
 * if `what` == 'r'|'R', right eigenvectors are computed,
 * if `what` == 'b'|'B', both eigenvectors are computed;
 * Returns `info` and pushes eigenvalue vector and eigenvector matrices. */
/* TODO [wishlist]: `nobalance` option, using Schur: gehrd, [orghr], hseqr,
 * [trevc] */
static int nl_eig (lua_State *L, nl_Matrix *m, char what, int hermitian) {
  int n = m->dim[0]; /* == m->dim[1] */
  int info, lwork = -1;
  nl_Buffer *work, *a;
  nl_Complex qwork;
  char jobvl = (what == 'l' || what == 'L' || what == 'a' || what == 'A')
    ? 'V' : 'N';
  char jobvr = (what == 'r' || what == 'R' || what == 'a' || what == 'A')
    ? 'V' : 'N';
  nl_Matrix *vl = NULL, *vr = NULL;
  nl_Matrix *w = pushmatrix(L, !hermitian, 1, &n, 1, n, NULL, NULL);
  a = nl_getbuffer(L, (m->iscomplex) ? 2 * m->size : m->size);
  setdatatovector(m, 1, 0, a->data.bnum); /* a := m */
  if (m->iscomplex) {
    if (hermitian) {
      char jobz = (jobvl == 'V' || jobvr == 'V') ? 'V' : 'N';
      char uplo = 'U';
      nl_Buffer *rwork = nl_getbuffer(L, 3 * n - 2);
      /* query lwork */
      ZHEEV(&jobz, &uplo, &n, CPX(a->data.bnum), &n, w->data,
          &qwork, &lwork, rwork->data.bnum, &info, 1, 1);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, 2 * lwork);
      /* compute eigen decomposition */
      ZHEEV(&jobz, &uplo, &n, CPX(a->data.bnum), &n, w->data,
          CPX(work->data.bnum), &lwork, rwork->data.bnum, &info, 1, 1);
      nl_freebuffer(rwork);
      if (info == 0 && jobz == 'V') { /* setup eigenvectors? */
        nl_Matrix *v = pushmatrix(L, 1, 2, m->dim, 1, m->size, NULL, NULL);
        ZCOPY(&m->size, CPX(a->data.bnum), &one, CPX(v->data), &one);
      }
    }
    else { /* general case */
      nl_Buffer *rwork = nl_getbuffer(L, 2 * n);
      if (jobvl == 'V')
        vl = pushmatrix(L, 1, 2, m->dim, 1, m->size, NULL, NULL);
      if (jobvr == 'V')
        vr = pushmatrix(L, 1, 2, m->dim, 1, m->size, NULL, NULL);
      /* query lwork */
      ZGEEV(&jobvl, &jobvr, &n, CPX(a->data.bnum), &n, CPX(w->data),
          vl ? CPX(vl->data) : NULL, &n, vr ? CPX(vr->data) : NULL, &n,
          &qwork, &lwork, rwork->data.bnum, &info, 1, 1);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, 2 * lwork);
      /* compute eigen decomposition */
      ZGEEV(&jobvl, &jobvr, &n, CPX(a->data.bnum), &n, CPX(w->data),
          vl ? CPX(vl->data) : NULL, &n, vr ? CPX(vr->data) : NULL, &n,
          CPX(work->data.bnum), &lwork, rwork->data.bnum, &info, 1, 1);
      nl_freebuffer(rwork);
    }
  }
  else {
    if (hermitian) {
      char jobz = (jobvl == 'V' || jobvr == 'V') ? 'V' : 'N';
      char uplo = 'U';
      /* query lwork */
      DSYEV(&jobz, &uplo, &n, a->data.bnum, &n, w->data,
          (lua_Number *) &qwork, &lwork, &info, 1, 1);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, lwork);
      /* compute eigen decomposition */
      DSYEV(&jobz, &uplo, &n, a->data.bnum, &n, w->data,
          work->data.bnum, &lwork, &info, 1, 1);
      if (info == 0 && jobz == 'V') { /* setup eigenvectors? */
        nl_Matrix *v = pushmatrix(L, 0, 2, m->dim, 1, m->size, NULL, NULL);
        DCOPY(&m->size, a->data.bnum, &one, v->data, &one);
      }
    }
    else { /* general case */
      nl_Buffer *wr = nl_getbuffer(L, n);
      nl_Buffer *wi = nl_getbuffer(L, n);
      if (jobvl == 'V')
        vl = pushmatrix(L, 0, 2, m->dim, 1, m->size, NULL, NULL);
      if (jobvr == 'V')
        vr = pushmatrix(L, 0, 2, m->dim, 1, m->size, NULL, NULL);
      /* query lwork */
      DGEEV(&jobvl, &jobvr, &n, a->data.bnum, &n,
          wr->data.bnum, wi->data.bnum,
          vl ? vl->data : NULL, &n, vr ? vr->data : NULL, &n,
          (lua_Number *) &qwork, &lwork, &info, 1, 1);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, lwork);
      /* compute eigen decomposition */
      DGEEV(&jobvl, &jobvr, &n, a->data.bnum, &n,
          wr->data.bnum, wi->data.bnum,
          vl ? vl->data : NULL, &n, vr ? vr->data : NULL, &n,
          work->data.bnum, &lwork, &info, 1, 1);
      if (info == 0) { /* setup eigenvectors? */
        DCOPY(&n, wr->data.bnum, &one, w->data, &two);
        DCOPY(&n, wi->data.bnum, &one, w->data + 1, &two);
      }
      nl_freebuffer(wr); nl_freebuffer(wi);
    }
  }
  nl_freebuffer(a); nl_freebuffer(work);
  return info;
}


/* nl_balance: balances `a` in-place, pushing a *similarity* diagonal matrix
 * `t` such that balanced(a) = inv(t) * a * t to compensate eigenvectors
 * -- note that `a` keeps the same eigenvalues -- computed from balanced(a),
 * do va = t * vba / || t * vba ||.
 * `what` specifies balancing operation:
 * * what = 'N': none
 * * what = 'P': permute only
 * * what = 'S': scale only
 * * what = 'B': both permute and scale [default]
 */
static int nl_balance (lua_State *L, nl_Matrix *a, char what) {
  int n = a->dim[0]; /* == a->dim[1] */
  int lda = LD(a, 0);
  int ilo, ihi, info;
  nl_Buffer *scale = nl_getbuffer(L, n);
  /* balance */
  if (a->iscomplex)
    ZGEBAL(&what, &n, CPX(a->data), &lda, &ilo, &ihi,
        scale->data.bnum, &info, 1);
  else
    DGEBAL(&what, &n, a->data, &lda, &ilo, &ihi,
        scale->data.bnum, &info, 1);
  if (info == 0) { /* compute similarity matrix? */
    nl_Matrix *t = pushmatrix(L, a->iscomplex, 2, NULL, 1, n * n, NULL, NULL); 
    int offset = ihi - ilo + 1;
    int stride = t->iscomplex ? 2 * (n + 1) : (n + 1);
    int i, j;
    lua_Number s;
    t->dim[0] = t->dim[1] = n;
    setdatatoscalar(t->iscomplex, t->size, 0, 1, 0, t->data);
    ilo--; ihi--; /* fix for zero offset */
    /* fill internal positions */
    DCOPY(&offset, scale->data.bnum, &one, t->data, &stride);
    /* swap columns */
    if (t->iscomplex) {
      for (i = n - 1; i > ihi; i--) {
        s = scale->data.bnum[i] - 1;
        lua_number2int(j, s);
        if (i != j)
          ZSWAP(&n, CPX(t->data) + i * n, &one, CPX(t->data) + j * n, &one);
      }
    }
    else {
      for (i = n - 1; i > ihi; i--) {
        s = scale->data.bnum[i] - 1;
        lua_number2int(j, s);
        if (i != j)
          DSWAP(&n, t->data + i * n, &one, t->data + j * n, &one);
      }
    }
  }
  nl_freebuffer(scale);
  return info;
}


/* nl_ls: compute lss ||a * x - b||; if s != NULL, use SVD decomposition of
 * `a` to compute, and assume that s->size = min(a->dim[0], a->dim[1]),
 * otherwise use QR decomposition. The effective rank of `a` relative to `tol`
 * is returned in `rank`. If `inplace` is true, assume that a->dim[0] >
 * a->dim[1] and return the solution in-place in `b`; otherwise push solution
 * in stack */
static int nl_ls (lua_State *L, nl_Matrix *a, nl_Matrix *b, nl_Matrix *s,
    lua_Number tol, int inplace, int *rank) {
  int m = a->dim[0];
  int n = a->dim[1];
  int nrhs = (b->ndims == 1) ? 1 : b->dim[1];
  int lda = LD(a, 0);
  int ldb = (m > n) ? m : n;
  int info, lwork = -1;
  nl_Buffer *work, *bufa, *bufb = NULL;
  nl_Complex qwork;
  /* setup buffers */
  bufa = nl_getbuffer(L, (a->iscomplex) ? 2 * a->size : a->size);
  setdatatovector(a, 1, 0, bufa->data.bnum); /* bufa := a */
  if (!inplace) {
    bufb = nl_getbuffer(L, (b->iscomplex) ? 2 * ldb * nrhs : ldb * nrhs);
    setdatatovector(b, 1, 0, bufb->data.bnum); /* bufb := b */
  }
  /* solve */
  if (a->iscomplex) {
    nl_Buffer *rwork;
    if (s != NULL) { /* svd? */
      rwork = nl_getbuffer(L, 5 * s->size); /* 5 * min(m, n) */
      /* query lwork */
      ZGELSS(&m, &n, &nrhs, CPX(bufa->data.bnum), &lda,
          inplace ? CPX(b->data) : CPX(bufb->data.bnum), &ldb,
          s->data, &tol, rank, &qwork, &lwork, rwork->data.bnum, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, 2 * lwork);
      /* compute ls */
      ZGELSS(&m, &n, &nrhs, CPX(bufa->data.bnum), &lda,
          inplace ? CPX(b->data) : CPX(bufb->data.bnum), &ldb,
          s->data, &tol, rank, CPX(work->data.bnum), &lwork,
          rwork->data.bnum, &info);
      nl_freebuffer(rwork);
    }
    else { /* qr */
      nl_Buffer *jpvt = nl_getbuffer(L, n);
      int i;
      for (i = 0; i < n; i++)
        jpvt->data.bint[i] = 0; /* free columns */
      rwork = nl_getbuffer(L, 2 * n);
      /* query lwork */
      ZGELSY(&m, &n, &nrhs, CPX(bufa->data.bnum), &lda,
          inplace ? CPX(b->data) : CPX(bufb->data.bnum), &ldb,
          jpvt->data.bint, &tol, rank, &qwork, &lwork,
          rwork->data.bnum, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, 2 * lwork);
      /* compute ls */
      ZGELSY(&m, &n, &nrhs, CPX(bufa->data.bnum), &lda,
          inplace ? CPX(b->data) : CPX(bufb->data.bnum), &ldb,
          jpvt->data.bint, &tol, rank, CPX(work->data.bnum), &lwork,
          rwork->data.bnum, &info);
      nl_freebuffer(rwork); nl_freebuffer(jpvt);
    }
  }
  else {
    if (s != NULL) { /* svd? */
      /* query lwork */
      DGELSS(&m, &n, &nrhs, bufa->data.bnum, &lda,
          inplace ? b->data : bufb->data.bnum, &ldb,
          s->data, &tol, rank, (lua_Number *) &qwork, &lwork, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, lwork);
      /* compute ls */
      DGELSS(&m, &n, &nrhs, bufa->data.bnum, &lda,
          inplace ? b->data : bufb->data.bnum, &ldb,
          s->data, &tol, rank, work->data.bnum, &lwork, &info);
    }
    else { /* qr */
      nl_Buffer *jpvt = nl_getbuffer(L, n);
      int i;
      for (i = 0; i < n; i++)
        jpvt->data.bint[i] = 0; /* free columns */
      /* query lwork */
      DGELSY(&m, &n, &nrhs, bufa->data.bnum, &lda,
          inplace ? b->data : bufb->data.bnum, &ldb,
          jpvt->data.bint, &tol, rank, (lua_Number *) &qwork, &lwork, &info);
      lua_number2int(lwork, creal(qwork));
      work = nl_getbuffer(L, lwork);
      /* compute ls */
      DGELSY(&m, &n, &nrhs, bufa->data.bnum, &lda,
          inplace ? b->data : bufb->data.bnum, &ldb,
          jpvt->data.bint, &tol, rank, work->data.bnum, &lwork, &info);
      nl_freebuffer(jpvt);
    }
  }
  if (!inplace && info == 0) { /* setup result? */
    nl_Matrix *x;
    int i;
    if (nrhs == 1)
      x = pushmatrix(L, b->iscomplex, 1, &n, 1, n, NULL, NULL);
    else {
      x = pushmatrix(L, b->iscomplex, 2, NULL, 1, nrhs * n, NULL, NULL);
      x->dim[0] = nrhs; x->dim[1] = n;
    }
    if (b->iscomplex) {
      nl_Complex *pb = CPX(bufb->data.bnum);
      nl_Complex *px = CPX(x->data);
      for (i = 0; i < nrhs; i++, pb += m, px += n) /* each column */
        zcopy_(&n, pb, &one, px, &one);
    }
    else {
      lua_Number *pb = bufb->data.bnum;
      lua_Number *px = x->data;
      for (i = 0; i < nrhs; i++, pb += m, px += n) /* each column */
        dcopy_(&n, pb, &one, px, &one);
    }
  }
  nl_freebuffer(work); nl_freebuffer(bufa);
  if (!inplace) nl_freebuffer(bufb);
  return info;
}



/* {=====================================================================
 *    Metamethods
 * ======================================================================} */

static int matrix_get (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  if (lua_isnumber(L, 2)) {
    int i, n = lua_gettop(L) - 1;
    int shift = 0, stride = m->stride, size = m->size;
    if (n > m->ndims) n = m->ndims;
    for (i = 0; i < n; i++) {
      int k = lua_tointeger(L, i + 2);
      luaL_argcheck(L, k != 0, i + 2, "null index");
      k = CIRC(k, m->dim[i]);
      shift += (k - 1) * STEP(m, i) * stride;
      stride *= LD(m, i);
      size /= m->dim[i];
    }
    if (n == m->ndims) { /* entry? */
      if (m->iscomplex)
        nl_pushcomplex(L, CPX(m->data)[shift]);
      else
        lua_pushnumber(L, m->data[shift]);
    }
    else { /* submatrix */
      lua_pushvalue(L, 1);
      lua_rawget(L, lua_upvalueindex(1)); /* push data */
      pushmatrix(L, m->iscomplex, m->ndims - n, m->dim + n,
          stride, size, (m->section) ? m->section + n : NULL,
          (m->iscomplex) ? (lua_Number *)(CPX(m->data) + shift)
            : m->data + shift);
    }
  }
  else if (lua_type(L, 2) == LUA_TUSERDATA) { /* vector? */
    int i;
    nl_Matrix *r, *v = checkmatrix(L, 2);
    lua_Number *f = v->data;
    checkrealvector(L, v, 2);
    r = pushmatrix(L, m->iscomplex, 1, &v->size, 1, v->size, NULL, NULL);
    if (m->iscomplex) {
      if (m->section) {
        for (i = 0; i < v->size; i++, f += v->stride) {
          int e = (int) (*f);
          if (e == 0)
            luaL_error(L, "null index");
          e = CIRC(e, m->size) - 1; /* zero based */
          CPX(r->data)[i] = *CSHIFT(m, e);
        }
      }
      else {
        for (i = 0; i < v->size; i++, f += v->stride) {
          int e = (int) (*f);
          if (e == 0)
            luaL_error(L, "null index");
          e = CIRC(e, m->size) - 1; /* zero based */
          CPX(r->data)[i] = CPX(m->data)[e * m->stride];
        }
      }
    }
    else {
      if (m->section) {
        for (i = 0; i < v->size; i++, f += v->stride) {
          int e = (int) (*f);
          if (e == 0)
            luaL_error(L, "null index");
          e = CIRC(e, m->size) - 1; /* zero based */
          r->data[i] = *DSHIFT(m, e);
        }
      }
      else {
        for (i = 0; i < v->size; i++, f += v->stride) {
          int e = (int) (*f);
          if (e == 0)
            luaL_error(L, "null index");
          e = CIRC(e, m->size) - 1; /* zero based */
          r->data[i] = m->data[e * m->stride];
        }
      }
    }
  }
  else { /* meta lookup? */
    lua_pushvalue(L, 2);
    lua_rawget(L, lua_upvalueindex(2)); /* class */
  }
  return 1;
}

static int matrix_set (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i;
  if (lua_gettop(L) == 2) { /* default index (full assign)? */
    lua_pushliteral(L, "_");
    lua_insert(L, -2);
  }
  /* numerical index? */
  if (lua_isnumber(L, 2)) {
    int n = lua_gettop(L) - 2;
    int shift = 0, stride = m->stride, size = m->size;
    /* setup submatrix parameters: shift, stride, size */
    if (n > m->ndims) {
      n = m->ndims;
      lua_settop(L, n + 2);
    }
    for (i = 0; i < n; i++) {
      int k = lua_tointeger(L, i + 2);
      luaL_argcheck(L, k != 0, i + 2, "null index");
      k = CIRC(k, m->dim[i]);
      shift += (k - 1) * STEP(m, i) * stride;
      stride *= LD(m, i);
      size /= m->dim[i];
    }
    /* set submatrix data */
    settoarg(L, m, n, stride, size, shift, n + 2);
  }
  else if (lua_type(L, 2) == LUA_TUSERDATA) { /* vector? */
    nl_Matrix *k = checkmatrix(L, 2);
    int ic;
    nl_Complex vc = nl_tocomplex(L, 3, &ic);
    lua_Number *f = k->data;
    checkrealvector(L, k, 2);
    /* number/complex? */ 
    if (ic) {
      if (m->iscomplex) {
        if (m->section) {
          for (i = 0; i < k->size; i++, f += k->stride) {
            int e = (int) (*f);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            *CSHIFT(m, e) = vc;
          }
        }
        else {
          for (i = 0; i < k->size; i++, f += k->stride) {
            int e = (int) (*f);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            CPX(m->data)[e * m->stride] = vc;
          }
        }
      }
      else {
        lua_Number rvc = creal(vc);
        if (m->section) {
          for (i = 0; i < k->size; i++, f += k->stride) {
            int e = (int) (*f);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            *DSHIFT(m, e) = rvc;
          }
        }
        else {
          for (i = 0; i < k->size; i++, f += k->stride) {
            int e = (int) (*f);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            m->data[e * m->stride] = rvc;
          }
        }
      }
    }
    /* vector? */
    else {
      nl_Matrix *v = checkmatrix(L, 3);
      checksection(L, v, 3);
      luaL_argcheck(L, k->size == v->size, 3, CONF_ERROR);
      if (m->iscomplex) {
        if (m->section) {
          for (i = 0; i < k->size; i++) {
            int e = (int) (k->data[i * k->stride]);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            if (v->iscomplex)
              *CSHIFT(m, e) = CPX(v->data)[i * v->stride];
            else
              *CSHIFT(m, e) = v->data[i * v->stride];
          }
        }
        else {
          for (i = 0; i < k->size; i++) {
            int e = (int) (k->data[i * k->stride]);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            if (v->iscomplex)
              CPX(m->data)[e * m->stride] = CPX(v->data)[i * v->stride];
            else
              CPX(m->data)[e * m->stride] = v->data[i * v->stride];
          }
        }
      }
      else {
        if (m->section) {
          for (i = 0; i < k->size; i++) {
            int e = (int) (k->data[i * k->stride]);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            if (v->iscomplex)
              *DSHIFT(m, e) = creal(CPX(v->data)[i * v->stride]);
            else
              *DSHIFT(m, e) = v->data[i * v->stride];
          }
        }
        else {
          for (i = 0; i < k->size; i++) {
            int e = (int) (k->data[i * k->stride]);
            luaL_argcheck(L, e != 0, 2, "null index");
            e = CIRC(e, m->size) - 1; /* zero based */
            if (v->iscomplex)
              m->data[e * m->stride] = creal(CPX(v->data)[i * v->stride]);
            else
              m->data[e * m->stride] = v->data[i * v->stride];
          }
        }
      }
    }
  }
  else if (lua_type(L, 2) == LUA_TSTRING) {
    const char what = *lua_tostring(L, 2);
    if (what == '_') { /* all entries? */
      if (lua_type(L, 3) != LUA_TUSERDATA
          || m != lua_touserdata(L, 3)) { /* not redundant? */ 
        settoarg(L, m, 0, m->stride, m->size, 0, 3);
      }
    }
    else {
      checkarray(L, m, 1);
      checksection(L, m, 1);
      if (what == 'd' || what == 'D') { /* diagonal? */
        settoarg(L, m, 0, m->stride * (m->dim[0] + 1),
            (m->dim[0] < m->dim[1]) ? m->dim[0] : m->dim[1], 0, 3);
      }
      else if (what == 'l' || what == 'L'
          || what == 'u' || what == 'U') { /* triangle? */
        int iscomplex;
        nl_Complex vc = nl_tocomplex(L, 3, &iscomplex);
        if (iscomplex) /* number/complex? */
          settriangtoscalar(vc, what, m);
        else { /* vector */
          nl_Matrix *v = checkmatrix(L, 3);
          checkarray(L, v, 3);
          checksection(L, v, 3);
          luaL_argcheck(L, m->dim[0] == v->dim[0] && m->dim[1] == v->dim[1]
              && m->iscomplex == v->iscomplex, 3, CONF_ERROR);
          settriangtovector(v, what, m);
        }
      }
      else
        luaL_error(L, "unknown option: %c", what);
    }
  }
  else
    luaL_error(L, "unexpected argument");
  lua_settop(L, 1); /* return self */
  return 1;
}


#define checktriplet(L,n,f,l,s) \
  do { \
    f = (f == 0) ? 1 : f; \
    l = (l == 0) ? n : l; \
    s = (s == 0) ? 1 : s; \
    f = CIRC(f, n); \
    l = CIRC(l, n); \
    if ((l > f && s < 0) || (l < f && s > 0)) \
      luaL_error(L, "inconsistent step argument"); \
  } while (0)

static int matrix_slice (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Matrix *s;
  int first = luaL_optinteger(L, 2, 1);
  int last = luaL_optinteger(L, 3, m->dim[0]);
  int step = luaL_optinteger(L, 4, 1);
  checktriplet(L, m->dim[0], first, last, step);
  /* create slice */
  lua_pushvalue(L, 1);
  lua_rawget(L, lua_upvalueindex(1)); /* push data */
  s = pushmatrix(L, m->iscomplex, m->ndims, m->dim, m->stride,
      m->size / m->dim[0],
      (m->ndims == 1) ? NULL : SECTION_STUB,
      DATA_STUB);
  s->size *= s->dim[0] = (last - first) / step + 1;
  first = (first - 1) * m->stride; /* shift */
  if (m->ndims == 1) /* vector? */
    s->stride *= step;
  else { /* set section */
    int i;
    for (i = 0; i < m->ndims; i++) {
      s->section[i].ld = LD(m, i);
      s->section[i].step = STEP(m, i);
    }
    s->section[0].step *= step;
    first *= STEP(m, 0);
  }
  s->data = (m->iscomplex) ? (lua_Number *)(CPX(m->data) + first)
    : m->data + first;
  return 1;
}

static int matrix_section (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Matrix *s;
  int i, first, last, step, shift, stride;
  lua_settop(L, 2);
  luaL_argcheck(L, lua_type(L, 2) == LUA_TTABLE, 2, "section table expected");
  lua_pushvalue(L, 1);
  lua_rawget(L, lua_upvalueindex(1)); /* push data */
  s = pushmatrix(L, m->iscomplex, m->ndims, NULL, m->stride, 1,
      (m->ndims == 1) ? NULL : SECTION_STUB, DATA_STUB); /* stub */
  /* setup section */
  shift = 0;
  stride = m->stride;
  for (i = 0; i < m->ndims; i++) {
    first = 1; last = m->dim[i]; step = 1;
    lua_rawgeti(L, 2, i + 1); /* t[i] */
    if (lua_type(L, -1) == LUA_TTABLE) {
      lua_rawgeti(L, -1, 1);
      first = luaL_optinteger(L, -1, first); /* t[i][1] */
      lua_rawgeti(L, -2, 2);
      last = luaL_optinteger(L, -1, last); /* t[i][2] */
      lua_rawgeti(L, -3, 3);
      step = luaL_optinteger(L, -1, step); /* t[i][3] */
      lua_pop(L, 3);
    }
    checktriplet(L, m->dim[i], first, last, step);
    s->dim[i] = (last - first) / step + 1;
    s->size *= s->dim[i]; 
    if (m->ndims == 1) { /* vector? */
      s->stride *= step;
      shift += (first - 1) * stride;
    }
    else { /* set section */
      s->section[i].ld = LD(m, i);
      s->section[i].step = STEP(m, i);
      shift += (first - 1) * s->section[i].step * stride;
      stride *= s->section[i].ld;
      s->section[i].step *= step;
    }
    lua_pop(L, 1); /* t[i] */
  }
  s->data = (m->iscomplex) ? (lua_Number *)(CPX(m->data) + shift)
    : m->data + shift;
  return 1;
}


static int matrix__tostring (lua_State *L) {
  lua_pushfstring(L, MATRIX_LIBNAME ": %p", lua_touserdata(L, 1));
  return 1;
}


static int matrix__unm (lua_State *L) {
  nl_Matrix *m = (nl_Matrix *) lua_touserdata(L, 1);
  if (nl_opmode) /* inplace? */
    lua_settop(L, 1);
  else { /* copy m */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  if (m->section) {
    int i;
    if (m->iscomplex) {
      for (i = 0; i < m->size; i++)
        *CSHIFT(m, i) *= -1;
    }
    else {
      for (i = 0; i < m->size; i++)
        *DSHIFT(m, i) *= -1;
    }
  }
  else {
    if (m->iscomplex)
      ZDSCAL(&m->size, &minusone, CPX(m->data), &m->stride);
    else
      DSCAL(&m->size, &minusone, m->data, &m->stride);
  }
  return 1;
}


static int matrix_add (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  int i, iscomplex;
  nl_Complex s = nl_tocomplex(L, 2, &iscomplex);
  int iarg = iscomplex ? 3 : 4; /* in-place arg */
  int inplace = nl_inplace(L, iarg);
  if (iscomplex) { /* just shift? */
    if (inplace)
      lua_settop(L, 1);
    else { /* copy */
      a = pushframe(L, a);
      settoarg(L, a, 0, 1, a->size, 0, 1);
    }
    if (a->iscomplex) {
      nl_Complex *e = CPX(a->data);
      if (a->section) {
        for (i = 0; i < a->size; i++, e = CSHIFT(a, i))
          *e += s;
      }
      else {
        for (i = 0; i < a->size; i++, e += a->stride)
          *e += s;
      }
    }
    else {
      lua_Number rs = creal(s);
      lua_Number *e = a->data;
      if (a->section) {
        for (i = 0; i < a->size; i++, e = DSHIFT(a, i))
          *e += rs;
      }
      else {
        for (i = 0; i < a->size; i++, e += a->stride)
          *e += rs;
      }
    }
  }
  else {
    nl_Matrix *b = checkmatrix(L, 2);
    nl_Complex alpha = nl_optcomplex(L, 3, 1);
    luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
    if (inplace)
      lua_settop(L, 2);
    else { /* copy */
      a = pushframe(L, a);
      settoarg(L, a, 0, 1, a->size, 0, 1);
    }
    if (a->section || b->section) {
      if (a->iscomplex) {
        nl_Complex *e, *f;
        for (i = 0; i < a->size; i++) {
          e = CPX(a->data) + nl_mshift(a, i);
          f = CPX(b->data) + nl_mshift(b, i);
          *e += alpha * (*f);
        }
      }
      else {
        lua_Number ra = creal(alpha);
        lua_Number *e, *f;
        for (i = 0; i < a->size; i++) {
          e = a->data + nl_mshift(a, i);
          f = b->data + nl_mshift(b, i);
          *e += ra * (*f);
        }
      }
    }
    else {
      if (a->iscomplex)
        ZAXPY(&b->size, &alpha, CPX(b->data), &b->stride,
            CPX(a->data), &a->stride);
      else {
        lua_Number ra = creal(alpha);
        DAXPY(&b->size, &ra, b->data, &b->stride, a->data, &a->stride);
      }
    }
    if (inplace) lua_pop(L, 1); /* b */
  }
  return 1; /* a */
}


static int matrix_mul (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  int i, iscomplex;
  nl_Complex s = nl_tocomplex(L, 2, &iscomplex);
  int inplace = nl_inplace(L, 3);
  if (!inplace) { /* copy? */
    a = pushframe(L, a);
    settoarg(L, a, 0, 1, a->size, 0, 1);
  }
  if (iscomplex) { /* just scale? */
    if (inplace) lua_settop(L, 1);
    if (a->section) {
      if (a->iscomplex) {
        for (i = 0; i < a->size; i++)
          *CSHIFT(a, i) *= s;
      }
      else {
        lua_Number rs = creal(s);
        for (i = 0; i < a->size; i++) {
          *DSHIFT(a, i) *= rs;
        }
      }
    }
    else {
      if (a->iscomplex)
        ZSCAL(&a->size, &s, CPX(a->data), &a->stride);
      else {
        lua_Number r = creal(s);
        DSCAL(&a->size, &r, a->data, &a->stride);
      }
    }
  }
  else {
    nl_Matrix *b = checkmatrix(L, 2);
    luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
    if (inplace) lua_settop(L, 2);
    if (a->section || b->section) {
      if (a->iscomplex) {
        nl_Complex *e, *f;
        for (i = 0; i < a->size; i++) {
          e = CPX(a->data) + nl_mshift(a, i);
          f = CPX(b->data) + nl_mshift(b, i);
          *e *= *f;
        }
      }
      else {
        lua_Number *e, *f;
        for (i = 0; i < a->size; i++) {
          e = a->data + nl_mshift(a, i);
          f = b->data + nl_mshift(b, i);
          *e *= *f;
        }
      }
    }
    else {
      if (a->iscomplex) {
        nl_Complex *e = CPX(a->data), *f = CPX(b->data);
        for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
          *e *= *f;
      }
      else {
        lua_Number *e = a->data, *f = b->data;
        for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
          *e *= *f;
      }
    }
    if (inplace) lua_pop(L, 1); /* b */
  }
  return 1; /* a */
}


static int matrix_div (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  int i, iscomplex;
  nl_Complex s = nl_tocomplex(L, 2, &iscomplex);
  int right = lua_toboolean(L, 3);
  int inplace = nl_inplace(L, 4);
  if (!inplace) { /* copy? */
    a = pushframe(L, a);
    settoarg(L, a, 0, 1, a->size, 0, 1);
  }
  if (iscomplex) { /* just scale? */
    if (inplace) lua_settop(L, 1);
    if (!right) { /* left division? */
      if (a->section) {
        if (a->iscomplex) {
          for (i = 0; i < a->size; i++)
            *CSHIFT(a, i) /= s;
        }
        else {
          lua_Number rs = creal(s);
          for (i = 0; i < a->size; i++)
            *DSHIFT(a, i) /= rs;
        }
      }
      else {
        s = 1 / s;
        if (a->iscomplex)
          ZSCAL(&a->size, &s, CPX(a->data), &a->stride);
        else {
          lua_Number r = creal(s);
          DSCAL(&a->size, &r, a->data, &a->stride);
        }
      }
    }
    else { /* right division */
      if (a->iscomplex) {
        nl_Complex *e = CPX(a->data);
        if (a->section) {
          for (i = 0; i < a->size; i++, e = CSHIFT(a, i))
            *e = s / *e;
        }
        else {
          for (i = 0; i < a->size; i++, e += a->stride)
            *e = s / *e;
        }
      }
      else {
        lua_Number rs = creal(s);
        lua_Number *e = a->data;
        if (a->section) {
          for (i = 0; i < a->size; i++, e = DSHIFT(a, i))
            *e = rs / *e;
        }
        else {
          for (i = 0; i < a->size; i++, e += a->stride)
            *e = rs / *e;
        }
      }
    }
  }
  else {
    nl_Matrix *b = checkmatrix(L, 2);
    luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
    if (inplace) lua_settop(L, 2);
    if (a->section || b->section) {
      if (a->iscomplex) {
        nl_Complex *e, *f;
        if (!right) {
          for (i = 0; i < a->size; i++) {
            e = CPX(a->data) + nl_mshift(a, i);
            f = CPX(b->data) + nl_mshift(b, i);
            *e /= *f;
          }
        }
        else {
          for (i = 0; i < a->size; i++) {
            e = CPX(a->data) + nl_mshift(a, i);
            f = CPX(b->data) + nl_mshift(b, i);
            *e = *f / (*e);
          }
        }
      }
      else {
        lua_Number *e, *f;
        if (!right) {
          for (i = 0; i < a->size; i++) {
            e = a->data + nl_mshift(a, i);
            f = b->data + nl_mshift(b, i);
            *e /= *f;
          }
        }
        else {
          for (i = 0; i < a->size; i++) {
            e = a->data + nl_mshift(a, i);
            f = b->data + nl_mshift(b, i);
            *e = *f / (*e);
          }
        }
      }
    }
    else {
      if (a->iscomplex) {
        nl_Complex *e = CPX(a->data), *f = CPX(b->data);
        if (!right) {
          for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
            *e /= *f;
        }
        else {
          for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
            *e = *f / *e;
        }
      }
      else {
        lua_Number *e = a->data, *f = b->data;
        if (!right) {
          for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
            *e /= *f;
        }
        else {
          for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
            *e = *f / *e;
        }
      }
    }
    if (inplace) lua_pop(L, 1); /* b */
  }
  return 1; /* a */
}


static int matrix_pow (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  int i, iscomplex;
  nl_Complex s = nl_tocomplex(L, 2, &iscomplex);
  int inplace = nl_inplace(L, 3);
  if (!inplace) { /* copy? */
    a = pushframe(L, a);
    settoarg(L, a, 0, 1, a->size, 0, 1);
  }
  if (iscomplex) { /* scalar power? */
    if (inplace) lua_settop(L, 1);
    if (a->iscomplex) {
      nl_Complex *e = CPX(a->data);
      if (a->section) {
        for (i = 0; i < a->size; i++, e = CSHIFT(a, i))
          *e = cpow(*e, s);
      }
      else {
        for (i = 0; i < a->size; i++, e += a->stride)
          *e = cpow(*e, s);
      }
    }
    else {
      lua_Number rs = creal(s);
      lua_Number *e = a->data;
      if (a->section) {
        for (i = 0; i < a->size; i++, e = DSHIFT(a, i))
          *e = pow(*e, rs);
      }
      else {
        for (i = 0; i < a->size; i++, e += a->stride)
          *e = pow(*e, rs);
      }
    }
  }
  else {
    nl_Matrix *b = checkmatrix(L, 2);
    luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
    if (inplace) lua_settop(L, 2);
    if (a->section || b->section) {
      if (a->iscomplex) {
        nl_Complex *e, *f;
        for (i = 0; i < a->size; i++) {
          e = CPX(a->data) + nl_mshift(a, i);
          f = CPX(b->data) + nl_mshift(b, i);
          *e = cpow(*e, *f);
        }
      }
      else {
        lua_Number *e, *f;
        for (i = 0; i < a->size; i++) {
          e = a->data + nl_mshift(a, i);
          f = b->data + nl_mshift(b, i);
          *e = pow(*e, *f);
        }
      }
    }
    else {
      if (a->iscomplex) {
        nl_Complex *e = CPX(a->data), *f = CPX(b->data);
        for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
          *e = cpow(*e, *f);
      }
      else {
        lua_Number *e = a->data, *f = b->data;
        for (i = 0; i < a->size; i++, e += a->stride, f += b->stride)
          *e = pow(*e, *f);
      }
    }
    if (inplace) lua_pop(L, 1); /* b */
  }
  return 1; /* a */
}


/* {=====================================================================
 *    Methods
 * ======================================================================} */

static int matrix_new (lua_State *L) {
  nl_Matrix *m;
  int i, e, size = 1;
  int iscomplex = 0;
  int ndims = lua_gettop(L);
  if (ndims == 0)
    luaL_error(L, "no dimensions given");
  if (lua_isboolean(L, ndims) || lua_isnil(L, ndims))
    iscomplex = lua_toboolean(L, ndims--);
  for (i = 0; i < ndims; i++) {
    e = lua_tointeger(L, i + 1);
    luaL_argcheck(L, e > 0, i + 1, "invalid dimension");
    size *= e;
  }
  m = pushmatrix(L, iscomplex, ndims, NULL, 1, size, NULL, NULL);
  for (i = 0; i < ndims; i++)
    m->dim[i] = lua_tointeger(L, i + 1);
  return 1;
}


/* element order position from index;
 * stack has index with `m->ndims` integers that get popped */
static int eorderaux (lua_State *L, nl_Matrix *m) {
  int i;
  int shift = 1; /* one based */
  int stride = 1;
  for (i = 0; i < m->ndims; i++) {
    int d = lua_tointeger(L, i - m->ndims) - 1;
    shift += d * stride;
    stride *= m->dim[i];
  }
  lua_pop(L, m->ndims);
  return shift;
}

static int matrix_eorder (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  if (lua_gettop(L) != m->ndims + 1)
    luaL_error(L, "wrong number of indices: %d expected", m->ndims);
  lua_pushinteger(L, eorderaux(L, m));
  return 1;
}


/* index from element order (zero-based): returns 0 if successful;
 * pushes indexes in stack */
static int eindexaux (lua_State *L, nl_Matrix *m, int pos) {
  int i;
  for (i = 0; i < m->ndims; i++) {
    int d = pos % m->dim[i];
    pos = (pos - d) / m->dim[i];
    lua_pushinteger(L, d + 1);
  }
  return pos;
}

static int matrix_eindex (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  eindexaux(L, m, luaL_checkinteger(L, 2) - 1);
  return m->ndims;
}


static int entriesiter (lua_State *L) {
  nl_Matrix *m = (nl_Matrix *) lua_touserdata(L, 1);
  int pos = lua_tointeger(L, 2);
  if (pos >= m->size) return 0;
  lua_pushinteger(L, pos + 1);
  if (m->iscomplex)
    nl_pushcomplex(L, CPX(m->data)[nl_mshift(m, pos)]);
  else
    lua_pushnumber(L, m->data[nl_mshift(m, pos)]);
  return 2;
}

static int entriesaux (lua_State *L) {
  nl_Matrix *m = (nl_Matrix *) lua_touserdata(L, lua_upvalueindex(1));
  int pos = lua_tointeger(L, lua_upvalueindex(2));
  if (pos >= m->size) return 0;
  eindexaux(L, m, pos);
  if (m->iscomplex)
    nl_pushcomplex(L, CPX(m->data)[nl_mshift(m, pos)]);
  else
    lua_pushnumber(L, m->data[nl_mshift(m, pos)]);
  lua_pushinteger(L, pos + 1);
  lua_replace(L, lua_upvalueindex(2));
  return m->ndims + 1; /* indices and entry */
}

static int matrix_entries (lua_State *L) {
  checkmatrix(L, 1);
  int eorder = lua_toboolean(L, 2);
  if (eorder) {
    lua_pushcfunction(L, entriesiter);
    lua_pushvalue(L, 1);
    lua_pushinteger(L, 0); /* eorder */
    return 3;
  }
  else {
    lua_settop(L, 1);
    lua_pushinteger(L, 0); /* eorder */
    lua_pushcclosure(L, entriesaux, 2);
    return 1;
  }
}


static int matrix_size (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  if (lua_type(L, 2) == LUA_TSTRING) {
    const char *s = lua_tostring(L, 2);
    if (*s == '#')
      lua_pushinteger(L, m->ndims);
    else if (*s == '*')
      lua_pushinteger(L, m->size);
    else
      lua_pushnil(L);
  }
  else {
    int i = (lua_touserdata(L, 2) == m) ? 1 /* __len? */
      : luaL_optinteger(L, 2, 1);
    if (i < 1 || i > m->ndims)
      lua_pushnil(L);
    else
      lua_pushinteger(L, m->dim[i - 1]);
  }
  return 1;
}


static int matrix_shape (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, d = luaL_optinteger(L, 2, 1);
  int iscomplex = lua_toboolean(L, 3);
  if (d < 1 || d > m->ndims) {
    lua_pushnil(L);
    return 1;
  }
  for (i = --d; i < m->ndims; i++)
    lua_pushinteger(L, m->dim[i]);
  if (iscomplex)
    lua_pushboolean(L, m->iscomplex);
  return m->ndims - d + iscomplex;
}


static int matrix_iscomplex (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_pushboolean(L, m->iscomplex);
  return 1;
}


static int matrix_real (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_settop(L, 1);
  if (m->iscomplex) {
    lua_rawget(L, lua_upvalueindex(1)); /* push data */
    pushmatrix(L, 0, m->ndims, m->dim, 2 * m->stride, m->size,
        m->section, m->data);
  }
  /* else leave self at top of stack */
  return 1;
}


static int matrix_imag (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_settop(L, 1);
  if (m->iscomplex) {
    lua_rawget(L, lua_upvalueindex(1)); /* push data */
    pushmatrix(L, 0, m->ndims, m->dim, 2 * m->stride, m->size,
        m->section, m->data + 1);
  }
  else {
    int i;
    lua_Number *data = lua_newuserdata(L, m->size * sizeof(lua_Number));
    for (i = 0; i < m->size; i++) /* unit stride */
      m->data[i] = 0;
    pushmatrix(L, 0, m->ndims, m->dim, 1, m->size, NULL, data);
  }
  return 1;
}


static int matrix_complex (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = (lua_isnoneornil(L, 2)) ? NULL : checkmatrix(L, 2);
  lua_Number *data;
  int i, stride;
  if (b)
    luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
  data = (lua_Number *) lua_newuserdata(L, a->size * sizeof(nl_Complex));
  /* copy real part */
  if (a->section) {
    nl_Complex *f = CPX(data);
    if (a->iscomplex)
      for (i = 0; i < a->size; i++, f++)
        *f = creal(*CSHIFT(a, i)) + cimag(*f) * I;
    else
      for (i = 0; i < a->size; i++, f++)
        *f = creal(*DSHIFT(a, i)) + cimag(*f) * I;
  }
  else {
    stride = (a->iscomplex) ? 2 * a->stride : a->stride;
    DCOPY(&a->size, a->data, &stride, data, &two);
  }
  /* copy imag part */
  if (!b) { /* single argument? */
    nl_Complex *e = CPX(data);
    for (i = 0; i < a->size; i++, e++)
      *e = creal(*e); /* cimag(*e) = 0 */
  }
  else {
    if (b->section) {
      nl_Complex *f = CPX(data);
      if (b->iscomplex)
        for (i = 0; i < b->size; i++, f++)
          *f = creal(*f) + cimag(*CSHIFT(b, i)) * I;
      else
        for (i = 0; i < b->size; i++, f++)
          *f = creal(*f) + cimag(*DSHIFT(b, i)) * I;
    }
    else {
      stride = (b->iscomplex) ? 2 * b->stride : b->stride;
      DCOPY(&b->size, b->data, &stride, data + 1, &two);
    }
  }
  pushmatrix(L, 1, a->ndims, a->dim, 1, a->size, NULL, data);
  return 1;
}


static int matrix_conj (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int inplace = nl_inplace(L, 2);
  if (inplace)
    lua_settop(L, 1);
  else { /* copy */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  if (m->iscomplex) {
    if (m->section) {
      int i;
      nl_Complex *c = CPX(m->data);
      for (i = 0; i < m->size; i++, c = CSHIFT(m, i))
        *c = conj(*c);
    }
    else {
      int inc = 2 * m->stride;
      DSCAL(&m->size, &minusone, m->data + 1, &inc);
    }
  }
  return 1;
}


static int matrix_copy (lua_State *L) {
  nl_Matrix *v = checkmatrix(L, 1);
  int n = lua_gettop(L);
  nl_Matrix *m = pushframe(L, v);
  if (n > 1) { /* set to x? */
    int iscomplex;
    nl_Complex x = nl_tocomplex(L, 2, &iscomplex);
    luaL_argcheck(L, iscomplex, 2, "number or complex expected");
    setdatatoscalar(m->iscomplex, m->size, x, 1, 0, m->data);
  }
  else
    setdatatovector(v, 1, 0, m->data);
  return 1;
}


static int matrix_reshape (lua_State *L) {
  nl_Matrix *m, *v = checkmatrix(L, 1);
  int i, e, size = 1;
  int ndims = lua_gettop(L) - 1;
  if (ndims == 0)
    luaL_error(L, "no dimensions given");
  checksection(L, v, 1);
  for (i = 0; i < ndims; i++) {
    e = lua_tointeger(L, i + 2);
    luaL_argcheck(L, e > 0, i + 2, "invalid dimension");
    size *= e;
  }
  if (size != v->size)
    luaL_error(L, "sizes are not consistent");
  lua_pushvalue(L, 1);
  lua_rawget(L, lua_upvalueindex(1)); /* push data */
  m = pushmatrix(L, v->iscomplex, ndims, NULL, 1, size, NULL, v->data);
  for (i = 0; i < ndims; i++)
    m->dim[i] = lua_tointeger(L, i + 2);
  return 1;
}


static int matrix_spread (lua_State *L) {
  nl_Matrix *m, *v = checkmatrix(L, 1);
  int dim = luaL_optinteger(L, 2, 1);
  int count = luaL_optinteger(L, 3, 1);
  int i, c, bsize, step, size;
  int sv = 0, sm = 0; /* shifts */
  checksection(L, v, 1);
  luaL_argcheck(L, dim >= 1 && dim <= v->ndims + 1, 2,
      "inconsistent dimension");
  luaL_argcheck(L, count > 0, 3, "positive count expected");
  size = v->size * count;
  m = pushmatrix(L, v->iscomplex, v->ndims + 1, NULL, 1, size, NULL, NULL);
  /* setup dims */
  bsize = 1;
  for (i = 0; i < dim - 1; i++)
    bsize *= m->dim[i] = v->dim[i];
  m->dim[dim - 1] = count;
  for (i = dim - 1; i < v->ndims; i++)
    m->dim[i + 1] = v->dim[i];
  /* fill matrix */
  step = v->size / bsize;
  if (v->iscomplex) {
    for (i = 0; i < step; i++, sv += bsize * v->stride)
      for (c = 0; c < count; c++, sm += bsize)
        ZCOPY(&bsize, CPX(v->data) + sv, &v->stride, CPX(m->data) + sm, &one);
  }
  else {
    for (i = 0; i < step; i++, sv += bsize * v->stride)
      for (c = 0; c < count; c++, sm += bsize)
        DCOPY(&bsize, v->data + sv, &v->stride, m->data + sm, &one);
  }
  return 1;
}


/* {=======   Logical   =======} */

static int matrix_find (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int reverse = lua_toboolean(L, 3);
  int i, found = 0;
  lua_settop(L, 2);
  /* predicate ? */
  if (lua_type(L, 2) == LUA_TFUNCTION) {
    int last;
    if (m->iscomplex) {
      nl_Complex *e = CPX(m->data);
      if (m->section) {
        for (i = 0; i < m->size && !found; i++, e = CSHIFT(m, i)) {
          lua_pushvalue(L, 2);
          nl_pushcomplex(L, *e);
          lua_call(L, 1, 1);
          last = lua_toboolean(L, -1);
          if ((!reverse && last) || (reverse && !last))
            found = i + 1;
          lua_pop(L, 1);
        }
      }
      else {
        for (i = 0; i < m->size && !found; i++, e += m->stride) {
          lua_pushvalue(L, 2);
          nl_pushcomplex(L, *e);
          lua_call(L, 1, 1);
          last = lua_toboolean(L, -1);
          if ((!reverse && last) || (reverse && !last))
            found = i + 1;
          lua_pop(L, 1);
        }
      }
    }
    else {
      lua_Number *e = m->data;
      if (m->section) {
        for (i = 0; i < m->size && !found; i++, e = DSHIFT(m, i)) {
          lua_pushvalue(L, 2);
          lua_pushnumber(L, *e);
          lua_call(L, 1, 1);
          last = lua_toboolean(L, -1);
          if ((!reverse && last) || (reverse && !last))
            found = i + 1;
          lua_pop(L, 1);
        }
      }
      else {
        for (i = 0; i < m->size && !found; i++, e += m->stride) {
          lua_pushvalue(L, 2);
          lua_pushnumber(L, *e);
          lua_call(L, 1, 1);
          last = lua_toboolean(L, -1);
          if ((!reverse && last) || (reverse && !last))
            found = i + 1;
          lua_pop(L, 1);
        }
      }
    }
  }
  else {
    nl_Matrix *v;
    int iscomplex;
    nl_Complex p = nl_tocomplex(L, 2, &iscomplex);
    /* number? */
    if (iscomplex) {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        if (m->section) {
          for (i = 0; i < m->size && !found; i++, e = CSHIFT(m, i)) {
            if ((!reverse && *e == p) || (reverse && *e != p))
              found = i + 1;
          }
        }
        else {
          for (i = 0; i < m->size && !found; i++, e += m->stride) {
            if ((!reverse && *e == p) || (reverse && *e != p))
              found = i + 1;
          }
        }
      }
      else {
        lua_Number rp = creal(p);
        lua_Number *e = m->data;
        if (m->section) {
          for (i = 0; i < m->size && !found; i++, e = DSHIFT(m, i)) {
            if ((!reverse && *e == rp) || (reverse && *e != rp))
              found = i + 1;
          }
        }
        else {
          for (i = 0; i < m->size && !found; i++, e += m->stride) {
            if ((!reverse && *e == rp) || (reverse && *e != rp))
              found = i + 1;
          }
        }
      }
    }
    /* vector? */
    else if ((v = tomatrix(L, 2)) != NULL) {
      luaL_argcheck(L, conformable(m, v), 2, CONF_ERROR);
      if (m->iscomplex) {
        nl_Complex *e, *f;
        if (m->section || v->section) {
          for (i = 0; i < m->size && !found; i++) {
            e = CPX(m->data) + nl_mshift(m, i);
            f = CPX(v->data) + nl_mshift(v, i);
            if ((!reverse && *e == *f) || (reverse && *e != *f))
              found = i + 1;
          }
        }
        else {
          e = CPX(m->data); f = CPX(v->data);
          for (i = 0; i < m->size && !found; i++) {
            if ((!reverse && *e == *f) || (reverse && *e != *f))
              found = i + 1;
            e += m->stride;
            f += v->stride;
          }
        }
      }
      else {
        lua_Number *e, *f;
        if (m->section || v->section) {
          for (i = 0; i < m->size && !found; i++) {
            e = m->data + nl_mshift(m, i);
            f = v->data + nl_mshift(v, i);
            if ((!reverse && *e == *f) || (reverse && *e != *f))
              found = i + 1;
          }
        }
        else {
          e = m->data; f = v->data;
          for (i = 0; i < m->size && !found; i++) {
            if ((!reverse && *e == *f) || (reverse && *e != *f))
              found = i + 1;
            e += m->stride;
            f += v->stride;
          }
        }
      }
    }
    /* default cond */
    else {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        if (m->section) {
          for (i = 0; i < m->size && !found; i++, e = CSHIFT(m, i)) {
            if ((!reverse && *e != 0) || (reverse && *e == 0))
              found = i + 1;
          }
        }
        else {
          for (i = 0; i < m->size && !found; i++, e += m->stride) {
            if ((!reverse && *e != 0) || (reverse && *e == 0))
              found = i + 1;
          }
        }
      }
      else {
        lua_Number *e = m->data;
        if (m->section) {
          for (i = 0; i < m->size && !found; i++, e = DSHIFT(m, i)) {
            if ((!reverse && *e != 0) || (reverse && *e == 0))
              found = i + 1;
          }
        }
        else {
          for (i = 0; i < m->size && !found; i++, e += m->stride) {
            if ((!reverse && *e != 0) || (reverse && *e == 0))
              found = i + 1;
          }
        }
      }
    }
  }
  if (found)
    lua_pushinteger(L, found);
  else
    lua_pushnil(L);
  return 1;
}


static int matrix_ifelse (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Matrix *x, *y = NULL;
  int i, iscomplex;
  nl_Complex c;
  /* setup args: x */
  c = nl_tocomplex(L, 3, &iscomplex);
  if (iscomplex) {
    x = pushframe(L, m);
    setdatatoscalar(x->iscomplex, x->size, c, 1, 0, x->data);
  }
  else {
    x = checkmatrix(L, 3);
    luaL_argcheck(L, conformable(m, x), 3, CONF_ERROR);
  }
  /* setup args: y */
  if (!lua_isnil(L, 4)) {
    c = nl_tocomplex(L, 4, &iscomplex);
    if (iscomplex) {
      y = pushframe(L, m);
      setdatatoscalar(y->iscomplex, y->size, c, 1, 0, y->data);
    }
    else {
      y = checkmatrix(L, 4);
      luaL_argcheck(L, conformable(m, y), 4, CONF_ERROR);
    }
  }
  /* predicate? */
  if (lua_type(L, 2) == LUA_TFUNCTION) {
    if (m->iscomplex) {
      nl_Complex *e = CPX(m->data);
      for (i = 0; i < m->size; ) {
        lua_pushvalue(L, 2);
        nl_pushcomplex(L, *e);
        lua_call(L, 1, 1);
        if (lua_toboolean(L, -1))
          *e = CPX(x->data)[nl_mshift(x, i)];
        else if (y)
          *e = CPX(y->data)[nl_mshift(y, i)];
        i++;
        e = (m->section) ? CSHIFT(m, i) : e + m->stride;
        lua_pop(L, 1);
      }
    }
    else {
      lua_Number *e = m->data;
      for (i = 0; i < m->size; ) {
        lua_pushvalue(L, 2);
        lua_pushnumber(L, *e);
        lua_call(L, 1, 1);
        if (lua_toboolean(L, -1))
          *e = x->data[nl_mshift(x, i)];
        else if (y)
          *e = y->data[nl_mshift(y, i)];
        i++;
        e = (m->section) ? DSHIFT(m, i) : e + m->stride;
        lua_pop(L, 1);
      }
    }
  }
  else {
    nl_Matrix *v;
    int iscomplex;
    nl_Complex p = nl_tocomplex(L, 2, &iscomplex);
    /* number? */
    if (iscomplex) {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        for (i = 0; i < m->size; ) {
          if (*e == p)
            *e = CPX(x->data)[nl_mshift(x, i)];
          else if (y)
            *e = CPX(y->data)[nl_mshift(y, i)];
          i++;
          e = (m->section) ? CSHIFT(m, i) : e + m->stride;
        }
      }
      else {
        lua_Number rp = creal(p);
        lua_Number *e = m->data;
        for (i = 0; i < m->size; ) {
          if (*e == rp)
            *e = x->data[nl_mshift(x, i)];
          else if (y)
            *e = y->data[nl_mshift(y, i)];
          i++;
          e = (m->section) ? DSHIFT(m, i) : e + m->stride;
        }
      }
    }
    /* vector? */
    else if ((v = tomatrix(L, 2)) != NULL) {
      luaL_argcheck(L, conformable(m, v), 2, CONF_ERROR);
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data), *f = CPX(v->data);
        if (m->section || v->section) {
          for (i = 0; i < m->size; i++) {
            e = CPX(m->data) + nl_mshift(m, i);
            f = CPX(v->data) + nl_mshift(v, i);
            if (*e == *f)
              *e = CPX(x->data)[nl_mshift(x, i)];
            else if (y)
              *e = CPX(y->data)[nl_mshift(y, i)];
          }
        }
        else {
          for (i = 0; i < m->size; i++) {
            if (*e == *f)
              *e = CPX(x->data)[nl_mshift(x, i)];
            else if (y)
              *e = CPX(y->data)[nl_mshift(y, i)];
            e += m->stride;
            f += v->stride;
          }
        }
      }
      else {
        lua_Number *e = m->data, *f = v->data;
        if (m->section || v->section) {
          for (i = 0; i < m->size; i++) {
            e = m->data + nl_mshift(m, i);
            f = v->data + nl_mshift(v, i);
            if (*e == *f)
              *e = x->data[nl_mshift(x, i)];
            else if (y)
              *e = y->data[nl_mshift(y, i)];
          }
        }
        else {
          for (i = 0; i < m->size; i++) {
            if (*e == *f)
              *e = x->data[nl_mshift(x, i)];
            else if (y)
              *e = y->data[nl_mshift(y, i)];
            e += m->stride;
            f += v->stride;
          }
        }
      }
    }
    else { /* default cond */
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        for (i = 0; i < m->size; ) {
          if (*e != 0)
            *e = CPX(x->data)[nl_mshift(x, i)];
          else if (y)
            *e = CPX(y->data)[nl_mshift(y, i)];
          i++;
          e = (m->section) ? CSHIFT(m, i) : e + m->stride;
        }
      }
      else {
        lua_Number *e = m->data;
        for (i = 0; i < m->size; ) {
          if (*e != 0)
            *e = x->data[nl_mshift(x, i)];
          else if (y)
            *e = y->data[nl_mshift(y, i)];
          i++;
          e = (m->section) ? DSHIFT(m, i) : e + m->stride;
        }
      }
    }
  }
  lua_settop(L, 1);
  return 1; /* self */
}


static int matrix_which (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, n = 0;
  const char what = (lua_type(L, 3) == LUA_TSTRING)
    ? *lua_tostring(L, 3) : 'k';
  nl_Buffer *buf = NULL;
  luaL_argcheck(L, what == 'k' || what == 'v' || what == '#', 3,
      "unknown option");
  if (what != '#')
    buf = nl_getbuffer(L, (m->iscomplex) ? 2 * m->size : m->size);
  /* predicate? */
  if (lua_type(L, 2) == LUA_TFUNCTION) {
    if (m->iscomplex) {
      nl_Complex *e = CPX(m->data);
      if (m->section) {
        for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
          lua_pushvalue(L, 2);
          nl_pushcomplex(L, *e);
          lua_call(L, 1, 1);
          if (lua_toboolean(L, -1)) {
            if (what == 'k')
              buf->data.bnum[n] = (lua_Number) (i + 1);
            else if (what == 'v')
              CPX(buf->data.bnum)[n] = *e;
            n++;
          }
          lua_pop(L, 1);
        }
      }
      else {
        for (i = 0; i < m->size; i++, e += m->stride) {
          lua_pushvalue(L, 2);
          nl_pushcomplex(L, *e);
          lua_call(L, 1, 1);
          if (lua_toboolean(L, -1)) {
            if (what == 'k')
              buf->data.bnum[n] = (lua_Number) (i + 1);
            else if (what == 'v')
              CPX(buf->data.bnum)[n] = *e;
            n++;
          }
          lua_pop(L, 1);
        }
      }
    }
    else {
      lua_Number *e = m->data;
      if (m->section) {
        for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
          lua_pushvalue(L, 2);
          lua_pushnumber(L, *e);
          lua_call(L, 1, 1);
          if (lua_toboolean(L, -1)) {
            if (what == 'k')
              buf->data.bnum[n] = (lua_Number) (i + 1);
            else if (what == 'v')
              buf->data.bnum[n] = *e;
            n++;
          }
          lua_pop(L, 1);
        }
      }
      else {
        for (i = 0; i < m->size; i++, e += m->stride) {
          lua_pushvalue(L, 2);
          lua_pushnumber(L, *e);
          lua_call(L, 1, 1);
          if (lua_toboolean(L, -1)) {
            if (what == 'k')
              buf->data.bnum[n] = (lua_Number) (i + 1);
            else if (what == 'v')
              buf->data.bnum[n] = *e;
            n++;
          }
          lua_pop(L, 1);
        }
      }
    }
  }
  else {
    nl_Matrix *v;
    int iscomplex;
    nl_Complex p = nl_tocomplex(L, 2, &iscomplex);
    /* number? */
    if (iscomplex) {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        if (m->section) {
          for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
            if (*e == p) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                CPX(buf->data.bnum)[n] = *e;
              n++;
            }
          }
        }
        else {
          for (i = 0; i < m->size; i++, e += m->stride) {
            if (*e == p) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                CPX(buf->data.bnum)[n] = *e;
              n++;
            }
          }
        }
      }
      else {
        lua_Number rp = creal(p);
        lua_Number *e = m->data;
        if (m->section) {
          for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
            if (*e == rp) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                buf->data.bnum[n] = *e;
              n++;
            }
          }
        }
        else {
          for (i = 0; i < m->size; i++, e += m->stride) {
            if (*e == rp) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                buf->data.bnum[n] = *e;
              n++;
            }
          }
        }
      }
    }
    /* vector? */
    else if ((v = tomatrix(L, 2)) != NULL) {
      luaL_argcheck(L, conformable(m, v), 2, CONF_ERROR);
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data), *f = CPX(v->data);
        if (m->section || v->section) {
          for (i = 0; i < m->size; i++) {
            e = CPX(m->data) + nl_mshift(m, i);
            f = CPX(v->data) + nl_mshift(v, i);
            if (*e == *f) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                CPX(buf->data.bnum)[n] = *e;
              n++;
            }
          }
        }
        else {
          for (i = 0; i < m->size; i++) {
            if (*e == *f) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                CPX(buf->data.bnum)[n] = *e;
              n++;
            }
            e += m->stride;
            f += v->stride;
          }
        }
      }
      else {
        lua_Number *e = m->data, *f = v->data;
        if (m->section || v->section) {
          for (i = 0; i < m->size; i++) {
            e = m->data + nl_mshift(m, i);
            f = v->data + nl_mshift(v, i);
            if (*e == *f) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                buf->data.bnum[n] = *e;
              n++;
            }
          }
        }
        else {
          for (i = 0; i < m->size; i++) {
            if (*e == *f) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                buf->data.bnum[n] = *e;
              n++;
            }
            e += m->stride;
            f += v->stride;
          }
        }
      }
    }
    /* default cond */
    else {
      if (m->iscomplex) {
        nl_Complex *e = CPX(m->data);
        if (m->section) {
          for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
            if (*e != 0) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                CPX(buf->data.bnum)[n] = *e;
              n++;
            }
          }
        }
        else {
          for (i = 0; i < m->size; i++, e += m->stride) {
            if (*e != 0) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                CPX(buf->data.bnum)[n] = *e;
              n++;
            }
          }
        }
      }
      else {
        lua_Number *e = m->data;
        if (m->section) {
          for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
            if (*e != 0) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                buf->data.bnum[n] = *e;
              n++;
            }
          }
        }
        else {
          for (i = 0; i < m->size; i++, e += m->stride) {
            if (*e != 0) {
              if (what == 'k')
                buf->data.bnum[n] = (lua_Number) (i + 1);
              else if (what == 'v')
                buf->data.bnum[n] = *e;
              n++;
            }
          }
        }
      }
    }
  }
  if (what == '#')
    lua_pushinteger(L, n);
  else {
    if (n > 0) {
      int iscomplex = what == 'v' && m->iscomplex;
      nl_Matrix *v = pushmatrix(L, iscomplex, 1, &n, 1, n, NULL, NULL);
      if (iscomplex)
        ZCOPY(&n, CPX(buf->data.bnum), &one, CPX(v->data), &one);
      else
        DCOPY(&n, buf->data.bnum, &one, v->data, &one);
    }
    else
      lua_pushnil(L);
    nl_freebuffer(buf);
  }
  return 1;
}


/* {=======   Functional   =======} */

static int matrix_apply (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, iscomplex, eorder = lua_toboolean(L, 3);
  int nargs = eorder ? 2 : m->ndims + 1;
  luaL_argcheck(L, lua_type(L, 2) == LUA_TFUNCTION, 2, "function expected");
  lua_settop(L, 2);
  if (m->iscomplex) {
    nl_Complex x, *e = CPX(m->data);
    if (m->section) {
      for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
        lua_pushvalue(L, 2); /* function */
        if (eorder)
          lua_pushinteger(L, i + 1);
        else
          eindexaux(L, m, i);
        nl_pushcomplex(L, *e);
        lua_call(L, nargs, 1);
        x = nl_tocomplex(L, 3, &iscomplex);
        if (iscomplex) *e = x;
        lua_pop(L, 1);
      }
    }
    else {
      for (i = 0; i < m->size; i++, e += m->stride) {
        lua_pushvalue(L, 2); /* function */
        if (eorder)
          lua_pushinteger(L, i + 1);
        else
          eindexaux(L, m, i);
        nl_pushcomplex(L, *e);
        lua_call(L, nargs, 1);
        x = nl_tocomplex(L, 3, &iscomplex);
        if (iscomplex) *e = x;
        lua_pop(L, 1);
      }
    }
  }
  else {
    lua_Number *e = m->data;
    if (m->section) {
      for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
        lua_pushvalue(L, 2); /* function */
        if (eorder)
          lua_pushinteger(L, i + 1);
        else
          eindexaux(L, m, i);
        lua_pushnumber(L, *e);
        lua_call(L, nargs, 1);
        if (lua_isnumber(L, 3)) *e = lua_tonumber(L, 3);
        lua_pop(L, 1);
      }
    }
    else {
      for (i = 0; i < m->size; i++, e += m->stride) {
        lua_pushvalue(L, 2); /* function */
        if (eorder)
          lua_pushinteger(L, i + 1);
        else
          eindexaux(L, m, i);
        lua_pushnumber(L, *e);
        lua_call(L, nargs, 1);
        if (lua_isnumber(L, 3)) *e = lua_tonumber(L, 3);
        lua_pop(L, 1);
      }
    }
  }
  lua_pop(L, 1); /* function */
  return 1; /* return self */
}

static int matrix_fold (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i;
  luaL_checktype(L, 2, LUA_TFUNCTION); /* fold(acc, element) */
  lua_settop(L, 3);
  if (m->iscomplex) {
    nl_Complex *e = CPX(m->data);
    if (m->section) {
      for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
        lua_pushvalue(L, 2); /* function */
        lua_insert(L, -2);
        nl_pushcomplex(L, *e); /* push element */
        lua_call(L, 2, 1);
      }
    }
    else {
      for (i = 0; i < m->size; i++, e += m->stride) {
        lua_pushvalue(L, 2); /* function */
        lua_insert(L, -2);
        nl_pushcomplex(L, *e); /* push element */
        lua_call(L, 2, 1);
      }
    }
  }
  else {
    lua_Number *e = m->data;
    if (m->section) {
      for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
        lua_pushvalue(L, 2); /* function */
        lua_insert(L, -2);
        lua_pushnumber(L, *e); /* push element */
        lua_call(L, 2, 1);
      }
    }
    else {
      for (i = 0; i < m->size; i++, e += m->stride) {
        lua_pushvalue(L, 2); /* function */
        lua_insert(L, -2);
        lua_pushnumber(L, *e); /* push element */
        lua_call(L, 2, 1);
      }
    }
  }
  return 1;
} 


static int matrix_map (lua_State *L) {
  int i, n = lua_gettop(L);
  nl_Matrix *m = checkmatrix(L, 1);
  luaL_argcheck(L, lua_type(L, n), n, "function expected");
  if (n == 2) { /* single matrix? */
    if (m->iscomplex) {
      int iscomplex;
      nl_Complex r, *e = CPX(m->data);
      if (m->section) {
        for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) {
          lua_pushvalue(L, 2); /* function */
          nl_pushcomplex(L, *e);
          lua_call(L, 1, 1);
          r = nl_tocomplex(L, 3, &iscomplex);
          if (iscomplex) *e = r;
          lua_pop(L, 1);
        }
      }
      else {
        for (i = 0; i < m->size; i++, e += m->stride) {
          lua_pushvalue(L, 2); /* function */
          nl_pushcomplex(L, *e);
          lua_call(L, 1, 1);
          r = nl_tocomplex(L, 3, &iscomplex);
          if (iscomplex) *e = r;
          lua_pop(L, 1);
        }
      }
    }
    else {
      lua_Number *e = m->data;
      if (m->section) {
        for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) {
          lua_pushvalue(L, 2); /* function */
          lua_pushnumber(L, *e);
          lua_call(L, 1, 1);
          if (lua_isnumber(L, -1)) *e = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
      }
      else {
        for (i = 0; i < m->size; i++, e += m->stride) {
          lua_pushvalue(L, 2); /* function */
          lua_pushnumber(L, *e);
          lua_call(L, 1, 1);
          if (lua_isnumber(L, -1)) *e = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
      }
    }
  }
  else {
    int j, anysections = 0;
    nl_Matrix *v;
    /* check args */
    for (j = 2; j < n; j++) {
      luaL_argcheck(L, conformable(m, checkmatrix(L, j)), j, CONF_ERROR);
      anysections = anysections || (m->section != NULL);
    }
    if (m->iscomplex) {
      int iscomplex;
      nl_Complex r, *e = CPX(m->data);
      if (anysections) {
        for (i = 0; i < m->size; i++) {
          lua_pushvalue(L, n); /* function */
          e = CPX(m->data) + nl_mshift(m, i);
          nl_pushcomplex(L, *e);
          for (j = 2; j < n; j++) {
            v = (nl_Matrix *) lua_touserdata(L, j);
            nl_pushcomplex(L, CPX(v->data)[nl_mshift(v, i)]);
          }
          lua_call(L, n - 1, 1);
          r = nl_tocomplex(L, n + 1, &iscomplex);
          if (iscomplex) *e = r;
          lua_pop(L, 1);
        }
      }
      else { /* no sections */
        for (i = 0; i < m->size; i++, e += m->stride) {
          lua_pushvalue(L, n); /* function */
          nl_pushcomplex(L, *e);
          for (j = 2; j < n; j++) {
            v = (nl_Matrix *) lua_touserdata(L, j);
            nl_pushcomplex(L, CPX(v->data)[i * v->stride]);
          }
          lua_call(L, n - 1, 1);
          r = nl_tocomplex(L, n + 1, &iscomplex);
          if (iscomplex) *e = r;
          lua_pop(L, 1);
        }
      }
    }
    else {
      lua_Number *e = m->data;
      if (anysections) {
        for (i = 0; i < m->size; i++) {
          lua_pushvalue(L, n); /* function */
          e = m->data + nl_mshift(m, i);
          lua_pushnumber(L, *e);
          for (j = 2; j < n; j++) {
            v = (nl_Matrix *) lua_touserdata(L, j);
            lua_pushnumber(L, v->data[nl_mshift(v, i)]);
          }
          lua_call(L, n - 1, 1);
          if (lua_isnumber(L, -1)) *e = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
      }
      else { /* no sections */
        for (i = 0; i < m->size; i++, e += m->stride) {
          lua_pushvalue(L, n); /* function */
          lua_pushnumber(L, *e);
          for (j = 2; j < n; j++) {
            v = (nl_Matrix *) lua_touserdata(L, j);
            lua_pushnumber(L, v->data[i * v->stride]);
          }
          lua_call(L, n - 1, 1);
          if (lua_isnumber(L, -1)) *e = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
      }
    }
  }
  lua_settop(L, 1); /* self */
  return 1;
}


/* {=======   Statistical   =======} */

static int matrix_sum (lua_State *L) { /* linear fold */
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Complex alpha = nl_optcomplex(L, 2, 1);
  nl_Complex x = nl_optcomplex(L, 3, 0);
  int i;
  if (m->iscomplex) {
    nl_Complex *e = CPX(m->data);
    if (m->section) {
      if (alpha == 1.0) { /* plain sum? */
        for (i = 0; i < m->size; i++, e = CSHIFT(m, i))
          x += *e;
      }
      else {
        for (i = 0; i < m->size; i++, e = CSHIFT(m, i))
          x = x * alpha + (*e);
      }
    }
    else {
      if (alpha == 1.0) { /* plain sum? */
        for (i = 0; i < m->size; i++, e += m->stride)
          x += *e;
      }
      else {
        for (i = 0; i < m->size; i++, e += m->stride)
          x = x * alpha + (*e);
      }
    }
    nl_pushcomplex(L, x);
  }
  else {
    lua_Number ra = creal(alpha), rx = creal(x);
    lua_Number *e = m->data;
    if (m->section) {
      if (ra == 1.0) { /* plain sum? */
        for (i = 0; i < m->size; i++, e = DSHIFT(m, i))
          rx += *e;
      }
      else {
        for (i = 0; i < m->size; i++, e = DSHIFT(m, i))
          rx = rx * ra + (*e);
      }
    }
    else {
      if (ra == 1.0) { /* plain sum? */
        for (i = 0; i < m->size; i++, e += m->stride)
          rx += *e;
      }
      else {
        for (i = 0; i < m->size; i++, e += m->stride)
          rx = rx * ra + (*e);
      }
    }
    lua_pushnumber(L, rx);
  }
  return 1;
}


static int matrix_min (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, im = 0;
  if (m->iscomplex) {
    nl_Complex x, mm = CPX(m->data)[0];
    if (m->section) {
      for (i = 1; i < m->size; i++) {
        x = *CSHIFT(m, i);
        if (clt(x, mm)) { /* new min? */
          mm = x; im = i;
        }
      }
    }
    else {
      for (i = 1; i < m->size; i++) {
        x = CPX(m->data)[i * m->stride];
        if (clt(x, mm)) { /* new min? */
          mm = x; im = i;
        }
      }
    }
    nl_pushcomplex(L, mm);
  }
  else {
    lua_Number x, mm = m->data[0];
    if (m->section) {
      for (i = 1; i < m->size; i++) {
        x = *DSHIFT(m, i);
        if (x < mm) { /* new min? */
          mm = x; im = i;
        }
      }
    }
    else {
      for (i = 1; i < m->size; i++) {
        x = m->data[i * m->stride];
        if (x < mm) { /* new min? */
          mm = x; im = i;
        }
      }
    }
    lua_pushnumber(L, mm);
  }
  lua_pushinteger(L, im + 1);
  return 2;
}


static int matrix_max (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, im = 0;
  if (m->iscomplex) {
    nl_Complex x, mm = CPX(m->data)[0];
    if (m->section) {
      for (i = 1; i < m->size; i++) {
        x = *CSHIFT(m, i);
        if (cgt(x, mm)) { /* new max? */
          mm = x; im = i;
        }
      }
    }
    else {
      for (i = 1; i < m->size; i++) {
        x = CPX(m->data)[i * m->stride];
        if (cgt(x, mm)) { /* new max? */
          mm = x; im = i;
        }
      }
    }
    nl_pushcomplex(L, mm);
  }
  else {
    lua_Number x, mm = m->data[0];
    if (m->section) {
      for (i = 1; i < m->size; i++) {
        x = *DSHIFT(m, i);
        if (x > mm) { /* new max? */
          mm = x; im = i;
        }
      }
    }
    else {
      for (i = 1; i < m->size; i++) {
        x = m->data[i * m->stride];
        if (x > mm) { /* new max? */
          mm = x; im = i;
        }
      }
    }
    lua_pushnumber(L, mm);
  }
  lua_pushinteger(L, im + 1);
  return 2;
}

int sort1c (nl_Matrix *mx);
int sort1d (nl_Matrix *mx);
int sort2c (nl_Matrix *mx, nl_Matrix *my);
int sort2d (nl_Matrix *mx, nl_Matrix *my);
static int matrix_sort (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int decreasing = lua_toboolean(L, 2);
  int retindex = lua_toboolean(L, 3); /* return index? */
  lua_settop(L, 1);
  checksection(L, m, 1);
  if (m->iscomplex) {
    if (decreasing) /* adjust m? */
      ZDSCAL(&m->size, &minusone, CPX(m->data), &m->stride);
    if (retindex) {
      int i;
      nl_Matrix *index;
      index = pushmatrix(L, 0, 1, &m->size, 1, m->size, NULL, NULL);
      for (i = 0; i < index->size; i++)
        index->data[i] = i + 1;
      sort2c(m, index);
    }
    else sort1c(m);
    if (decreasing) /* fix m? */
      ZDSCAL(&m->size, &minusone, CPX(m->data), &m->stride);
  }
  else {
    if (decreasing) /* adjust m? */
      DSCAL(&m->size, &minusone, m->data, &m->stride);
    if (retindex) {
      int i;
      nl_Matrix *index;
      index = pushmatrix(L, 0, 1, &m->size, 1, m->size, NULL, NULL);
      for (i = 0; i < index->size; i++)
        index->data[i] = i + 1;
      sort2d(m, index);
    }
    else sort1d(m);
    if (decreasing) /* fix m? */
      DSCAL(&m->size, &minusone, m->data, &m->stride);
  }
  return 1;
}


/* {=======   Math Functions   =======} */

#define F(name) \
  static int matrix_ ## name (lua_State *L) { \
    nl_Matrix *m = checkmatrix(L, 1); \
    int i, inplace = nl_inplace(L, 2); \
    if (inplace) \
      lua_settop(L, 1); \
    else { /* copy */ \
      m = pushframe(L, m); \
      settoarg(L, m, 0, 1, m->size, 0, 1); \
    } \
    if (m->iscomplex) { \
      nl_Complex *e = CPX(m->data); \
      if (m->section) { \
        for (i = 0; i < m->size; i++, e = CSHIFT(m, i)) \
          *e = c ## name(*e); \
      } \
      else { \
        for (i = 0; i < m->size; i++, e += m->stride) \
          *e = c ## name(*e); \
      } \
    } \
    else { \
      lua_Number *e = m->data; \
      if (m->section) { \
        for (i = 0; i < m->size; i++, e = DSHIFT(m, i)) \
          *e = name(*e); \
      } \
      else { \
        for (i = 0; i < m->size; i++, e += m->stride) \
          *e = name(*e); \
      } \
    } \
    return 1; /* self */ \
  }

#define abs fabs
F(abs)
#undef abs
F(acos)
F(acosh)
F(asin)
F(asinh)
F(atan)
F(atanh)
F(cos)
F(cosh)
F(exp)
F(log)
F(sin)
F(sinh)
F(sqrt)
F(tan)
F(tanh)


static int matrix_linspace (lua_State *L) {
  nl_Complex a = nl_checkcomplex(L, 1);
  nl_Complex b = nl_checkcomplex(L, 2);
  int i, n, iscomplex;
  iscomplex = (cimag(a) != 0 || cimag(b) != 0);
  if (iscomplex) {
    nl_Complex *data;
    nl_Complex s = b - a;
    n = luaL_optinteger(L, 3, cabs(s) + 1);
    luaL_argcheck(L, n > 0, 3, "number of steps is non-positive");
    lua_settop(L, 0);
    data = lua_newuserdata(L, n * sizeof(nl_Complex));
    s /= n - 1;
    data[0] = a;
    for (i = 1; i < n; i++)
      data[i] = data[i - 1] + s;
    pushmatrix(L, 1, 1, &n, 1, n, NULL, (lua_Number *) data);
  }
  else {
    lua_Number *data;
    lua_Number s = creal(b) - creal(a);
    n = luaL_optinteger(L, 3, fabs(s) + 1);
    luaL_argcheck(L, n > 0, 3, "number of steps is non-positive");
    lua_settop(L, 0);
    data = lua_newuserdata(L, n * sizeof(lua_Number));
    s /= n - 1;
    data[0] = creal(a);
    for (i = 1; i < n; i++)
      data[i] = data[i - 1] + s;
    pushmatrix(L, 0, 1, &n, 1, n, NULL, data);
  }
  return 1;
}


static int matrix_dot (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  int trans = lua_toboolean(L, 3);
  luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
  if (a->section || b->section) {
    int i;
    if (a->iscomplex) {
      nl_Complex t;
      nl_Complex *dot = nl_pushcomplex(L, 0);
      for (i = 0; i < a->size; i++) {
        t = CPX(a->data)[nl_mshift(a, i)];
        if (trans) t = conj(t);
        *dot += t * CPX(b->data)[nl_mshift(b, i)];
      }
    }
    else {
      lua_Number dot = 0;
      for (i = 0; i < a->size; i++)
        dot += a->data[nl_mshift(a, i)] * b->data[nl_mshift(b, i)];
      lua_pushnumber(L, dot);
    }
  }
  else {
    if (a->iscomplex) {
      nl_Complex *dot = nl_newcomplex(L);
      if (trans)
        ZDOTU(dot, &a->size, CPX(a->data), &a->stride,
            CPX(b->data), &b->stride);
      else
        ZDOTC(dot, &a->size, CPX(a->data), &a->stride,
            CPX(b->data), &b->stride);
    }
    else
      lua_pushnumber(L, DDOT(&a->size, a->data, &a->stride,
            b->data, &b->stride));
  }
  return 1;
}

static int matrix_cross (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  int inplace = nl_inplace(L, 3);
  luaL_argcheck(L, a->ndims == 1 && a->dim[0] == 3, 1,
      "ternary vector expected");
  luaL_argcheck(L, b->ndims == 1 && b->dim[0] == 3, 2,
      "ternary vector expected");
  if (a->iscomplex != b->iscomplex)
    luaL_error(L, CONF_ERROR);
  if (inplace)
    lua_settop(L, 2);
  else { /* copy */
    a = pushframe(L, a);
    settoarg(L, a, 0, 1, a->size, 0, 1);
  }
  if (a->iscomplex) {
    nl_Complex *a1 = CPX(a->data);
    nl_Complex *a2 = a1 + nl_mshift(a, 1);
    nl_Complex *a3 = a1 + nl_mshift(a, 2);
    nl_Complex *b1 = CPX(b->data);
    nl_Complex *b2 = b1 + nl_mshift(b, 1);
    nl_Complex *b3 = b1 + nl_mshift(b, 2);
    nl_Complex u, v, t;
    u = (*a2) * (*b3) - (*a3) * (*b2);
    v = (*a3) * (*b1) - (*a1) * (*b3);
    t = (*a1) * (*b2) - (*a2) * (*b1);
    *a1 = u; *a2 = v; *a3 = t;
  }
  else {
    lua_Number *a1 = a->data;
    lua_Number *a2 = a1 + nl_mshift(a, 1);
    lua_Number *a3 = a1 + nl_mshift(a, 2);
    lua_Number *b1 = b->data;
    lua_Number *b2 = b1 + nl_mshift(b, 1);
    lua_Number *b3 = b1 + nl_mshift(b, 2);
    lua_Number u, v, t;
    u = (*a2) * (*b3) - (*a3) * (*b2);
    v = (*a3) * (*b1) - (*a1) * (*b3);
    t = (*a1) * (*b2) - (*a2) * (*b1);
    *a1 = u; *a2 = v; *a3 = t;
  }
  if (inplace) lua_pop(L, 1); /* b */
  return 1; /* a */
}


/* {=======   Two-dimensional matrices (arrays)   =======} */

static int matrix_col (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int p, k = luaL_checkinteger(L, 2);
  checkarray(L, m, 1);
  check2dsection(L, m, 1);
  luaL_argcheck(L, k != 0, 2, "null index");
  k = CIRC(k, m->dim[1]);
  p = (k - 1) * (m->section ? LD(m, 0) : m->dim[0] * m->stride); /* offset */
  lua_pushvalue(L, 1); /* matrix */
  lua_rawget(L, lua_upvalueindex(1)); /* push data */
  /* push new column */
  pushmatrix(L, m->iscomplex, 1, &m->dim[0], m->stride, m->dim[0],
      NULL, m->iscomplex ? (lua_Number *) (CPX(m->data) + p) : m->data + p);
  return 1;
}


static int matrix_transpose (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int hermitian = lua_toboolean(L, 2);
  check2dmatrix(L, m, 1);
  nl_transpose(L, m, hermitian);
  return 1;
}


static int matrix_diag (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int dshift = luaL_optinteger(L, 2, 0); /* diagonal shift */
  check2dmatrix(L, m, 1);
  if (m->ndims == 1) { /* vector? */
    int n = m->size + abs(dshift);
    int shift = (dshift <= 0) ? -dshift : dshift * n;
    nl_Matrix *d = pushmatrix(L, m->iscomplex, 2, NULL, 1, n * n, NULL, NULL);
    d->dim[0] = d->dim[1] = n;
    setdatatoscalar(d->iscomplex, d->size, 0, 1, 0, d->data); /* d := 0 */
    setdatatovector(m, n + 1, shift, d->data); /* diag(d, dshift) := m */
  }
  else { /* matrix */
    int shift, stride;
    int n = (m->dim[0] < m->dim[1]) ? m->dim[0] : m->dim[1]; /* min(m, n) */
    n -= abs(dshift);
    luaL_argcheck(L, n > 0, 2, "diagonal shift is larger than min dimension");
    lua_pushvalue(L, 1);
    lua_rawget(L, lua_upvalueindex(1)); /* push data */
    if (m->section) {
      shift = (dshift <= 0) ? -dshift * m->section[0].step
        : dshift * m->section[1].step * m->section[0].ld;
      stride = m->section[0].step + m->section[1].step * m->section[0].ld;
    }
    else {
      shift = (dshift <= 0) ? -dshift : dshift * m->dim[0];
      stride = m->dim[0] + 1;
    }
    pushmatrix(L, m->iscomplex, 1, &n, m->stride * stride, n,
        NULL, m->data + shift);
  }
  return 1;
}


/* assumes there are `n` matrices on stack */
static int colcat (lua_State *L, int n) {
  nl_Matrix *m, *v;
  int i, s, c, dim[2];
  int iscomplex = 0;
  /* test if at most two-dimensional matrices */
  for (i = 0; i < n; i++) {
    v = checkmatrix(L, i + 1);
    check2dmatrix(L, v, i + 1);
    if (i == 0) {
      dim[0] = v->dim[0];
      dim[1] = 0;
      iscomplex = v->iscomplex;
    }
    else
      luaL_argcheck(L, dim[0] == v->dim[0] && iscomplex == v->iscomplex,
          i + 1, CONF_ERROR);
    dim[1] += (v->ndims == 1) ? 1 : v->dim[1];
  }
  /* fill new matrix */
  m = pushmatrix(L, iscomplex, 2, dim, 1, dim[0] * dim[1], NULL, NULL);
  s = 0; /* shift: partial # columns */
  for (i = 0; i < n; i++) {
    v = (nl_Matrix *) lua_touserdata(L, i + 1);
    c = (v->ndims == 1) ? 1 : v->dim[1]; /* # columns */
    setdatatovector(v, 1, dim[0] * s, m->data);
    s += c;
  }
  return 1;
}

/* assumes there are `n` matrices on stack */
static int rowcat (lua_State *L, int n) {
  nl_Matrix *m, *v;
  int i, s, dim[2];
  int iscomplex = 0;
  /* test if at most two-dimensional matrices */
  for (i = 0; i < n; i++) {
    v = checkmatrix(L, i + 1);
    check2dmatrix(L, v, i + 1);
    if (i == 0) {
      dim[0] = 0;
      dim[1] = (v->ndims == 1) ? v->dim[0] : v->dim[1];
      iscomplex = v->iscomplex;
    }
    else
      luaL_argcheck(L, (dim[1] == ((v->ndims == 1) ? v->dim[0] : v->dim[1]))
            && iscomplex == v->iscomplex,
          i + 1, CONF_ERROR);
    dim[0] += (v->ndims == 1) ? 1 : v->dim[0];
  }
  /* fill new matrix */
  m = pushmatrix(L, iscomplex, 2, dim, 1, dim[0] * dim[1], NULL, NULL);
  s = 0; /* shift: partial # rows */
  for (i = 0; i < n; i++) {
    v = (nl_Matrix *) lua_touserdata(L, i + 1);
    if (v->ndims == 1) /* vector? */
      setdatatovector(v, dim[0], s++, m->data); /* direct, row-wise */
    else { /* assign column-wise */
      int j;
      int stride = v->stride * STEP(v, 0);
      int shift = LD(v, 0) * stride;
      if (v->iscomplex) {
        nl_Complex *md = CPX(m->data) + s;
        nl_Complex *vd = CPX(v->data);
        for (j = 0; j < dim[1]; j++) { /* each column */
          ZCOPY(&v->dim[0], vd, &stride, md, &one);
          md += dim[0];
          vd += shift;
        }
      }
      else {
        lua_Number *md = m->data + s;
        lua_Number *vd = v->data;
        for (j = 0; j < dim[1]; j++) { /* each column */
          DCOPY(&v->dim[0], vd, &stride, md, &one);
          md += dim[0];
          vd += shift;
        }
      }
      s += v->dim[0];
    }
  }
  return 1;
}

static int matrix_concat (lua_State *L) {
  int n = lua_gettop(L);
  return (lua_isboolean(L, n) && lua_toboolean(L, n--)) /* column cat? */
    ? colcat(L, n) : rowcat(L, n);
}


static int matrix_c (lua_State *L) {
  int i, ic, iscomplex = 0, n = lua_gettop(L);
  int d = 0; /* new dimension */
  nl_Complex vc;
  nl_Matrix *v, *m;
  /* check args */
  for (i = 0; i < n; i++) {
    vc = nl_tocomplex(L, i + 1, &ic);
    if (ic) { /* number/complex? */
      if (i == 0)
        iscomplex = (cimag(vc) != 0);
      else
        luaL_argcheck(L, iscomplex || (!iscomplex && (cimag(vc) == 0)),
            i + 1, CONF_ERROR);
      d++;
    }
    else {
      v = checkmatrix(L, i + 1);
      checkvector(L, v, i + 1);
      if (i == 0)
        iscomplex = v->iscomplex;
      else
        luaL_argcheck(L, iscomplex == v->iscomplex, i + 1, CONF_ERROR);
      d += v->size;
    }
  }
  /* fill result */
  m = pushmatrix(L, iscomplex, 1, &d, 1, d, NULL, NULL);
  d = 0;
  for (i = 0; i < n; i++) {
    vc = nl_tocomplex(L, i + 1, &ic);
    if (ic) { /* number/complex? */
      if (m->iscomplex)
        CPX(m->data)[d++] = vc;
      else
        m->data[d++] = creal(vc);
    }
    else {
      v = (nl_Matrix *) lua_touserdata(L, i + 1);
      settoarg(L, m, 0, 1, v->size, d, i + 1);
      d += v->size;
    }
  }
  return 1;
}


/* {=======   Linear Algebra   =======} */

#define checkinfo(L,name) \
  do { \
    if (info < 0) { \
      lua_pushnil(L); \
      lua_pushfstring(L, "illegal argument to " name ": info = %d", info); \
      return 2; \
    } \
  } while (0)

static int matrix_swap (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  checkvector(L, a, 1);
  checkvector(L, b, 2);
  luaL_argcheck(L, conformable(a, b), 2, CONF_ERROR);
  if (a->iscomplex)
    ZSWAP(&a->size, CPX(a->data), &a->stride, CPX(b->data), &b->stride);
  else
    DSWAP(&a->size, a->data, &a->stride, b->data, &b->stride);
  return 0;
}

static int matrix_pivot (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *p = checkmatrix(L, 2);
  int colpivot = lua_toboolean(L, 3);
  int inplace = nl_inplace(L, 4);
  int n, l, stride, step; /* max index, length, stride, step */
  int i, j;
  lua_Number *e = p->data;
  /* check args */
  check2dmatrix(L, a, 1);
  checkvector(L, p, 2);
  if (inplace)
    check2dsection(L, a, 1);
  else { /* copy */
    a = pushframe(L, a);
    settoarg(L, a, 0, 1, a->size, 0, 1);
  }
  if (colpivot && a->ndims == 2) {
    n = a->dim[1];
    l = a->dim[0];
    stride = a->stride;
    step = LD(a, 0);
  }
  else {
    n = a->dim[0];
    l = a->dim[1];
    stride = LD(a, 0);
    step = a->stride;
  }
  /* swap */
  if (a->iscomplex) {
    for (i = 0; i < p->size; i++, e += p->stride) {
      lua_number2int(j, *e); j--;
      if (i != j && i < n && j < n) { /* swap? */
        if (a->ndims == 1) { /* vector? */
          nl_Complex c = CPX(a->data)[i * a->stride];
          CPX(a->data)[i * a->stride] = CPX(a->data)[j * a->stride];
          CPX(a->data)[j * a->stride] = c;
        }
        else
          ZSWAP(&l, CPX(a->data) + i * step, &stride,
              CPX(a->data) + j * step, &stride);
      }
    }
  }
  else {
    for (i = 0; i < p->size; i++, e += p->stride) {
      lua_number2int(j, *e); j--;
      if (i != j && i < n && j < n) { /* swap? */
        if (a->ndims == 1) { /* vector? */
          lua_Number c = a->data[i * a->stride];
          a->data[i * a->stride] = a->data[j * a->stride];
          a->data[j * a->stride] = c;
        }
        else
          DSWAP(&l, a->data + i * step, &stride,
              a->data + j * step, &stride);
      }
    }
  }
  if (inplace) lua_settop(L, 1); /* just self */
  return 1;
}


static int matrix_trmul (lua_State *L) {
  nl_Matrix *x = checkmatrix(L, 1);
  nl_Matrix *a = checkmatrix(L, 2);
  const char *u = luaL_optstring(L, 3, "L"); /* lower */
  int invert = lua_toboolean(L, 4);
  const char *t = luaL_optstring(L, 5, "N"); /* no transpose */
  const char *s = luaL_optstring(L, 6, "L"); /* left */
  char uplo = u[0];
  char trans = t[0];
  char side = s[0];
  /* check args */
  check2dmatrix(L, x, 1);
  check2dsection(L, x, 1);
  checksquare(L, a, 2);
  check2dsection(L, a, 2);
  luaL_argcheck(L, uplo == 'u' || uplo == 'U' || uplo == 'l' || uplo == 'L',
      3, "unknown triangle option");
  luaL_argcheck(L, trans == 'n' || trans == 'N'
      || trans == 't' || trans == 'T' || trans == 'c' || trans == 'C',
      5, "unknown transpose option");
  luaL_argcheck(L, side == 'l' || side == 'L' || side == 'r' || side == 'R',
      6, "unknown side option");
  if (x->iscomplex != a->iscomplex)
    luaL_error(L, DOMC_ERROR);
  if (x->ndims == 1) {
    if (x->size != a->dim[0])
      luaL_error(L, CONF_ERROR);
  }
  else {
    if ((side == 'r' || side == 'R') && x->dim[1] != a->dim[0])
      luaL_error(L, CONF_ERROR);
    if ((side == 'l' || side == 'L') && x->dim[0] != a->dim[0])
      luaL_error(L, CONF_ERROR);
  }
  nl_trmul(x, a, uplo, invert, trans, side);
  lua_settop(L, 1);
  return 1; /* self */
}


static int matrix_hemul (lua_State *L) {
  nl_Matrix *x = checkmatrix(L, 1);
  nl_Matrix *a = checkmatrix(L, 2);
  int inner = lua_toboolean(L, 3);
  const char *w = luaL_optstring(L, 4, "F");
  lua_Number alpha = luaL_optnumber(L, 5, 1.);
  char what = w[0];
  /* check args */
  checksquare(L, x, 1);
  check2dsection(L, x, 1);
  check2dmatrix(L, a, 2);
  check2dsection(L, a, 2);
  if (x->iscomplex != a->iscomplex)
    luaL_error(L, DOMC_ERROR);
  if ((a->ndims == 1 && a->dim[0] != x->dim[0])
      || (a->ndims == 2 && ((inner && a->dim[1] != x->dim[0])
          || (!inner && a->dim[0] != x->dim[0]))))
      luaL_error(L, CONF_ERROR);
  luaL_argcheck(L, what == 'l' || what == 'L'
      || what == 'u' || what == 'U' || what == 'f' || what == 'F', 4,
      "unknown triangle option");
  nl_hemul(x, a, inner, what, alpha);
  lua_settop(L, 1);
  return 1; /* self */
}


static int matrix_mmul (lua_State *L) {
  nl_Matrix *c = checkmatrix(L, 1);
  nl_Matrix *a = checkmatrix(L, 2);
  nl_Matrix *b = checkmatrix(L, 3);
  const char *ta = luaL_optstring(L, 4, "N");
  const char *tb = luaL_optstring(L, 5, "N");
  nl_Complex alpha = nl_optcomplex(L, 6, onec);
  char transa = ta[0];
  char transb = tb[0];
  /* check args */
  check2dmatrix(L, c, 1);
  check2dsection(L, c, 1);
  check2dmatrix(L, a, 2);
  check2dsection(L, a, 2);
  check2dmatrix(L, b, 3);
  check2dsection(L, b, 3);
  luaL_argcheck(L, transa == 'n' || transa == 'N'
      || transa == 't' || transa == 'T'
      || transa == 'c' || transa == 'C', 4, "unknown transpose option");
  luaL_argcheck(L, transb == 'n' || transb == 'N'
      || transb == 't' || transb == 'T'
      || transb == 'c' || transb == 'C', 5, "unknown transpose option");
  if (c->iscomplex != a->iscomplex || c->iscomplex != b->iscomplex)
    luaL_error(L, DOMC_ERROR);
  if (b->ndims == 1) {
    if (a->ndims == 1) {
      if (c->ndims == 1 || c->dim[0] != a->size || c->dim[1] != b->size)
        luaL_error(L, CONF_ERROR);
    }
    else {
      if (c->ndims == 2
          || ((transa == 'n' || transa == 'N')
            && (c->size != a->dim[0] || b->size != a->dim[1]))
          || ((transa != 'n' && transa != 'N')
            && (c->size != a->dim[1] || b->size != a->dim[0])))
        luaL_error(L, CONF_ERROR);
    }
  }
  else {
    if (c->ndims == 1 || a->ndims == 1)
      luaL_error(L, CONF_ERROR);
    if (transa == 'n' || transa == 'N')  {
      if ((transb == 'n' || transb == 'N')
          && (c->dim[0] != a->dim[0] || c->dim[1] != b->dim[1]
            || a->dim[1] != b->dim[0]))
        luaL_error(L, CONF_ERROR);
      if ((transb != 'n' && transb != 'N')
          && (c->dim[0] != a->dim[0] || c->dim[1] != b->dim[0]
            || a->dim[1] != b->dim[1]))
        luaL_error(L, CONF_ERROR);
    }
    else { /* transa != 'n'/'N' */
      if ((transb == 'n' || transb == 'N')
          && (c->dim[0] != a->dim[1] || c->dim[1] != b->dim[1]
            || a->dim[0] != b->dim[0]))
        luaL_error(L, CONF_ERROR);
      if ((transb != 'n' && transb != 'N')
          && (c->dim[0] != a->dim[1] || c->dim[1] != b->dim[0]
            || a->dim[0] != b->dim[1]))
        luaL_error(L, CONF_ERROR);
    }
  }
  nl_mmul(c, a, b, transa, transb, alpha);
  lua_settop(L, 1);
  return 1; /* self */
}


static int matrix_norm (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  char what = (lua_isnumber(L, 2)) ? 0 : *luaL_optstring(L, 2, "F");
  lua_Number p = 0;
  int argm;
  check2dmatrix(L, m, 1);
  check2dsection(L, m, 1);
  /* adjust args */
  if (what == 0) {
    p = lua_tonumber(L, 2);
    if (p == 1) what = 'O';
  }
  else if (what == 'F' || what == 'f' || what == 'E' || what == 'e') {
    p = 2;
    what = 0;
  }
  luaL_argcheck(L, what == 'O' || what == 'o' || what == 'I' || what == 'i'
      || what == 'M' || what == 'm' || what == 0, 2,
      "unknown norm option");
  lua_pushnumber(L, nl_norm(m, what, p, &argm));
  if (what == 'M' || what == 'm' || what == 'I' || what == 'i') {
    lua_pushinteger(L, argm);
    return 2;
  }
  return 1;
}


static int matrix_chol (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  const char *s = luaL_optstring(L, 2, "L");
  int inplace = nl_inplace(L, 3);
  char what = s[0];
  int info;
  checksquare(L, m, 1);
  if (what == 'l') what = 'L';
  if (what == 'u') what = 'U';
  luaL_argcheck(L, what == 'L' || what == 'U', 2, "unknown triangle option");
  lua_settop(L, 1);
  if (inplace)
    check2dsection(L, m, 1);
  else { /* copy */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  info = nl_chol(m, what);
  checkinfo(L, "chol");
  if (info > 0) {
    lua_pushboolean(L, 0);
    lua_pushfstring(L, "matrix is not positive definite:"
        " leading minor of order %d is not posdef", info);
    return 2;
  }
  settriangtoscalar(0, what == 'L' ? 'U' : 'L', m); /* m.(L or U) := 0 */
  return 1; /* self or copy */
}


static int matrix_lu (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int inplace = nl_inplace(L, 2);
  nl_Matrix *pvt;
  int i, mn, info;
  nl_Buffer *ipiv;
  check2dmatrix(L, m, 1);
  mn = (m->dim[0] < m->dim[1]) ? m->dim[0] : m->dim[1]; /* min(m,n) */
  ipiv = nl_getbuffer(L, mn);
  if (inplace) {
    check2dsection(L, m, 1);
    lua_settop(L, 1);
    info = nl_lu(m, ipiv);
  }
  else { /* full matrix return */
    nl_Matrix *l, *u;
    l = pushmatrix(L, m->iscomplex, 2, NULL, 1, m->dim[0] * mn, NULL, NULL);
    l->dim[0] = m->dim[0]; l->dim[1] = mn;
    u = pushmatrix(L, m->iscomplex, 2, NULL, 1, m->dim[1] * mn, NULL, NULL);
    u->dim[0] = mn; u->dim[1] = m->dim[1];
    if (m->dim[0] < m->dim[1]) { /* copy to u? */
      setdatatovector(m, u->stride, 0, u->data); /* u := m */
      info = nl_lu(u, ipiv);
      settriangtovector(u, 'L', l); /* l.L := u.L */
      settriangtoscalar(0, 'L', u); /* u.L := 0 */
      settriangtoscalar(0, 'U', l); /* l.U := 0 */
    }
    else { /* copy to l */
      int dsl = m->dim[0] + 1; /* diagonal stride for l */
      int dsu = mn + 1; /* diagonal stride for u */
      setdatatovector(m, l->stride, 0, l->data);
      info = nl_lu(l, ipiv);
      settriangtovector(l, 'U', u); /* u.U := l.U */
      settriangtoscalar(0, 'U', l); /* l.U := 0 */
      settriangtoscalar(0, 'L', u); /* u.L := 0 */
      /* u.D := l.D: */
      if (m->iscomplex)
        ZCOPY(&mn, CPX(l->data), &dsl, CPX(u->data), &dsu);
      else
        DCOPY(&mn, l->data, &dsl, u->data, &dsu);
    }
    /* set unit diagonal of l */
    setdatatoscalar(l->iscomplex, mn, 1, m->dim[0] + 1, 0, l->data);
  }
  pvt = pushmatrix(L, 0, 1, &mn, 1, mn, NULL, NULL);
  for (i = 0; i < mn; i++)
    pvt->data[i] = (lua_Number) ipiv->data.bint[i];
  nl_freebuffer(ipiv);
  checkinfo(L, "lu");
  if (inplace) return 2; /* m, pvt */
  return 3; /* l, u, pvt */
}


static int matrix_rcond (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  const char *s = luaL_optstring(L, 2, "g");
  int inplace = nl_inplace(L, 3);
  char what = s[0];
  int info;
  lua_Number r;
  checksquare(L, m, 1);
  luaL_argcheck(L, what == 'd' || what == 'D'
      || what == 'l' || what == 'L' || what == 'u' || what == 'U'
      || what == 'p' || what == 'P' || what == 'g' || what == 'G', 2,
      "unknown matrix option");
  if (inplace)
    check2dsection(L, m, 1);
  else { /* copy */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  if (what == 'g' || what == 'G') {
    nl_Buffer *ipiv = nl_getbuffer(L, m->dim[0]);
    r = nl_rcond(L, m, what, ipiv->data.bint, &info);
    nl_freebuffer(ipiv);
  }
  else
    r = nl_rcond(L, m, what, NULL, &info);
  checkinfo(L, "rcond");
  lua_pushnumber(L, r);
  return 1;
}


static int matrix_inv (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  const char *s = luaL_optstring(L, 2, "g");
  int inplace = nl_inplace(L, 3);
  int norcond = lua_toboolean(L, 4);
  char what = s[0];
  int info;
  lua_Number rcond;
  checksquare(L, m, 1);
  luaL_argcheck(L, what == 'd' || what == 'D'
      || what == 'l' || what == 'L' || what == 'u' || what == 'U'
      || what == 'p' || what == 'P' || what == 'g' || what == 'G', 2,
      "unknown matrix option");
  if (inplace) {
    check2dsection(L, m, 1);
    lua_settop(L, 1);
  }
  else {
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  info = nl_inv(L, m, what, norcond ? NULL : &rcond);
  checkinfo(L, "inv");
  if (info > 0) {
    lua_pushboolean(L, 0);
    lua_pushfstring(L, "matrix is singular: (%d, %d) element is zero",
        info, info);
    return 2;
  }
  if (norcond) return 1;
  lua_pushnumber(L, rcond);
  return 2;
}


static int matrix_svd (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  const char *s = luaL_optstring(L, 2, "a");
  char what = s[0];
  int info;
  check2dmatrix(L, m, 1);
  luaL_argcheck(L, what == 'n' || what == 'N'
      || what == 'l' || what == 'L' || what == 'r' || what == 'R'
      || what == 'a' || what == 'A', 2, "unknown job option");
  info = nl_svd(L, m, what);
  checkinfo(L, "svd");
  if (info > 0) {
    lua_pushboolean(L, 0);
    lua_pushfstring(L, "failed to converge: info = %d", info);
    return 2;
  }
  if (what != 'a' && what != 'A') return 1; /* only S */
  return 3;
}


static int matrix_qr (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int permute = lua_toboolean(L, 2);
  int inplace = nl_inplace(L, 3);
  int info;
  nl_Buffer *pvt = NULL;
  check2dmatrix(L, m, 1);
  if (permute) {
    int i;
    pvt = nl_getbuffer(L, m->dim[1]);
    for (i = 0; i < m->dim[1]; i++)
      pvt->data.bint[i] = 0; /* free columns */
  }
  if (inplace) {
    check2dsection(L, m, 1);
    lua_settop(L, 1);
  }
  else { /* copy m */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  info = nl_qr(L, m, pvt);
  checkinfo(L, "qr");
  lua_insert(L, -2); /* q before r */
  if (permute) { /* build pivot vector? */
    int i, j, n = m->dim[1];
    nl_Matrix *jpvt = pushmatrix(L, 0, 1, &n, 1, n, NULL, NULL);
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        if (pvt->data.bint[j] - 1 == i) { /* swap i and j? */
          jpvt->data[i] = j + 1;
          pvt->data.bint[j] = pvt->data.bint[i];
          break;
        }
      }
    }
    nl_freebuffer(pvt);
    return 3; /* q, r, pvt */
  }
  return 2; /* q, r */
}


static int matrix_eig (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  const char *s = luaL_optstring(L, 2, "R");
  int hermitian = lua_toboolean(L, 3);
  char what = s[0];
  int info;
  checksquare(L, m, 1);
  luaL_argcheck(L, what == 'a' || what == 'A'
      || what == 'l' || what == 'L' || what == 'r' || what == 'R'
      || what == 'n' || what == 'N', 2, "unknown job option");
  info = nl_eig(L, m, what, hermitian);
  checkinfo(L, "eig");
  if (info > 0) {
    lua_pushboolean(L, 0);
    lua_pushfstring(L, "failed to converge: info = %d", info);
    return 2;
  }
  if (what == 'n' || what == 'N') return 1; /* just eigenvalues */
  if (what == 'a' || what == 'A') return 3; /* eigenvalues & eigenvectors */
  return 2; /* eigenvalues and one-side eigenvectors */
}


static int matrix_balance (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  const char *s = luaL_optstring(L, 2, "B");
  char what = s[0];
  int info;
  lua_settop(L, 1);
  checksquare(L, a, 1);
  if (what == 'n') what = 'N';
  if (what == 's') what = 'S';
  if (what == 'p') what = 'P';
  if (what == 'b') what = 'B';
  luaL_argcheck(L, what == 'N' || what == 'S' || what == 'P' || what == 'B',
      2, "unknown balance option");
  info = nl_balance(L, a, what);
  checkinfo(L, "balance");
  return 2;
}


static int matrix_ls (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  int svd = lua_toboolean(L, 3);
  lua_Number tol = luaL_optnumber(L, 4, 0);
  int inplace = nl_inplace(L, 5);
  int info, rank;
  nl_Matrix *s = NULL;
  checkarray(L, a, 1);
  check2dsection(L, a, 1);
  check2dmatrix(L, b, 2);
  if (inplace) {
    checksection(L, b, 2);
    luaL_argcheck(L, a->dim[0] > a->dim[1], 1,
        "not enough space to store solution in-place");
  }
  if (a->iscomplex != b->iscomplex)
    luaL_error(L, DOMC_ERROR);
  lua_settop(L, 2);
  if (tol <= 0) /* set default tolerance? */
    tol = (a->dim[0] > a->dim[1] ? a->dim[0] : a->dim[1]) * DBL_EPSILON;
  if (svd) {
    int mn = (a->dim[0] > a->dim[1]) ? a->dim[1] : a->dim[0];
    s = pushmatrix(L, 0, 1, &mn, 1, mn, NULL, NULL);
  }
  info = nl_ls(L, a, b, s, tol, inplace, &rank);
  checkinfo(L, "ls");
  if (info > 0) {
    lua_pushboolean(L, 0);
    lua_pushfstring(L, "failed to converge: info = %d", info);
    return 2;
  }
  lua_pushinteger(L, rank);
  if (svd) {
    if (inplace)
      lua_insert(L, -2); /* rank before s */
    else
      lua_pushvalue(L, -3); /* s at the top */
  }
  if (svd) return 3; /* sol, rank, s */
  return 2; /* sol, rank */
}


/* {=======   Fourier Transforms   =======} */

fftw_plan nl_createplan (lua_State *L,
    nl_Matrix *m, int inverse, unsigned flags, lua_Number *scale);

static int matrix_fft (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int inverse = lua_toboolean(L, 2);
  int wisdomonly = lua_toboolean(L, 3);
  int inplace = nl_inplace(L, 4);
  lua_Number scale;
  fftw_plan plan;
  luaL_argcheck(L, m->iscomplex, 1, "complex matrix expected");
  if (inplace) {
    checksection(L, m, 1);
    lua_settop(L, 1);
  }
  else { /* copy */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  plan = nl_createplan(L, m, inverse, 
      wisdomonly ? FFTW_WISDOM_ONLY : FFTW_ESTIMATE, &scale);
  if (plan == NULL)
    luaL_error(L, "cannot create plan");
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  if (inverse) /* scale? */
    ZDSCAL(&m->size, &scale, CPX(m->data), &m->stride);
  return 1;
}

static int matrix_fct (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int inverse = lua_toboolean(L, 2);
  int wisdomonly = lua_toboolean(L, 3);
  int inplace = nl_inplace(L, 4);
  lua_Number scale;
  fftw_plan plan;
  luaL_argcheck(L, !m->iscomplex, 1, "real matrix expected");
  if (inplace) {
    checksection(L, m, 1);
    lua_settop(L, 1);
  }
  else { /* copy */
    m = pushframe(L, m);
    settoarg(L, m, 0, 1, m->size, 0, 1);
  }
  plan = nl_createplan(L, m, inverse, 
      wisdomonly ? FFTW_WISDOM_ONLY : FFTW_ESTIMATE, &scale);
  if (plan == NULL)
    luaL_error(L, "cannot create plan");
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  if (inverse) /* scale? */
    DSCAL(&m->size, &scale, m->data, &m->stride);
  return 1;
}


/* {=======   Serialization (HDF5)   =======} */

typedef struct { double re, im; } TComplex;

static int matrix_save (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  const char *path = luaL_checkstring(L, 2);
  hid_t file_id, type_id, dataspace_id, dataset_id;
  herr_t status;
  hsize_t *dim;
  int i;
  if (m->section || m->stride > 1) { /* copy to buffer? */
    nl_Matrix *x = pushframe(L, m);
    settoarg(L, x, 0, 1, x->size, 0, 1);
    m = x;
  }
  file_id = H5Fcreate(path, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  type_id = H5T_IEEE_F64LE;
  if (m->iscomplex) {
    type_id = H5Tcreate(H5T_COMPOUND, sizeof(TComplex));
    H5Tinsert(type_id, "real", HOFFSET(TComplex, re), H5T_IEEE_F64LE);
    H5Tinsert(type_id, "imag", HOFFSET(TComplex, im), H5T_IEEE_F64LE);
  }
  dim = (hsize_t *) lua_newuserdata(L, m->ndims * sizeof(hsize_t));
  for (i = 0; i < m->ndims; i++)
    dim[i] = m->dim[i];
  dataspace_id = H5Screate_simple(m->ndims, dim, NULL);
  dataset_id = H5Dcreate(file_id, "/matrix", type_id, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      m->data);
  if (m->iscomplex)
    status = H5Tclose(type_id);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  return 0;
}

static int matrix_load (lua_State *L) {
  nl_Matrix *m;
  const char *path = luaL_checkstring(L, 1);
  hid_t file_id, type_id, dataspace_id, dataset_id;
  herr_t status;
  int i, size, ndims, iscomplex;
  hsize_t *dim;
  /* read matrix properties */
  file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, "/matrix", H5P_DEFAULT);
  type_id = H5Dget_type(dataset_id);
  iscomplex = H5Tget_class(type_id) == H5T_COMPOUND; /* assume complex */
  dataspace_id = H5Dget_space(dataset_id);
  ndims = H5Sget_simple_extent_ndims(dataspace_id);
  dim = (hsize_t *) lua_newuserdata(L, ndims * sizeof(hsize_t));
  H5Sget_simple_extent_dims(dataspace_id, dim, dim);
  size = 1;
  for (i = 0; i < ndims; i++)
    size *= dim[i];
  m = pushmatrix(L, iscomplex, ndims, NULL, 1, size, NULL, NULL);
  for (i = 0; i < ndims; i++)
    m->dim[i] = dim[i];
  /* write data from disk */
  status = H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      m->data);
  status = H5Sclose(dataspace_id);
  if (iscomplex)
    status = H5Tclose(type_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  return 1;
}


/* {=====================================================================
 *    Interface
 * ======================================================================} */

static const luaL_Reg lmatrix_func[] = {
  /* metamethods */
  {"set", matrix_set},
  {"slice", matrix_slice},
  {"section", matrix_section},
  {"add", matrix_add},
  {"mul", matrix_mul},
  {"div", matrix_div},
  {"pow", matrix_pow},
  /* methods */
  {"new", matrix_new},
  {"eorder", matrix_eorder},
  {"eindex", matrix_eindex},
  {"entries", matrix_entries},
  {"size", matrix_size},
  {"shape", matrix_shape},
  {"iscomplex", matrix_iscomplex},
  {"real", matrix_real},
  {"imag", matrix_imag},
  {"complex", matrix_complex},
  {"conj", matrix_conj},
  {"copy", matrix_copy},
  {"reshape", matrix_reshape},
  {"spread", matrix_spread},
  /* logical */
  {"find", matrix_find},
  {"ifelse", matrix_ifelse},
  {"which", matrix_which},
  /* functional */
  {"apply", matrix_apply},
  {"fold", matrix_fold},
  {"map", matrix_map},
  /* statistical */
  {"sum", matrix_sum},
  {"min", matrix_min},
  {"max", matrix_max},
  {"sort", matrix_sort},
  /* functions */
  {"abs", matrix_abs},
  {"acos", matrix_acos},
  {"acosh", matrix_acosh},
  {"asin", matrix_asin},
  {"asinh", matrix_asinh},
  {"atan", matrix_atan},
  {"atanh", matrix_atanh},
  {"cos", matrix_cos},
  {"cosh", matrix_cosh},
  {"exp", matrix_exp},
  {"log", matrix_log},
  {"sin", matrix_sin},
  {"sinh", matrix_sinh},
  {"sqrt", matrix_sqrt},
  {"tan", matrix_tan},
  {"tanh", matrix_tanh},
  {"linspace", matrix_linspace},
  {"dot", matrix_dot},
  {"cross", matrix_cross},
  /* 2d */
  {"col", matrix_col},
  {"transpose", matrix_transpose},
  {"diag", matrix_diag},
  {"concat", matrix_concat},
  {"c", matrix_c},
  /* linear algebra */
  {"swap", matrix_swap},
  {"pivot", matrix_pivot},
  {"trmul", matrix_trmul},
  {"hemul", matrix_hemul},
  {"mmul", matrix_mmul},
  {"norm", matrix_norm},
  {"chol", matrix_chol},
  {"lu", matrix_lu},
  {"rcond", matrix_rcond},
  {"inv", matrix_inv},
  {"svd", matrix_svd},
  {"qr", matrix_qr},
  {"eig", matrix_eig},
  {"balance", matrix_balance},
  {"ls", matrix_ls},
  /* FFT */
  {"fft", matrix_fft},
  {"fct", matrix_fct},
  /* HDF5 */
  {"save", matrix_save},
  {"load", matrix_load},
  {NULL, NULL}
};

static const luaL_Reg lmatrix_mt[] = {
  {"__newindex", matrix_set},
  {"__len", matrix_size},
  {"__concat", matrix_concat},
  {"__tostring", matrix__tostring},
  {"__unm", matrix__unm},
  {"__pow", matrix_pow},
  {NULL, NULL}
};

int luaopen_numlua_lmatrix (lua_State *L) {
  /* load lmatrix */
  lua_newtable(L); /* new metatable */
  lua_pushliteral(L, MATRIX_LIBNAME);
  lua_setfield(L, -2, "__type");
  lua_createtable(L, 0, 1); /* meta metatable */
  lua_pushliteral(L, "k");
  lua_setfield(L, -2, "__mode");
  lua_setmetatable(L, -2); /* set metatable as weak-keyed table */
  lua_pushlightuserdata(L, MATRIX_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX);
  /* push metamethods */
  lua_pushvalue(L, -1);
  nl_register(L, lmatrix_mt, 1);
  luaL_newlibtable(L, lmatrix_func);
  lua_pushvalue(L, -2); /* metatable */
  nl_register(L, lmatrix_func, 1);
  /* push __index */
  lua_pushvalue(L, -2); /* metatable */
  lua_pushvalue(L, -2); /* class */
  lua_pushcclosure(L, matrix_get, 2); /* mt and class as upvalues */
  lua_setfield(L, -3, "__index");
  /* push `get` */
  lua_pushvalue(L, -2); /* metatable */
  lua_pushvalue(L, -2); /* class */
  lua_pushcclosure(L, matrix_get, 2); /* mt and class as upvalues */
  lua_setfield(L, -2, "get");
  return 1;
}

