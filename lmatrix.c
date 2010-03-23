#include <lua.h>
#include <lauxlib.h>
#include "numlua.h"
#include "cblas.h"
#include "clapack.h"

static int matrix_mt_ = 0;
#define MATRIX_MT ((void *)&matrix_mt_)

#define DCOPY cblas_dcopy
#define ZCOPY cblas_zcopy
#define DSCAL cblas_dscal
#define ZSCAL cblas_zscal
#define DAXPY cblas_daxpy
#define ZAXPY cblas_zaxpy
#define DDOT cblas_ddot
#define ZDOTU cblas_zdotu
#define ZDOTC cblas_zdotc

static nl_Complex C0 = {0, 0};
static nl_Complex C1 = {1, 0};

#define CPX(x) ((nl_Complex *) (x))

/* Internal: meant for use on functions that have MT as env */
toudata(matrix, MATRIX_LIBNAME, nl_Matrix)

/* assumes data block is at top of stack */
static nl_Matrix *pushmatrix (lua_State *L, int iscomplex, int ndims, int *dim,
    int stride, int size, lua_Number *data) {
  int i;
  nl_Matrix *m = lua_newuserdata(L, sizeof(nl_Matrix)
      + ndims * sizeof(int));
  lua_pushvalue(L, LUA_ENVIRONINDEX); /* env */
  lua_pushvalue(L, -2);
  lua_pushvalue(L, -4); /* data block */
  lua_rawset(L, -3); /* env[matrix] = data */
  m->ndims = ndims;
  m->stride = stride;
  m->size = size;
  m->iscomplex = iscomplex;
  m->data = data;
  for (i = 0; i < ndims; i++) m->dim[i] = dim[i];
  lua_setmetatable(L, -2);
  if (ndims > 1) { /* not a vector? */
    lua_newtable(L); /* new matrix environment */
    lua_setfenv(L, -2);
  }
  lua_replace(L, -2); /* pop data block */
  return m;
}

static int checkconsistency (nl_Matrix *a, nl_Matrix *b) {
  int i;
  if (a->ndims != b->ndims || a->size != b->size) return 0;
  for (i = 0; i < a->ndims; i++)
    if (a->dim[i] != b->dim[i]) return 0;
  return 1;
}


/* API */
libtoudata(matrix, MATRIX_LIBNAME, MATRIX_MT, nl_Matrix)

/* assumes data block is at top of stack */
NUMLUA_API nl_Matrix *nl_pushmatrix (lua_State *L, int iscomplex, int ndims,
    int *dim, int stride, int size, lua_Number *data) {
  int i;
  nl_Matrix *m = lua_newuserdata(L, sizeof(nl_Matrix)
      + ndims * sizeof(int));
  lua_pushlightuserdata(L, MATRIX_MT);
  lua_rawget(L, LUA_REGISTRYINDEX); /* env */
  lua_pushvalue(L, -2);
  lua_pushvalue(L, -4); /* data block */
  lua_rawset(L, -3); /* env[matrix] = data */
  m->ndims = ndims;
  m->stride = stride;
  m->size = size;
  m->iscomplex = iscomplex;
  m->data = data;
  for (i = 0; i < ndims; i++) m->dim[i] = dim[i];
  lua_setmetatable(L, -2);
  if (ndims > 1) { /* not a vector? */
    lua_newtable(L); /* new matrix environment */
    lua_setfenv(L, -2);
  }
  lua_replace(L, -2); /* pop data block */
  return m;
}


/* Metamethods */

static int matrix__index (lua_State *L) {
  nl_Matrix *m = (nl_Matrix *) lua_touserdata(L, 1);
  if (lua_isnumber(L, 2)) {
    int p, k = lua_tointeger(L, 2);
    if (k < 1 || k > m->dim[0]) luaL_error(L, "matrix index out of range");
    p = (k - 1) * m->stride; /* data offset */
    if (m->ndims == 1) { /* vector? */
      if (m->iscomplex) nl_pushcomplex(L, CPX(m->data)[p]);
      else lua_pushnumber(L, m->data[p]);
    }
    else {
      lua_getfenv(L, 1); /* matrix env */
      lua_rawgeti(L, -1, k);
      if (lua_isnil(L, -1)) { /* isn't row k interned? */
        lua_pushvalue(L, 1); /* matrix */
        lua_rawget(L, LUA_ENVIRONINDEX); /* push data */
        /* push new row */
        if (m->iscomplex)
          pushmatrix(L, m->iscomplex, m->ndims - 1, m->dim + 1,
              m->dim[0] * m->stride, m->size / m->dim[0],
              (lua_Number *) (CPX(m->data) + p));
        else
          pushmatrix(L, m->iscomplex, m->ndims - 1, m->dim + 1,
              m->dim[0] * m->stride, m->size / m->dim[0],
              m->data + p);
        lua_pushvalue(L, -1);
        lua_rawseti(L, 3, k); /* matrix env[k] = new row */
      }
    }
  }
  else { /* meta lookup? */
    lua_pushvalue(L, lua_upvalueindex(1)); /* class */
    lua_pushvalue(L, 2);
    lua_rawget(L, -2);
  }
  return 1;
}

static int matrix__newindex (lua_State *L) {
  nl_Matrix *m = (nl_Matrix *) lua_touserdata(L, 1);
  int p, k, iscomplex;
  nl_Complex v;
  if (m->ndims > 1) luaL_error(L, "cannot assign to matrix row");
  v = nl_tocomplex(L, 3, &iscomplex);
  if (!lua_isnumber(L, 2) || !iscomplex)
    luaL_error(L, "wrong type to matrix assignment");
  k = lua_tointeger(L, 2);
  p = (k - 1) * m->stride; /* data offset */
  /*if (m->iscomplex) p *= 2;*/
  if (k < 1 || k > m->dim[0]) luaL_error(L, "matrix index out of range");
  if (m->iscomplex) CPX(m->data)[p] = v;
  else m->data[p] = v.re;
  return 0;
}

static int matrix__len (lua_State *L) {
  nl_Matrix *m = lua_touserdata(L, 1);
  lua_pushinteger(L, m->ndims);
  return 1;
}

static int matrix__tostring (lua_State *L) {
  lua_pushfstring(L, MATRIX_LIBNAME ": %p", lua_touserdata(L, 1));
  return 1;
}


/* Methods */

/* if last arg is true, matrix is complex */
static int matrix_zeros (lua_State *L) {
  nBuffer *buf;
  int i, size = 1;
  int iscomplex = 0;
  int ndims = lua_gettop(L);
  lua_Number *data;
  if (ndims == 0) luaL_error(L, "no dimensions given");
  if (lua_type(L, ndims) == LUA_TBOOLEAN)
    iscomplex = lua_toboolean(L, ndims--);
  buf = nl_getbuffer(L, ndims);
  for (i = 0; i < ndims; i++) {
    int e;
    e = lua_tointeger(L, i + 1);
    luaL_argcheck(L, e > 0, i + 1, "invalid dimension");
    buf->data.bint[i] = e;
    size *= e;
  }
  data = lua_newuserdata(L, size
      * (iscomplex ? sizeof(nl_Complex) : sizeof(lua_Number)));
  if (iscomplex)
    for (i = 0; i < size; i++) CPX(data)[i] = C0;
  else
    for (i = 0; i < size; i++) data[i] = 0;
  pushmatrix(L, iscomplex, ndims, buf->data.bint, 1, size, data);
  nl_freebuffer(buf);
  return 1;
}

static int matrix_size (lua_State *L) {
  int i;
  nl_Matrix *m = checkmatrix(L, 1);
  lua_settop(L, 2);
  if (lua_isnil(L, 2)) { /* return all dims? */
    for (i = 0; i < m->ndims; i++) lua_pushinteger(L, m->dim[i]);
    return m->ndims;
  }
  i = lua_tointeger(L, 2);
  if (i < 1 || i > m->ndims) lua_pushnil(L);
  else lua_pushinteger(L, m->dim[i - 1]);
  return 1;
}

static int matrix_iscomplex (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_pushboolean(L, m->iscomplex);
  return 1;
}

/* matrix.real(a) returns a matrix m such that m[.] = real(a[i]) if a is
 * complex or a if a is real. */
static int matrix_real (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_Number *data;
  lua_settop(L, 1);
  if (!m->iscomplex) return 1;
  data = lua_newuserdata(L, m->size * sizeof(lua_Number));
  pushmatrix(L, 0, m->ndims, m->dim, 1, m->size, data);
  DCOPY(m->size, m->data, 2 * m->stride, data, 1);
  return 1;
}

/* matrix.imag(a) returns a matrix m such that m[.] = imag(a[i]) if a is
 * complex or fills a with zeros if a is real. */
static int matrix_imag (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_settop(L, 1);
  if (m->iscomplex) {
    lua_Number *data = lua_newuserdata(L, m->size * sizeof(lua_Number));
    pushmatrix(L, 0, m->ndims, m->dim, 1, m->size, data);
    DCOPY(m->size, m->data + 1, 2 * m->stride, data, 1);
  }
  else {
    int i;
    for (i = 0; i < m->size; i++) m->data[i] = 0;
  }
  return 1;
}

/* matrix.complex(a [, b]) returns @a if a is complex or a complex matrix m
 * such that real(m[.]) = a[.] and imag(m[.]) = b[.] if a is real. */
static int matrix_complex (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_Number *data;
  int i;
  lua_settop(L, 1);
  if (m->iscomplex) return 1;
  data = lua_newuserdata(L, m->size * sizeof(nl_Complex));
  pushmatrix(L, 1, m->ndims, m->dim, 1, m->size, data);
  DCOPY(m->size, m->data, m->stride, data, 2); /* real part */
  data++;
  for (i = 0; i < 2 * m->size; i += 2) data[i] = 0;
  return 1;
}

/* matrix.conj(a) sets each entry in a to its conjugate, a[.] = conj(a[.]). */
static int matrix_conj (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  lua_settop(L, 1);
  if (m->iscomplex) DSCAL(m->size, -1, m->data + 1, 2);
  return 1;
}

/* matrix.copy(m [, x]): returns a copy of m or, if x is a consistent matrix
 * with m, m is copied to x. */
static int matrix_copy (lua_State *L) {
  nl_Matrix *x, *m = checkmatrix(L, 1);
  lua_Number *data;
  int inc;
  lua_settop(L, 2);
  if (!lua_isnil(L, 2)) {
    x = checkmatrix(L, 2);
    if (!checkconsistency(m, x))
      luaL_error(L, "matrices are not consistent");
    data = x->data;
    inc = x->stride;
  }
  else {
    data = lua_newuserdata(L, m->size
        * ((m->iscomplex) ? sizeof(nl_Complex) : sizeof(lua_Number)));
    inc = 1;
    pushmatrix(L, m->iscomplex, m->ndims, m->dim, inc, m->size, data);
  }
  /* copy */
  if (m->iscomplex)
    ZCOPY(m->size, CPX(m->data), m->stride, CPX(data), inc);
  else
    DCOPY(m->size, m->data, m->stride, data, inc);
  return 1;
}

/* matrix.trcopy(m [, lower, [, x]]): returns a copy of the entries in the
 * upper triangle of two-dimensional matrix m, that is, above and including
 * the diagonal, if @lower is false or a copy of the entries in the lower
 * triangle of m, below and including the diagonal, if @lower is true.
 * Similarly to matrix.copy, if x is a consistent matrix with m then x serves
 * as copy destination. */
static int matrix_trcopy (lua_State *L) {
  nl_Matrix *x, *m = checkmatrix(L, 1);
  int lower = lua_toboolean(L, 2);
  lua_Number *data;
  int i, n, inc;
  luaL_argcheck(L, m->ndims == 2, 1, "two-dimensional matrix expected");
  lua_settop(L, 3);
  if (!lua_isnil(L, 3)) {
    x = checkmatrix(L, 3);
    luaL_argcheck(L, checkconsistency(m, x), 3,
        "matrices are not consistent");
    data = x->data;
    inc = x->stride;
  }
  else {
    data = lua_newuserdata(L, m->size
        * ((m->iscomplex) ? sizeof(nl_Complex) : sizeof(lua_Number)));
    inc = 1;
    pushmatrix(L, m->iscomplex, m->ndims, m->dim, inc, m->size, data);
  }
  /* triangular copy */
  n = m->dim[0];
  if (m->iscomplex) {
    if (lower)
      for (i = 0; i < m->dim[1]; i++) {
        ZCOPY(n, CPX(m->data) + i * (m->dim[0] + 1) * m->stride, m->stride,
            CPX(data) + i * (m->dim[0] + 1) * inc, inc);
        n--;
      }
    else /* upper */
      for (i = m->dim[1] - 1; i >= 0; i--) {
        ZCOPY(n, CPX(m->data) + i * m->dim[0] * m->stride, m->stride,
            CPX(data) + i * m->dim[0] * inc, inc);
        if (i < n) n--;
      }
  }
  else { /* real */
    if (lower)
      for (i = 0; i < m->dim[1]; i++) {
        DCOPY(n, m->data + i * (m->dim[0] + 1) * m->stride, m->stride,
            data + i * (m->dim[0] + 1) * inc, inc);
        n--;
      }
    else /* upper */
      for (i = m->dim[1] - 1; i >= 0; i--) {
        DCOPY(n, m->data + i * m->dim[0] * m->stride, m->stride,
            data + i * m->dim[0] * inc, inc);
        if (i < n) n--;
      }
  }
  return 1;
}

static int matrix_linspace (lua_State *L) {
  nl_Complex a = nl_checkcomplex(L, 1);
  nl_Complex b = nl_checkcomplex(L, 2);
  int i, n, iscomplex;
  lua_Number *data;
  iscomplex = (a.im != 0 || b.im != 0);
  if (iscomplex) {
    nl_Complex s;
    s.re = b.re - a.re;
    s.im = b.im - a.im;
    n = luaL_optinteger(L, 3, nl_abscomplex(s) + 1);
    luaL_argcheck(L, n > 0, 3, "number of steps is non-positive");
    lua_settop(L, 0);
    data = lua_newuserdata(L, n * sizeof(nl_Complex));
    s.re /= n - 1;
    s.im /= n - 1;
    CPX(data)[0] = a;
    for (i = 1; i < n; i++) {
      CPX(data)[i].re = CPX(data)[i - 1].re + s.re;
      CPX(data)[i].im = CPX(data)[i - 1].im + s.im;
    }
    pushmatrix(L, 1, 1, &n, 1, n, data);
  }
  else {
    lua_Number s;
    n = luaL_optinteger(L, 3, fabs(b.re - a.re) + 1);
    luaL_argcheck(L, n > 0, 3, "number of steps is non-positive");
    lua_settop(L, 0);
    data = lua_newuserdata(L, n * sizeof(lua_Number));
    s = (b.re - a.re) / (n - 1);
    data[0] = a.re;
    for (i = 1; i < n; i++) data[i] = data[i - 1] + s;
    pushmatrix(L, 0, 1, &n, 1, n, data);
  }
  return 1;
}

static int matrix_scale (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Complex s = nl_checkcomplex(L, 2);
  lua_settop(L, 1);
  if (m->iscomplex)
    ZSCAL(m->size, &s, CPX(m->data), m->stride);
  else
    DSCAL(m->size, s.re, m->data, m->stride);
  return 1;
}

static int matrix_fill (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Complex s = nl_checkcomplex(L, 2);
  int i, n = m->size * m->stride;
  lua_settop(L, 1);
  if (m->iscomplex)
    for (i = 0; i < n; i += m->stride) CPX(m->data)[i] = s;
  else
    for (i = 0; i < n; i += m->stride) m->data[i] = s.re;
  return 1;
}

static int matrix_shift (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Complex s = nl_checkcomplex(L, 2);
  int i, n = m->size * m->stride;
  lua_settop(L, 1);
  if (m->iscomplex)
    for (i = 0; i < n; i += m->stride) {
      CPX(m->data)[i].re += s.re;
      CPX(m->data)[i].im += s.im;
    }
  else
    for (i = 0; i < n; i += m->stride) m->data[i] += s.re;
  return 1;
}

static int matrix_abs (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, n = m->size * m->stride;
  lua_settop(L, 1);
  if (m->iscomplex) {
    for (i = 0; i < n; i += m->stride) {
      CPX(m->data)[i].re = nl_abscomplex(CPX(m->data)[i]);
      CPX(m->data)[i].im = 0;
    }
  }
  else {
    for (i = 0; i < n; i += m->stride)
      m->data[i] = fabs(m->data[i]);
  }
  return 1;
}

static int matrix_pow (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, n = m->size * m->stride;
  nl_Complex p = nl_checkcomplex(L, 2);
  if (!(m->iscomplex) && p.im != 0.0)
    luaL_error(L, "cannot take complex power of a real matrix");
  if (m->iscomplex) {
    if (p.im == 0.0 && p.re == floor(p.re)) { /* integer pow? */
      lua_Number x = p.re;
      int e;
      lua_number2int(e, x);
      for (i = 0; i < n; i += m->stride)
        nl_ipowcomplex(CPX(m->data) + i, CPX(m->data)[i], e);
    }
    else {
      for (i = 0; i < n; i += m->stride)
        nl_powcomplex(CPX(m->data) + i, CPX(m->data)[i], p);
    }
  }
  else {
    for (i = 0; i < n; i += m->stride)
      m->data[i] = pow(m->data[i], p.re);
  }
  lua_settop(L, 1);
  return 1;
}

#define matrixfunc(name) \
  static int matrix_ ## name (lua_State *L) { \
    nl_Matrix *m = checkmatrix(L, 1); \
    int i, n = m->size * m->stride; \
    lua_settop(L, 1); \
    if (m->iscomplex) \
      for (i = 0; i < n; i += m->stride) \
        nl_ ## name ## complex(CPX(m->data) + i, CPX(m->data)[i]); \
    else \
      for (i = 0; i < n; i += m->stride) \
        m->data[i] = name(m->data[i]); \
    return 1; \
  }

matrixfunc(exp)
matrixfunc(log)
matrixfunc(sin)
matrixfunc(cos)
matrixfunc(sqrt)


/* functional */

static int matrix_linfold (lua_State *L) { /* linear fold */
  nl_Matrix *m = checkmatrix(L, 1);
  nl_Complex alpha = nl_optcomplex(L, 2, C1);
  nl_Complex x = nl_optcomplex(L, 3, C0);
  int i, n = m->size * m->stride;
  if (m->iscomplex) {
    nl_Complex *d = (nl_Complex *) m->data;
    for (i = 0; i < n; i += m->stride) {
      lua_Number t = x.re * alpha.re - x.im * alpha.im + d[i].re;
      x.im = x.re * alpha.im + x.im * alpha.re + d[i].im;
      x.re = t;
    }
    nl_pushcomplex(L, x);
  }
  else {
    for (i = 0; i < n; i += m->stride)
      x.re = x.re * alpha.re + m->data[i];
    lua_pushnumber(L, x.re);
  }
  return 1;
}

static int matrix_fold (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, n = m->size * m->stride;
  luaL_checktype(L, 2, LUA_TFUNCTION); /* fold(acc, entry) */
  lua_settop(L, 3);
  for (i = 0; i < n; i += m->stride) {
    lua_pushvalue(L, 2); /* function */
    lua_insert(L, -2);
    /* push entry */
    if (m->iscomplex) nl_pushcomplex(L, CPX(m->data)[i]);
    else lua_pushnumber(L, m->data[i]);
    lua_call(L, 2, 1);
  }
  return 1;
} 

static int matrix_map (lua_State *L) {
  nl_Matrix *m = checkmatrix(L, 1);
  int i, iscomplex, n = m->size * m->stride;
  nl_Complex x;
  luaL_checktype(L, 2, LUA_TFUNCTION);
  lua_settop(L, 2);
  for (i = 0; i < n; i += m->stride) {
    lua_pushvalue(L, 2); /* function */
    /* push entry */
    if (m->iscomplex) nl_pushcomplex(L, CPX(m->data)[i]);
    else lua_pushnumber(L, m->data[i]);
    lua_call(L, 1, 1);
    x = nl_tocomplex(L, 3, &iscomplex);
    if (iscomplex) {
      if (m->iscomplex) CPX(m->data)[i] = x;
      else m->data[i] = x.re;
    }
    lua_pop(L, 1);
  }
  lua_pop(L, 1); /* function */
  return 1; /* matrix */
}


static int matrix_ewmul (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  int i;
  if (a->size != b->size || a->iscomplex != b->iscomplex)
    luaL_error(L, "matrices are not consistent");
  lua_settop(L, 2);
  if (a->iscomplex) {
    for (i = 0; i < a->size; i++) 
      nl_mulcomplex(CPX(a->data) + i * a->stride,
              CPX(a->data)[i * a->stride], CPX(b->data)[i * b->stride]);
  }
  else {
    for (i = 0; i < a->size; i++)
      a->data[i * a->stride] *= b->data[i * b->stride];
  }
  lua_pop(L, 1); /* b */
  return 1; /* a */
}

/* FIXME: complex, check strides for multi-D */
static int matrix_add (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  lua_Number alpha = luaL_optnumber(L, 3, 1.0);
  if (a->size != b->size)
    luaL_error(L, "matrices are not consistent");
  lua_settop(L, 2);
  DAXPY(b->size, alpha, b->data, b->stride, a->data, a->stride);
  lua_pop(L, 1); /* b */
  return 1; /* a */
}

/* FIXME: complex, check strides for multi-D */
static int matrix_dot (lua_State *L) {
  nl_Matrix *a = checkmatrix(L, 1);
  nl_Matrix *b = checkmatrix(L, 2);
  if (a->size != b->size)
    luaL_error(L, "matrices are not consistent");
  lua_pushnumber(L, DDOT(a->size, a->data, a->stride, b->data, b->stride));
  return 1;
}



/* TODO: ewdiv, ewldiv, ewmul, ewpow, min, max, range (__call),
 * add, dot, diag, ger, syr, mul, trmul, transpose, ctranspose,
 * rowcat, colcat, inner, outer, find,
 * eigs, eig, ls, lsq, chol,
 * fct, ifct, fst, ifst, fft, ifft, (fftw)
 * load, save (hdf5) */


/* library setup */

static const luaL_reg lmatrix_func[] = {
  /* info */
  {"size", matrix_size},
  {"iscomplex", matrix_iscomplex},
  /* level 1 */
  {"zeros", matrix_zeros},
  {"real", matrix_real},
  {"imag", matrix_imag},
  {"complex", matrix_complex},
  {"conj", matrix_conj},
  {"copy", matrix_copy},
  {"trcopy", matrix_trcopy},
  {"linspace", matrix_linspace},
  {"scale", matrix_scale},
  {"fill", matrix_fill},
  {"shift", matrix_shift},
  {"abs", matrix_abs},
  {"pow", matrix_pow},
  {"exp", matrix_exp},
  {"log", matrix_log},
  {"sin", matrix_sin},
  {"cos", matrix_cos},
  {"sqrt", matrix_sqrt},
  /* functional */
  {"linfold", matrix_linfold},
  {"fold", matrix_fold},
  {"map", matrix_map},
  /* FIXME */
  {"ewmul", matrix_ewmul},
  {"add", matrix_add},
  {"dot", matrix_dot},
  {NULL, NULL}
};

static const luaL_reg lmatrix_mt[] = {
  {"__newindex", matrix__newindex},
  {"__len", matrix__len},
  {"__tostring", matrix__tostring},
  {NULL, NULL}
};

static const luaL_Reg lmatrix_deps[] = {
  {NUMLUA_LIBNAME,  luaopen_numlua_base},
  {COMPLEX_LIBNAME, luaopen_numlua_complex},
  {FFT_LIBNAME, luaopen_numlua_fft},
  {NULL, NULL}
};

int luaopen_numlua_lmatrix (lua_State *L) {
  if (nl_requiredeps(L, MATRIX_LIBNAME, lmatrix_deps)) return 1;
  /* load lmatrix */
  lua_newtable(L); /* new environment/metatable */
  lua_createtable(L, 0, 1); /* env metatable */
  lua_pushliteral(L, "k");
  lua_setfield(L, -2, "__mode");
  lua_setmetatable(L, -2); /* set environment as weak-keyed table */
  lua_pushlightuserdata(L, MATRIX_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX);
  /* set as current environment */
  lua_pushvalue(L, -1);
  lua_replace(L, LUA_ENVIRONINDEX);
  /* push metamethods */
  luaL_register(L, NULL, lmatrix_mt);
  lua_pushliteral(L, MATRIX_LIBNAME);
  lua_setfield(L, -2, "__type");
  luaL_register(L, MATRIX_LIBNAME, lmatrix_func); /* fill class table */
  /* push class table as upvalue to __index */
  lua_pushvalue(L, -1);
  lua_pushcclosure(L, matrix__index, 1); /* class as upvalue for __index */
  lua_setfield(L, -3, "__index");
  return 1;
}

