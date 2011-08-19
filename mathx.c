/* {=================================================================
 *
 * mathx.c
 * Extra math routines for NumericLua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lua.h>
#include <lauxlib.h>
#include <float.h>
#include <math.h>
#include "numlua.h"
#include "cdflib.h"

/* C99 math functions: based on lhf's lmathx */

#define B1(name) \
  static int mathx_ ## name (lua_State *L) { \
    lua_pushboolean(L, name(luaL_checknumber(L, 1))); \
    return 1; \
  }

#define N1(name) \
  static int mathx_ ## name (lua_State *L) { \
    lua_pushnumber(L, name(luaL_checknumber(L, 1))); \
    return 1; \
  }

#define N2(name) \
  static int mathx_ ## name (lua_State *L) { \
    lua_pushnumber(L, name(luaL_checknumber(L, 1), luaL_checknumber(L, 2))); \
    return 1; \
  }

B1(isfinite)
B1(isinf)
B1(isnan)
B1(isnormal)
B1(signbit)
N1(acosh)
N1(asinh)
N1(atanh)
N1(cbrt)
N1(erf)
N1(erfc)
N1(exp2)
N1(expm1)
N1(lgamma)
N1(log1p)
N1(log2)
N1(logb)
N1(nearbyint)
N1(rint)
N1(round)
N1(tgamma)
N1(trunc)
N2(copysign)
N2(fdim)
N2(fmax)
N2(fmin)
N2(hypot)
N2(nextafter)
N2(remainder)

static int mathx_fpclassify (lua_State *L) {
  switch (fpclassify(luaL_checknumber(L, 1))) {
    case FP_INFINITE: lua_pushliteral(L, "inf"); break;
    case FP_NAN: lua_pushliteral(L, "nan"); break;
    case FP_NORMAL: lua_pushliteral(L, "normal"); break;
    case FP_SUBNORMAL: lua_pushliteral(L, "subnormal"); break;
    case FP_ZERO: lua_pushliteral(L, "zero"); break;
  }
  return 1;
}

static int mathx_fma (lua_State *L) {
  lua_pushnumber(L, fma(luaL_checknumber(L, 1),
        luaL_checknumber(L, 2), luaL_checknumber(L, 3)));
  return 1;
}

static int mathx_scalbn (lua_State *L) {
  lua_pushnumber(L, scalbn(luaL_checknumber(L, 1), luaL_checkinteger(L, 2)));
  return 1;
}


/* [[ Airy and Bessel functions ]] */

/* Adapted from Netlib's Amos library */
int zairy_(double *zr, double *zi, int *id, int *kode,
    double *air, double *aii, int *nz, int *ierr);
int zbiry_(double *zr, double *zi, int *id, int *kode,
    double *bir, double *bii, int *ierr);
int zbesh_(double *zr, double *zi, double *fnu, 
	int *kode, int *m, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesi_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesj_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesk_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesy_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz,
  double *cwrkr, double *cwrki, int *ierr);

static int mathx_airya (lua_State *L) {
  nl_Complex z = nl_checkcomplex(L, 1);
  int id = lua_toboolean(L, 2);
  int kode = lua_toboolean(L, 3) + 1;
  lua_Number zr = creal(z), zi = cimag(z);
  lua_Number air, aii;
  int nz, ierr;
  zairy_(&zr, &zi, &id, &kode, &air, &aii, &nz, &ierr);
  if (nz == 0 && (ierr == 0 || ierr == 3)) { /* push result? */
    nl_pushcomplex(L, air + aii * I);
    if (ierr == 0) lua_pushboolean(L, 1);
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    if (nz == 1) lua_pushliteral(L, "underflow");
    else switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


static int mathx_airyb (lua_State *L) {
  nl_Complex z = nl_checkcomplex(L, 1);
  int id = lua_toboolean(L, 2);
  int kode = lua_toboolean(L, 3) + 1;
  lua_Number zr = creal(z), zi = cimag(z);
  lua_Number bir, bii;
  int ierr;
  zbiry_(&zr, &zi, &id, &kode, &bir, &bii, &ierr);
  if (ierr == 0 || ierr == 3) { /* push result? */
    nl_pushcomplex(L, bir + bii * I);
    if (ierr == 0) lua_pushboolean(L, 1);
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


static int mathx_besselh (lua_State *L) {
  lua_Number nu = luaL_checknumber(L, 1);
  nl_Complex z = nl_checkcomplex(L, 2);
  int m = lua_toboolean(L, 3) + 1;
  int kode = lua_toboolean(L, 4) + 1;
  int n = luaL_optinteger(L, 5, 1);
  lua_Number zr = creal(z), zi = cimag(z);
  int nz, ierr;
  luaL_argcheck(L, nu >= 0, 1, "initial order must be non-negative");
  luaL_argcheck(L, cabs(z) > 0, 2, "argument must be different than zero");
  luaL_argcheck(L, n > 0, 5, "number of members must be positive");
  if (n == 1) { /* output single member? */
    lua_Number cyr, cyi;
    zbesh_(&zr, &zi, &nu, &kode, &m, &n, &cyr, &cyi, &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) /* push result? */
      nl_pushcomplex(L, cyr + cyi * I);
  }
  else {
    int one = 1, two = 2;
    nl_Buffer *buf = nl_getbuffer(L, 2 * n);
    zbesh_(&zr, &zi, &nu, &kode, &m, &n, buf->data.bnum, buf->data.bnum + n,
        &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) { /* normal return? */
      nl_Matrix *cy = nl_pushmatrix(L, 1, 1, &n, 1, n, NULL);
      DCOPY(&n, buf->data.bnum, &one, cy->data, &two);
      DCOPY(&n, buf->data.bnum + n, &one, cy->data + 1, &two);
    }
    nl_freebuffer(buf);
  }
  /* push status */
  if (nz == 0 && (ierr == 0 || ierr == 3)) {
    if (ierr == 0) lua_pushboolean(L, 1); /* no error */
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    if (nz > 0)
      lua_pushfstring(L, "underflow: first %d component(s) set to zero", nz);
    else switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


static int mathx_besseli (lua_State *L) {
  lua_Number nu = luaL_checknumber(L, 1);
  nl_Complex z = nl_checkcomplex(L, 2);
  int kode = lua_toboolean(L, 3) + 1;
  int n = luaL_optinteger(L, 4, 1);
  lua_Number zr = creal(z), zi = cimag(z);
  int nz, ierr;
  luaL_argcheck(L, nu >= 0, 1, "initial order must be non-negative");
  luaL_argcheck(L, n > 0, 4, "number of members must be positive");
  if (n == 1) { /* output single member? */
    lua_Number cyr, cyi;
    zbesi_(&zr, &zi, &nu, &kode, &n, &cyr, &cyi, &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) /* push result? */
      nl_pushcomplex(L, cyr + cyi * I);
  }
  else {
    int one = 1, two = 2;
    nl_Buffer *buf = nl_getbuffer(L, 2 * n);
    zbesi_(&zr, &zi, &nu, &kode, &n, buf->data.bnum, buf->data.bnum + n,
        &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) { /* push result? */
      nl_Matrix *cy = nl_pushmatrix(L, 1, 1, &n, 1, n, NULL);
      DCOPY(&n, buf->data.bnum, &one, cy->data, &two);
      DCOPY(&n, buf->data.bnum + n, &one, cy->data + 1, &two);
    }
    nl_freebuffer(buf);
  }
  /* push status */
  if (nz == 0 && (ierr == 0 || ierr == 3)) {
    if (ierr == 0) lua_pushboolean(L, 1); /* no error */
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    if (nz > 0)
      lua_pushfstring(L, "underflow: last %d component(s) set to zero", nz);
    else switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


static int mathx_besselj (lua_State *L) {
  lua_Number nu = luaL_checknumber(L, 1);
  nl_Complex z = nl_checkcomplex(L, 2);
  int kode = lua_toboolean(L, 3) + 1;
  int n = luaL_optinteger(L, 4, 1);
  lua_Number zr = creal(z), zi = cimag(z);
  int nz, ierr;
  luaL_argcheck(L, nu >= 0, 1, "initial order must be non-negative");
  luaL_argcheck(L, n > 0, 4, "number of members must be positive");
  if (n == 1) { /* output single member? */
    lua_Number cyr, cyi;
    zbesj_(&zr, &zi, &nu, &kode, &n, &cyr, &cyi, &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) /* push result? */
      nl_pushcomplex(L, cyr + cyi * I);
  }
  else {
    int one = 1, two = 2;
    nl_Buffer *buf = nl_getbuffer(L, 2 * n);
    zbesj_(&zr, &zi, &nu, &kode, &n, buf->data.bnum, buf->data.bnum + n,
        &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) { /* push result? */
      nl_Matrix *cy = nl_pushmatrix(L, 1, 1, &n, 1, n, NULL);
      DCOPY(&n, buf->data.bnum, &one, cy->data, &two);
      DCOPY(&n, buf->data.bnum + n, &one, cy->data + 1, &two);
    }
    nl_freebuffer(buf);
  }
  /* push status */
  if (nz == 0 && (ierr == 0 || ierr == 3)) {
    if (ierr == 0) lua_pushboolean(L, 1); /* no error */
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    if (nz > 0)
      lua_pushfstring(L, "underflow: last %d component(s) set to zero", nz);
    else switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


static int mathx_besselk (lua_State *L) {
  lua_Number nu = luaL_checknumber(L, 1);
  nl_Complex z = nl_checkcomplex(L, 2);
  int kode = lua_toboolean(L, 3) + 1;
  int n = luaL_optinteger(L, 4, 1);
  lua_Number zr = creal(z), zi = cimag(z);
  int nz, ierr;
  luaL_argcheck(L, nu >= 0, 1, "initial order must be non-negative");
  luaL_argcheck(L, cabs(z) > 0, 2, "argument cannot be zero");
  luaL_argcheck(L, n > 0, 4, "number of members must be positive");
  if (n == 1) { /* output single member? */
    lua_Number cyr, cyi;
    zbesk_(&zr, &zi, &nu, &kode, &n, &cyr, &cyi, &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) /* push result? */
      nl_pushcomplex(L, cyr + cyi * I);
  }
  else {
    int one = 1, two = 2;
    nl_Buffer *buf = nl_getbuffer(L, 2 * n);
    zbesk_(&zr, &zi, &nu, &kode, &n, buf->data.bnum, buf->data.bnum + n,
        &nz, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) { /* push result? */
      nl_Matrix *cy = nl_pushmatrix(L, 1, 1, &n, 1, n, NULL);
      DCOPY(&n, buf->data.bnum, &one, cy->data, &two);
      DCOPY(&n, buf->data.bnum + n, &one, cy->data + 1, &two);
    }
    nl_freebuffer(buf);
  }
  /* push status */
  if (nz == 0 && (ierr == 0 || ierr == 3)) {
    if (ierr == 0) lua_pushboolean(L, 1); /* no error */
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    if (nz > 0)
      lua_pushfstring(L, "underflow: first %d component(s) set to zero", nz);
    else switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


static int mathx_bessely (lua_State *L) {
  lua_Number nu = luaL_checknumber(L, 1);
  nl_Complex z = nl_checkcomplex(L, 2);
  int kode = lua_toboolean(L, 3) + 1;
  int n = luaL_optinteger(L, 4, 1);
  lua_Number zr = creal(z), zi = cimag(z);
  int nz, ierr;
  luaL_argcheck(L, nu >= 0, 1, "initial order must be non-negative");
  luaL_argcheck(L, n > 0, 4, "number of members must be positive");
  if (n == 1) { /* output single member? */
    lua_Number cyr, cyi;
    lua_Number cwrkr, cwrki;
    zbesy_(&zr, &zi, &nu, &kode, &n, &cyr, &cyi, &nz,
        &cwrkr, &cwrki, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) /* push result? */
      nl_pushcomplex(L, cyr + cyi * I);
  }
  else {
    int one = 1, two = 2;
    nl_Buffer *buf = nl_getbuffer(L, 2 * n);
    nl_Buffer *cwrk = nl_getbuffer(L, 2 * n);
    zbesy_(&zr, &zi, &nu, &kode, &n, buf->data.bnum, buf->data.bnum + n,
        &nz, cwrk->data.bnum, cwrk->data.bnum + n, &ierr);
    if (nz == 0 && (ierr == 0 || ierr == 3)) { /* push result? */
      nl_Matrix *cy = nl_pushmatrix(L, 1, 1, &n, 1, n, NULL);
      DCOPY(&n, buf->data.bnum, &one, cy->data, &two);
      DCOPY(&n, buf->data.bnum + n, &one, cy->data + 1, &two);
    }
    nl_freebuffer(buf);
    nl_freebuffer(cwrk);
  }
  /* push status */
  if (nz == 0 && (ierr == 0 || ierr == 3)) {
    if (ierr == 0) lua_pushboolean(L, 1); /* no error */
    else /* ierr == 3 */
      lua_pushliteral(L, "abs(z) large: loss of machine accuracy");
  }
  else {
    lua_pushnil(L);
    if (nz > 0)
      lua_pushfstring(L, "underflow: %d component(s) set to zero", nz);
    else switch (ierr) {
      case 1: lua_pushliteral(L, "input error"); break;
      case 2: lua_pushliteral(L, "overflow"); break;
      case 4:
        lua_pushliteral(L, "abs(z) too large: complete loss of accuracy");
        break;
      case 5: lua_pushliteral(L, "failed to converge"); break;
    }
  }
  return 2;
}


/* [[ extra ]] */

static int mathx_feq (lua_State *L) {
  lua_Number x1 = luaL_checknumber(L, 1);
  lua_Number x2 = luaL_checknumber(L, 2);
  lua_Number epsilon = luaL_optnumber(L, 3, DBL_EPSILON);
  int exp; /* frexp(max(|x1|,|x2|)) */
  frexp(fabs(x1) > fabs(x2) ? x1 : x2, &exp);
  lua_pushboolean(L, fabs(x1 - x2) <= ldexp(epsilon, exp));
  return 1;
}


#define LOGEPS (-36.043653389117) /* double precision */

/* computes log(1 + exp(x)) */
static int mathx_log1pe (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number d = (x > 0) ? -x : x;
  d = (d < LOGEPS) ? 0 : log1p(exp(d)); /* avoid calls if possible */
  if (x > 0) d += x;
  lua_pushnumber(L, d);
  return 1;
}

NUMLUA_API lua_Number nl_lse (lua_Number w1, lua_Number w2) {
  lua_Number d, w;
  if (!isfinite(w1)) return w2;
  if (!isfinite(w2)) return w1;
  if (w1 > w2) {
    d = w2 - w1; w = w1;
  }
  else {
    d = w1 - w2; w = w2;
  }
  if (d < LOGEPS) return w;
  return w + log1p(exp(d));
}

static int mathx_lse (lua_State *L) {
  lua_Number w1 = luaL_optnumber(L, 1, -HUGE_VAL);
  lua_Number w2 = luaL_optnumber(L, 2, -HUGE_VAL);
  lua_pushnumber(L, nl_lse(w1, w2));
  return 1;
}


/* {=====================================================================
 *    Statistical
 * ======================================================================} */

static int mathx_lbeta (lua_State *L) {
  lua_Number a = luaL_checknumber(L, 1);
  lua_Number b = luaL_checknumber(L, 2);
  lua_pushnumber(L, dlnbet(&a, &b));
  return 1;
}

static int mathx_beta (lua_State *L) {
  lua_Number a = luaL_checknumber(L, 1);
  lua_Number b = luaL_checknumber(L, 2);
  lua_pushnumber(L, exp(dlnbet(&a, &b)));
  return 1;
}


static int mathx_digamma (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  lua_pushnumber(L, psi(&x));
  return 1;
}

static lua_Number fchoose (lua_Number n, lua_Number k) {
  lua_Number a = n - k + 1;
  lua_Number b = k + 1;
  return -dlnbet(&a, &b) - log(n + 1);
}

static int mathx_lchoose (lua_State *L) {
  lua_Number n = luaL_checknumber(L, 1);
  lua_Number k = luaL_checknumber(L, 2);
  lua_Number c;
  if (k < 0) c = -HUGE_VAL;
  else if (k == 0) c = 0;
  /* k > 0 */
  else if (n < 0) c = fchoose(k - n - 1, k);
  else if (n < k) c = -HUGE_VAL;
  /* k <= n */
  else c = fchoose(n, k);
  lua_pushnumber(L, c);
  return 1;
}

static int mathx_choose (lua_State *L) {
  lua_Number n = luaL_checknumber(L, 1);
  lua_Number k = luaL_checknumber(L, 2);
  lua_Number c;
  if (k < 0) c = 0;
  else if (k == 0) c = 1;
  /* k > 0 */
  else if (n < 0) c = exp(fchoose(k - n - 1, k));
  else if (n < k) c = 0;
  /* k <= n */
  else c = exp(fchoose(n, k));
  lua_pushnumber(L, c);
  return 1;
}




static const luaL_Reg mathx_lib[] = {
  {"isfinite", mathx_isfinite},
  {"isinf", mathx_isinf},
  {"isnan", mathx_isnan},
  {"isnormal", mathx_isnormal},
  {"signbit", mathx_signbit},
  {"acosh", mathx_acosh},
  {"asinh", mathx_asinh},
  {"atanh", mathx_atanh},
  {"cbrt", mathx_cbrt},
  {"erf", mathx_erf},
  {"erfc", mathx_erfc},
  {"exp2", mathx_exp2},
  {"expm1", mathx_expm1},
  {"lgamma", mathx_lgamma},
  {"log1p", mathx_log1p},
  {"log2", mathx_log2},
  {"logb", mathx_logb},
  {"nearbyint", mathx_nearbyint},
  {"rint", mathx_rint},
  {"round", mathx_round},
  {"gamma", mathx_tgamma},
  {"trunc", mathx_trunc},
  {"copysign", mathx_copysign},
  {"fdim", mathx_fdim},
  {"fmax", mathx_fmax},
  {"fmin", mathx_fmin},
  {"hypot", mathx_hypot},
  {"nextafter", mathx_nextafter},
  {"remainder", mathx_remainder},
  {"fpclassify", mathx_fpclassify},
  {"fma", mathx_fma},
  {"scalbn", mathx_scalbn},
  {"airya", mathx_airya},
  {"airyb", mathx_airyb},
  {"besselh", mathx_besselh},
  {"besseli", mathx_besseli},
  {"besselj", mathx_besselj},
  {"besselk", mathx_besselk},
  {"bessely", mathx_bessely},
  {"feq", mathx_feq},
  {"log1pe", mathx_log1pe},
  {"lse", mathx_lse},
  {"lbeta", mathx_lbeta},
  {"beta", mathx_beta},
  {"digamma", mathx_digamma},
  {"lchoose", mathx_lchoose},
  {"choose", mathx_choose},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_mathx (lua_State *L) {
  luaL_newlib(L, mathx_lib);
  lua_pushnumber(L, NAN);
  lua_setfield(L, -2, "nan");
  lua_pushnumber(L, INFINITY);
  lua_setfield(L, -2, "inf");
  lua_pushnumber(L, DBL_EPSILON);
  lua_setfield(L, -2, "eps");
  lua_pushnumber(L, DBL_MAX);
  lua_setfield(L, -2, "realmax");
  lua_pushnumber(L, DBL_MIN);
  lua_setfield(L, -2, "realmin");
  lua_pushnumber(L, INT_MAX);
  lua_setfield(L, -2, "intmax");
  lua_pushnumber(L, INT_MIN);
  lua_setfield(L, -2, "intmin");
  return 1;
}

