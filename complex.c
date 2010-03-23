#include <lua.h>
#include <lauxlib.h>
#include <string.h> /* strncmp */
#include "numlua.h"

#define COMPLEX_LBL_REAL "real"
#define COMPLEX_LBL_IMAG "imag"
#define COMPLEX_MAX_LBL  (4)

static int complex_mt_ = 0;
#define COMPLEX_MT ((void *)&complex_mt_)

/* Internal: meant for use on functions that have MT as env */

static nl_Complex *tocomplexP (lua_State *L, int narg) {
  nl_Complex *p = NULL;
  if (lua_getmetatable(L, narg)) { /* has metatable? */
    if (lua_rawequal(L, -1, LUA_ENVIRONINDEX)) /* MT == env? */
      p = lua_touserdata(L, narg);
    lua_pop(L, 1); /* MT */
  }
  return p;
}

static nl_Complex tocomplex (lua_State *L, int narg, int *iscomplex) {
  nl_Complex c, *p = tocomplexP(L, narg);
  *iscomplex = 0;
  if (p == NULL) { /* not complex? */
    c.re = lua_tonumber(L, narg); c.im = 0;
    *iscomplex = (c.re != 0 || lua_isnumber(L, narg));
  }
  else {
    c = *p;
    *iscomplex = 1;
  }
  return c;
}

static nl_Complex checkcomplex (lua_State *L, int narg) {
  int iscomplex;
  nl_Complex c = tocomplex(L, narg, &iscomplex);
  if (!iscomplex) luaL_typerror(L, narg, "number or complex");
  return c;
}

static nl_Complex *newcomplex (lua_State *L) {
  nl_Complex *z = (nl_Complex *) lua_newuserdata(L, sizeof(nl_Complex));
  lua_pushvalue(L, LUA_ENVIRONINDEX);
  lua_setmetatable(L, -2);
  return z;
}

static nl_Complex *pushcomplex (lua_State *L, nl_Complex c) {
  nl_Complex *z = newcomplex(L);
  *z = c;
  return z;
}


/* API: get MT from registry */

nl_Complex nl_tocomplex (lua_State *L, int narg, int *iscomplex) {
  nl_Complex c, *p = NULL;
  *iscomplex = 0;
  if (lua_getmetatable(L, narg)) { /* has metatable? */
    lua_pushlightuserdata(L, COMPLEX_MT);
    lua_rawget(L, LUA_REGISTRYINDEX);
    if (lua_rawequal(L, -1, -2)) /* right MT? */
      p = lua_touserdata(L, narg);
    lua_pop(L, 2); /* MTs */
  }
  if (p == NULL) { /* not complex? */
    c.re = lua_tonumber(L, narg); c.im = 0;
    *iscomplex = (c.re != 0 || lua_isnumber(L, narg));
  }
  else {
    c = *p;
    *iscomplex = 1;
  }
  return c;
}

nl_Complex nl_checkcomplex (lua_State *L, int narg) {
  int iscomplex;
  nl_Complex c = nl_tocomplex(L, narg, &iscomplex);
  if (!iscomplex) luaL_typerror(L, narg, "number or complex");
  return c;
}

nl_Complex nl_optcomplex (lua_State *L, int narg, nl_Complex def) {
  return luaL_opt(L, nl_checkcomplex, narg, def);
}

nl_Complex *nl_newcomplex (lua_State *L) {
  nl_Complex *z = (nl_Complex *) lua_newuserdata(L, sizeof(nl_Complex));
  lua_pushlightuserdata(L, COMPLEX_MT);
  lua_rawget(L, LUA_REGISTRYINDEX);
  lua_setmetatable(L, -2);
  return z;
}

nl_Complex *nl_pushcomplex (lua_State *L, nl_Complex c) {
  nl_Complex *z = nl_newcomplex(L);
  *z = c;
  return z;
}


/* Auxiliar */

NUMLUA_API lua_Number nl_abscomplex (nl_Complex c) {
  lua_Number r = fabs(c.re);
  lua_Number i = fabs(c.im);
  lua_Number t;
  if ((r + i) == r) return r;
  if (i > r) { t = r; r = i; i = t; }
  t = i / r;
  return r * sqrt(1.0 + t * t);
}

NUMLUA_API void nl_expcomplex (nl_Complex *c, nl_Complex a) {
  lua_Number t = exp(a.re);
  c->re = t * cos(a.im);
  c->im = t * sin(a.im);
}

NUMLUA_API void nl_logcomplex (nl_Complex *c, nl_Complex a) {
  lua_Number t = atan2(a.im, a.re);
  c->re = log(nl_abscomplex(a));
  c->im = t;
}

NUMLUA_API void nl_sincomplex (nl_Complex *c, nl_Complex a) {
  lua_Number t = cos(a.re) * sinh(a.im);
  c->re = sin(a.re) * cosh(a.im);
  c->im = t;
}

NUMLUA_API void nl_coscomplex (nl_Complex *c, nl_Complex a) {
  lua_Number t = -sin(a.re) * sinh(a.im);
  c->re = cos(a.re) * cosh(a.im);
  c->im = t;
}

NUMLUA_API void nl_sqrtcomplex (nl_Complex *c, nl_Complex a) {
  lua_Number t = nl_abscomplex(a);
  if (t == 0.0) c->re = c->im = 0.0;
  else if (a.re > 0.0) {
    c->re = sqrt(0.5 * (t + a.re));
    c->im = a.im / c->re / 2.0;
  }
  else {
    t = sqrt(0.5 * (t - a.re));
    if (a.im < 0.0) t = -t;
    c->re = a.im / t / 2.0;
    c->im = t;
  }
}

NUMLUA_API void nl_mulcomplex (nl_Complex *c, nl_Complex a, nl_Complex b) {
  lua_Number re = a.re;
  c->re = re * b.re - a.im * b.im;
  c->im = re * b.im + a.im * b.re;
}

NUMLUA_API void nl_divcomplex (nl_Complex *c, nl_Complex a, nl_Complex b) {
  /* adapted from libF77 */
  lua_Number ratio, den, abr, abi;
  if ((abr = b.re) < 0.0) abr = -abr;
  if ((abi = b.im) < 0.0) abi = -abi;
  if (abr <= abi) {
    if (abi == 0) {
      if (a.im != 0 || a.re != 0) abi = 1.0;
      c->im = c->re = abi / abr;
      return;
    }
    ratio = b.re / b.im;
    den = b.im * (1 + ratio * ratio);
    c->re = (a.re * ratio + a.im) / den;
    c->im = (a.im * ratio - a.re) / den;
  }
  else {
    ratio = b.im / b.re;
    den = b.re * (1 + ratio * ratio);
    c->re = (a.re + a.im * ratio) / den;
    c->im = (a.im - a.re * ratio) / den;
  }
}

NUMLUA_API void nl_ipowcomplex (nl_Complex *f, nl_Complex c, int e) {
  /* adapted from libF77 */
  unsigned long u;
  lua_Number t;
  nl_Complex x, q = {1.0, 0.0};
  if (e == 0) {
    *f = q;
    return;
  }
  if (e < 0) {
    e = -e;
    nl_divcomplex(&x, q, c); /* q = (1,0), invert c to x */
  }
  else x = c;
  for (u = e; ; ) {
    if (u & 01) {
      t = q.re * x.re - q.im * x.im;
      q.im = q.re * x.im + q.im * x.re;
      q.re = t;
    }
    if (u >>= 1) {
      t = x.re * x.re - x.im * x.im;
      x.im = 2 * x.re * x.im;
      x.re = t;
    }
    else break;
  }
  *f = q;
}

NUMLUA_API void nl_powcomplex (nl_Complex *f, nl_Complex a, nl_Complex b) {
  lua_Number logr = log(nl_abscomplex(a));
  lua_Number logi = atan2(a.im, a.re);
  lua_Number x = exp(logr * b.re - logi * b.im);
  lua_Number y = logr * b.im + logi * b.re;
  f->re = x * cos(y);
  f->im = x * sin(y);
}


/* Metamethods */

static int complex__index (lua_State *L) {
  nl_Complex *c = (nl_Complex *) lua_touserdata(L, 1);
  const char *k = lua_tostring(L, 2);
  if (k == NULL) luaL_error(L, "complex index must be a string");
  if (strncmp(k, COMPLEX_LBL_REAL, COMPLEX_MAX_LBL) == 0)
    lua_pushnumber(L, c->re);
  else if (strncmp(k, COMPLEX_LBL_IMAG, COMPLEX_MAX_LBL) == 0)
    lua_pushnumber(L, c->im);
  else
    lua_rawget(L, LUA_ENVIRONINDEX); /* class lookup */
  return 1;
}

static int complex__newindex (lua_State *L) {
  nl_Complex *c = (nl_Complex *) lua_touserdata(L, 1);
  const char *k = lua_tostring(L, 2);
  lua_Number x = lua_tonumber(L, 3);
  if (k == NULL) luaL_error(L, "complex index must be a string");
  if (x == 0 && !lua_isnumber(L, 3))
    luaL_error(L, "number expected for assignment, got %s", luaL_typename(L, 3));
  if (strncmp(k, COMPLEX_LBL_REAL, COMPLEX_MAX_LBL) == 0)
    c->re = x;
  else if (strncmp(k, COMPLEX_LBL_IMAG, COMPLEX_MAX_LBL) == 0)
    c->im = x;
  else
    luaL_error(L, "complex has no `%s' property", k);
  return 0;
}

static int complex__tostring (lua_State *L) {
  nl_Complex *c = (nl_Complex *) lua_touserdata(L, 1);
  if (c->im < 0)
    lua_pushfstring(L, "%f - %fi", c->re, -c->im);
  else
    lua_pushfstring(L, "%f + %fi", c->re, c->im);
  return 1;
}

/* absolute value */
static int complex__len (lua_State *L) {
  nl_Complex *c = (nl_Complex *) lua_touserdata(L, 1);
  lua_pushnumber(L, nl_abscomplex(*c));
  return 1;
}

static int complex__add (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex b = checkcomplex(L, 2);
  nl_Complex *c = newcomplex(L);
  c->re = a.re + b.re;
  c->im = a.im + b.im;
  return 1;
}

static int complex__sub (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex b = checkcomplex(L, 2);
  nl_Complex *c = newcomplex(L);
  c->re = a.re - b.re;
  c->im = a.im - b.im;
  return 1;
}

static int complex__unm (lua_State *L) {
  nl_Complex *a = (nl_Complex *) lua_touserdata(L, 1);
  nl_Complex *c = newcomplex(L);
  c->re = -a->re;
  c->im = -a->im;
  return 1;
}



static int complex__mul (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex b = checkcomplex(L, 2);
  nl_mulcomplex(newcomplex(L), a, b);
  return 1;
}

static int complex__div (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex b = checkcomplex(L, 2);
  nl_divcomplex(newcomplex(L), a, b);
  return 1;
}

static int complex__pow (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex b = checkcomplex(L, 2);
  nl_Complex *c = newcomplex(L);
  if (b.im == 0.0 && b.re == floor(b.re)) { /* integer pow? */
    lua_Number x = b.re;
    int p;
    lua_number2int(p, x);
    nl_ipowcomplex(c, a, p);
  }
  else nl_powcomplex(c, a, b);
  return 1;
}

/* Methods */

static int complex_real (lua_State *L) {
  nl_Complex c = checkcomplex(L, 1);
  lua_pushnumber(L, c.re);
  return 1;
}

static int complex_imag (lua_State *L) {
  nl_Complex c = checkcomplex(L, 1);
  lua_pushnumber(L, c.im);
  return 1;
}

static int complex_conj (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex *c = newcomplex(L);
  c->re = a.re; c->im = -a.im;
  return 1;
}

static int complex_abs (lua_State *L) {
  lua_pushnumber(L, nl_abscomplex(checkcomplex(L, 1)));
  return 1;
}

static int complex_arg (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  lua_pushnumber(L, atan2(a.im, a.re));
  return 1;
}

static int complex_exp (lua_State *L) {
  nl_expcomplex(newcomplex(L), checkcomplex(L, 1));
  return 1;
}

static int complex_log (lua_State *L) {
  nl_logcomplex(newcomplex(L), checkcomplex(L, 1));
  return 1;
}

static int complex_sin (lua_State *L) {
  nl_sincomplex(newcomplex(L), checkcomplex(L, 1));
  return 1;
}

static int complex_cos (lua_State *L) {
  nl_coscomplex(newcomplex(L), checkcomplex(L, 1));
  return 1;
}

static int complex_sqrt (lua_State *L) {
  nl_sqrtcomplex(newcomplex(L), checkcomplex(L, 1));
  return 1;
}

static int complex_copy (lua_State *L) {
  nl_Complex a = checkcomplex(L, 1);
  nl_Complex *c = tocomplexP(L, 2);
  if (c == NULL) /* no destination? */
    c = pushcomplex(L, a);
  else
    *c = a;
  return 1;
}


/* interface */

static int complexMT__call (lua_State *L) {
  /* complex, or one or two numbers on stack */
  int iscomplex;
  nl_Complex *c = pushcomplex(L, tocomplex(L, 2, &iscomplex));
  if (lua_isnumber(L, 3)) c->im = lua_tonumber(L, 3);
  return 1;
}

static const luaL_reg complex_func[] = {
  {"real", complex_real},
  {"imag", complex_imag},
  {"conj", complex_conj},
  {"abs", complex_abs},
  {"arg", complex_arg},
  {"exp", complex_exp},
  {"log", complex_log},
  {"sin", complex_sin},
  {"cos", complex_cos},
  {"sqrt", complex_sqrt},
  {"copy", complex_copy},
  {NULL, NULL}
};

static const luaL_reg complex_mt[] = {
  {"__newindex", complex__newindex},
  {"__len", complex__len},
  {"__tostring", complex__tostring},
  {"__add", complex__add},
  {"__sub", complex__sub},
  {"__unm", complex__unm},
  {"__mul", complex__mul},
  {"__div", complex__div},
  {"__pow", complex__pow},
  {NULL, NULL}
};

static const luaL_Reg complex_deps[] = {
  {NUMLUA_LIBNAME, luaopen_numlua_base},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_complex (lua_State *L) {
  nl_Complex i = {0, 1};
  if (nl_requiredeps(L, COMPLEX_LIBNAME, complex_deps)) return 1;
  /* complex MT */
  lua_newtable(L);
  lua_pushlightuserdata(L, COMPLEX_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX);
  /* set as current environment */
  lua_pushvalue(L, -1);
  lua_replace(L, LUA_ENVIRONINDEX);
  /* push metamethods */
  nl_register(L, complex_mt, 0);
  lua_pushliteral(L, COMPLEX_LIBNAME);
  lua_setfield(L, -2, "__type");
  /* class table */
  luaL_register(L, COMPLEX_LIBNAME, complex_func);
  pushcomplex(L, i);
  lua_pushvalue(L, -1);
  lua_setfield(L, -3, "i");
  lua_setfield(L, -2, "j");
  /* push class table as environment to __index */
  lua_pushcfunction(L, complex__index);
  lua_pushvalue(L, -2);
  lua_setfenv(L, -2);
  lua_setfield(L, -3, "__index");
  /* push class metatable */
  lua_createtable(L, 0, 2);
  lua_pushcfunction(L, complexMT__call);
  lua_setfield(L, -2, "__call");
  lua_pushvalue(L, -2); /* class */
  lua_setfield(L, -2, "__metatable"); /* protect */
  lua_setmetatable(L, -2);
  return 1;
}

