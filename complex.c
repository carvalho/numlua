/* {=================================================================
 *
 * complex.c
 * C99 complex numbers for NumericLua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lua.h>
#include <lauxlib.h>
#include "numlua.h"

static int complex_mt_ = 0;
#define COMPLEX_MT ((void *)&complex_mt_)


/* {=====================================================================
 *    API
 * ======================================================================} */

/* Internal: meant for use on functions that have MT as env */

static nl_Complex *tocomplexP (lua_State *L, int narg) {
  nl_Complex *p = NULL;
  if (lua_type(L, narg) == LUA_TUSERDATA /* userdata? */
      && lua_getmetatable(L, narg)) { /* has metatable? */
    if (lua_rawequal(L, -1, lua_upvalueindex(1))) /* MT == upvalue? */
      p = lua_touserdata(L, narg);
    lua_pop(L, 1); /* MT */
  }
  return p;
}

static nl_Complex tocomplex (lua_State *L, int narg, int *iscomplex) {
  nl_Complex c, *p = tocomplexP(L, narg);
  if (p == NULL) { /* not complex? */
    c = lua_tonumber(L, narg);
    if (iscomplex != NULL) *iscomplex = lua_isnumber(L, narg);
  }
  else {
    c = *p;
    if (iscomplex != NULL) *iscomplex = 1;
  }
  return c;
}

static nl_Complex checkcomplex (lua_State *L, int narg) {
  int iscomplex;
  nl_Complex c = tocomplex(L, narg, &iscomplex);
  if (!iscomplex) nl_typeerror(L, narg, "number or complex");
  return c;
}

static nl_Complex *newcomplex (lua_State *L) {
  nl_Complex *z = (nl_Complex *) lua_newuserdata(L, sizeof(nl_Complex));
  lua_pushvalue(L, lua_upvalueindex(1));
  lua_setmetatable(L, -2);
  return z;
}

static nl_Complex *pushcomplex (lua_State *L, nl_Complex c) {
  nl_Complex *z = newcomplex(L);
  *z = c;
  return z;
}


/* API: get MT from registry */

NUMLUA_API nl_Complex nl_tocomplex (lua_State *L, int narg, int *iscomplex) {
  nl_Complex c, *p = NULL;
  *iscomplex = 0;
  if (lua_type(L, narg) == LUA_TUSERDATA /* userdata? */
      && lua_getmetatable(L, narg)) { /* has metatable? */
    nl_getmetatable(L, COMPLEX_MT);
    if (lua_rawequal(L, -1, -2)) /* right MT? */
      p = lua_touserdata(L, narg);
    lua_pop(L, 2); /* MTs */
  }
  if (p == NULL) { /* not complex? */
    c = lua_tonumber(L, narg);
    *iscomplex = (creal(c) != 0 || lua_isnumber(L, narg));
  }
  else {
    c = *p;
    *iscomplex = 1;
  }
  return c;
}

NUMLUA_API nl_Complex nl_checkcomplex (lua_State *L, int narg) {
  int iscomplex;
  nl_Complex c = nl_tocomplex(L, narg, &iscomplex);
  if (!iscomplex) nl_typeerror(L, narg, "number or complex");
  return c;
}

NUMLUA_API nl_Complex nl_optcomplex (lua_State *L, int narg, nl_Complex def) {
  return luaL_opt(L, nl_checkcomplex, narg, def);
}

NUMLUA_API nl_Complex *nl_newcomplex (lua_State *L) {
  nl_Complex *z = (nl_Complex *) lua_newuserdata(L, sizeof(nl_Complex));
  nl_getmetatable(L, COMPLEX_MT);
  lua_setmetatable(L, -2);
  return z;
}

NUMLUA_API nl_Complex *nl_pushcomplex (lua_State *L, nl_Complex c) {
  nl_Complex *z = nl_newcomplex(L);
  *z = c;
  return z;
}


NUMLUA_API lua_Number clogabs (nl_Complex c) {
  lua_Number r = fabs(creal(c));
  lua_Number i = fabs(cimag(c));
  lua_Number t;
  if ((r + i) == r) return log(r);
  if (i > r) { t = r; r = i; i = t; }
  t = i / r;
  return log(r) + 0.5 * log1p(t * t);
}



/* {=====================================================================
 *    Metamethods
 * ======================================================================} */

/* complex.newindex: there are three possible assignments:
 *  o c._ = x, in-place copy
 *  o c.r = x, creal(c) = x
 *  o c.i = x, cimag(c) = x
 */
static int complex_newindex (lua_State *L) {
  nl_Complex *c = (nl_Complex *) lua_touserdata(L, 1);
  size_t l;
  const char *k = lua_tolstring(L, 2, &l);
  int iscomplex;
  nl_Complex x = tocomplex(L, 3, &iscomplex);
  if (k == NULL) luaL_error(L, "complex index must be a string");
  if (l > 1) luaL_error(L, "invalid index: %s", k);
  if (!iscomplex)
    luaL_error(L, "number expected for assignment, got %s",
        luaL_typename(L, 3));
  switch (k[0]) {
    case '_': *c = x; break; /* assign */
    case 'r': *c = creal(x) + cimag(*c) * I; break; /* real assign */
    case 'i': *c = creal(*c) + creal(x) * I; break; /* imag assign */
    default: luaL_error(L, "invalid index: %s", k);
  }
  return 0;
}

static int complex_tostring (lua_State *L) {
  nl_Complex *c = (nl_Complex *) lua_touserdata(L, 1);
  if (!signbit(cimag(*c)))
    lua_pushfstring(L, "%f+%fi", creal(*c), cimag(*c));
  else
    lua_pushfstring(L, "%f%fi", creal(*c), cimag(*c));
  return 1;
}


/* {=====================================================================
 *    Functions
 * ======================================================================} */

#define cconj conj

#define F0(f) \
  static int complex_##f (lua_State *L) { \
    nl_Complex a = checkcomplex(L, 1); \
    lua_pushnumber(L, c##f(a)); \
    return 1;\
  }

#define F1(f) \
  static int complex_##f (lua_State *L) { \
    if (nl_inplace(L, 2)) { \
      nl_Complex *a = tocomplexP(L, 1); \
      if (a == NULL) nl_typeerror(L, 1, "complex"); \
      *a = c##f(*a); \
      lua_settop(L, 1); \
    } \
    else { \
      nl_Complex a = checkcomplex(L, 1); \
      nl_Complex *c = newcomplex(L); \
      *c = c##f(a); \
    } \
    return 1; \
  }

#define F2(f) \
  static int complex_##f (lua_State *L) { \
    nl_Complex b = checkcomplex(L, 2); \
    if (nl_inplace(L, 3)) { \
      nl_Complex *a = tocomplexP(L, 1); \
      if (a == NULL) nl_typeerror(L, 1, "complex"); \
      *a = c##f(*a,b); \
      lua_settop(L, 1); \
    } \
    else { \
      nl_Complex a = checkcomplex(L, 1); \
      nl_Complex *c = newcomplex(L); \
      *c = c##f(a,b); \
    } \
    return 1;\
  }

F1(unm)
F2(add)
F2(sub)
F2(mul)
F2(div)
F2(pow)

/* Methods */
F0(abs)
F0(logabs)
F0(arg)
F0(real)
F0(imag)
F1(acos)
F1(acosh)
F1(asin)
F1(asinh)
F1(atan)
F1(atanh)
F1(conj)
F1(cos)
F1(cosh)
F1(exp)
F1(log)
F1(proj)
F1(sin)
F1(sinh)
F1(sqrt)
F1(tan)
F1(tanh)


/* {=====================================================================
 *    Interface
 * ======================================================================} */

static int complexMT__call (lua_State *L) {
  /* complex, or one or two numbers on stack */
  nl_Complex *c;
  lua_remove(L, 1); /* class table */
  c = pushcomplex(L, checkcomplex(L, 1));
  if (lua_isnumber(L, 2)) /* set imag? */
    *c = creal(*c) + lua_tonumber(L, 2) * I;
  return 1;
}

static const luaL_Reg complex_func[] = {
  {"abs", complex_abs},
  {"logabs", complex_logabs},
  {"arg", complex_arg},
  {"real", complex_real},
  {"imag", complex_imag},
  {"acos", complex_acos},
  {"acosh", complex_acosh},
  {"asin", complex_asin},
  {"asinh", complex_asinh},
  {"atan", complex_atan},
  {"atanh", complex_atanh},
  {"conj", complex_conj},
  {"cos", complex_cos},
  {"cosh", complex_cosh},
  {"exp", complex_exp},
  {"log", complex_log},
  {"proj", complex_proj},
  {"sin", complex_sin},
  {"sinh", complex_sinh},
  {"sqrt", complex_sqrt},
  {"tan", complex_tan},
  {"tanh", complex_tanh},
  {"add", complex_add},
  {"sub", complex_sub},
  {"mul", complex_mul},
  {"div", complex_div},
  {"pow", complex_pow},
  {NULL, NULL}
};

static const luaL_Reg complex_mt[] = {
  {"__newindex", complex_newindex},
  {"__len", complex_abs},
  {"__tostring", complex_tostring},
  {"__unm", complex_unm},
  {"__add", complex_add},
  {"__sub", complex_sub},
  {"__mul", complex_mul},
  {"__div", complex_div},
  {"__pow", complex_pow},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_complex (lua_State *L) {
  nl_Complex i = 1 * I;
  /* complex metatable */
  lua_newtable(L);
  lua_pushliteral(L, COMPLEX_LIBNAME);
  lua_setfield(L, -2, "__type");
  lua_pushlightuserdata(L, COMPLEX_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX);
  /* push metamethods */
  lua_pushvalue(L, -1);
  nl_register(L, complex_mt, 1);
  /* class table */
  lua_newtable(L);
  lua_pushvalue(L, -2); /* metatable */
  nl_register(L, complex_func, 1);
  pushcomplex(L, i);
  /* fix metatable for `i` (no upvalues here) */
  lua_pushvalue(L, -3); /* metatable */
  lua_setmetatable(L, -2);
  lua_pushvalue(L, -1);
  lua_setfield(L, -3, "i");
  lua_setfield(L, -2, "j");
  /* push class table for __index */
  lua_pushvalue(L, -1);
  lua_setfield(L, -3, "__index");
  /* push class metatable */
  lua_createtable(L, 0, 2);
  lua_pushvalue(L, -3); /* metatable */
  lua_pushcclosure(L, complexMT__call, 1);
  lua_setfield(L, -2, "__call");
  lua_pushvalue(L, -2); /* class */
  lua_setfield(L, -2, "__metatable"); /* protect */
  lua_setmetatable(L, -2);
  return 1;
}

