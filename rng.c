/* {=================================================================
 *
 * luarng.c
 * Random number generator (RNG) library for Lua
 * Luis Carvalho (lexcarvalho @ gmail.com)
 * See Copyright Notice in luarng.h
 * $Id: luarng.c,v 1.5 2006-04-29 17:26:51 carvalho Exp $
 *
 * ==================================================================} */

#include <lauxlib.h>
#include <math.h>
#include "numlua.h"
#include "rng.h"

static int rng_mt_ = 0;
#define RNG_MT ((void *)&rng_mt_)
#define getrng(L) ((nl_RNG *) lua_touserdata((L), lua_upvalueindex(1)))

toudata(rng, RNG_LIBNAME, nl_RNG)
libtoudata(rng, RNG_LIBNAME, RNG_MT, nl_RNG)

/* {=================================================================
 *      Class metamethods
 * ==================================================================} */

static int rng__tostring (lua_State *L) {
  lua_pushfstring(L, RNG_LIBNAME ": %p", lua_touserdata(L, 1));
  return 1;
}

/* {=================================================================
 *      Class methods
 * ==================================================================} */

static int new_rng (lua_State *L) {
  unsigned long s = luaL_optlong(L, 1, RNG_SEED);
  nl_RNG *r = lua_newuserdata(L, sizeof(nl_RNG));
  init_genrand(r, s);
  lua_pushvalue(L, LUA_ENVIRONINDEX);
  lua_setmetatable(L, -2);
  return 1;
}

static int copy_rng (lua_State *L) {
  nl_RNG *c = checkrng(L, 1);
  nl_RNG *r = (nl_RNG *) lua_newuserdata(L, sizeof(nl_RNG));
  lua_pushvalue(L, LUA_ENVIRONINDEX);
  lua_setmetatable(L, -2);
  *r = *c;
  return 1;
}

static int seed_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  init_genrand(r, luaL_optlong(L, 2, RNG_SEED));
  lua_settop(L, 1);
  return 1;
}

static int rseed_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  init_genrand(r, luaL_optlong(L, 1, RNG_SEED));
  return 0;
}


static int seedarray_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  unsigned long init_key[RNG_MAXSTATES];
  int pos = 0;
  luaL_checktype(L, 2, LUA_TTABLE);
  /* traverse table and store keys in array */
  lua_pushnil(L);
  while (lua_next(L, -2) != 0) {
    init_key[pos++] = luaL_optlong(L, -1, RNG_SEED);
    lua_pop(L, 1);
  }
  init_by_array(r, init_key, pos);
  lua_settop(L, 1);
  return 1;
}

/* {=================================================================
 *      Main routines
 * ==================================================================} */

#define checknp(L,x) \
  if ((x) <= 0) luaL_error((L), "nonpositive parameter: %f", (x))
#define checkneg(L,x) \
  if ((x) < 0) luaL_error((L), "negative parameter: %f", (x))

static int beta_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number a = luaL_checknumber(L, 2);
  lua_Number b = luaL_checknumber(L, 3);
  checknp(L, a);
  checknp(L, b);
  lua_pushnumber(L, genbet(r, a, b));
  return 1;
}

static int rbeta_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number a = luaL_checknumber(L, 1);
  lua_Number b = luaL_checknumber(L, 2);
  checknp(L, a);
  checknp(L, b);
  lua_pushnumber(L, genbet(r, a, b));
  return 1;
}

static int chisq_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number xnonc = luaL_optnumber(L, 3, 0);
  checknp(L, df);
  checkneg(L, xnonc);
  lua_pushnumber(L, (xnonc == 0) ? genchi(r, df) : gennch(r, df, xnonc));
  return 1;
}

static int rchisq_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number df = luaL_checknumber(L, 1);
  lua_Number xnonc = luaL_optnumber(L, 2, 0);
  checknp(L, df);
  checkneg(L, xnonc);
  lua_pushnumber(L, (xnonc == 0) ? genchi(r, df) : gennch(r, df, xnonc));
  return 1;
}

static int exp_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number av = luaL_optnumber(L, 2, 1);
  lua_pushnumber(L, genexp(r, av));
  return 1;
}

static int rexp_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number av = luaL_optnumber(L, 1, 1);
  lua_pushnumber(L, genexp(r, av));
  return 1;
}

static int f_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number dfn = luaL_checknumber(L, 2);
  lua_Number dfd = luaL_checknumber(L, 3);
  lua_Number xnonc = luaL_optnumber(L, 4, 0);
  checknp(L, dfn);
  checknp(L, dfd);
  checkneg(L, xnonc);
  lua_pushnumber(L, (xnonc == 0) ? genf(r, dfn, dfd)
      : gennf(r, dfn, dfd, xnonc));
  return 1;
}

static int rf_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number dfn = luaL_checknumber(L, 1);
  lua_Number dfd = luaL_checknumber(L, 2);
  lua_Number xnonc = luaL_optnumber(L, 3, 0);
  checknp(L, dfn);
  checknp(L, dfd);
  checkneg(L, xnonc);
  lua_pushnumber(L, (xnonc == 0) ? genf(r, dfn, dfd)
      : gennf(r, dfn, dfd, xnonc));
  return 1;
}

static int gamma_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number a = luaL_checknumber(L, 2);
  lua_Number s = luaL_optnumber(L, 3, 1);
  lua_pushnumber(L, gengam(r, s, a));
  return 1;
}

static int rgamma_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number a = luaL_checknumber(L, 1);
  lua_Number s = luaL_optnumber(L, 2, 1);
  lua_pushnumber(L, gengam(r, s, a));
  return 1;
}

static int norm_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number av = luaL_optnumber(L, 2, 0);
  lua_Number v = luaL_optnumber(L, 3, 1);
  checknp(L, v);
  lua_pushnumber(L, gennor(r, av, v));
  return 1;
}

static int rnorm_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number av = luaL_optnumber(L, 1, 0);
  lua_Number v = luaL_optnumber(L, 2, 1);
  checknp(L, v);
  lua_pushnumber(L, gennor(r, av, v));
  return 1;
}

static int unif_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number low = luaL_optnumber(L, 2, 0);
  lua_Number high = luaL_optnumber(L, 3, 1);
  if (low > high)
    luaL_error(L, "inconsistent parameters: %f > %f", low, high);
  lua_pushnumber(L, genunf(r, low, high));
  return 1;
}

static int runif_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number low = luaL_optnumber(L, 1, 0);
  lua_Number high = luaL_optnumber(L, 2, 1);
  if (low > high)
    luaL_error(L, "inconsistent parameters: %f > %f", low, high);
  lua_pushnumber(L, genunf(r, low, high));
  return 1;
}

static int binom_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  int n = luaL_checkinteger(L, 2);
  lua_Number p = luaL_checknumber(L, 3);
  checkneg(L, n);
  if (p <= 0 || p >= 1) luaL_error(L, "parameter is out of range: %f", p);
  lua_pushinteger(L, ignbin(r, n, p));
  return 1;
}

static int rbinom_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  int n = luaL_checkinteger(L, 1);
  lua_Number p = luaL_checknumber(L, 2);
  checkneg(L, n);
  if (p <= 0 || p >= 1) luaL_error(L, "parameter is out of range: %f", p);
  lua_pushinteger(L, ignbin(r, n, p));
  return 1;
}

static int nbinom_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  int n = luaL_checkinteger(L, 2);
  lua_Number p = luaL_checknumber(L, 3);
  checkneg(L, n);
  if (p <= 0 || p >= 1) luaL_error(L, "parameter is out of range: %f", p);
  lua_pushinteger(L, ignnbn(r, n, p));
  return 1;
}

static int rnbinom_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  int n = luaL_checkinteger(L, 1);
  lua_Number p = luaL_checknumber(L, 2);
  checkneg(L, n);
  if (p <= 0 || p >= 1) luaL_error(L, "parameter is out of range: %f", p);
  lua_pushinteger(L, ignnbn(r, n, p));
  return 1;
}

static int pois_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  lua_Number mu = luaL_checknumber(L, 2);
  checkneg(L, mu);
  lua_pushinteger(L, ignpoi(r, mu));
  return 1;
}

static int rpois_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number mu = luaL_checknumber(L, 1);
  checkneg(L, mu);
  lua_pushinteger(L, ignpoi(r, mu));
  return 1;
}

static int unifint_rng (lua_State *L) {
  nl_RNG *r = checkrng(L, 1);
  int low = luaL_optinteger(L, 2, 0);
  int high = luaL_optinteger(L, 3, 0x7ffffffeUL);
  if (low > high)
    luaL_error(L, "inconsistent parameters: %d > %d", low, high);
  lua_pushinteger(L, ignuin(r, low, high));
  return 1;
}

static int runifint_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  int low = luaL_optinteger(L, 1, 0);
  int high = luaL_optinteger(L, 2, 0x7ffffffeUL);
  if (low > high)
    luaL_error(L, "inconsistent parameters: %d > %d", low, high);
  lua_pushinteger(L, ignuin(r, low, high));
  return 1;
}

/* {=================================================================
 *      Interface
 * ==================================================================} */

static const luaL_reg rng_lib[] = {
  {"new", new_rng},
  {"rseed", rseed_rng},
  /* deviates (common rng) */
  {"rbeta", rbeta_rng},
  {"rchisq", rchisq_rng},
  {"rexp", rexp_rng},
  {"rf", rf_rng},
  {"rgamma", rgamma_rng},
  {"rnorm", rnorm_rng},
  {"runif", runif_rng},
  {"rbinom", rbinom_rng},
  {"rnbinom", rnbinom_rng},
  {"rpois", rpois_rng},
  {"runifint", runifint_rng},
  {NULL, NULL}
};

static const luaL_reg rng_index[] = {
  {"copy", copy_rng},
  {"seed", seed_rng},
  {"seedarray", seedarray_rng},
  /* deviates */
  {"beta", beta_rng},
  {"chisq", chisq_rng},
  {"exp", exp_rng},
  {"f", f_rng},
  {"gamma", gamma_rng},
  {"norm", norm_rng},
  {"unif", unif_rng},
  {"binom", binom_rng},
  {"nbinom", nbinom_rng},
  {"pois", pois_rng},
  {"unifint", unifint_rng},
  {NULL, NULL}
};

static const luaL_Reg rng_deps[] = {
  {NUMLUA_LIBNAME,  luaopen_numlua_base},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_rng (lua_State *L) {
  if (nl_requiredeps(L, RNG_LIBNAME, rng_deps)) return 1;
  /* new environment/metatable */
  lua_newtable(L);
  lua_pushlightuserdata(L, RNG_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX);
  /* set as current environment */
  lua_pushvalue(L, -1);
  lua_replace(L, LUA_ENVIRONINDEX);
  /* push metamethods */
  lua_newtable(L); /* __index */
  nl_register(L, rng_index, 0);
  lua_setfield(L, -2, "__index");
  lua_pushcfunction(L, rng__tostring);
  lua_setfield(L, -2, "__tostring");
  lua_pushliteral(L, RNG_LIBNAME);
  lua_setfield(L, -2, "__type");
  /* class table */
  luaL_register(L, RNG_LIBNAME, rng_lib);
  init_genrand((nl_RNG *) lua_newuserdata(L, sizeof(nl_RNG)), RNG_SEED);
  nl_register(L, rng_lib, 1);
  return 1;
}

