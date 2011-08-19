/* {=================================================================
 *
 * rng.c
 * Random number generator (RNG) library for NumericLua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lauxlib.h>
#include <math.h>
#include "numlua.h"
#include "rng.h"

/* TODO [wishlist]: Wishart, multinomial */

static int rng_mt_ = 0;
#define RNG_MT ((void *)&rng_mt_)
#define getrng(L) ((nl_RNG *) lua_touserdata((L), lua_upvalueindex(1)))

#define checkrvector(L,m,n) \
  luaL_argcheck(L,!(m)->section && !(m)->iscomplex,(n),"real vector expected")

/* pushes new rng object in stack */
static nl_RNG *newrng (lua_State *L); /* forward: needs rng_lib */


/* {=====================================================================
 *    Metamethods
 * ======================================================================} */

static int rng__len (lua_State *L) {
  lua_getuservalue(L, 1);
  return 1;
}

static int rng__tostring (lua_State *L) {
  lua_pushfstring(L, RNG_LIBNAME ": %p", lua_touserdata(L, 1));
  return 1;
}


/* {=====================================================================
 *    Core
 * ======================================================================} */

static int new_rng (lua_State *L) {
  init_genrand(newrng(L), luaL_optlong(L, 1, RNG_SEED));
  return 1;
}


static int copy_rng (lua_State *L) {
  nl_RNG *c = getrng(L);
  nl_RNG *r = newrng(L);
  *r = *c;
  return 1;
}

static int seed_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  if (lua_isnoneornil(L, 1))
    init_genrand(r, RNG_SEED);
  else if (lua_isnumber(L, 1)) /* seed? */
    init_genrand(r, lua_tointeger(L, 1));
  else { /* vector */
    unsigned long initkey[RNG_MAXSTATES];
    int i, k;
    lua_Number *e;
    nl_Matrix *m = nl_checkmatrix(L, 1);
    checkrvector(L, m, 1);
    for (i = 0, e = m->data; i < m->size; i++, e += m->stride) {
      lua_number2int(k, *e);
      initkey[i] = (unsigned long) k;
    }
    init_by_array(r, initkey, m->size);
  }
  return 0;
}


/* {=====================================================================
 *    Deviate methods
 * ======================================================================} */

#define checknp(L,x) \
  if ((x) <= 0) luaL_error((L), "nonpositive parameter: %f", (x))
#define checkneg(L,x) \
  if ((x) < 0) luaL_error((L), "negative parameter: %f", (x))

#define setdeviate(type,exp,n) \
  do { \
    lua_settop(L, n); \
    if (lua_isnil(L, n)) /* no destination */ \
      lua_push ## type(L, exp); \
    else { \
      int i; lua_Number *e; \
      nl_Matrix *v = nl_checkmatrix(L, n); \
      checkrvector(L, v, n); \
      for (i=0, e=v->data; i<v->size; i++, e+=v->stride) \
        *e = exp; \
    } \
  } while (0) 


static int rbeta_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number a = luaL_checknumber(L, 1);
  lua_Number b = luaL_checknumber(L, 2);
  checknp(L, a);
  checknp(L, b);
  setdeviate(number, genbet(r, a, b), 3);
  return 1;
}


static int rchisq_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number df = luaL_checknumber(L, 1);
  lua_Number xnonc = luaL_optnumber(L, 2, 0);
  checknp(L, df);
  checkneg(L, xnonc);
  setdeviate(number,
      (xnonc == 0) ? genchi(r, df) : gennch(r, df, xnonc), 3);
  return 1;
}


static int rexp_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number av = luaL_optnumber(L, 1, 1);
  setdeviate(number, genexp(r, av), 2);
  return 1;
}


static int rf_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number dfn = luaL_checknumber(L, 1);
  lua_Number dfd = luaL_checknumber(L, 2);
  lua_Number xnc = luaL_optnumber(L, 3, 0);
  checknp(L, dfn);
  checknp(L, dfd);
  checkneg(L, xnc);
  setdeviate(number,
      (xnc == 0) ? genf(r, dfn, dfd) : gennf(r, dfn, dfd, xnc), 4);
  return 1;
}


static int rgamma_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number a = luaL_checknumber(L, 1);
  lua_Number s = luaL_optnumber(L, 2, 1);
  setdeviate(number, gengam(r, s, a), 3);
  return 1;
}


static int rnorm_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number mean = luaL_optnumber(L, 1, 0);
  lua_Number sd = luaL_optnumber(L, 2, 1);
  checknp(L, sd);
  setdeviate(number, gennor(r, mean, sd), 3);
  return 1;
}


static int rmvnorm_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  nl_Matrix *m = nl_checkmatrix(L, 1);
  nl_Matrix *S = nl_checkmatrix(L, 2);
  nl_Matrix *u;
  int i, n = m->size;
  lua_Number *em, *ev, *eu;
  /* check args */
  checkrvector(L, m, 1);
  luaL_argcheck(L, !S->iscomplex, 2, "real matrix expected");
  if (S->ndims == 1) {
    luaL_argcheck(L, S->size == n, 2, "arguments are not conformable");
    for (i = 0, ev = S->data; i < n; i++, ev += S->stride)
      luaL_argcheck(L, *ev > 0, 2, "variance is not positive");
  }
  else
    luaL_argcheck(L, S->ndims == 2 && S->dim[0] == n && S->dim[1] == n, 2,
        "arguments are not conformable");
  /* setup destination */
  lua_settop(L, 3);
  if (lua_isnil(L, 3))
    u = nl_pushmatrix(L, 0, 1, &n, 1, n,
        lua_newuserdata(L, n * sizeof(lua_Number)));
  else {
    u = nl_checkmatrix(L, 3);
    checkrvector(L, u, 3);
    luaL_argcheck(L, u->size == n, 3, "arguments are not conformable");
  }
  /* sample */
  if (S->ndims == 1) {
    em = m->data; ev = S->data; eu = u->data;
    for (i = 0; i < n; i++) {
      *eu = gennor(r, *em, *ev);
      em += m->stride; ev += S->stride; eu += u->stride;
    }
  }
  else {
    char uplo = 'L', trans = 'N', diag = 'N';
    lua_Number one = 1.0;
    /* u ~ N(0, I_n) */
    eu = u->data;
    for (i = 0; i < n; i++, eu += u->stride)
      *eu = gennor(r, 0, 1);
    /* u = S * u */
    if (S->stride != 1 /* non-unitary stride? */
        || (S->section != NULL /* non-block section? */
          && (S->section[0].step != 1 || S->section[1].step != 1))) {
      nl_Buffer *buf = nl_getbuffer(L, n * n);
      /* copy S to buffer */
      for (i = 0; i < S->size; i++)
        buf->data.bnum[i] = S->data[nl_mshift(S, i)];
      DTRMV(&uplo, &trans, &diag, &n, buf->data.bnum, &n,
          u->data, &u->stride, 1, 1, 1);
      nl_freebuffer(buf);
    }
    else {
      int ld = S->section ? S->section[0].ld : S->dim[0];
      DTRMV(&uplo, &trans, &diag, &n, S->data, &ld,
          u->data, &u->stride, 1, 1, 1);
    }
    /* u = u + m */
    DAXPY(&n, &one, m->data, &m->stride, u->data, &u->stride);
  }
  return 1;
}


static int runif_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number low = luaL_optnumber(L, 1, 0);
  lua_Number high = luaL_optnumber(L, 2, 1);
  if (low > high)
    luaL_error(L, "inconsistent parameters: %f > %f", low, high);
  if (low == 0 && high == 1)
    setdeviate(number, ranf(r), 3);
  else
    setdeviate(number, genunf(r, low, high), 3);
  return 1;
}


static int runifx_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number low = luaL_optnumber(L, 1, 0);
  lua_Number high = luaL_optnumber(L, 2, 1);
  lua_Number range = high - low;
  if (range < 0)
    luaL_error(L, "inconsistent parameters: %f > %f", low, high);
  if (low == 0 && high == 1)
    setdeviate(number, genrand_res53(r), 3); /* note: [0, 1) */
  else
    setdeviate(number, low + range * genrand_res53(r), 3);
  return 1;
}

static int rdirichlet_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  nl_Matrix *v, *alpha = nl_checkmatrix(L, 1);
  lua_Number *ea, *ev, s;
  int i;
  checkrvector(L, alpha, 1);
  for (i = 0, ea = alpha->data; i < alpha->size; i++, ea += alpha->stride)
    luaL_argcheck(L, *ea > 0, 1, "nonpositive entry");
  lua_settop(L, 2);
  if (lua_isnil(L, 2))
    v = nl_pushmatrix(L, 0, 1, alpha->dim, 1, alpha->size,
        lua_newuserdata(L, alpha->size * sizeof(lua_Number)));
  else {
    v = nl_checkmatrix(L, 2);
    checkrvector(L, v, 2);
    luaL_argcheck(L, alpha->size == v->size, 2, "vector sizes differ");
  }
  /* sample gammas */
  ea = alpha->data;
  ev = v->data;
  s = 0;
  for (i = 0; i < v->size; i++) {
    s += *ev = gengam(r, *ea, 1);
    ev += v->stride;
    ea += alpha->stride;
  }
  /* normalize */
  for (i = 0, ev = v->data; i < v->size; i++, ev += v->stride)
    *ev /= s;
  return 1;
}


static int rbinom_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  int n = luaL_checkinteger(L, 1);
  lua_Number p = luaL_checknumber(L, 2);
  checkneg(L, n);
  if (p <= 0 || p >= 1)
    luaL_error(L, "parameter is out of range: %f", p);
  setdeviate(integer, ignbin(r, n, p), 3);
  return 1;
}

static int rnbinom_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  int n = luaL_checkinteger(L, 1);
  lua_Number p = luaL_checknumber(L, 2);
  checkneg(L, n);
  if (p <= 0 || p >= 1)
    luaL_error(L, "parameter is out of range: %f", p);
  setdeviate(integer, ignnbn(r, n, p), 3);
  return 1;
}


static int rpois_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  lua_Number mu = luaL_checknumber(L, 1);
  checkneg(L, mu);
  setdeviate(integer, ignpoi(r, mu), 2);
  return 1;
}


static int runifint_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  int low = luaL_optinteger(L, 1, 0);
  int high = luaL_optinteger(L, 2, 0x7ffffffeUL);
  if (low > high)
    luaL_error(L, "inconsistent parameters: %d > %d", low, high);
  setdeviate(integer, ignuin(r, low, high), 3);
  return 1;
}


static int sample_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  nl_Matrix *w = nl_checkmatrix(L, 1);
  int i, normalized = lua_toboolean(L, 2);
  lua_Number *e, c, u, z = 0;
  checkrvector(L, w, 1);
  /* compute normalizing constant */
  if (!normalized) {
    e = w->data; z = *e;
    for (i = 1; i < w->size; i++) {
      e += w->stride;
      z += *e;
    }
  }
  /* sample */
  c = 0; u = ranf(r); e = w->data;
  for (i = 0; i < w->size && u >= c; i++) {
    c += *e / z;
    e += w->stride;
  }
  lua_pushinteger(L, i);
  return 1;
}


static int lsample_rng (lua_State *L) {
  nl_RNG *r = getrng(L);
  nl_Matrix *w = nl_checkmatrix(L, 1);
  int i, normalized = lua_toboolean(L, 2);
  lua_Number *e, c, u, z = 0;
  checkrvector(L, w, 1);
  /* compute normalizing constant */
  if (!normalized) {
    e = w->data; z = *e;
    for (i = 1; i < w->size; i++) {
      e += w->stride;
      z = nl_lse(z, *e);
    }
  }
  /* sample */
  c = 0; u = ranf(r); e = w->data;
  for (i = 0; i < w->size && u >= c; i++) {
    if (!isinf(*e)) c += exp(*e - z);
    e += w->stride;
  }
  lua_pushinteger(L, i);
  return 1;
}



/* {=====================================================================
 *    Interface
 * ======================================================================} */

static const luaL_Reg rng_lib[] = {
  /* core */
  {"new", new_rng},
  {"copy", copy_rng},
  {"seed", seed_rng},
  /* deviates */
  {"rbeta", rbeta_rng},
  {"rchisq", rchisq_rng},
  {"rexp", rexp_rng},
  {"rf", rf_rng},
  {"rgamma", rgamma_rng},
  {"rnorm", rnorm_rng},
  {"rmvnorm", rmvnorm_rng},
  {"runif", runif_rng},
  {"runifx", runifx_rng},
  {"rdirichlet", rdirichlet_rng},
  {"rbinom", rbinom_rng},
  {"rnbinom", rnbinom_rng},
  {"rpois", rpois_rng},
  {"runifint", runifint_rng},
  /* sample */
  {"sample", sample_rng},
  {"lsample", lsample_rng},
  {NULL, NULL}
};

static nl_RNG *newrng (lua_State *L) {
  /* create rng object */
  nl_RNG *r = (nl_RNG *) lua_newuserdata(L, sizeof(nl_RNG));
  nl_getmetatable(L, RNG_MT);
  lua_setmetatable(L, -2);
  /* set user value */
  luaL_newlibtable(L, rng_lib);
  lua_pushvalue(L, -2); /* rng object */
  nl_register(L, rng_lib, 1);
  lua_setuservalue(L, -2);
  return r;
}

NUMLUA_API int luaopen_numlua_rng (lua_State *L) {
  luaL_newlibtable(L, rng_lib); /* "class" table */
  init_genrand((nl_RNG *) lua_newuserdata(L, sizeof(nl_RNG)), RNG_SEED);
  /* metatable */
  lua_newtable(L);
  lua_pushcfunction(L, rng__tostring);
  lua_setfield(L, -2, "__tostring");
  lua_pushcfunction(L, rng__len);
  lua_setfield(L, -2, "__len");
  lua_pushliteral(L, RNG_LIBNAME);
  lua_setfield(L, -2, "__type");
  lua_pushlightuserdata(L, RNG_MT);
  lua_pushvalue(L, -2); /* MT */
  lua_rawset(L, LUA_REGISTRYINDEX);
  /* register in lib rng object */
  lua_setmetatable(L, -2);
  nl_register(L, rng_lib, 1);
  return 1;
}

