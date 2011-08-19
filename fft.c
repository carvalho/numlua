/* {=================================================================
 *
 * fft.c
 * Fast Fourier transforms from FFTW for NumericLua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lua.h>
#include <lauxlib.h>
#include "numlua.h"

#define PLAN_LIBNAME "plan"

typedef struct {
  fftw_plan plan;
  nl_Matrix *cache;
  int inverse; /* inverse transform? */
  lua_Number scale; /* for inversion */
} nl_Plan;


/* Guidelines: do not cross domains: real-real, complex-complex; operations
 * are in-place */

/* {=====================================================================
 *    Plans
 * ======================================================================} */

fftw_plan nl_createplan (lua_State *L,
    nl_Matrix *m, int inverse, unsigned flags, lua_Number *scale) {
  fftw_plan plan;
  int i;
  nl_Buffer *dim = nl_getbuffer(L, m->ndims);
  for (i = 0; i < m->ndims; i++) /* reverse dims */
    dim->data.bint[i] = m->dim[m->ndims - 1 - i];
  *scale = 1.0 / m->size;
  if (m->iscomplex) { /* fft plan? */
    /* in-place, howmany == 1, dist ignored, nembed == n */
    plan = fftw_plan_many_dft(m->ndims, (const int *) dim->data.bint, 1,
        (fftw_complex *) m->data, NULL, m->stride, 0,
        (fftw_complex *) m->data, NULL, m->stride, 0,
        inverse ? FFTW_BACKWARD : FFTW_FORWARD, flags);
  }
  else { /* fct plan? */
    nl_Buffer *kind = nl_getbuffer(L, m->ndims);
    if (inverse) {
      for (i = 0; i < m->ndims; i++) {
        kind->data.bint[i] = FFTW_REDFT01;
        *scale *= 0.5;
      }
    }
    else {
      for (i = 0; i < m->ndims; i++)
        kind->data.bint[i] = FFTW_REDFT10;
    }
    /* in-place, howmany == 1, dist ignored, nembed == n */
    plan = fftw_plan_many_r2r(m->ndims, (const int *) dim->data.bint, 1,
        m->data, NULL, m->stride, 0,
        m->data, NULL, m->stride, 0,
        (const fftw_r2r_kind *) kind->data.bint, flags);
    nl_freebuffer(kind);
  }
  nl_freebuffer(dim);
  return plan;
}

static int fft_plan (lua_State *L) {
  nl_Matrix *m = nl_checkmatrix(L, 1);
  int inverse = lua_toboolean(L, 2);
  unsigned flags = (unsigned) luaL_optinteger(L, 3, FFTW_MEASURE);
  nl_Plan *p = (nl_Plan *) lua_newuserdata(L, sizeof(nl_Plan));
  p->cache = m;
  p->inverse = inverse;
  p->plan = nl_createplan(L, m, inverse, flags, &p->scale);
  if (p->plan == NULL) {
    lua_pushnil(L);
    lua_pushliteral(L, "cannot create plan");
    return 2;
  }
  lua_pushvalue(L, lua_upvalueindex(1)); /* metatable */
  lua_pushlightuserdata(L, (void *) p);
  lua_pushvalue(L, 1); /* cache (matrix) */
  lua_rawset(L, -3); /* mt[light(p)] = cache */
  lua_setmetatable(L, -2); /* mt(p) = upvalue */
  return 1;
}


/* get cache */
static int fft_plan__len (lua_State *L) {
  lua_pushlightuserdata(L, lua_touserdata(L, 1));
  lua_rawget(L, lua_upvalueindex(1)); /* mt[light(plan)] = cache */
  return 1;
}

/* execute plan */
static int fft_plan__call (lua_State *L) {
  nl_Plan *p = (nl_Plan *) lua_touserdata(L, 1);
  fftw_execute(p->plan);
  if (p->inverse) {
    if (p->cache->iscomplex) /* fft? */
      ZDSCAL(&p->cache->size, &p->scale, (nl_Complex *) p->cache->data,
          &p->cache->stride);
    else
      DSCAL(&p->cache->size, &p->scale, p->cache->data, &p->cache->stride);
  }
  return 0;
}

static int fft_plan__gc (lua_State *L) {
  nl_Plan *p = (nl_Plan *) lua_touserdata(L, 1);
  fftw_destroy_plan(p->plan);
  lua_pushvalue(L, lua_upvalueindex(1));
  lua_pushlightuserdata(L, (void *) p);
  lua_pushnil(L);
  lua_rawset(L, -3); /* mt[light(plan)] = nil, free ref to cache */
  return 0;
}

static int fft_plan__tostring (lua_State *L) {
  lua_pushfstring(L, PLAN_LIBNAME ": %p", lua_touserdata(L, 1));
  return 1;
}


/* {=====================================================================
 *    Wisdom
 * ======================================================================} */

static void write_char (char c, void *v) {
  luaL_Buffer *b = (luaL_Buffer *) v;
  luaL_addchar(b, c);
}

static int fft_wisdom (lua_State *L) {
  if (lua_isnoneornil(L, 1)) { /* return wisdom? */
    luaL_Buffer b;
    luaL_buffinit(L, &b);
    fftw_export_wisdom(write_char, (void *)&b);
    luaL_pushresult(&b);
    return 1;
  }
  if (lua_toboolean(L, 1)) /* forget wisdom? */
    fftw_forget_wisdom();
  else /* get wisdom */
    fftw_import_wisdom_from_string(luaL_checkstring(L, 1));
  return 0;
}

/* {=====================================================================
 *    Interface
 * ======================================================================} */


static const luaL_Reg fft_func[] = {
  {"wisdom", fft_wisdom},
  {NULL, NULL}
};

static const luaL_Reg plan_mt[] = {
  {"__len", fft_plan__len},
  {"__call", fft_plan__call},
  {"__gc", fft_plan__gc},
  {"__tostring", fft_plan__tostring},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_fft (lua_State *L) {
  fftw_import_system_wisdom();
  /* plan */
  lua_newtable(L); /* metatable */
  lua_pushvalue(L, -1);
  nl_register(L, plan_mt, 1);
  lua_pushliteral(L, PLAN_LIBNAME);
  lua_setfield(L, -2, "__type");
  /* fft */
  luaL_newlib(L, fft_func);
  lua_pushvalue(L, -2); /* plan metatable */
  lua_pushcclosure(L, fft_plan, 1);
  lua_setfield(L, -2, "plan"); /* fft.plan */
  lua_newtable(L); /* planner flags */
  lua_pushinteger(L, FFTW_ESTIMATE);
  lua_setfield(L, -2, "estimate");
  lua_pushinteger(L, FFTW_MEASURE);
  lua_setfield(L, -2, "measure");
  lua_pushinteger(L, FFTW_PATIENT);
  lua_setfield(L, -2, "patient");
  lua_pushinteger(L, FFTW_EXHAUSTIVE);
  lua_setfield(L, -2, "exhaustive");
  lua_pushinteger(L, FFTW_WISDOM_ONLY);
  lua_setfield(L, -2, "wisdomonly");
  lua_setfield(L, -2, "flag"); /* fft.flag */
  return 1;
}

