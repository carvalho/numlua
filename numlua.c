/* {=================================================================
 *
 * numlua.c
 * Numeric Lua (base routines)
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lua.h>
#include <lauxlib.h>
#include "numlua.h"


/* {=====================================================================
 *    API
 * ======================================================================} */

NUMLUA_API int nl_typeerror (lua_State *L, int narg, const char *tname) {
  const char *msg = lua_pushfstring(L, "%s expected, got %s",
      tname, luaL_typename(L, narg));
  return luaL_argerror(L, narg, msg);
}

#if LUA_VERSION_NUM <= 501
/* adapted from Lua 5.2's lauxlib.c */
NUMLUA_API void nl_require (lua_State *L, const char *modname,
    lua_CFunction openf, int glb) {
  lua_pushcfunction(L, openf);
  lua_pushstring(L, modname);
  lua_call(L, 1, 1); /* open module */
  luaL_findtable(L, LUA_REGISTRYINDEX, "_LOADED", 1);
  lua_pushvalue(L, -2);
  lua_setfield(L, -2, modname); /* _LOADED[modname] = module */
  lua_pop(L, 1); /* _LOADED table */
  if (glb) {
    lua_pushvalue(L, -1); /* module */
    lua_setfield(L, LUA_GLOBALSINDEX, modname); /* _G[modname] = module */
  }
}
#endif

/* Numeric buffer */
static int nl_buffer_ = 0;
#define NL_BUFFER ((void *)&nl_buffer_)
static int nl_buffer_mt_ = 0;
#define NL_BUFFER_MT ((void *)&nl_buffer_mt_)

static nl_Buffer *nbuffer_new (lua_State *L, int size) {
  nl_Buffer *nb = (nl_Buffer *) lua_newuserdata(L, sizeof(nl_Buffer)
      + size * sizeof(buf_Number));
  nb->size = size;
  nb->busy = 0; /* free */
  return nb;
}

NUMLUA_API nl_Buffer *nl_getbuffer (lua_State *L, int size) {
  nl_Buffer *nb = NULL;
  int found = 0;
  int i, n;
  lua_pushlightuserdata(L, NL_BUFFER);
  lua_rawget(L, LUA_REGISTRYINDEX); /* buffer table */
  n = (int) lua_rawlen(L, -1);
  /* search best nl_Buffer:
   * either max size or the first one s.t. b->size >= size */
  for (i = 1; i <= n && !found; i++) {
    nl_Buffer *curr_nb;
    lua_rawgeti(L, -1, i);
    curr_nb = (nl_Buffer *) lua_touserdata(L, -1);
    if (!curr_nb->busy) {
      if (curr_nb->size >= size) { /* found */
        nb = curr_nb;
        found = 1;
      }
      else if (!nb || curr_nb->size > nb->size)
        nb = curr_nb;
    }
    lua_pop(L, 1);
  }
  if (!nb || !found) { /* need new buffer? */
    nb = nbuffer_new(L, size);
    lua_rawseti(L, -2, n + 1);
  }
  nb->busy = 1; /* busy */
  lua_pop(L, 1); /* buffer table */
  return nb;
}

NUMLUA_API void nl_getbuffertable (lua_State *L) {
  lua_pushlightuserdata(L, NL_BUFFER);
  lua_rawget(L, LUA_REGISTRYINDEX);
}

NUMLUA_API int nl_releasebuffer (lua_State *L, int thold) {
  int i, l, n = 0;
  lua_pushlightuserdata(L, NL_BUFFER);
  lua_rawget(L, LUA_REGISTRYINDEX); /* buffer table */
  l = lua_rawlen(L, -1);
  lua_createtable(L, l, 0); /* store at most l buffers */
  lua_pushlightuserdata(L, NL_BUFFER_MT);
  lua_rawget(L, LUA_REGISTRYINDEX); /* weak-valued MT */
  lua_setmetatable(L, -2); /* weak-valued MT */
  for (i = 1; i <= l; i++) {
    nl_Buffer *nb;
    lua_rawgeti(L, -2, i);
    nb = (nl_Buffer *) lua_touserdata(L, -1);
    if (nb->busy || nb->size < thold)
      lua_rawseti(L, -2, ++n);
    else lua_pop(L, 1);
  }
  lua_pushlightuserdata(L, NL_BUFFER);
  lua_insert(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX); /* reg[NL_BUFFER] = new buffer table */
  lua_pop(L, 1); /* old buffer table */
  return l - n; /* #released buffers */
}


/* {=====================================================================
 *    Interface
 * ======================================================================} */

static int numlua_buffer (lua_State *L) {
  nl_Buffer *nb;
  static const char *const opts[] = {"release", "status", NULL};
  int opt = luaL_checkoption(L, 1, "status", opts);
  switch (opt) {
    case 0: {
      lua_pushinteger(L, nl_releasebuffer(L, luaL_optinteger(L, 2, 0)));
      return 1;
    }
    case 1: {
      int status = lua_toboolean(L, 2);
      int i, l, n, size;
      nl_getbuffertable(L);
      l = lua_rawlen(L, -1);
      n = size = 0;
      for (i = 1; i <= l; i++) {
        lua_rawgeti(L, -1, i);
        nb = (nl_Buffer *) lua_touserdata(L, -1);
        if (nb->busy == status) {
          n++;
          size += nb->size;
        }
        lua_pop(L, 1);
      }
      lua_pushinteger(L, size);
      lua_pushinteger(L, n);
      return 2;
    }
  }
  return 0; /* avoid compiler warnings */
}


static int numlua_type (lua_State *L) {
  luaL_checkany(L, 1);
  if (lua_type(L, 1) == LUA_TUSERDATA) {
    if (lua_getmetatable(L, 1)) {
      lua_getfield(L, -1, "__type");
      if (!lua_isnil(L, -1)) return 1;
    }
  }
  lua_pushstring(L, luaL_typename(L, 1));
  return 1;
}

int nl_opmode = 0; /* in-place? */
static int numlua_opmode (lua_State *L) {
  if (lua_gettop(L) > 0) { /* set mode? */
    lua_pushboolean(L, nl_opmode);
    nl_opmode = lua_toboolean(L, 1);
  }
  else lua_pushboolean(L, nl_opmode);
  return 1;
}


static const luaL_Reg numlua_func[] = {
  {"buffer", numlua_buffer},
  {"type", numlua_type},
  {"opmode", numlua_opmode},
  {NULL, NULL}
};




NUMLUA_API int luaopen_numlua_base (lua_State *L) {
  /* init buffer table */
  lua_pushlightuserdata(L, NL_BUFFER);
  lua_newtable(L);
  lua_createtable(L, 0, 1);
  lua_pushliteral(L, "v");
  lua_setfield(L, -2, "__mode");
  lua_pushlightuserdata(L, NL_BUFFER_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX); /* reg[NL_BUFFER_MT] = mt */
  lua_setmetatable(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX); /* reg[NL_BUFFER] = buffer table */
  /* register base lib */
  luaL_newlib(L, numlua_func);
  lua_pushliteral(L, NUMLUA_VERSION);
  lua_setfield(L, -2, "_VERSION");
  return 1;
}

static const luaL_Reg numlua_modules[] = {
  {NUMLUA_LIBNAME,  luaopen_numlua_base},
  {COMPLEX_LIBNAME, luaopen_numlua_complex},
  {FFT_LIBNAME, luaopen_numlua_fft},
  {MATRIX_LIBNAME, luaopen_numlua_lmatrix},
  {RNG_LIBNAME, luaopen_numlua_rng},
  {STAT_LIBNAME, luaopen_numlua_stat},
  {MATHX_LIBNAME, luaopen_numlua_mathx},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua (lua_State *L) {
  const luaL_Reg *mod;
  for (mod = numlua_modules; mod->func; mod++) {
    nl_require(L, mod->name, mod->func, 1);
    lua_pop(L, 1);
  }
  lua_pushglobaltable(L);
  lua_getfield(L, -1, "require");
  lua_pushliteral(L, "numlua.matrix");
  lua_call(L, 1, 0);
  return 0;
}

