#include <lua.h>
#include <lauxlib.h>
#include "numlua.h"

/* =======   API   ======= */

NUMLUA_API int nl_isloaded (lua_State *L, const char *name) {
  lua_getfield(L, LUA_REGISTRYINDEX, "_LOADED");
  lua_getfield(L, -1, name);
  lua_replace(L, -2);
  return (lua_type(L, -1) == LUA_TTABLE);
}

/* returns 1 if lname is already loaded */
NUMLUA_API int nl_requiredeps (lua_State *L, const char *lname,
    const luaL_Reg *dep) {
  if (lname != NULL && nl_isloaded(L, lname)) return 1;
  if (dep == NULL) return 0;
  for (; dep->name; dep++) {
    lua_pushcfunction(L, dep->func);
    lua_pushstring(L, dep->name);
    if (lua_pcall(L, 1, 0, 0))
      luaL_error(L, "error loading module `%s': %s",
          dep->name, lua_tostring(L, -1));
  }
  return 0;
}

/* assumes there is table and upvalues on top of stack */
NUMLUA_API void nl_register (lua_State *L, const luaL_Reg *l, int nup) {
  if (l == NULL) return;  /* nothing to register? */
  for (; l->name; l++) {  /* else fill the table with given functions */
    int i;
    for (i=0; i<nup; i++) /* copy upvalues to the top */
      lua_pushvalue(L, -nup);
    lua_pushcclosure(L, l->func, nup);
    lua_setfield(L, -(nup+2), l->name);
  }
  lua_pop(L, nup); /* remove upvalues */
}

/* Numeric buffer */
static int nbuffer_ = 0;
#define NBUFFER ((void *)&nbuffer_)
static int nbuffer_mt_ = 0;
#define NBUFFER_MT ((void *)&nbuffer_mt_)

static nBuffer *nbuffer_new (lua_State *L, int size) {
  nBuffer *nb = (nBuffer *) lua_newuserdata(L, sizeof(nBuffer)
      + size * sizeof(buf_Number));
  nb->size = size;
  nb->busy = 0; /* free */
  return nb;
}

NUMLUA_API nBuffer *nl_getbuffer (lua_State *L, int size) {
  nBuffer *nb = NULL;
  int found = 0;
  int i, n;
  lua_pushlightuserdata(L, NBUFFER);
  lua_rawget(L, LUA_REGISTRYINDEX); /* buffer table */
  n = (int) lua_objlen(L, -1);
  /* search best nBuffer:
   * either max size or the first one s.t. b->size >= size */
  for (i = 1; i <= n && !found; i++) {
    nBuffer *curr_nb;
    lua_rawgeti(L, -1, i);
    curr_nb = (nBuffer *) lua_touserdata(L, -1);
    if (!curr_nb->busy) {
      if (curr_nb->size >= size) { /* found */
        nb = curr_nb;
        found = 1;
      }
      else if (!nb || curr_nb->size > nb->size) nb = curr_nb;
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
  lua_pushlightuserdata(L, NBUFFER);
  lua_rawget(L, LUA_REGISTRYINDEX);
}

NUMLUA_API int nl_releasebuffer (lua_State *L, int thold) {
  int i, l, n = 0;
  lua_pushlightuserdata(L, NBUFFER);
  lua_rawget(L, LUA_REGISTRYINDEX); /* buffer table */
  l = luaL_getn(L, -1);
  lua_createtable(L, l, 0); /* store at most l buffers */
  lua_pushlightuserdata(L, NBUFFER_MT);
  lua_rawget(L, LUA_REGISTRYINDEX); /* weak-valued MT */
  lua_setmetatable(L, -2); /* weak-valued MT */
  for (i = 1; i <= l; i++) {
    nBuffer *nb;
    lua_rawgeti(L, -2, i);
    nb = (nBuffer *) lua_touserdata(L, -1);
    if (nb->busy || nb->size < thold)
      lua_rawseti(L, -2, ++n);
    else lua_pop(L, 1);
  }
  lua_pushlightuserdata(L, NBUFFER);
  lua_insert(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX); /* reg[NBUFFER] = new buffer table */
  lua_pop(L, 1); /* old buffer table */
  return l - n; /* #released buffers */
}


/* =======   Interface   ======= */

/*
 * numlua.buffer(opt [, arg])
 *   o "release": release all buffers with size at least @arg. @arg defaults
 *      to 0.
 *   o "status": returns the total size and number of buffers that are "busy"
 *      if @arg is true or "free" otherwise.
*/
static int numlua_buffer (lua_State *L) {
  nBuffer *nb;
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
      l = luaL_getn(L, -1);
      n = size = 0;
      for (i = 1; i <= l; i++) {
        lua_rawgeti(L, -1, i);
        nb = (nBuffer *) lua_touserdata(L, -1);
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


/* numlua.type(obj) */
static int numlua_type (lua_State *L) {
  luaL_checkany(L, 1);
  if (lua_type(L, 1) == LUA_TUSERDATA) {
    if (lua_getmetatable(L, 1)) {
      lua_getfield(L, -1, "__type");
      if (lua_isstring(L, -1)) return 1;
    }
  }
  lua_pushstring(L, luaL_typename(L, 1));
  return 1;
}

static const luaL_Reg numlua_func[] = {
  {"buffer", numlua_buffer},
  {"type", numlua_type},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_base (lua_State *L) {
  if (nl_isloaded(L, NUMLUA_LIBNAME)) return 1;
  /* init buffer table */
  lua_pushlightuserdata(L, NBUFFER);
  lua_newtable(L);
  lua_createtable(L, 0, 1);
  lua_pushliteral(L, "v");
  lua_setfield(L, -2, "__mode");
  lua_pushlightuserdata(L, NBUFFER_MT);
  lua_pushvalue(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX); /* reg[NBUFFER_MT] = mt */
  lua_setmetatable(L, -2);
  lua_rawset(L, LUA_REGISTRYINDEX); /* reg[NBUFFER] = buffer table */
  /* register base lib */
  luaL_register(L, NUMLUA_LIBNAME, numlua_func);
  lua_pushliteral(L, NUMLUA_VERSION);
  lua_setfield(L, -2, "_VERSION");
  return 1;
}

static const luaL_Reg numlua_modules[] = {
  {NUMLUA_LIBNAME,  luaopen_numlua_base},
  {COMPLEX_LIBNAME, luaopen_numlua_complex},
  {MATRIX_LIBNAME, luaopen_numlua_lmatrix},
  {FFT_LIBNAME, luaopen_numlua_fft},
  {RNG_LIBNAME, luaopen_numlua_rng},
  {STAT_LIBNAME, luaopen_numlua_stat},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua (lua_State *L) {
  nl_requiredeps(L, NULL, numlua_modules); /* call submodules */
  return 0;
}

