/* {=================================================================
 *
 * numlua.h
 * Numeric Lua
 * Luis Carvalho (lecarval @ math.bu.edu)
 * See Copyright Notice at the bottom of this file
 * $Id: luarng.h,v 1.4 2006/02/18 02:06:52 carvalho Exp $
 *
 * ==================================================================} */

#ifndef numlua_h
#define numlua_h

#include <lua.h>
#include <lauxlib.h>
#include <math.h>
#include <fftw3.h>

#define NUMLUA_VERSION "NumericLua 0.3"
#define NUMLUA_API extern

#define NUMLUA_LIBNAME "numlua"
#define COMPLEX_LIBNAME "complex"
#define MATRIX_LIBNAME "matrix"
#define FFT_LIBNAME "fft"
#define RNG_LIBNAME "rng"
#define STAT_LIBNAME "stat"

#define toudata(name,libname,type) \
  static type *to ## name (lua_State *L, int narg) { \
    type *p = NULL; \
    if (lua_getmetatable(L, narg)) {  /* does it have a metatable? */ \
      if (lua_rawequal(L, -1, LUA_ENVIRONINDEX)) /* MT == env? */ \
        p = (type *) lua_touserdata(L, narg); \
      lua_pop(L, 1); /* remove metatable */ \
    } \
    return p; \
  } \
  \
  static type *check ## name (lua_State *L, int narg) { \
    type *p = to ## name(L, narg); \
    if (p == NULL) luaL_typerror(L, narg, libname); \
    return p; \
  }

#define libtoudata(name,libname,mt,type) \
  NUMLUA_API type *nl_to ## name (lua_State *L, int narg) { \
    type *p = NULL; \
    if (lua_getmetatable(L, narg)) { /* does it have a metatable? */ \
      lua_pushlightuserdata(L, mt); \
      lua_rawget(L, LUA_REGISTRYINDEX); \
      if (lua_rawequal(L, -1, -2)) /* right MT? */ \
        p = lua_touserdata(L, narg); \
      lua_pop(L, 2); /* MTs */ \
    } \
    return p; \
  } \
  \
  NUMLUA_API type *nl_check ## name (lua_State *L, int narg) { \
    type *p = nl_to ## name(L, narg); \
    if (p == NULL) luaL_typerror(L, narg, libname); \
    return p; \
  }

NUMLUA_API int nl_isloaded (lua_State *L, const char *name);
NUMLUA_API int nl_requiredeps (lua_State *L, const char *lname,
    const luaL_Reg *dep);
NUMLUA_API void nl_register (lua_State *L, const luaL_Reg *l, int nup);

/* nBuffer */
typedef union {
  int bint[1];
  lua_Number bnum[1];
} buf_Number;

typedef struct {
  int size;
  int busy; /* status */
  buf_Number data;
} nBuffer;

#define nl_freebuffer(nb) ((nb)->busy = 0)
NUMLUA_API nBuffer *nl_getbuffer (lua_State *L, int size);
NUMLUA_API void nl_getbuffertable (lua_State *L);
NUMLUA_API int nl_releasebuffer (lua_State *L, int threshold);
NUMLUA_API int luaopen_numlua_base (lua_State *L);
NUMLUA_API int luaopen_numlua (lua_State *L);

/* complex */
typedef struct {
  lua_Number re;
  lua_Number im;
} nl_Complex;

NUMLUA_API lua_Number nl_abscomplex (nl_Complex c);
NUMLUA_API void nl_expcomplex (nl_Complex *c, nl_Complex a);
NUMLUA_API void nl_logcomplex (nl_Complex *c, nl_Complex a);
NUMLUA_API void nl_sincomplex (nl_Complex *c, nl_Complex a);
NUMLUA_API void nl_coscomplex (nl_Complex *c, nl_Complex a);
NUMLUA_API void nl_sqrtcomplex (nl_Complex *c, nl_Complex a);
NUMLUA_API void nl_mulcomplex (nl_Complex *c, nl_Complex a, nl_Complex b);
NUMLUA_API void nl_divcomplex (nl_Complex *c, nl_Complex a, nl_Complex b);
NUMLUA_API void nl_ipowcomplex (nl_Complex *f, nl_Complex c, int e);
NUMLUA_API void nl_powcomplex (nl_Complex *f, nl_Complex a, nl_Complex b);

NUMLUA_API nl_Complex nl_tocomplex (lua_State *L, int narg,
    int *iscomplex);
NUMLUA_API nl_Complex nl_checkcomplex (lua_State *L, int narg);
NUMLUA_API nl_Complex nl_optcomplex (lua_State *L, int narg,
    nl_Complex def);
NUMLUA_API nl_Complex *nl_newcomplex (lua_State *L);
NUMLUA_API nl_Complex *nl_pushcomplex (lua_State *L, nl_Complex c);
NUMLUA_API int luaopen_numlua_complex (lua_State *L);

/* matrix */
typedef struct {
  int iscomplex;
  int ndims;
  int stride;
  int size;
  lua_Number *data;
  int dim[1];
} nl_Matrix;

/* TODO: alignment */

NUMLUA_API nl_Matrix *nl_tomatrix (lua_State *L, int narg);
NUMLUA_API nl_Matrix *nl_checkmatrix (lua_State *L, int narg);
NUMLUA_API nl_Matrix *nl_pushmatrix (lua_State *L, int iscomplex,
    int ndims, int *dim, int stride, int size, lua_Number *data);
int luaopen_numlua_lmatrix (lua_State *L);

/* rng */

/* Struct to hold state vector and current position */
#define RNG_MAXSTATES (624)  /* _N_ in mt19937ar.c */
typedef struct {
  unsigned long v[RNG_MAXSTATES];  /* the array for the state vector */
  int i;  /* i==N+1 means v[N] is not initialized */
} nl_RNG;

NUMLUA_API nl_RNG *nl_checkrng (lua_State *L, int pos);
NUMLUA_API nl_RNG *nl_torng (lua_State *L, int pos);
NUMLUA_API int luaopen_numlua_rng (lua_State *L);

/* TODO: fft, plot, hdf5, stat */

NUMLUA_API int luaopen_numlua_fft (lua_State *L);

NUMLUA_API int luaopen_numlua_stat (lua_State *L);

/* {=================================================================
*
* Copyright (c) 2005-2010 Luis Carvalho
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation files
* (the "Software"), to deal in the Software without restriction,
* including without limitation the rights to use, copy, modify,
* merge, publish, distribute, sublicense, and/or sell copies of the
* Software, and to permit persons to whom the Software is furnished
* to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
* BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
* ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
* ==================================================================} */

#endif
