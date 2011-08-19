/* {=================================================================
 *
 * numlua.h
 * Numeric Lua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice at the bottom of this file
 *
 * ==================================================================} */

#ifndef numlua_h
#define numlua_h

#define H5_NO_DEPRECATED_SYMBOLS /* use HDF5 1.8 */

#include <lua.h>
#include <lauxlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <hdf5.h>
#include "blas.h"
#include "lapack.h"

#define NUMLUA_VERSION "NumericLua 0.3"
#define NUMLUA_API extern

#if LUA_VERSION_NUM <= 501
#define lua_rawlen lua_objlen
#define nl_register(L,l,n) luaL_openlib(L,NULL,l,n)
#define luaL_newlibtable(L,l) \
  lua_createtable(L, 0, sizeof(l)/sizeof((l)[0]) - 1)
#define luaL_newlib(L,l) (luaL_newlibtable(L,l), nl_register(L,l,0))
#define lua_setuservalue lua_setfenv
#define lua_getuservalue lua_getfenv
NUMLUA_API void nl_require (lua_State *L, const char *modname,
    lua_CFunction openf, int glb);
#define lua_pushglobaltable(L) lua_pushvalue(L, LUA_GLOBALSINDEX)
#else
#define lua_number2int(i,n) ((i)=(int)(n))
#define nl_register luaL_setfuncs
#define nl_require luaL_requiref
#endif
NUMLUA_API int nl_typeerror (lua_State *L, int narg, const char *tname);

NUMLUA_API int nl_opmode;
#define nl_inplace(L,n) \
  (lua_isnoneornil(L,(n)) ? nl_opmode : lua_toboolean(L,(n)))

#define NUMLUA_LIBNAME "numlua"
#define COMPLEX_LIBNAME "complex"
#define MATRIX_LIBNAME "matrix"
#define FFT_LIBNAME "fft"
#define RNG_LIBNAME "rng"
#define STAT_LIBNAME "stat"
#define MATHX_LIBNAME "mathx"

#define nl_getmetatable(L,mt) (lua_pushlightuserdata(L, (mt)), \
    lua_rawget(L, LUA_REGISTRYINDEX))

typedef struct {
  int ld; /* leading dimension */
  int step; /* per dimension stride */
} nl_Section;

/* nl_Buffer */
typedef union {
  int bint[1];
  lua_Number bnum[1];
} buf_Number;

typedef struct {
  int size;
  int busy; /* status */
  buf_Number data;
} nl_Buffer;

#define nl_freebuffer(nb) ((nb)->busy = 0)
NUMLUA_API nl_Buffer *nl_getbuffer (lua_State *L, int size);
NUMLUA_API void nl_getbuffertable (lua_State *L);
NUMLUA_API int nl_releasebuffer (lua_State *L, int threshold);

NUMLUA_API int luaopen_numlua_base (lua_State *L);
NUMLUA_API int luaopen_numlua (lua_State *L);


/* complex */
typedef double complex nl_Complex;
#define CPX(x) ((nl_Complex *) (x))

NUMLUA_API lua_Number clogabs (nl_Complex c);
#define cunm(c) (-(c))
#define cadd(a,b) ((a) + (b))
#define csub(a,b) ((a) - (b))
#define cmul(a,b) ((a) * (b))
#define cdiv(a,b) ((a) / (b))
#define ceq(a,b) ((a) == (b))
#define clt(a,b) \
  (creal(a)<creal(b) || ((creal(a)==creal(b))&&(cimag(a)<cimag(b))))
#define cgt(a,b) \
  (creal(a)>creal(b) || ((creal(a)==creal(b))&&(cimag(a)>cimag(b))))

NUMLUA_API nl_Complex nl_tocomplex (lua_State *L, int narg, int *iscomplex);
NUMLUA_API nl_Complex nl_checkcomplex (lua_State *L, int narg);
NUMLUA_API nl_Complex nl_optcomplex (lua_State *L, int narg, nl_Complex def);
NUMLUA_API nl_Complex *nl_newcomplex (lua_State *L);
NUMLUA_API nl_Complex *nl_pushcomplex (lua_State *L, nl_Complex c);
NUMLUA_API int luaopen_numlua_complex (lua_State *L);

/* matrix */
typedef struct {
  int iscomplex;
  int ndims;
  int stride; /* internal for indexing */
  int size; /* overall size (product of dims) */
  nl_Section *section;
  lua_Number *data;
  int dim[1];
} nl_Matrix;

NUMLUA_API nl_Matrix *nl_tomatrix (lua_State *L, int narg);
NUMLUA_API nl_Matrix *nl_checkmatrix (lua_State *L, int narg);
NUMLUA_API nl_Matrix *nl_pushmatrix (lua_State *L, int iscomplex,
    int ndims, int *dim, int stride, int size, lua_Number *data);

NUMLUA_API int nl_msshift (nl_Matrix *m, int eo);
#define nl_mshift(m,i) \
  (((m)->section) ? nl_msshift((m),(i)) : ((i) * (m)->stride))

int luaopen_numlua_lmatrix (lua_State *L);

/* rng */
NUMLUA_API int luaopen_numlua_rng (lua_State *L);

/* fft */
NUMLUA_API int luaopen_numlua_fft (lua_State *L);

/* stat */
NUMLUA_API int luaopen_numlua_stat (lua_State *L);

/* C99 mathx */
NUMLUA_API lua_Number nl_lse (lua_Number w1, lua_Number w2);
NUMLUA_API int luaopen_numlua_mathx (lua_State *L);

/* BLAS level 1 */
#define DSWAP dswap_
#define ZSWAP zswap_
#define DCOPY dcopy_
#define ZCOPY zcopy_
#define DSCAL dscal_
#define ZSCAL zscal_
#define ZDSCAL zdscal_
#define DAXPY daxpy_
#define ZAXPY zaxpy_
#define DDOT ddot_
#define ZDOTC zdotc_
#define ZDOTU zdotu_
#define DNRM2 dnrm2_
#define DZNRM2 dznrm2_
#define DASUM dasum_
#define DZASUM dzasum_
#define IDAMAX idamax_
#define IZAMAX izamax_
/* BLAS level 2 */
#define DTRMV dtrmv_
#define DTRSV dtrsv_
#define ZTRMV ztrmv_
#define ZTRSV ztrsv_
#define ZHER zher_
#define DSYR dsyr_
#define DGER dger_
#define ZGERC zgerc_
#define DGEMV dgemv_
#define ZGEMV zgemv_
/* BLAS level 3 */
#define DTRMM dtrmm_
#define DTRSM dtrsm_
#define ZTRMM ztrmm_
#define ZTRSM ztrsm_
#define ZHERK zherk_
#define DSYRK dsyrk_
#define DGEMM dgemm_
#define ZGEMM zgemm_

/* LAPACK */
/* chol */
#define DPOTRF dpotrf_
#define ZPOTRF zpotrf_
/* lu */
#define DGETRF dgetrf_
#define ZGETRF zgetrf_
/* inv */
#define DTRTRI dtrtri_
#define ZTRTRI ztrtri_
#define DPOTRI dpotri_
#define ZPOTRI zpotri_
#define DGETRI dgetri_
#define ZGETRI zgetri_
/* rcond */
#define DTRCON dtrcon_
#define ZTRCON ztrcon_
#define DPOCON dpocon_
#define ZPOCON zpocon_
#define DGECON dgecon_
#define ZGECON zgecon_
/* svd */
#define ZGESVD zgesvd_
#define DGESVD dgesvd_
/* qr */
#define ZGEQRF zgeqrf_
#define DGEQRF dgeqrf_
#define ZGEQP3 zgeqp3_
#define DGEQP3 dgeqp3_
#define ZUNGQR zungqr_
#define DORGQR dorgqr_
/* eig */
#define ZHEEV zheev_
#define DSYEV dsyev_
#define ZGEEV zgeev_
#define DGEEV dgeev_
/* balance */
#define ZGEBAL zgebal_
#define DGEBAL dgebal_
/* ls */
#define ZGELSY zgelsy_
#define DGELSY dgelsy_
#define ZGELSS zgelss_
#define DGELSS dgelss_


/* {=================================================================
*
* Copyright (c) 2005 Luis Carvalho
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
