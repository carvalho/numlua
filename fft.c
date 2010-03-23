#include <lua.h>
#include <lauxlib.h>
#include "numlua.h"

static void write_char (char c, void *v) {
  luaL_Buffer *b = (luaL_Buffer *) v;
  luaL_addchar(b, c);
}

static int fft_wisdom (lua_State *L) {
  if (lua_gettop(L) == 0) { /* return wisdom? */
    luaL_Buffer b;
    luaL_buffinit(L, &b);
    fftw_export_wisdom(write_char, (void *)&b);
    luaL_pushresult(&b);
    return 1;
  }
  /* get wisdom */
  fftw_import_wisdom_from_string(luaL_checkstring(L, 1));
  return 0;
}

static const luaL_reg fft_func[] = {
  {"wisdom", fft_wisdom},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_fft (lua_State *L) {
  if (nl_isloaded(L, FFT_LIBNAME)) return 1;
  fftw_import_system_wisdom();
  luaL_register(L, FFT_LIBNAME, fft_func);
  return 1;
}

