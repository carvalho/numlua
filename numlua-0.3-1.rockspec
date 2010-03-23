package = "NumLua"
version = "0.3-1"

source = {
  url = "http://...",
  md5 = "asda",
  dir = "numlua-0.3"
}

description = {
  summary = "Numeric routines for Lua",
  detailed = [[
    NumericLua provides numeric routines...
  ]],
  homepage = "http://numlua.luaforge.net",
  license = "MIT"
}

dependencies = {
  "lua >= 5.1"
}
-- TODO: luahelp?

external_dependencies = {
  BLAS = { header = "cblas.h" },
  LAPACK = { header = "clapack.h" },
  FFTW3 = { header = "fftw3.h" },
}

build = {
  type = "builtin",
  modules = {
    numlua = {
      sources = {"numlua.c", "complex.c", "fft.c", "lmatrix.c",
        "mt.c", "ranlib.c", "rng.c", -- rng
        "dcdflib.c", "ipmpar.c", "stat.c", -- stat
      },
      libraries = {"fftw3", "lapack", "cblas", "atlas"},
    }
  }
}

