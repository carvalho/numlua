package = "NumLua"
version = "0.3-1"

source = {
  url = "https://github.com/downloads/carvalho/numlua/numlua-0.3-1.tar.gz",
  md5 = "abc29945daadd1499c7596ebbfb0f88a",
  dir = "numlua-0.3"
}

description = {
  summary = "Numerical routines for Lua",
  detailed = [[
    Numeric Lua is a numerical package for the Lua programming language. It
    includes support for complex numbers, multidimensional matrices, random
    number generation, fast Fourier transforms, and special functions. Most of
    the routines are simple wrappers for well known numerical libraries:
    complex numbers and part of the extended math modules come from C99; other
    special functions, including statistical functions, are adapted from
    Netlib's SLATEC and DCDFLIB; random number generation is based on Takuji
    Nishimura and Makoto Matsumoto's Mersenne Twister generator as the
    "engine" (uniform deviates) and Netlib's RANLIB for the remaining
    deviates; fast Fourier transforms are implemented from FFTW; and the
    matrix package draws most of its numeric intensive routines from Netlib's
    ubiquitous BLAS and LAPACK packages.
  ]],
  homepage = "https://github.com/carvalho/numlua",
  license = "MIT"
}

dependencies = {
  "lua >= 5.1"
}

external_dependencies = {
  FFTW3 = { header = "fftw3.h" },
  HDF5 = { header = "hdf5.h" },
}

build = {
  type = "builtin",
  modules = {
    numlua = {
      sources = {"numlua.c", "complex.c", "fft.c",
        "msort.c", "lmatrix.c", -- matrix
        "mt.c", "ranlib.c", "rng.c", -- rng
        "dcdflib.c", "ipmpar.c", "stat.c", -- stat
        "amos.c", "mathx.c", -- C99 math
      },
      incdirs = {"$(FFTW3_INCDIR)", "$(HDF5_INCDIR)"},
      -- assume blas is in same libdir as lapack
      libdirs = {"$(FFTW3_LIBDIR)", "$(HDF5_LIBDIR)", "$(LAPACK_LIBDIR)"},
      --libraries = {"hdf5", "fftw3", "lapack", "f77blas", "atlas"},
      libraries = {"hdf5", "fftw3", "lapack", "blas"}, -- f77 blas
    },
    ["numlua.matrix"] = "matrix.lua",
    ["numlua.seeall"] = "seeall.lua",
  },
  platforms = {
    mingw32 = { -- assumes that libs are linked with "-l%s" format
      modules = {
        numlua = { -- needs gfortran for blas and lapack
          libraries = {"hdf5", "fftw3", "lapack", "blas", "gfortran"}
        }
      }
    }
  }
}

-- lhp files
local baselhp = {
  "lhp/complex.lhp", "lhp/mathx.lhp", "lhp/matrix.lhp", "lhp/rng.lhp",
}

local liblhp = {
  -- factor
  "lhp/factor/design.lhp", "lhp/factor/fold.lhp", "lhp/factor/partition.lhp",
  -- fft
  "lhp/fft/plan.lhp", "lhp/fft/wisdom.lhp",
  -- mathx
  "lhp/mathx/airya.lhp", "lhp/mathx/airyb.lhp", "lhp/mathx/besselh.lhp",
  "lhp/mathx/besseli.lhp", "lhp/mathx/besselj.lhp", "lhp/mathx/besselk.lhp",
  "lhp/mathx/bessely.lhp", "lhp/mathx/beta.lhp", "lhp/mathx/choose.lhp",
  "lhp/mathx/digamma.lhp", "lhp/mathx/feq.lhp", "lhp/mathx/lbeta.lhp",
  "lhp/mathx/lchoose.lhp", "lhp/mathx/log1pe.lhp", "lhp/mathx/lse.lhp",
  -- matrix
  "lhp/matrix/add.lhp", "lhp/matrix/apply.lhp", "lhp/matrix/balance.lhp",
  "lhp/matrix/chol.lhp", "lhp/matrix/c.lhp", "lhp/matrix/col.lhp",
  "lhp/matrix/complex.lhp", "lhp/matrix/concat.lhp", "lhp/matrix/conj.lhp",
  "lhp/matrix/copy.lhp", "lhp/matrix/cross.lhp", "lhp/matrix/diag.lhp",
  "lhp/matrix/div.lhp", "lhp/matrix/dot.lhp", "lhp/matrix/eig.lhp",
  "lhp/matrix/eindex.lhp", "lhp/matrix/entries.lhp", "lhp/matrix/eorder.lhp",
  "lhp/matrix/fct.lhp", "lhp/matrix/fft.lhp", "lhp/matrix/find.lhp",
  "lhp/matrix/fold.lhp", "lhp/matrix/get.lhp", "lhp/matrix/hemul.lhp",
  "lhp/matrix/ifelse.lhp", "lhp/matrix/imag.lhp", "lhp/matrix/inv.lhp",
  "lhp/matrix/iscomplex.lhp", "lhp/matrix/linspace.lhp", "lhp/matrix/load.lhp",
  "lhp/matrix/ls.lhp", "lhp/matrix/lu.lhp", "lhp/matrix/map.lhp",
  "lhp/matrix/max.lhp", "lhp/matrix/min.lhp", "lhp/matrix/mmul.lhp",
  "lhp/matrix/mul.lhp", "lhp/matrix/new.lhp", "lhp/matrix/norm.lhp",
  "lhp/matrix/ones.lhp", "lhp/matrix/pivot.lhp", "lhp/matrix/pow.lhp",
  "lhp/matrix/qr.lhp", "lhp/matrix/rcond.lhp", "lhp/matrix/real.lhp",
  "lhp/matrix/reshape.lhp", "lhp/matrix/save.lhp", "lhp/matrix/section.lhp",
  "lhp/matrix/set.lhp", "lhp/matrix/shape.lhp", "lhp/matrix/size.lhp",
  "lhp/matrix/slice.lhp", "lhp/matrix/sort.lhp", "lhp/matrix/spread.lhp",
  "lhp/matrix/sum.lhp", "lhp/matrix/svd.lhp", "lhp/matrix/swap.lhp",
  "lhp/matrix/transpose.lhp", "lhp/matrix/trmul.lhp", "lhp/matrix/which.lhp",
  "lhp/matrix/zeros.lhp",
  -- numlua
  "lhp/numlua/buffer.lhp", "lhp/numlua/opmode.lhp", "lhp/numlua/type.lhp",
  -- rng
  "lhp/rng/copy.lhp", "lhp/rng/lsample.lhp", "lhp/rng/new.lhp",
  "lhp/rng/rbeta.lhp", "lhp/rng/rbinom.lhp", "lhp/rng/rchisq.lhp",
  "lhp/rng/rdirichlet.lhp", "lhp/rng/rexp.lhp", "lhp/rng/rf.lhp",
  "lhp/rng/rgamma.lhp", "lhp/rng/rmvnorm.lhp", "lhp/rng/rnbinom.lhp",
  "lhp/rng/rnorm.lhp", "lhp/rng/rpois.lhp", "lhp/rng/runifint.lhp",
  "lhp/rng/runif.lhp", "lhp/rng/runifx.lhp", "lhp/rng/sample.lhp",
  "lhp/rng/seed.lhp",
  -- stat
  "lhp/stat/dbeta.lhp", "lhp/stat/dbinom.lhp", "lhp/stat/dchisq.lhp",
  "lhp/stat/dexp.lhp", "lhp/stat/df.lhp", "lhp/stat/dgamma.lhp",
  "lhp/stat/dhyper.lhp", "lhp/stat/dnbinom.lhp", "lhp/stat/dnorm.lhp",
  "lhp/stat/dpois.lhp", "lhp/stat/dt.lhp", "lhp/stat/factor.lhp",
}

local lhp = {}
for i = 1, #baselhp do
  local f = baselhp[i]
  local t = {f:match("^(lhp)/([%w_]+)%.lhp$")}
  lhp[t[1] .. "." .. t[2]] = f
end
for i = 1, #liblhp do
  local f = liblhp[i]
  local t = {f:match("^(lhp)/(%w+)/([%w_]+)%.lhp$")}
  lhp[t[1] .. "." .. t[2] .. "." .. t[3]] = f
end
build.install = {lua = lhp}

-- vim: set syn=lua :
