Numeric Lua
===========

_Numeric Lua_ is a numerical package for the Lua programming language. It
includes support for complex numbers, multidimensional matrices, random number
generation, fast Fourier transforms, and special functions. Most of the
routines are simple wrappers for well known numerical libraries: complex
numbers and part of the extended math modules come from C99; other special
functions, including statistical functions, are adapted from Netlib's SLATEC
and DCDFLIB; random number generation is based on Takuji Nishimura and Makoto
Matsumoto's Mersenne Twister generator as the "engine" (uniform deviates) and
Netlib's RANLIB for the remaining deviates; fast Fourier transforms are
implemented from FFTW; and the matrix package draws most of its numeric
intensive routines from Netlib's ubiquitous BLAS and LAPACK packages.

Numeric Lua tries to maintain Lua's minimalist approach by providing bare-bone
wrappers to the numerical routines. The user can use the outputs for further
computations and is then fully responsible for the results. Other Lua features
are also available, such as OO simulation through metamethods and functional
facilities. A basic API is provided in order to promote extensibility. Also,
check `numlua.seeall` for a quick way to start using Numeric Lua.

Numeric Lua is licensed under the same license as Lua -- the MIT license --
and so can be freely used for academic and commercial purposes.


Documentation
-------------

The documentation for Numeric Lua can be found in the distribution in two
formats: _html_ in `docs` and "Lua help pages" (_lhp_). Lua help pages are
viewed in the interpreter and are thus more suitable for quick on-line
reference checks --- similar to man pages --- but they require `luahelp`.


Installation
------------

Numeric Lua depends on
[BLAS/LAPACK](http://www.netlib.org/lapack "BLAS/LAPACK"),
[FFTW](http://www.fftw.org "FFTW"), and
[HDF5](http://www.hdfgroup.org/HDF5 "HDF5"). These external dependencies are
more easily and conveniently handled by package managers in Linux (`apt` and
`yum`, for example) or Mac OS X (`fink`, `macports`); in Windows there is more
work to be done, especially if tailored, optimized versions of these libraries
are desired. In any case, Numeric Lua can be installed or built with
[`luarocks`](http://luarocks.org "luarocks").

Building Numeric Lua in `luarocks` should be straightforward and mostly
requires specifying include and lib dirs for the dependencies. In
Debian/Ubuntu systems, for example, you would typically just need to issue

      $ sudo apt-get install libblas-dev liblapack-dev libfftw3-dev libhdf5-serial-dev
      $ luarocks make numlua-0.3-1.rockspec

since the dirs happen to be system defaults. Similarly, on Mac OS X using
`fink`,

        $ fink install atlas fftw3 hdf5
        $ luarocks HDF5_DIR=/sw FFTW3_DIR=/sw LAPACK=/sw make numlua-0.3-1.rockspec

Note that ATLAS optimized libraries are being used here, instead of the
vanilla BLAS/LAPACK (a similar package exists in `apt`.) In this case, you
need to change the rockspec to match your libraries.


### Building on Windows ###

Building on Windows is trickier, so we provide more details. We will be using
[MinGW](http://www.mingw.org "Minimalist GNU for Windows") in the following
steps:

  1. Install MinGW and utils with `mingw-get install mingw-utils`

  2. Install the `luarocks` binaries for Windows:
    * Download and install with: `install.bat /MW` (use MinGW)
    * Change `add_flags(extras, "%s.lib", libraries)` to
      `add_flags(extras, "-l%s", libraries)` in `compile_library` at
      `c:\LuaRocks\2.0\lua\luarocks\build\builtin.lua`

  3. Download BLAS/LAPACK (note: _unoptimized_!) from Netlib and build:
    * Copy `make.inc.example` to `make.inc`, set `PLAT=`, and remove the `-g`
      flags and add `-O2` in `OPTS`
    * Now build: `cd install && make`, `cd blas\src && make`, `cd src && make`
      (for LAPACK)
    * Rename `blas.a` and `lapack.a` to `libblas.a` and `liblapack.a`

  4. Download `fftw3` and compile it (check instructions in FFTW
     documentation.)

  5. Download the binary distribution of HDF5.
    * Patch `include/H5public.h` by adding:

          #undef H5_SIZEOF_SSIZE_T
          #define H5_SIZEOF_SSIZE_T H5_SIZEOF_LONG

    * Wrap the dll for MinGW:

          pexports hdf5dll.dll > hdf5dll.def
          dlltool -d hdf5dll.def -l libhdf5.a

    * Put `hdf5dll.dll`, `zlib1.dll`, and `szip.dll` in your `PATH`

  6. Finally, run `luarocks make`:

        luarocks HDF5_INCDIR=\path\to\hdf5\include HDF5_LIBDIR=\path\tohdf5\dll
        FFTW3_INCDIR=\path\to\fftw3\include FFTW3_LIBDIR=\path\to\fftw3\lib
        LAPACK_LIBDIR=\path\to\lapack\lib make numlua-0.3-1.rockspec

      For example, `\path\to\fftw3` is usually `c:\MinGW\msys\1.0\local`.

