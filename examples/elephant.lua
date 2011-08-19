--[[
Implementation for:
"Drawing an elephant with four complex parameters"
Jurgen Mayer, Khaled Khairy, and Jonathon Howard,
Am. J. Phys. 78, 648 (2010), doi:10.1119/1.3254017
--]]

require "numlua"

local zeros, c, linspace = matrix.zeros, matrix.c, matrix.linspace
local sin, cos, set = matrix.sin, matrix.cos, matrix.set
local j, creal, cimag = complex.j, complex.real, complex.imag

-- expand series on `t` with coefficients `A` and `B`
local function fourier0 (t, A, B) -- naive: many temp objects
  local f = zeros(#t)
  for k = 1, #A do
    f = f + A[k] * cos((k - 1) * t) + B[k] * sin((k - 1) * t)
  end
  return f
end

local function fourier (t, A, B)
  local op = numlua.opmode(true) -- in-place
  local f, w = zeros(#t), zeros(#t) -- `w` is workspace
  for k = 1, #A do
    f = f + A[k] * cos((k - 1) * set(w, t))
    f = f + B[k] * sin((k - 1) * set(w, t))
  end
  numlua.opmode(op) -- restore opmode
  return f
end

function elephant (p, nsamples) -- `p` are shape parameters
  local Ax = c(0, 0, 0, creal(p[3]), 0, creal(p[4]))
  local Bx = c(0, creal(p[1]), creal(p[2]), 0, 0, 0)
  local Ay = c(0, cimag(p[4]), 0, 0, 0, 0)
  local By = c(0, cimag(p[1]), cimag(p[2]), cimag(p[3]), 0, 0)
  -- generate curve
  local t = linspace(0, 2 * math.pi, nsamples or 100)
  local x, y = fourier(t, Ax, Bx), fourier(t, Ay, By)
  return y, -x -- reflect and rotate
end

-- parameters from Table I in above reference (no "wiggling"):
local param = c(50 - 30*j, 18 + 8*j, 12 - 10*j, -14 - 60*j)
local x, y = elephant(param) 
for i = 1, #x do print(x[i], y[i]) end

