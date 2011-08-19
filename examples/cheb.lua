-- Based on "A Primer of Scientific Computing in Lua"

require "numlua"

local seq, dot, cos = matrix.seq, matrix.dot, matrix.cos
local zeros, map, set, fft = matrix.zeros, matrix.map, matrix.set, matrix.fft
local real, complex, c = matrix.real, matrix.complex, matrix.c
local setmetatable, opmode, pi = setmetatable, numlua.opmode, math.pi

local function dct (x)
  local n = #x - 1
  local f = zeros(2 * n, true)
  local fr = f:real()
  local g = fr(1, n + 1)
  g._ = x
  fr(n + 2, 2 * n)._ = x(n, 2, -1)
  f = fft(f, false, false, true) -- in-place
  return g
end

function clenshawcurtis (f, n)
  local x = cos(pi * seq(0, n) / n)
  local w = 0 * x; w(1, n + 1, 2)._ = 2 / (1 - seq(0, n, 2) ^ 2)
  local fx = map(x, f) / (2 * n)
  local g = dct(fx) -- real(fft(complex(c(fx, fx(n, 2, -1)))))
  g(2, n)._ = 2 * g(2, n)
  return dot(w, g)
end


function new0 (f, n) -- naive
  local x = cos(pi * seq(0, n) / n) -- Chebyshev points
  local w = 0 * x; w(1, n + 1, 2)._ = 2 / (1 - seq(0, n, 2) ^ 2) -- weights
  return setmetatable({f=f, n=n}, {
    __call = function (_, a, b) -- integral by linear transform
      local a, b = a or -1, b or 1
      local s = (b - a) / 2
      local fx = map(s * x + (b + a) / 2, f) / (2 * n)
      local g = real(fft(complex(c(fx, fx(n, 2, -1))))) -- dct
      g(2, n)._ = 2 * g(2, n)
      return s * dot(w, g(1, n + 1))
    end
  })
end

function new (f, n)
  local x = cos(pi * seq(0, n) / n) -- Chebyshev points
  local w = 0 * x; w(1, n + 1, 2)._ = 2 / (1 - seq(0, n, 2) ^ 2) -- weights
  local d = zeros(2 * n, true) -- dct buffer
  local dr = d:real() -- real ref
  local g, r = dr(1, n + 1), dr(n + 2, 2 * n) -- real halfs
  local g1 = g(2, n) -- inner ref
  -- object
  return setmetatable({f=f, n=n}, {
    __call = function (_, a, b) -- integral by linear transform
      local op = opmode(true) -- in-place
      local a, b = a or -1, b or 1
      local s = (b - a) / 2
      g = map(s * set(g, x) + (b + a) / 2, f) / (2 * n)
      r._ = g(n, 2, -1); d = fft(d) -- dct
      g1 = 2 * g1 -- g(2, n)._ = 2 * g(2, n)
      opmode(op) -- restore opmode
      return s * dot(w, g)
    end
  })
end

