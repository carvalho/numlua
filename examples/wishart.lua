-- Sample from a Wishart distribution
require "numlua"

local zeros, chol = matrix.zeros, matrix.chol
local add, trmul, hemul = matrix.add, matrix.trmul, matrix.hemul
local rnorm = rng.rnorm

function mvnorm (mu, S)
  local L = assert(chol(S, "L")) -- Cholesky factor: S = L * L^T
  return function (dest)
    local dest = dest or zeros(#mu)
    local s = rnorm(0, 1, dest) -- s ~ N(0, I_n)
    s = trmul(s, L) -- s = L * s, s ~ N(0, S)
    s = add(s, mu, true) -- s = s + mu, s ~ N(mu, S)
    return s
  end
end

function wishart (S, n)
  local m = #S
  local L = assert(chol(S, "L")) -- Cholesky factor: S = L * L^T
  local c = zeros(m) -- cache
  return function ()
    local w = zeros(m, m)
    for i = 1, n - 1 do
      c = trmul(rnorm(0, 1, c), L) -- c ~ N(0, S)
      w = hemul(w, c, false, "L") -- w = w + c * c^T
    end
    return hemul(w, trmul(rnorm(0, 1, c), L)) -- n-th run: full hemul
  end
end

