-- Levenberg-Marquadt

-- solve (A + mu * I) * h = -g for h
function solveh (A, mu, g, h)
  A:diag():add(mu, true) -- A = A + mu * I
  chol(A, true) -- in-place
  h._ = g; h:mul(-1, true) -- h = -g
  trmul(h, A, "l", true) -- h = inv(A) * h
  trmul(h, A, "l", true, "t") -- h = inv(A)' * h
end

-- f : R^n -> R^m
local TAU, EPS1, EPS2, KMAX = 1e-3, 1e-4, 1e-8, 100
function marquadt (n, m, ff, Jf, x0, tau, tolg, tolx, maxeval)
  local tau, eps1, eps2, kmax = TAU, EPS1, EPS2, KMAX
  local x, J, f = new(n), new(m, n), new(m) -- arg, Jacobian, function
  local xnew, C = new(n), new(n, n) -- caches
  local h = new(n) -- step
  -- initialize
  local k, nu = 0, 2
  x:set(x0);  J:apply(Jf(x));  f:apply(ff(x))
  local A = zeros(n, n):hemul(J, true) -- A = J' * J
  local g = zeros(n):mmul(J, f, "t") -- g = J' * f
  local F = dot(f, f) -- objective
  local found = g:norm"inf" <= eps1
  local mu = tau * A:diag():norm"inf"
  -- iterate
  while not found and k < kmax do
    k = k + 1
    solveh(set(C, A), mu, g, h) -- solve (A + mu * I) * h = -g for h
    if h:norm() <= eps2 * (eps2 + x:norm()) then
      found = true
    else
      xnew:set(x):add(h, true) -- xnew = x + h
      f:apply(ff(xnew))
      local Fnew = dot(f, f)
      local rho = (F - Fnew) / (mu * dot(h, h) - dot(h, g)) -- gain ratio
      if rho > 0 then -- step acceptable?
        x:set(xnew);  J:apply(Jf(xnew))
        A:set(0):hemul(J, true) -- A = J' * J
        g:set(0):mmul(J, f, "t") -- g = J' * f
        F = Fnew
        found = g:norm"inf" <= eps1
        mu = mu * max(1 / 3, 1 - (2 * rho - 1) ^ 3)
        vu = 2
      else
        mu = mu * vu
        vu = vu * 2
      end
    end
  end
  return x
end

-- Example 1.1
local t = {...}
local y = {...}
local exp = math.exp

local f = function (x)
  local x1, x2, x3, x4 = x[1], x[2], x[3], x[4]
  return function (i) -- each f_i(x)
    local ti = t[i]
    return y[i] - x3 * exp(x1 * ti) - x4 * exp(x2 * ti)
  end
end

local J = function (x) -- Jacobian
  local x1, x2, x3, x4 = x[1], x[2], x[3], x[4]
  local df = {
    [1] = function (i) local ti = t[i] return -x3 * ti * exp(x1 * ti) end,
    [2] = function (i) local ti = t[i] return -x4 * ti * exp(x2 * ti) end,
    [3] = function (i) return -exp(x1 * ti) end,
    [4] = function (i) return -exp(x2 * ti) end
  }
  return function (i, j) return df[j](i) end -- each df_i(x)/dx_j
end

