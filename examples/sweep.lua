-- sweep operator from cols i to j (defaults to i):
-- when solving y = X * b using least-squares,
--
-- [ X' * X | X' * y ]   SWEEP   [ (X' * X) ^ (-1) |  b^ ]
-- [-----------------] --------> [-----------------------]
-- [ y' * X | y' * y ]  X' * X   [       -b^'      | RSS ]

require "numlua"

local function sweepaux (a, k) -- single column sweep
  local ck = a:col(k)
  local b = ck:copy()
  local ak, d = a[k], b[k]
  ak:div(d, false, true) -- a[k] = a[k] / b[k], in-place
  for i = 1, #a do
    if i ~= k then
      a[i]:add(ak, -b[i], true) -- a[i] = a[i] - b[i] * a[k], in-place
    end
  end
  ck._ = b:div(-d, false, true) -- ck = b, b = -b / d in-place
  ck[k] = 1 / d
  return a
end

function sweep (a, i, j)
  for k = i, (j or i) do sweepaux(a, k) end
  return a
end

-- example from [...]
local concat, ones, c = matrix.concat, matrix.ones, matrix.c
local seq, t = matrix.seq, matrix.transpose
local inv, norm = matrix.inv, matrix.norm
local X = concat(ones(3), seq(1, 3), ones(3), true)
      .. concat(ones(3), seq(1, 3), -ones(3), true)
local y = matrix.c(1, 3, 3, 2, 2, 1)
local xtx, xty = t(X) * X, t(X) * y
local a = concat(xtx, xty, true) .. c(xty, norm(y) ^ 2)
print("inv(X' * X):")
print(matrix.pretty(inv(xtx)))
print("LS solution:", matrix.pretty(X % y))
print("RSS:", norm(y - X * (X % y)) ^ 2)
print("Sweep:")
print(matrix.pretty(sweep(a, 1, 3)))

