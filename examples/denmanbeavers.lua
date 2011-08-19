require "numlua"

-- Square-root of a matrix A by Denman-Beavers algorithm:
local TOL, MAXITERS, inf = mathx.eps, 1000, mathx.inf
local inv, norm, abs = matrix.inv, matrix.norm, math.abs
function dbsqrtm (A, tol, maxiters)
  local opmode = numlua.opmode(true) -- set in-place operations
  local tol, maxiters = tol or TOL, maxiters or MAXITERS
  local iY, iZ = matrix.new(A:shape()), matrix.new(A:shape()) -- buffers
  local s, n = inf, 1 -- norm(Y), #iterations
  local Y, Z = A:copy(), matrix.eye(#A) -- Y_0, Z_0 = A, I
  while true do
    -- Y_{k+1} = 1/2 * (Y_k + inv(Z_k))
    -- Z_{k+1} = 1/2 * (Z_k + inv(Y_k))
    iY._, iZ._ = Y, Z
    Y, Z = (Y + inv(iZ)) / 2, (Z + inv(iY)) / 2
    -- check termination
    local f = norm(Y)
    if abs(f - s) <= tol or n > maxiters then break end
    s, n = f, n + 1
  end
  numlua.opmode(opmode) -- restore previous mode
  return Y, Z, s
end

