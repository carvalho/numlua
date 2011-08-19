require "numlua"

-- Gauss-Seidel linear solver for A * x = b, where A is a square matrix
local zeros, swap = matrix.zeros, matrix.swap
function gaussseidel (A, b, x0, tol)
  local tol = tol or 1e-6
  local n = A:size(1)
  assert(A:size"#" == 2 and A:size(2) == n, "square matrix expected")
  assert(b:size"#" == 1 and b:size(1) == n, "consistent vector expected")
  local S = zeros(n, n)
  S.U = A -- S has upper triangle of A
  S:mul(-1, true) -- S = -S (inplace)
  S:trmul(A, "l", true) -- S = A.L^(-1) * S
  local u = b:copy():trmul(A, "l", true) -- u = A.L^(-1) * b
  local x = x0:copy()
  local xc = zeros(n) -- current solution
  repeat
    xc._ = u
    xc:mmul(S, x) -- xc = xc + S * x
    swap(x, xc)
    xc:add(x, -1, true) -- xc = xc - x (inplace)
  until xc:norm() < tol
  return x
end

