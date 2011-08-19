require "numlua"

local zeros, dot, sqrt = matrix.zeros, matrix.dot, math.sqrt

-- Cholesky decomposition: given positive-definite matrix A, returns U such
-- that A = U' * U
function cholesky (A)
  local n = A:size(1)
  assert(A:size"#" == 2 and A:size(2) == n, "square matrix expected")
  local U = zeros(n, n)
  local C = {} -- columns of U
  for i = 1, n do C[i] = U:col(i) end
  -- factorization:
  -- i = 1
  local Ai, Ui = A[1], U[1]
  local uii = sqrt(Ai[1])
  Ui[1] = uii
  for j = 2, n do Ui[j] = Ai[j] / uii end
  -- i > 1
  for i = 2, n do
    Ai, Ui = A[i], U[i]
    local uci = C[i](1, i - 1) -- U(i, 1:(i-1))
    uii = sqrt(Ai[i] - dot(uci, uci))
    Ui[i] = uii
    for j = i + 1, n do
      local ucj = C[j](1, i - 1) -- U(j, 1:(i-1))
      Ui[j] = (Ai[j] - dot(uci, ucj)) / uii
    end
  end
  return U
end

-- test: cholesky(pascal(n))
function pascal (n)
  local A = zeros(n, n)
  A[1]._, A:col(1)._ = 1, 1 -- base
  for i = 2, n do
    for j = 2, n do
      A[i][j] = A[i - 1][j] + A[i][j - 1]
    end
  end
  return A
end

