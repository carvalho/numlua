require "numlua"

local zeros, set = matrix.zeros, matrix.set
local eig, mul, dot = matrix.eig, matrix.mul, matrix.dot

-- if F[n] is the n-th Fibonacci number,
-- {{F[n-1], F[n]}, {F[n], F[n+1]}} = {{0, 1}, {1, 1}} ^ n
local s, V = eig(matrix{{0, 1}, {1, 1}}, 'r', true)
local v = mul(V[1], V[2]) -- element-wise
local w = zeros(2) -- workspace
function fib (n)
  return dot(v, set(w, s) ^ n) -- V * diag(s ^ n) * t(V)
end

for i = 1, 100 do print(i, fib(i)) end

