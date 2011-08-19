-- a version of linear system solver
-- modified from "Fortran 90/95 explained"
-- by Metcalf and Reid (p. 131-132)
-- to show array sections

require "numlua"

-- a is nxn matrix, b is a vector
-- the solution is returned in b
local swap = matrix.swap
function solve (a, b)
  local n = #b
  assert(a:size"#" == 2 and a:size(1) == a:size(2), "a must be square")
  assert(#a == n, "a and b sizes must be equal")

  for j = 1, n do -- update elements in column j
    for i = 1, j - 1 do
      a{{i+1,n}, {j,j}}._ = a{{i+1,n}, {j,j}} - a{{i+1,n}, {i,i}} * a[i][j]
      --a{{i+1,n}, {j,j}}:add(a{{i+1,n}, {i,i}}, -a[i][j], true) -- in-place
    end

    if j < n then
      -- find pivot and check its size
      local _, i = a{{j,n},{j,j}}:norm"m" -- argmax(abs(a[[j:n, j]]))
      i = i + j - 1

      -- if necessary apply row interchange
      if i ~= j then
        swap(a[i], a[j])
        b[j], b[i] = b[i], b[j]
      end

      a{{j+1,n},{j,j}}._ = a{{j+1,n},{j,j}} / a[j][j]
      --a{{j+1,n},{j,j}}:div(a[j][j], false, true) -- in-place
    end
  end

  -- forward substitution
  for i = 1, n - 1 do
    b{{i+1,n}}._ = b{{i+1,n}} - a{{i+1,n},{i,i}} * b[i]
    --b{{i+1,n}}:add(a{{i+1,n},{i,i}}, -b[i], true) -- in-place
  end

  -- back substitution
  for j = n, 2, -1 do
    b[j] = b[j] / a[j][j]
    b{{1,j-1}}._ = b{{1,j-1}} - a{{1,j-1},{j,j}} * b[j]
    --b{{1,j-1}}:add(a{{1,j-1},{j,j}}, -b[j], true) -- in-place
  end
  b[1] = b[1] / a[1][1] -- last case

  return b
end

-- main
--[
local iscomplex = true
local n = 200
a = matrix.zeros(n, n, iscomplex)
a:apply(function (i,j) return i + j * (0.5 - rng.runif()) end)

b = matrix.zeros(n, iscomplex)
b:apply(function (i) return i * rng.runif() end)
--]]
--[[
a = matrix{{2,1,-1},{-3,-1,2},{-2,1,2}}
b = matrix{8, -11, -3}
--]]

--[[
print"a matrix"
list(a)
print"b vector"
list(b)
--]]

print"solution by matrix"
matrix.list(a:inv() * b)

print"present solution"
matrix.list(solve(a:copy(), b:copy()))

