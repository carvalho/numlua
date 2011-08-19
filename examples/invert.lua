-- inverting a matrix
--
-- adapted from p. 124 of Redwine:
-- "Updating to Fortran 90"
-- to show array sections in numlua

require "numlua"

-- a is nxn matrix
local swap, eye = matrix.swap, matrix.eye
local function invert (a)
  assert(a:size"#" == 2 and a:size(1) == a:size(2),
      "square matrix expected")
  local n = #a

  local a = a:concat(eye(n, a:iscomplex()), true) -- colcat

  for j = 1, n do
    -- find pivot and check its size
    local _, i = a:section{{j,n},{j,j}}:norm"m" -- argmax(abs(a[[j:n,j]]))
    i = i + j - 1

    -- if necessary apply row interchange
    if i ~= j then swap(a[i], a[j]) end

    --a[j]._ = a[j] / a[j][j]
    a[j]:div(a[j][j], false, true) -- in-place

    -- loop over each row of matrix
    for i = 1, n do
      if i ~= j then
        --a[i]._ = a[i] - a[j] * a[i][j]
        a[i]:add(a[j], -a[i][j], true) -- inplace
      end
    end
  end

  return a:section{{}, {n+1, 2*n}}:copy()
end

-- main

local iscomplex = true
local n = 5
a = matrix.zeros(n, n, iscomplex)
a:apply(function (i, j) return i + j * (0.5 - math.random()) end)

print"a matrix"
a:list()

print"solution by matrix"
a:inv():list()

print"present solution"
invert(a):list()

