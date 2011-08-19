--[[
-- matrix.lua
-- Multidimensional matrix library for NumericLua
-- Luis Carvalho (lexcarvalho@gmail.com)
-- See Copyright Notice in numlua.h
--]]

-- Methods
local new, copy = matrix.new, matrix.copy
local sum, diag = matrix.sum, matrix.diag
local get, set, size = matrix.get, matrix.set, matrix.size
local concat, shape = matrix.concat, matrix.shape
local type, eps = numlua.type, mathx.eps

local unpack, assert, ipairs = unpack, assert, ipairs
local setmetatable, select = setmetatable, select
local floor, max = math.floor, math.max

local function checkmatrix (m)
  local t = type(m)
  assert(t == "matrix", "matrix expected, got " .. t)
  return m
end

local transpose = matrix.transpose
matrix.t = transpose -- handy alias
function matrix.ctranspose (m) return transpose(m, true) end

local zeros = function (...) return set(new(...), 0) end
function matrix.ones (...) return set(new(...), 1) end
function matrix.eye (n, c) return set(set(new(n, n, c), 0), "D", 1) end
matrix.zeros = zeros

-- [ Metamethods ]
local mt = getmetatable(new(1))
local add, mul, mmul = matrix.add, matrix.mul, matrix.mmul
local div, ls = matrix.div, matrix.ls
local section, slice = matrix.section, matrix.slice

mt.__add = function (a, b)
  if type(a) == "number" or type(a) == "complex" then
    return add(b, a)
  end
  return add(a, b)
end

mt.__sub = function (a, b)
  if type(a) == "number" or type(a) == "complex" then
    return add(-b, a)
  end
  if type(b) == "number" or type(b) == "complex" then
    return add(a, -b)
  end
  return add(a, b, -1)
end

mt.__mul = function (a, b)
  if type(a) == "number" or type(a) == "complex" then
    return mul(b, a)
  end
  if type(b) == "number" or type(b) == "complex" then
    return mul(a, b)
  end
  local da, db = size(a, "#"), size(b, "#")
  local n, m = size(a, 1), size(b, 2)
  local iscomplex = a:iscomplex() or b:iscomplex()
  if da == 1 and db == 1 then -- outer product?
    return mmul(zeros(n, n, iscomplex), a, b)
  end
  if da == 1 then -- v * A?
    return mmul(zeros(m, iscomplex), b, a, "T")
  end
  if db == 1 then -- A * v?
    return mmul(zeros(n, iscomplex), a, b)
  end
  -- da = db = 2:
  return mmul(zeros(n, m, iscomplex), a, b)
end

mt.__mod = ls

mt.__div = function (a, b)
  if type(a) == "number" or type(a) == "complex" then
    return div(b, a, true)
  end
  if type(b) == "number" or type(b) == "complex" then
    return div(a, b)
  end
  local x = ls(transpose(b), transpose(a))
  return size(x, "#") == 2 and transpose(x) or x
end

-- TODO: __call using string triplets (based on section)
-- [ _section_(m,  "f1:l1:s1, f2:l2:s2, ...") <=> m[[f1:l1:s1,...]] ]
--                 ^----'triplet string'---^
mt.__call = function (a, ...)
  local t = select(1, ...)
  if type(t) == "table" then return section(a, t) end
  return slice(a, ...)
end

local cabs, linspace = complex.abs, matrix.linspace
function matrix.seq (a, b, step)
  local s = step or 1
  local n = floor(cabs((b - a) / s + 1))
  return linspace(a, b, n)
end

function matrix.trace (m)
  local r, c = shape(checkmatrix(m))
  assert(size(m, "#") == 2 and r == c, "square matrix expected")
  return sum(diag(m))
end


-- [ Logical ]
local find, ifelse, which = matrix.find, matrix.ifelse, matrix.which

function matrix.any (m, cond) return find(m, cond) ~= nil end
function matrix.all (m, cond) return find(m, cond, true) == nil end

-- count(m, cond) <=> fold(m, \a,e(a + (cond(e) and 1 or 0)), 0)
--                <=> sum(ifelse(copy(m), cond, 1, 0))
function matrix.count (m, cond) return which(m, cond, "#") end
function matrix.merge (x, y, mask) return ifelse(copy(mask), 1, x, y) end
function matrix.pack (m, mask) return which(m, mask, "v") end
function matrix.unpack (v, mask, m) return set(m, which(m, mask), v) end


-- [ From/To table conversions ]

local function checkvector (t, iscomplex)
  local isvector, iscomplex = true, iscomplex or t.complex
  for i, v in ipairs(t) do
    if type(v) ~= "number" and type(v) ~= "complex" then
      isvector = false
      if type(v) == "matrix" and size(v, "#") == 1 then
        v = v[1]
        t[i] = v
        isvector = true
      end
    end
    if isvector then
      iscomplex = iscomplex or type(v) == "complex"
    else
      break
    end
  end
  return isvector, iscomplex
end

local function fromtable (t, iscomplex)
  assert(type(t) == "table", "table expected")
  local isvector, iscomplex = checkvector(t, iscomplex)
  if isvector then -- base case?
    local v = new(#t, iscomplex)
    for i, e in ipairs(t) do v[i] = e end
    return v
  end
  -- recursion
  for i, v in ipairs(t) do
    if type(v) == "table" then -- recurse?
      t[i] = fromtable(v, iscomplex)
    end
  end
  -- fix if complex
  iscomplex = false
  for _, v in ipairs(t) do
    iscomplex = iscomplex or v:iscomplex()
  end
  if iscomplex then
    for i, v in ipairs(t) do
      if not v:iscomplex() then
        t[i] = v:complex()
      end
    end
  end
  return concat(unpack(t))
end
matrix.fromtable = fromtable

local function totable (m)
  assert(type(m) == "matrix", "matrix expected")
  local d, t = size(m, "#"), {}
  for i = 1, #m do
    t[i] = d == 1 and m[i] or totable(m[i])
  end
  return t
end
matrix.totable = totable

function matrix.list (m)
  checkmatrix(m)
  for i, e in m:entries(true) do
    local t = {m:eindex(i)}
    t[#t + 1] = e
    print(unpack(t))
  end
end


local function formatnumber (x, d)
  local fmt = d and ("%." .. d .. "f") or "%g"
  return fmt:format(x)
end

local signbit = mathx.signbit
local function formatcomplex (c, d)
  local re, im = c:real(), c:imag()
  local fmt = signbit(im) and "%s%si" or "%s+%si"
  return fmt:format(formatnumber(re, d), formatnumber(im, d))
end

local function getmaxlen (fmt, d)
  return function (l, e) return max(l, #fmt(e, d)) end
end

local tconcat = table.concat
local function prettyaux (v, ml, fmt, d) -- print vector with max length ml
  local t = {}
  for i = 1, #v do
    local vi = fmt(v[i], d)
    t[i] = (" "):rep(3 + ml - #vi) .. vi
  end
  return tconcat(t)
end

function matrix.pretty (m, d) -- `d` is number of decimal places
  assert(size(checkmatrix(m), "#") <= 2, "two-dimensional matrix expected")
  local fmt = m:iscomplex() and formatcomplex or formatnumber
  local ml = m:fold(getmaxlen(fmt, d), 0) -- max length
  if size(m, "#") == 1 then
    return prettyaux(m, ml, fmt, d)
  else -- m:size"#" == 2
    local t = {}
    for i = 1, #m do t[i] = prettyaux(m[i], ml, fmt, d) end
    return tconcat(t, "\n")
  end
end

-- set metatable for class
matrix = setmetatable(matrix, {
  __call = function(_, ...)
    return type(select(1, ...)) == "table" and fromtable(...) or new(...)
  end
})


-- [ Aggregators ]

local function opfold (f, init)
  local c
  return function (i, e)
    if i == 1 then c = init end
    c = f(c, v)
    return c
  end
end

local sum2 = function(x, y) return x + y end
function matrix.cumsum (m)
  return m:apply(opfold(sum2, 0), true)
end

local prod2 = function(x, y) return x * y end
function matrix.cumprod (m)
  return m:apply(opfold(prod2, 1), true)
end
local prod = function (m) return m:fold(prod2, 1) end
matrix.prod = prod


-- [ Linear algebra ]
local chol, lu, svd = matrix.chol, matrix.lu, matrix.svd

function matrix.kronecker (a, b)
  assert(size(a, "#") == 2 and size(b, "#") == 2,
    "two-dimensional matrix expected")
  local ra, ca, ica = shape(a, 1, true)
  local rb, cb, icb = shape(b, 1, true)
  local iscomplex = ica or icb
  if iscomplex then
    if not ica then a = a:complex() end
    if not icb then b = b:complex() end
  end
  local x = new(ra * rb, ca * cb, iscomplex)
  local indexr, indexc = {}, {}
  local index = {indexr, indexc}
  for i = 1, ra do
    local ai = a[i]
    indexr[1], indexr[2] = (i - 1) * rb + 1, i * rb
    for j = 1, ca do
      indexc[1], indexc[2] = (j - 1) * cb + 1, j * cb
      mul(set(section(x, index), b), ai[j], true) -- x[index] = a[i][j] * b
    end
  end
  return x
end

function matrix.isposdef (m)
  local c, msg = chol(checkmatrix(m))
  if c == nil then error(msg) end
  return not c == false
end

function matrix.det (m)
  local c = assert(lu(copy(checkmatrix(m)), true))
  return prod(diag(c))
end

function matrix.cond (m)
  local s = assert(svd(checkmatrix(m), "n")) -- just singular values
  return s[1] / s[#m]
end

-- effective rank from singular values `s`, max dim `m`, tolerance `tol`
local lt = function (x) return function(e) return e < x end end
local function srank (s, m, tol)
  local tol = tol or 0
  if tol <= 0 then -- set default tolerance?
    tol = m * eps * s[1]
  end
  local r = s:find(lt(tol))
  return r and r - 1 or #s
end

function matrix.rank (m, tol)
  local s = assert(svd(checkmatrix(m), "n")) -- just singular values
  return srank(s, max(shape(m)), tol)
end

function matrix.null (m, tol)
  local u, s, vh = assert(svd(checkmatrix(m)))
  local nr, nc = shape(m)
  local rank = srank(s, max(nr, nc), tol)
  return rank < nc and slice(vh, rank + 1) or nil
end

function matrix.orth (m, tol)
  local u = copy(checkmatrix(m))
  local s = assert(svd(u, "l"))
  local rank = srank(s, max(shape(m)), tol)
  return u{{}, {1, rank}} -- columns from 1 to rank
end

-- pseudo-inverse
function matrix.pinv (m, tol)
  local u, s, vh = assert(svd(checkmatrix(m)))
  local nr, nc = shape(m)
  local rank = srank(s, max(nr, nc), tol)
  local v = slice(vh, 1, rank)
  for i = 1, rank do -- inv(s) * vh
    v[i]:div(s[i], false, true) -- v[i,:] = v[i,:] / s[i], in-place
  end
  return zeros(nc, nr):mmul(v, u{{}, {1, rank}}, "c", "c")
end


-- basic LS linear model fitting
function matrix.lm (a, b, svd)
  local m, n = shape(checkmatrix(a))
  assert(m >= n, "system is underdetermined")
  assert(checkmatrix(b):size"#" == 1, "single RHS expected")
  local x, rank = ls(a, b, svd)
  -- report summary statistics
  local coef = slice(x, 1, n)
  local rss = (b - a % coef):norm() ^ 2
  local rss0 = (b - b:sum() / m):norm() ^ 2
  local df = m - rank
  local F = df / (rank - 1) * (rss0 / rss - 1)
  local pvalue = 1 - stat.pf(F, rank - 1, df)
  return {coef = coef, rss = rss, df = df, F = F, pvalue =  pvalue}
end

