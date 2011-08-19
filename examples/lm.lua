-- Linear models

require "numlua"

local zeros, eye, copy = matrix.zeros, matrix.eye, matrix.copy
local ones, concat, shape = matrix.ones, matrix.concat, matrix.shape
local mul, mmul, trmul = matrix.mul, matrix.mmul, matrix.trmul
local diag, sqrt, norm = matrix.diag, matrix.sqrt, matrix.norm
local linspace, sum, qr = matrix.linspace, matrix.sum, matrix.qr
local min, max, abs = math.min, math.max, math.abs
local eps, pt, pf = mathx.eps, stat.pt, stat.pf
local type = numlua.type

local function fmt (x, nd)
  local nd = nd or 4 -- #digits
  return string.format("%." .. nd .. "f", x)
end
local function label (n)
  n = n - 1
  return n == 0 and "(Intercept)" or ("X" .. n)
end
local function printtable (t, ns)
  local ns = ns or 1 -- # spaces
  local ml = 0
  for _, r in ipairs(t) do
    for i, v in ipairs(r) do
      local s = type(v) == "string" and v or fmt(v)
      r[i] = s; s = #s
      if ml < s then ml = s end
    end
  end
  for _, r in ipairs(t) do
    for i, v in ipairs(r) do
      r[i] = (" "):rep(ns + ml - #v) .. v
    end
    print(table.concat(r))
  end
end

local function summary (m)
  local coef, pvt, rank = m.coef, m.pivot, m.rank
  local res, fit = m.residuals, m.fitted
  local sigma = sqrt(diag(m.covmatrix))
  -- statistics
  local n = #res
  local df = n - rank
  local dd = rank - 1 -- df0 - df (with intercept)
  local mss = norm(fit - sum(fit) / #fit) ^ 2
  local rss = norm(res) ^ 2
  -- report
  local order = linspace(1, #pvt):pivot(pvt)
  print("Coefficients:")
  local t = {}
  -- header
  t[1] = {"Predictor", "Estimate", "Std. Error", "t value", "Pr(>|t|)"}
  -- coefficients
  for i = 1, #pvt do
    local oi = order[i]
    if oi <= rank then
      local bi = coef[oi]
      local si = sigma[oi]
      local ti = bi / si
      local pi = 2 * pt(-abs(ti), df) -- p-value
      t[#t + 1] = {label(i), fmt(bi), fmt(si), fmt(ti), fmt(pi, 5)}
    end
  end
  printtable(t); print()
  -- regression
  local se = rss / df
  print("Residual standard error: " .. fmt(math.sqrt(se))
    .. " on " .. df .. " degrees of freedom")
  local r2 = mss / (mss + rss)
  local adjr2 = 1 - (1 - r2) * (n - 1) / df
  print("Multiple R-squared: " .. fmt(r2)
    .. ",  Adjusted R-squared: " .. fmt(adjr2))
  local f = mss / dd / se
  local pv = 1 - pf(f, dd, df)
  print("F-statistics: " .. fmt(f) .. " on " .. dd .. " and " .. df
    .. " DFs,  p-value: " .. fmt(pv, 5))
end


-- find rank from right triangular matrix, output from _,r,_=qr(a,true)
local function getrank (r, tol)
  local n = min(shape(r))
  local c = r[n][n] ^ 2
  while c < tol and n > 1 do
    n = n - 1
    local rn = r[n]
    c = max(c + rn[n + 1] ^ 2, rn[n] ^ 2)
  end
  return n
end


local mt = {__index = {summary = summary}}
function lm (y, ...)
  -- read variables and build design matrix
  assert(type(y) == "matrix" and y:size"#" == 1,
    "vector expected for response variable")
  local m = #y
  local t = {ones(m)}
  for i = 1, select("#", ...) do
    local x = select(i, ...)
    assert(type(x) == "matrix" and x:size"#" == 1,
      "vector expected for predictor")
    assert(#x == m, "inconsistent dimension for predictor")
    t[#t + 1] = x
  end
  local n = #t
  t[n + 1] = true
  local X = concat(unpack(t))
  local tol = max(m, n) * eps
  -- fit model
  local q, r, p = qr(X, true)
  local rank = getrank(r, tol)
  local qt = q{{}, {1, rank}}
  local rt = r{{1, rank}, {1, rank}}
  local yt = mmul(zeros(rank), qt, y, 't') -- bt = Q'(1:rank,:) * b
  local z = trmul(yt, rt, 'u', true) -- rt * z = bt
  local fit = qt * trmul(copy(z), rt, 'u') -- fitted values: qt * rt * z
  local res = y - fit -- residuals
  local xtxi = trmul(eye(rank), rt, 'u', true, 't')
  trmul(xtxi, rt, 'u', true) -- xtxi = (r' * r)^(-1)
  mul(xtxi, norm(res) ^ 2 / (m - rank), true) -- xtxi = xtxi * sigmahat^2
  return setmetatable({coef = z, pivot = p, rank = rank,
    residuals = res, fitted = fit, covmatrix = xtxi}, mt)
end

