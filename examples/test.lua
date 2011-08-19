require "numlua.seeall"

a = matrix{{1,2,3},{4,5,6},{7,8,0}}
b = matrix{7, 13, -8}

-- solve equation using LU
-- a:pivot(p) = l * u or, since a:pivot(p) = eye(#a):pivot(p) * a,
-- a = t(eye(#a):pivot(p)) * l * u
function lslu (a, b)
  local l, u, p = lu(a)
  local x = b:pivot(p) -- x = eye(#a):pivot(p) * b
  x:trmul(l, 'l', true) -- x = inv(l) * x
  x:trmul(u, 'u', true) -- x = inv(u) * x
  return x
end

print("LS by LU")
list(lslu(a,b))

-- solve equation using QR
-- q,r,p = qr(a,true) -- a = q * r * t(p)
function lsqr (a, b)
  local q,r = qr(a)
  local x = zeros(#b):mmul(a, b, 't') -- x = t(a) * b
  x:trmul(r, 'u', true, 't') -- x = t(inv(r)) * x
  x:trmul(r, 'u', true) -- x = inv(r) * x
  return x
end

print("LS by QR")
list(lsqr(a,b))

--[[
-- iterative refinement
h = b - a * x
e = zeros(#b):mmul(a, h, 't') -- e = t(a) * r
e:trmul(r, 'u', true, 't') -- e = t(inv(r)) * e
e:trmul(r, 'u', true) -- e = inv(r) * e
x:add(e, 1, true) -- in-place
list(x)
--]]

-- solve equation using LS
print("LS by LS [QR]")
list(ls(a, b)) -- using QR

print("LS by LS [SVD]")
list(ls(a, b, true)) -- using SVD

