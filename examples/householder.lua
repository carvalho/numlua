require "numlua"

-- Householder reflections:
-- householder returns a function that multiplies an orthonormal symmetric
-- matrix P by its argument, where P = I - u * u' / c, c = 1/2 * u' * u
local dot, size, add = matrix.dot, matrix.size, matrix.add
function householder (u)
  assert(size(u, "#") == 1, "vector expected")
  local n = size(u)
  local c = 0.5 * dot(u, u)
  return function (x)
    assert(size(x, "#") == 1, "vector expected")
    assert(size(x) == n, "inconsistent dimension")
    return add(x, u, -dot(u, x) / c)
  end
end

-- test: reflect1(x) returns -copysign(norm(x), x[1]) * e_1
local zeros, copysign = matrix.zeros, mathx.copysign
function reflect1 (x)
  local e1 = zeros(size(x))
  e1[1] = 1
  return householder(add(x, e1, copysign(norm(x), x[1])))(x)
end

