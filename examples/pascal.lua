require "numlua"

-- Pascal matrix of order `n`
local zeros, set = matrix.zeros, matrix.set
function pascal (n)
  local p = zeros(n, n)
  set(p[1], 1)
  set(p:col(1), 1)
  for i = 2, n do
    for j = 2, n do
      p[i][j] = p[i - 1][j] + p[i][j - 1]
    end
  end
  return p
end

