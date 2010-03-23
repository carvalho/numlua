require "numlua"

function pp1 (m)
  for i = 1, m:size(1) do print(i, m[i]) end
end

function pp2 (m)
  for i = 1, m:size(1) do
    local r = m[i]
    for j = 1, m:size(2) do
      print(i, j, r[j])
    end
  end
end

function pp3 (m)
  for i = 1, m:size(1) do
    local r = m[i]
    for j = 1, m:size(2) do
      local v = r[j]
      for k = 1, m:size(3) do
        print(i, j, k, v[k])
      end
    end
  end
end

m = matrix.zeros(2,3,4)
for i=1,m:size(1) do
  for j=1,m:size(2) do
    for k=1,m:size(3) do
      m[i][j][k] = math.random()
    end
  end
end

x = matrix.zeros(3,4,true)
for i=1,x:size(1) do
  for j=1,x:size(2) do
    x[i][j] = complex(math.random(), math.random())
  end
end

c = matrix.zeros(x:size(1), x:size(2), true)
