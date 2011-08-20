--[[
-- seeall.lua
-- Easy prototyping in NumericLua
-- Luis Carvalho (lexcarvalho@gmail.com)
-- See Copyright Notice in numlua.h
--]]

require "numlua"
pcall(require, "help")

local classes = {"math", "mathx", "complex", "matrix"}
local words = setmetatable({}, {__index = function() return 0 end})
for _, c in pairs(classes) do
  for k in pairs(_G[c]) do
    words[k] = words[k] + 1
  end
end

local nltype = numlua.type
local function registermath (name)
  _G[name] = function (x, ...)
    local xtype = nltype(x)
    if xtype == "number" then
      local f = math[name] or mathx[name] or complex[name]
      return f(x, ...)
    elseif xtype == "complex" then
      return complex[name](x, ...)
    else
      return matrix[name](x, ...)
    end
  end
end

local register = function (t, prefix)
  for k, v in pairs(t) do
    if not _G[k] then
      _G[k] = v
    elseif prefix then
      _G[prefix .. k] = v
    end
  end
end

-- register common math functions first
for k, v in pairs(words) do
  if v > 1 then -- conflict?
    registermath(k) -- register common function
  end
end
register(math, "x")
register(mathx, "x")
register(complex, "c")
register(matrix, "m")
register(stat, "s")
register(rng, "r") -- conflicts with matrix: new, copy

type = numlua.type -- override
opmode = numlua.opmode

-- set matrix __tostring to `pretty`: more convenient for interpreter
getmetatable(matrix(1)).__tostring = matrix.pretty

