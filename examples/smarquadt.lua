
-- Broyden's rank one update
local function broyden (B, dx, df)
  local u = mmul(df, B, dx, "n", "n", -1) -- u = df - B * dx
  return mmul(B, u, dx, "n", "t", dot(dx, dx)) -- B = B + u * dx' / u' * u
end

