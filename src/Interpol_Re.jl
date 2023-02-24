using Interpolations
using DataFrames,XLSX

Reinfo = DataFrame(XLSX.readtable("src/Data/Re_ch.xlsx", "Sheet1")...)

interp_linear = linear_interpolation([Reinfo.Z,Reinfo.Y], Reinfo.U)
unique(Reinfo.Z)
unique(Reinfo.Y)

xs = collect(1:0.2:5)
ys = collect(1:0.2:5)
A = @. log(xs) + log(ys')

LinearInterpolation((Reinfo.Z,Reinfo.Y),Reinfo.U)

Plots.scatter(xs,A)
Plots.scatter!(val,interp_linear(val))

# Data
x = collect(range(-2, 3, length=20))
y = collect(range(3, 4, length=10))
z = @. cos(x) + sin(y')
# Interpolatant object
itp = LinearInterpolation((x, y), z)