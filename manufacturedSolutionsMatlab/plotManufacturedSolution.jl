
#include("./manufacturedSolutionsMatlab/f_example81.jl")

include("./manufacturedSolutionsMatlab/f_example82.jl")

using LinearAlgebra
using Plots
x = range(0, stop = 2*pi, length = 100)
y = range(-1, stop = 1, length = 100)

surface(x, y, (x,y) -> fm(x,y))