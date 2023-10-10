
using GridapReynoldsEquation
using LinearAlgebra
using Plots

exampleNo = 1
u, f, resultFolder, ub, u₀ = selectExample(exampleNo)
x = range(0, stop = 2*pi, length = 100)
y = range(-1, stop = 1, length = 100)

pU = surface(x, y, (x,y)->u(VectorValue(x,y)))
pF = surface(x, y, (x,y)->f(VectorValue(x,y)))

display(plot(pU,pF,layout=(2,1),size=(600,600/1.618*2)))