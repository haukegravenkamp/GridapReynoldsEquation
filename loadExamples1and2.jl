## manufactured solutions used in Sections 8.1 and 8.2 in the paper

# The chosen solutions are shown below. The corresponding body loads have been computed using the Matlab script available at  https://github.com/haukegravenkamp/manufacturedSolutions and stored in the files f_example81.jl and f_example82.jl 

# The Matlab script returns the body load fm(x,y) as a function of x and y. For gridap, we need it as a function of x, with x[1] and x[2] being the two Cartesian coordinates. For this reason, we define below the function f(x).

if exampleNo == 1       # example 8.1
    # read body load from the file
    include("./manufacturedSolutionsMatlab/f_example81.jl")
    # choose folder to store results
    resultFolder = "./results/example81/"
    # chosen solution, Eq. (40) in paper 
    u(x) = 1 / 6 * (1 - cos(2x[1])) * sin(x[1]) * (1 + cos(π * x[2]))
elseif exampleNo == 2   # example 8.2
    # read body load from the file
    include("./manufacturedSolutionsMatlab/f_example82.jl")
    # choose folder to store results
    resultFolder = "./results/example82/"
    # chosen solution, Eq. (42) in paper 
    u(x) = 0.25 * ((1 - exp(100 * x[1] / (2 * π))) / (1 - exp(100)) - 1 + (cos(x[1] / 2) + 1) / 2) * (1 + cos(π * x[2]))
end

function f(x)
    # just because the load is stored as f(x,y) instead of f(x) with a vector x
    return fm(x[1], x[2])
end