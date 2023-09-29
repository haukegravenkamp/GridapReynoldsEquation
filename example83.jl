# load packages and modules
using Revise                    # helps tracking changes during development
using GridapReynoldsEquation    # use Reynolds module 
using Gridap                    # the FE packages Gridap           
using Plots                     # for plotting results


## INPUT
# where to store results for plotting in paraview
resultFolder = "./results/example83/"

# problem definition
paramProblem = Dict(
    :nx => 100,                 # number of elements in x
    :ny => 32,                  # number of elements in x
    :hRefinement => 0,          # refinement steps starting from nx, ny
    :order => 1,                # element order 
    :artDiff => ~true,          # use artificial diffusion 
    :stabilizationType => 2,    # 0: no stabilization, 1: ASGS, 2: OSGS
    :K => 1.0,                  # extra factor for artificial diffusion
    :ζ => 0.6,                  # constants ζ, xₐ in gap function, Eq. (2)
    :xₐ => 7 / 9 * pi,
    :u̅ => 0.98,                 # parameter in regularization, Eq. (3)
    :u => x -> 0.0,             # analytical solution, NA for this example
    :ub => x -> 0.0,            # boundary condition, function of x
    :f => x -> 0.0,             # body load as read from file 
    :u₀ => x -> 1.0,            # initial guess, function of x
    :isNL => true               # whether problem is nonlinear 
)

# solver parameters
paramSolver = Dict(
    :numIt => 500,              # max number of iterations
    :linType => 1,              # linearization type; 0: Picard, 1: Newton 
    :initialPicard => 3,        # in case of Newton, do this many Picard iterations first, then switch
    :backTracking => ~true      # use backtracking (only in case of Newton)
)

# solve problem using Reynolds module
uh, dΩ, hMin, h, solverCache, residuals = runReynolds(paramProblem, paramSolver, resultFolder);
# plot residual 
plotResiduals(residuals, paramProblem[:order])







