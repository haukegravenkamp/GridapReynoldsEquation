# load packages and modules
using Revise                    # helps tracking changes during development
using GridapReynoldsEquation    # use Reynolds module 
using Gridap                    # the FE packages Gridap           
using Plots                     # for plotting results

## INPUT
u, f, resultFolder, ub, u₀ = selectExample(3)

# problem definition
paramProblem = Dict(
    :nx => 100,                 # number of elements in x
    :ny => 32,                  # number of elements in y
    :hRefinement => 0,          # refinement steps starting from nx, ny
    :order => 1,                # element order 
    :artDiff => false,          # use artificial diffusion 
    :stabilizationType => 2,    # 0: no stabilization, 1: ASGS, 2: OSGS
    :K => 1.0,                  # extra factor for artificial diffusion
    :ζ => 0.6,                  # constants ζ, xₐ in gap function, Eq. (2)
    :xₐ => 7 / 9 * pi,
    :u̅ => 0.98,                 # parameter in regularization, Eq. (3)
    :u => u,                    # analytical solution, NA for this example
    :ub => ub,                  # boundary condition, function of x
    :f => f,                    # body load as read from file 
    :u₀ => u₀,                  # initial guess, function of x
    :isNL => true,              # whether problem is nonlinear 
    :γ => 5                     # slenderness ratio (only used in linear problem)
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
display(plotResiduals(residuals))







