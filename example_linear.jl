# load packages and modules
using Revise                    # helps tracking changes during development
using GridapReynoldsEquation    # use Reynolds module 
using Gridap                    # the FE packages Gridap           
using Plots                     # for plotting results


# simple example of the linear Reynolds equation, see Eq. (3) in
# Scaled boundary finite element method for hydrodynamic bearings in rotordynamic simulations
# Pfeil et al., International Journal of Mechanical Sciences 2021
#
# 
# This example should - except for slightly different boundary conditions - correspond to the one shown in Fig. 9a of the paper mentioned above. However, the linear version was only created for debugging purposes, and I haven't thoroughly verified this example.

## INPUT
# where to store results for plotting in paraview
resultFolder = "./results/example_linear/"

# problem definition
paramProblem = Dict(
    :nx => 100,                 # number of elements in x
    :ny => 32,                  # number of elements in x
    :hRefinement => 0:4,        # refinement steps starting from nx, ny
    :order => 1,                # element order 
    :artDiff => ~true,          # use artificial diffusion 
    :stabilizationType => 0,    # 0: no stabilization, 1: ASGS, 2: OSGS
    :K => 1.0,                  # extra factor for artificial diffusion
    :ζ => 0.5,                  # constants ζ, xₐ in gap function, Eq. (2)
    :xₐ => π,
    :u̅ => 0.98,                 # parameter in regularization, Eq. (3)
    :u => x -> 0.0,             # analytical solution, NA for this example 
    :ub => x -> 0.0,            # boundary condition, function of x
    :f => x -> 0.0,             # body load
    :u₀ => x -> 1.0,            # initial guess, function of x
    :isNL => false,             # whether problem is nonlinear 
    :γ => 5                     # slenderness ratio (only used in linear problem)
)

# solver parameters
# paramSolver = Dict(
#     :numIt => 50,               # max number of iterations
#     :linType => 1,              # linearization type; 0: Picard, 1: Newton 
#     :initialPicard => 4,        # in case of Newton, do this many Picard iterations first, then switch
#     :backTracking => true       # use backtracking (only in case of Newton)
# )

# solve problem using Reynolds module
uh, dΩ, hMin, h, solverCache, residuals = runReynolds(paramProblem, 0, resultFolder);
