# load packages and modules
using Revise                    # helps tracking changes during development
using GridapReynoldsEquation    # use Reynolds module 
using Gridap                    # the FE packages Gridap           
using Plots                     # for plotting results

## INPUT

# example number as in paper, Section 8.1 or 8.2
# (example 8.3 is in a different file)
exampleNo = 1

# read manufactured solution and corresponding body load, 
# as well as the name of the folder for storing the results
u, f, resultFolder, ub, u₀ = selectExample(exampleNo)

# whether to run convergence test or only one mesh
runConvergence = true

# problem definition
paramProblem = Dict(
    :nx => 12,                  # number of elements in x
    :ny => 4,                   # number of elements in y
    :hRefinement => 0:4,        # refinement steps starting from nx, ny
    :order => 1,                # element order 
    :artDiff => false,          # use artificial diffusion 
    :stabilizationType => 2,    # 0: no stabilization, 1: ASGS, 2: OSGS
    :K => 1.0,                  # extra factor for artificial diffusion
    :ζ => 0.5,                  # constants ζ, xₐ in gap function, Eq. (2)
    :xₐ => π,
    :u̅ => 0.98,                 # parameter in regularization, Eq. (3)
    :u => u,                    # analytical solution 
    :ub => ub,                  # boundary condition, function of x
    :f => f,                    # body load as read from file 
    :u₀ => u₀,                  # initial guess, function of x
    :isNL => true,              # whether problem is nonlinear 
    :γ => 5                     # slenderness ratio (only used in linear problem)
)

# solver parameters
paramSolver = Dict(
    :numIt => 50,               # max number of iterations
    :linType => 1,              # linearization type; 0: Picard, 1: Newton 
    :initialPicard => 4,        # in case of Newton, do this many Picard iterations first, then switch
    :backTracking => true,      # use backtracking (only in case of Newton)
    :showTrace => true          # show solver output
)

# solve using Reynolds module
if runConvergence               # perform convergence test
    el2, h, nEle, t, solverCache, residuals = convTest(paramProblem, paramSolver, resultFolder)
    # plot residual of finest mesh
    pR = plotResiduals(residuals[end])
    # plot errors for all meshes
    pE = plotErrors(h, el2)
    display(plot(pR,pE,layout=(2,1),size=(600,600/1.618*2)))
    # compute slope of convergence graph
    if length(h) > 1
        sl = slope(h[end-1:end], el2[end-1:end])
        println("convergence rate: ", sl)
    end

else                            # only compute base mesh
    uₕ, dΩ, hMin, h, solverCache, residuals = runReynolds(paramProblem, paramSolver, resultFolder)
    # plot residual 
    pR = plotResiduals(residuals)
    display(plot(pR, size=(600,600/1.618)))
    # compute error
    el2 = computeError(u, uₕ, dΩ)
    println("L2 error: ", el2)
end







