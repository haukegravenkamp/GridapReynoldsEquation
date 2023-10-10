
## INPUT
uT, fT, resultFolderT, ubT, u₀T = selectExample(3)

for iartDiff in [false, true]

# problem definition
paramProblem = Dict(
    :nx => 2,                   # number of elements in x
    :ny => 2,                   # number of elements in y
    :hRefinement => 0,          # refinement steps starting from nx, ny
    :order => 1,                # element order 
    :artDiff => iartDiff,       # use artificial diffusion 
    :stabilizationType => 2,    # 0: no stabilization, 1: ASGS, 2: OSGS
    :K => 1.0,                  # extra factor for artificial diffusion
    :ζ => 0.6,                  # constants ζ, xₐ in gap function, Eq. (2)
    :xₐ => 7 / 9 * pi,
    :u̅ => 0.98,                 # parameter in regularization, Eq. (3)
    :u => uT,                   # analytical solution, NA for this example
    :ub => ubT,                 # boundary condition, function of x
    :f => fT,                   # body load as read from file 
    :u₀ => u₀T,                 # initial guess, function of x
    :isNL => true,              # whether problem is nonlinear 
    :γ => 5                     # slenderness ratio (only used in linear problem)
)

# solver parameters
paramSolver = Dict(
    :numIt => 500,              # max number of iterations
    :linType => 1,              # linearization type; 0: Picard, 1: Newton 
    :initialPicard => 3,        # in case of Newton, do this many Picard iterations first, then switch
    :backTracking => false,     # use backtracking (only in case of Newton)
    :showTrace => false         # show solver output
)

# solve problem using Reynolds module
uh, dΩ, hMin, h, solverCache, residuals = runReynolds(paramProblem, paramSolver, resultFolderT)

end
