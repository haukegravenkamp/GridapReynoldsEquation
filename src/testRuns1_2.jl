## INPUT
# read manufactured solution and corresponding body load, 
# as well as the name of the folder for storing the results

println("precompiling typical function calls for faster start-up")
println("this may take a few minutes")


for iEx = 1:2
    for iartDiff in [false, true]
        for ibacktrack in [false, true]

                u1, f1, resultFolder1, ub1, u₀1 = selectExample(iEx)

                # problem definition
                paramProblem = Dict(
                    :nx => 2,                  # number of elements in x
                    :ny => 2,                   # number of elements in y
                    :hRefinement => 1,        # refinement steps starting from nx, ny
                    :order => 1,                # element order 
                    :artDiff => iartDiff,          # use artificial diffusion 
                    :stabilizationType => 2,    # 0: no stabilization, 1: ASGS, 2: OSGS
                    :K => 1.0,                  # extra factor for artificial diffusion
                    :ζ => 0.5,                  # constants ζ, xₐ in gap function, Eq. (2)
                    :xₐ => π,
                    :u̅ => 0.98,                 # parameter in regularization, Eq. (3)
                    :u => u1,                    # analytical solution 
                    :ub => ub1,                  # boundary condition, function of x
                    :f => f1,                    # body load as read from file 
                    :u₀ => u₀1,                  # initial guess, function of x
                    :isNL => true,              # whether problem is nonlinear 
                    :γ => 5                     # slenderness ratio (only used in linear problem)
                )

                # solver parameters
                paramSolver = Dict(
                    :numIt => 2,               # max number of iterations
                    :linType => 1,              # linearization type; 0: Picard, 1: Newton 
                    :initialPicard => 1,        # in case of Newton, do this many Picard iterations first, then switch
                    :backTracking => ibacktrack, # use backtracking (only in case of Newton)
                    :showTrace => false         # show solver output
                )

                # solve using Reynolds module
                el2, h, nEle, t, solverCache, residuals = convTest(paramProblem, paramSolver, resultFolder1)
                # plot residual of finest mesh
                pR = plotResiduals(residuals[end])
                # plot errors for all meshes

                # compute slope of convergence graph
                if length(h) > 1
                    sl = slope(h[end-1:end], el2[end-1:end])
                    # println("convergence rate: ", sl)
                end

                # uₕ, dΩ, hMin, h, solverCache, residuals = runReynolds(paramProblem, paramSolver, resultFolder1)
                # # compute error
                # el2 = computeError(u1, uₕ, dΩ)


        end
    end
end






