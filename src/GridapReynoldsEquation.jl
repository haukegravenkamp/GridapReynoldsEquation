module GridapReynoldsEquation
export domainDefinition, ReynoldsWrite, getFeSpaces, ReynoldsSolve, runReynolds, computeError, slope, convTest, plotResiduals, plotErrors
using Gridap
using LineSearches: BackTracking
using Plots
include("ReynoldsExampleLoads.jl")


# TODO: shock capturing
# TODO: defaults for optional arguments
# TODO: test higher-order elements

## DOMAIN
function domainDefinition(nx, ny, resultFolder)
  # minimum element edge length
  hMin = min(2 / ny, 2π / nx)
  hMax = max(2 / ny, 2π / nx)

  # domain and mesh definition
  domain = (0, 2π, -1, 1)
  partition = (nx, ny)
  model = CartesianDiscreteModel(domain, partition)

  # read number of elements (here just nx*ny)
  nCells = length(model.grid.cell_type)

  # write mesh for paraview
  filePathModel = resultFolder * "model"
  writevtk(model, filePathModel)

  return domain, model, partition, hMin, hMax, nCells
end

## FE SPACES
function getFeSpaces(model, order, ub)

  # space for u
  reffe = ReferenceFE(lagrangian, Float64, order)
  V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")

  # trial space with Dirichlet value for solution
  U = TrialFESpace(V, ub)
  # trial space without Dirichlet conditions for projections
  P = TrialFESpace(V)

  # combine to mixed space
  Y = MultiFieldFESpace([V, V])
  X = MultiFieldFESpace([U, P])

  # mesh 
  Ωₕ = Triangulation(model)
  # integration rule
  dΩ = Measure(Ωₕ, 2 * order)

  return U, V, X, Y, Ωₕ, dΩ
end


# complete definition of weak form and solver
function ReynoldsSolve(U, V, X, Y, paramProblem, dΩ, hMin, paramSolver)

  # read problem parameters
  ζ = paramProblem[:ζ]
  u̅ = paramProblem[:u̅]
  xₐ = paramProblem[:xₐ]
  K = paramProblem[:K]
  isNL = paramProblem[:isNL]
  artDiff = paramProblem[:artDiff]
  stabilizationType = paramProblem[:stabilizationType]
  f = paramProblem[:f]
  u₀ = paramProblem[:u₀]
  γ = paramProblem[:γ]

  # read solver parameters
  if isNL
    numIt = paramSolver[:numIt]
    linType = paramSolver[:linType]
    initialPicard = paramSolver[:initialPicard]
    backTracking = paramSolver[:backTracking]
    # in case of Picard, we make sure backTracking is disabled
    if linType == 0
      backTracking = false
    end
    # set a local flag for Picard iteration
    # used when switching automatically from Picard to Newton
    if linType == 0
      Picard = true
    else
      Picard = false
    end
  end

  ## FUNCTIONS DEFINING PHYSICAL PROBLEM 

  # gap function Eq. (2) and derivative wrt x
  H(x) = 1.0 - ζ * cos(x[1] - xₐ)
  Hₓ(x) = ζ * sin(x[1] - xₐ)

  # regularization, Eq. (3) and dg/du, d^2g/du^2
  argg(u) = u / (1 - u̅)
  g(u) = 1 / π * atan(argg(u)) + 0.5
  dg(u) = 1 / π * 1 / (argg(u)^2 + 1) / (1 - u̅)
  ddg(u) = -2 / π * argg(u) / (argg(u)^2 + 1)^2 * 1 / (1 - u̅)^2

  # right-hand side
  f̂(x) = f(x) - Hₓ(x)

  # shorthand for x- and y-derivative 
  # (as only the gradient of an FE function is defined in Gridap, I think)
  dx(∇u) = [1, 0] ⋅ ∇u
  dy(∇u) = [0, 1] ⋅ ∇u

  # helper function phi = g*u, Eq. (7)
  ϕ(u) = g(u) * u
  dϕ(u) = dg(u) * u + g(u)
  ddϕ(u) = ddg(u) * u + 2 * dg(u)
  # helper function psi = (g-1)u, Eq. (34)
  ψ(u) = (g(u) - 1) * u
  dψ(u) = dg(u) * u + (g(u) - 1)
  ddψ(u) = ddg(u) * u + 2 * dg(u)
  # helper function theta = aₓu, not in paper, simplifies artificial diffusion term
  dϑ(u, H) = dax(u, H) * u + ax(u, H)
  ddϑ(u, H) = ddax(u, H) * u + 2 * dax(u, H)

  # convection in x-direction, above Eq. (9)
  ax(u, H) = (g(u) - 1) ⋅ H
  dax(u, H) = dg(u) ⋅ H
  ddax(u, H) = ddg(u) ⋅ H
  # diffusivity, above Eq. (9)
  k(u, H) = H^3 ⋅ dϕ(u) / 12.0
  dk(u, H) = H^3 ⋅ ddϕ(u) / 12.0
  # reaction, above Eq. (9)
  s(u, ∇u, H, Hₓ) = dg(u) ⋅ dx(∇u) ⋅ H + (g(u) - 1) ⋅ Hₓ

  # stabilization parameter, Eq. (22)
  τ(u, ∇u, H, Hₓ) = (4 / hMin^2 * abs(k(u, H)) + 2 / hMin * abs(ax(u, H)) + abs(s(u, ∇u, H, Hₓ)))^(-1)

  ## FUNCTIONS DEFINING WEAK FORM
  # preceeding 'd' indicated their contribution to the Jacobian as in dC
  # subscript ₚ indicates Picard iteration, otherwise it's Newton-Raphson

  # convection, Eq. (9)
  C(u, v) = ∫(-(v ⋅ (dx ∘ ∇(u))) ⋅ (ax ∘ (u, H)))dΩ
  dC(u, du, v) = ∫(-(v ⋅ (dx ∘ ∇(du))) ⋅ H ⋅ (dψ ∘ u))dΩ
  dCₚ(u, du, v) = ∫(-(v ⋅ (dx ∘ ∇(du))) ⋅ (ax ∘ (u, H)))dΩ

  # diffusion, Eq. (9)
  D(u, v) = ∫(∇(v) ⋅ ∇(u) ⋅ (k ∘ (u, H)))dΩ
  dD(u, du, v) = ∫(∇(v) ⋅ ∇(u) ⋅ (dk ∘ (u, H)) ⋅ du)dΩ + ∫(∇(v) ⋅ ∇(du) ⋅ (k ∘ (u, H)))dΩ
  dDₚ(u, du, v) = ∫(∇(v) ⋅ ∇(du) ⋅ (k ∘ (u, H)))dΩ

  # reaction, Eq. (9)
  R(u, v) = ∫(-(v ⋅ u ⋅ (s ∘ (u, ∇(u), H, Hₓ))))dΩ
  dR(u, du, v) = ∫(-(v ⋅ du ⋅ H ⋅ (ddψ ∘ u) ⋅ (dx ∘ ∇(u))))dΩ + ∫(-(v ⋅ du ⋅ Hₓ ⋅ (dψ ∘ u)))dΩ
  dRₚ(u, du, v) = ∫(-(v ⋅ du ⋅ (s ∘ (u, ∇(u), H, Hₓ))))dΩ

  # artificial diffusion, Eq. (13)
  Da(u, v) = ∫(-((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(u))) ⋅ (hMin / 2 * K) ⋅ (dϑ ∘ (u, H)))dΩ
  dDa(u, du, v) = ∫(-((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(du))) ⋅ (hMin / 2 * K) ⋅ (dϑ ∘ (u, H)))dΩ +
                  ∫(-((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(u))) ⋅ (hMin / 2 * K) ⋅ (ddϑ ∘ (u, H)) ⋅ du)dΩ
  dDaₚ(u, du, v) = ∫(-((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(du))) ⋅ (hMin / 2 * K) ⋅ (dϑ ∘ (u, H)))dΩ

  # stabilization, Galerkin term, Eq. (21)
  S(u, v) = ∫(((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(u))) ⋅ (ax ∘ (u, H)) ⋅ (ax ∘ (u, H)) ⋅ (τ ∘ (u, ∇(u), H, Hₓ)))dΩ
  dSₚ(u, du, v) = ∫(((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(du))) ⋅ (ax ∘ (u, H)) ⋅ (ax ∘ (u, H)) ⋅ (τ ∘ (u, ∇(u), H, Hₓ)))dΩ

  # stabilization, contribution of projections, Eq. (29a)
  Pτ((u, ξ), (v, η)) = ∫(((dx ∘ ∇(v)) ⋅ ξ) ⋅ (ax ∘ (u, H)) ⋅ (τ ∘ (u, ∇(u), H, Hₓ)))dΩ
  dPτₚ((u, ξ), (du, dξ), (v, η)) = ∫(((dx ∘ ∇(v)) ⋅ dξ) ⋅ (ax ∘ (u, H)) ⋅ (τ ∘ (u, ∇(u), H, Hₓ)))dΩ

  # stabilization, compute projections, Eq. (29b) 
  P((u, ξ), (v, η)) = ∫(((dx ∘ ∇(u)) ⋅ η) ⋅ (ax ∘ (u, H)))dΩ
  dPₚ((u, ξ), (du, dξ), (v, η)) = ∫(((dx ∘ ∇(du)) ⋅ η) ⋅ (ax ∘ (u, H)))dΩ

  # mass matrix, Eq. (29b) 
  M((u, ξ), (v, η)) = ∫(η ⋅ ξ)dΩ
  dMₚ((u, ξ), (du, dξ), (v, η)) = ∫(η ⋅ dξ)dΩ

  # linear version
  # diffusion
  kₗ(H) = H^3
  Dₗ(u, v) = ∫((dx ∘ ∇(v)) ⋅ (dx ∘ ∇(u)) ⋅ (kₗ ∘ (H)) + (dy ∘ ∇(v)) ⋅ (dy ∘ ∇(u)) ⋅ (kₗ ∘ (H)) / 4 / (γ^2))dΩ

  # left-hand side, Galerkin terms
  function B(u, v)
    # combine C,D,R terms 
    if isNL
      Bt = C(u, v) + D(u, v) + R(u, v)
    else
      Bt = Dₗ(u, v)
    end
    # add artificial diffusion if requested
    if artDiff
      Bt = Bt + Da(u, v)
    end
    # add Galerkin term of stabilization
    if (stabilizationType in [1, 2])
      Bt = Bt + S(u, v)
    end
    return Bt
  end

  # left-hand side, projection terms
  B_OSGS((u, ξ), (v, η)) = P((u, ξ), (v, η)) - Pτ((u, ξ), (v, η)) - M((u, ξ), (v, η))

  # Jacobian
  function jac(u, du, v)
    # combine variation of C,D,R terms
    if isNL
      if Picard
        J = dCₚ(u, du, v) + dDₚ(u, du, v) + dRₚ(u, du, v)
      else # Newton
        J = dC(u, du, v) + dD(u, du, v) + dR(u, du, v)
      end
    else # linear case
      return 0
    end
    # add artificial diffusion
    if artDiff
      if Picard
        J = J + dDaₚ(u, du, v)
      else
        J = J + dDa(u, du, v)
      end
    end
    # add Galerkin term of stabilization
    if (stabilizationType in [1, 2])
      J = J + dSₚ(u, du, v)
    end
    return J
  end

  # Jacobian, projection terms
  jacOSS((u, ξ), (du, dξ), (v, η)) = jac(u, du, v) + dPₚ((u, ξ), (du, dξ), (v, η)) - dPτₚ((u, ξ), (du, dξ), (v, η)) - dMₚ((u, ξ), (du, dξ), (v, η))

  # right-hand side
  L(v) = ∫(v * f̂)dΩ

  # residual without projection
  res(u, v) = B(u, v) - L(v)

  # residual with projection
  resOSS((u, ξ), (v, η)) = B(u, v) + B_OSGS((u, ξ), (v, η)) - L(v)

  # read the residuals from solverCache, just for plotting
  function getResiduals(solverCache)

    # storage for residuals 
    residualNorms = Float64[]

    # number of results, including initial residual 
    nResults = solverCache.result.iterations + 1
    for i in 1:nResults
      resTemp = solverCache.result.trace[i].fnorm
      push!(residualNorms, resTemp)
    end
    return residualNorms, nResults
  end

  # define solver
  function NLsolverPhase(nIterations, useBacktracking)
    if useBacktracking
      solver = NLSolver(show_trace=true, store_trace=true, method=:newton, linesearch=BackTracking(), iterations=nIterations)
    else
      solver = NLSolver(show_trace=true, store_trace=true, method=:newton, iterations=nIterations)
    end
    # solution depending on stabilization 
    if stabilizationType == 2 # orthogonal subgrid scales (mixed elements with projection DOFs)
      op = FEOperator(resOSS, jacOSS, X, Y)
      uξₕ, solverCache = solve!(uξₕ, solver, op)
      uₕ = uξₕ[1]
    else # all other cases, only u as unknown
      op = FEOperator(res, jac, U, V)
      uₕ, solverCache = solve!(uₕ, solver, op)
    end
    # read residuals
    residualNorms, nCalls = getResiduals(solverCache)
    return uₕ, solverCache, residualNorms, nCalls
  end


  # do solution
  if isNL # nonlinear

    # initial guess
    # projections always initialized as 0
    ξ₀(x) = 0.0
    # combined FE-functions for u and projection
    uξₕ = interpolate_everywhere([u₀, ξ₀], X)
    # interpolate only u if projections are not used
    uₕ = interpolate_everywhere(u₀, U)

    #    Newton           initial Picard steps
    if (linType == 1) && (initialPicard > 0)
      # temporarily switch to Picard 
      Picard = true
      # solve initial Picard iterations without backtracking 
      println("\n" * "fixed-point iteration:")
      uₕ, solverCache1, residualNorms1 = NLsolverPhase(initialPicard, false)
      # switch to Newton 
      Picard = false
    else # pure Newton or pure Picard
      initialPicard = 0
      residualNorms1 = Float64[]
    end
    # do (remaining, in case of initialPicard) iterations 
    Picard && println("\n" * "fixed-point iteration:")
    ~Picard && println("\n" * "Newton's method:")
    uₕ, solverCache, residualNorms = NLsolverPhase(numIt - initialPicard, backTracking)
    residualNorms = [residualNorms1[1:end-1]; residualNorms]

  else # linear for debugging

    # linear operator
    op = AffineFEOperator(B, L, U, V)
    uₕ = solve(op)
    solverCache = []
    residualNorms = []

  end

  return uₕ, solverCache, residualNorms
end

# write to file
function ReynoldsWrite(uₕ, Ωₕ, resultFolder)
  writevtk(Ωₕ, resultFolder * "results", cellfields=["uₕ" => uₕ])
end

# perform all computations in order
function runReynolds(paramProblem, paramSolver, resultFolder)
  # this is the main file to call, unless an automatic convergence test is done
  nx = paramProblem[:nx]
  ny = paramProblem[:ny]
  order = paramProblem[:order]
  ub = paramProblem[:ub]
  tCPU = @elapsed begin
    domain, model, partition, hMin, hMax, nCells = domainDefinition(nx, ny, resultFolder)
    U, V, X, Y, Ωₕ, dΩ = getFeSpaces(model, order, ub)
    uₕ, solverCache, residualNorms = ReynoldsSolve(U, V, X, Y, paramProblem, dΩ, hMin, paramSolver)
  end
  println("finished run with $(nCells) elements in $(round(tCPU;digits = 2)) s")
  ReynoldsWrite(uₕ, Ωₕ, resultFolder)
  return uₕ, dΩ, hMin, hMax, solverCache, residualNorms, tCPU, nCells

end

# compute L2 error 
function computeError(u, uₕ, dΩ)
  e = u - uₕ
  return sqrt(sum(∫(e * e) * dΩ))
end

# slope of logarithmic data for computing convergence rate
function slope(hs, errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

# perform convergence test
function convTest(paramProblem, paramSolver, resultFolder)

  # initialize containers to store results
  el2s = Float64[]  # L2 errors
  hs = Float64[]    # element size 
  nEle = Float64[]  # number of elements
  ts = Float64[]    # cpu times 
  solverCache = []  # solver output data
  residuals = []    # residuals in nonlinear solver

  # read number of elements and use them for base mesh
  nx0 = paramProblem[:nx]
  ny0 = paramProblem[:ny]
  # refinement steps 
  ns = paramProblem[:hRefinement]
  # analytical solution
  u = paramProblem[:u]

  for n in ns # loop over refinement steps 
    # element numbers 
    paramProblem[:ny] = ny0 * 2^n
    paramProblem[:nx] = nx0 * 2^n

    # call main file and time complete computation using @elapsed macro
    uₕ, dΩ, hMin, hMax, sC, residualNorms, tCPU, nCells = runReynolds(paramProblem, paramSolver, resultFolder)

    # compute L2 norm
    el2 = computeError(u, uₕ, dΩ)

    # append results of current refinement step
    push!(el2s, el2)
    push!(hs, hMax)
    push!(nEle, nCells)
    push!(ts, tCPU)
    push!(solverCache, sC)
    push!(residuals, residualNorms)

  end

  return el2s, hs, nEle, ts, solverCache, residuals

end

# simple function for plotting errors of convergence study 
function plotErrors(h, errors)

  p = plot(h, errors, xaxis=:log, yaxis=:log, label="", xlabel="h", ylabel="error norm", linecolor=:black, marker=(:circle, 4, :white, stroke(1.5, :black)), fontfamily="Computer Modern", linewidth=1, framestyle=:box, grid=false )
  return p

end

# plotting residuals of result k if there are several
function plotResiduals(residuals)

  nCalls = length(residuals)
  p = plot(0:(nCalls-1), residuals, yaxis=:log, label="", xlabel="iteration", ylabel="residual norm", linecolor=:black, marker=(:circle, 4, :white, stroke(1.5, :black)),fontfamily="Computer Modern", linewidth=1, framestyle=:box, grid=false )
  return p
end

end
