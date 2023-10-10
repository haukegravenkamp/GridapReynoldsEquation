module GridapReynoldsEquation

using PrecompileTools: @setup_workload, @compile_workload

export domainDefinition, ReynoldsWrite, getFeSpaces, ReynoldsSolve, runReynolds, computeError, slope, convTest, plotResiduals, plotErrors
using Gridap
using LineSearches: BackTracking
using Plots

# manufactured solutions and loads for examples in paper
include("ReynoldsExampleLoads.jl")

# all other functions defining the Reynolds problem
include("Reynolds_functions.jl")

# precompilation for faster startup
# (increases run time at first execution)
@compile_workload begin
    # test runs using examples 1 and 2
    include("testRuns1_2.jl")
end
@compile_workload begin
    # test runs using example 3
    include("testRuns3.jl")
end

end

# TODO: shock capturing
# TODO: defaults for optional arguments
# TODO: test higher-order elements
