# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Solving LensOptimizationProblems with Optimization.jl. Loaded automatically when Optimization, OptimizationOptimJL and
# ForwardDiff are all loaded.
module OpticSimOptimizationExt

import OpticSim
using Optimization: OptimizationFunction, OptimizationProblem, AutoForwardDiff, solve
using OptimizationOptimJL: Optim
import ForwardDiff

function OpticSim.optimizelens(lens::OpticSim.LensOptimizationProblem; solver=Optim.Fminbox(Optim.LBFGS()), maxiters=300)
    lower, upper = first.(lens.bounds), last.(lens.bounds)
    unscale(s) = lower .+ s .* (upper .- lower)
    # one ForwardDiff pass carrying all the derivatives is faster than the default of several passes
    merit = OptimizationFunction((s, _) -> lens.merit(unscale(s)), AutoForwardDiff(chunksize=length(lens.x0)))
    s0 = (lens.x0 .- lower) ./ (upper .- lower)
    problem = OptimizationProblem(merit, s0; lb=zeros(length(s0)), ub=ones(length(s0)))
    solution = solve(problem, solver; maxiters)
    return (u=unscale(solution.u), retcode=solution.retcode, solution=solution)
end

end # module OpticSimOptimizationExt
