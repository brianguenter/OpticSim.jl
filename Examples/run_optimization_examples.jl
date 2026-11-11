# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Complete, runnable lens optimization examples: each creates an optimization problem from one of the
# OptimizationExamples in optimization_examples.jl, solves it with a bounded optimizer and prints the results.
#
# Run all of them from the repository root with:
#   julia --project=Examples Examples/run_optimization_examples.jl
#
# or run one interactively:
#   julia --project=Examples
#   julia> include("Examples/run_optimization_examples.jl")
#   julia> optimize_singlet()
#
# To also draw the designs before and after, load a Makie backend first (using GLMakie or using CairoMakie) and call
# drawresult(example, solution).

using OpticSim
using Optimization, OptimizationOptimJL # Optimization.jl and its wrapper for Optim.jl's bounded solvers
import ForwardDiff                      # gradients for Optimization.AutoForwardDiff()
using Printf

include(joinpath(@__DIR__, "Examples.jl"))

"""
    solveexample(example; solver = Optim.Fminbox(Optim.LBFGS()), maxiters = 300) -> (u, retcode, solution)

Create the optimization problem for `example`, with its starting values, bounds and merit function, and solve it.
The solver enforces the bounds. Gradients come from ForwardDiff. Returns the optimized variables `u`, the solver's
return code and its full solution.

The solver works with each variable scaled to run from 0 to 1 across its bounds. Lens variables have very different
scales: a curvature change of 0.001/mm affects the image about as much as a spacing change of millimeters, so without
scaling the gradient is thousands of times larger for curvatures, and gradient based solvers barely move spacings.
"""
function solveexample(example; solver=Optim.Fminbox(Optim.LBFGS()), maxiters=300)
    lower, upper = first.(example.bounds), last.(example.bounds)
    unscale(s) = lower .+ s .* (upper .- lower)
    # one ForwardDiff pass carrying all the derivatives is faster than the default of several passes
    merit = OptimizationFunction((s, _) -> example.merit(unscale(s)), Optimization.AutoForwardDiff(chunksize=length(example.x0)))
    s0 = (example.x0 .- lower) ./ (upper .- lower)
    problem = OptimizationProblem(merit, s0; lb=zeros(length(s0)), ub=ones(length(s0)))
    solution = solve(problem, solver; maxiters)
    return (u=unscale(solution.u), retcode=solution.retcode, solution=solution)
end

"""
    printresults(example, solution)

Print the variables before and after, the merit function, the RMS spot radius per field and the smallest separation
between each pair of neighboring lens surfaces.
"""
function printresults(example, solution)
    x0, x = example.x0, solution.u
    println("\n", example.name)
    println("  source: ", example.source)
    println("  solver status: ", solution.retcode)

    println("\n  variables                                         start         final")
    for (name, a, b, (lower, upper)) in zip(example.variables, x0, x, example.bounds)
        # solvers stop within a small fraction of the range from a bound they are pressing against
        tolerance = 1e-6 * (upper - lower)
        atbound = b - lower <= tolerance ? "  (at lower bound)" : upper - b <= tolerance ? "  (at upper bound)" : ""
        @printf("    %-44s %12.6g  %12.6g%s\n", name, a, b, atbound)
    end

    @printf("\n  merit (mean squared RMS spot radius, mm²): %.4e -> %.4e\n", example.merit(x0), example.merit(x))
    fields = join((replace(string(f), r"\.0$" => "") for f in example.fields), "°, ") * "°"
    println("  RMS spot radius (μm) for fields ", fields, ":")
    println("    start ", round.(Examples.rmsspotradii(example, x0), digits=1))
    println("    final ", round.(Examples.rmsspotradii(example, x), digits=1))

    println("\n  smallest separation between neighboring surfaces (negative = surfaces cross):")
    for s in Examples.surfaceseparations(example, x)
        @printf("    surfaces %d-%d (%-5s)  center %8.3f mm   smallest %8.3f mm at height %.2f mm%s\n",
            s.surfaces..., s.medium, s.center, s.minimum, s.height, s.minimum > 0 ? "" : "   <-- SURFACES CROSS")
    end
end

"""
    drawresult(example, solution)

Draw the design before and after optimization, one above the other. Needs a Makie backend: run `using GLMakie` or
`using CairoMakie` first.
"""
function drawresult(example, solution)
    rays = Examples.fieldrays(example)
    Vis.drawcomparison([example.build(example.x0), example.build(solution.u)]; raygenerator=rays,
        titles=["Before: " * example.name, "After"], trackallrays=true, test=true, colorbysourcenum=true)
end

#############################################################################
# One complete example per lens

"""Zemax's singlet tutorial: an F/4 N-BK7 singlet, focal length 100mm, fields 0, 3.5 and 5 degrees."""
function optimize_singlet()
    example = Examples.singlet_optimization()
    solution = solveexample(example)
    printresults(example, solution)
    return example, solution
end

"""The singlet with an even aspheric front surface: its conic and 4th and 6th order terms are free as well."""
function optimize_aspheric_singlet()
    example = Examples.aspheric_singlet_optimization()
    solution = solveexample(example)
    printresults(example, solution)
    return example, solution
end

"""Zemax's sample "Cooke 40 degree field": starts at Zemax's optimized design, so little should change."""
function optimize_cooke40()
    example = Examples.cooke40_optimization()
    solution = solveexample(example)
    printresults(example, solution)
    return example, solution
end

"""OpticSim's cooketriplet with its radii, air gaps and image distance free. Takes about a minute and a half."""
function optimize_cooke_triplet()
    example = Examples.cooke_triplet_optimization()
    solution = solveexample(example)
    printresults(example, solution)
    return example, solution
end

"""
Zeiss's f/1.4 Biotar (Merté, 1927), 5 fields and the F, d and C lines. First the patent's lens is refocused (only the
image distance free), then the curvatures, group spacings, stop position and glass thicknesses are optimized starting
from the refocused lens, and the two are compared. The optimization took 4.2 minutes with 8 threads on the performance
cores of a 24-core laptop processor, and about half an hour on one thread; start Julia with `julia -t auto` to use every
core.
"""
function optimize_biotar()
    patent = Examples.biotar_refocus_optimization()
    refocused = solveexample(patent)
    printresults(patent, refocused)
    example = Examples.biotar_optimization(image=refocused.u[1])
    solution = solveexample(example)
    printresults(example, solution)

    println("\nPatent lens refocused vs optimized:")
    for (name, ex, x) in (("patent, refocused", patent, refocused.u), ("optimized", example, solution.u))
        system = ex.build(x)
        rays = ex.aimedrays(x)
        # rays that pass through each lens surface once, in order, as the merit function counts them
        reached = sum(count(((r, _),) -> (lt = trace(system, r, test=true); lt !== nothing && OpticSim.nhits(lt) == Examples.lenssurfaces(system.prescription)), f) for f in rays)
        @printf("  %-18s merit %.4e mm²   rays reaching the image %d of %d\n", name, ex.merit(x), reached, sum(length, rays))
        println("      RMS spot radius (μm) for fields ", round.(ex.fields, digits=1), "°: ", round.(Examples.rmsspotradii(ex, x), digits=1))
    end
    return patent, refocused, example, solution
end

# run all of them when this file is run as a script
if abspath(PROGRAM_FILE) == @__FILE__
    optimize_singlet()
    optimize_aspheric_singlet()
    optimize_cooke40()
    optimize_cooke_triplet()
    optimize_biotar()
end
