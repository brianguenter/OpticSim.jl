# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Measures how long it takes to build a lens assembly with many elements and trace the first ray through it, including
# compilation, for the tuple-based LensAssembly and the old macro-generated MacroLensAssembly.
# Each measurement runs in a fresh Julia process so no assembly code has been compiled yet.
#
# Run from the repository root with:
#   julia --project=test test/Benchmarks/LensAssemblyCompileTime.jl

using Printf

const SIZES = (5, 20, 40, 80)
const KINDS = ("macro", "tuple")

function measure(kind, n)
    @eval using OpticSim, StaticArrays, AGFFileReader
    @eval using OpticSim: MacroLensAssembly
    Base.invokelatest(_measure, kind, n)
end

function _measure(kind, n)
    glass = AGFFileReader.Examples_N_BK7
    lenses = map(1:n) do i
        vertex = -6.0 * i
        if i % 3 == 0
            SphericalLens(glass, vertex, 30.0, -30.0, 3.0, 5.0)()
        elseif i % 3 == 1
            ConicLens(glass, vertex, 30.0, -0.5, -30.0, 0.0, 3.0, 5.0)()
        else
            SphericalLens(glass, vertex, Inf, -30.0, 3.0, 5.0)()
        end
    end
    detector = Rectangle(20.0, 20.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -6.0 * n - 20.0), interface=opaqueinterface())
    ray = OpticalRay(SVector(0.5, 0.3, 10.0), SVector(0.0, 0.0, -1.0), 1.0, 0.55)

    # warm up the per-lens code so only assembly-level code is compiled in the timings below
    for l in lenses[1:min(3, n)]
        surfaceintersection(l, OpticSim.ray(ray))
    end

    make = kind == "macro" ? MacroLensAssembly : LensAssembly
    tbuild = @elapsed assembly = make(lenses...)
    tsystem = @elapsed system = CSGOpticalSystem(assembly, detector)
    tfirst = @elapsed trace(system, ray, test=true)
    tsecond = @elapsed trace(system, ray, test=true)
    println(join((tbuild, tsystem, tfirst, tsecond), " "))
end

function driver()
    script = @__FILE__
    project = Base.active_project()
    @printf("%-6s %4s %10s %10s %12s %12s %10s\n", "kind", "n", "build s", "system s", "1st trace s", "2nd trace s", "total s")
    for n in SIZES, kind in KINDS
        out = read(`$(Base.julia_cmd()) --startup-file=no --project=$project $script $kind $n`, String)
        tbuild, tsystem, tfirst, tsecond = parse.(Float64, split(strip(last(split(strip(out), '\n')))))
        @printf("%-6s %4d %10.3f %10.3f %12.3f %12.6f %10.3f\n", kind, n, tbuild, tsystem, tfirst, tsecond, tbuild + tsystem + tfirst)
    end
end

if isempty(ARGS)
    driver()
else
    measure(ARGS[1], parse(Int, ARGS[2]))
end
