# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Compares the tuple-based LensAssembly (TupleLensAssembly) against the old macro-generated one (MacroLensAssembly)
# by tracing a ray through the same systems built with each.
#
# Run from the repository root with:
#   julia --project=test test/Benchmarks/LensAssemblyBenchmark.jl

using BenchmarkTools
using StaticArrays
using OpticSim
using OpticSim: MacroLensAssembly, TupleLensAssembly
using AGFFileReader
using Printf

include(joinpath(@__DIR__, "..", "TestData", "TestData.jl"))

csgsystem(sys::AxisymmetricOpticalSystem) = sys.system
csgsystem(sys::CSGOpticalSystem) = sys

# rebuild a system with the same elements in the given kind of lens assembly
function withassembly(makeassembly, sys)
    sys = csgsystem(sys)
    assembly = makeassembly(elements(sys.assembly)...; axis=OpticSim.axis(sys.assembly))
    return CSGOpticalSystem(assembly, sys.detector)
end

# a stack of n lenses of three different CSG types, to test assemblies with many elements of mixed types
function lensstack(n)
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
    return CSGOpticalSystem(LensAssembly(lenses...), detector)
end

const RAY = OpticalRay(SVector(0.5, 0.3, 10.0), SVector(0.0, 0.0, -1.0), 1.0, 0.55)

const SYSTEMS = [
    "double convex" => () -> TestData.doubleconvex(),
    "cooke triplet" => () -> TestData.cooketriplet(),
    "zernike lens" => () -> TestData.zernikesystem(),
    "aspheric lens" => () -> TestData.asphericsystem(),
    "chebyshev lens" => () -> TestData.chebyshevsystem(),
    "planar shapes" => () -> TestData.planarshapes(),
    "5 lens stack" => () -> lensstack(5),
    "20 lens stack" => () -> lensstack(20),
    "40 lens stack" => () -> lensstack(40),
    "80 lens stack" => () -> lensstack(80),
]

function run()
    @printf("%-16s %13s %13s %8s %13s %13s\n", "system", "macro min ns", "tuple min ns", "ratio", "macro allocs", "tuple allocs")
    for (name, makesystem) in SYSTEMS
        oldsys = withassembly(MacroLensAssembly, makesystem())
        newsys = withassembly(LensAssembly, makesystem())
        @assert newsys.assembly isa TupleLensAssembly

        # both versions must give the same answer
        oldres = trace(oldsys, RAY, test=true)
        newres = trace(newsys, RAY, test=true)
        @assert (oldres === nothing && newres === nothing) || point(oldres) == point(newres) "results differ for $name"

        old = @benchmark trace($oldsys, $RAY, test=true)
        new = @benchmark trace($newsys, $RAY, test=true)
        told = minimum(old).time
        tnew = minimum(new).time
        @printf("%-16s %13.0f %13.0f %8.2f %13d %13d\n", name, told, tnew, tnew / told, old.allocs, new.allocs)
    end
end

run()
