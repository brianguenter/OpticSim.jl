# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Test
using TestItems
using TestItemRunner

@testsnippet TestConstants begin
    const COMP_TOLERANCE = 25 * eps(Float64)
    const RTOLERANCE = 1e-10
    const ATOLERANCE = 10 * eps(Float64)
    const SEED = 12312487
end

@testsnippet Dependencies begin
    using Random
    using OpticSim.Emitters
    using OpticSim.Geometry
    using StaticArrays
    import StableRNGs
    using Colors
    using DataFrames
end

@testsnippet TEST_DATA begin
    include("TestData/TestData.jl")
end

"""Evaluate all functions not requiring arguments in a given module and test they don't throw anything"""
macro test_all_no_arg_functions(m)
    quote
        for n in names($m, all=true)
            # get all the valid functions
            if occursin("#", string(n)) || string(n) == string(nameof($m))
                continue
            end
            f = Core.eval($m, n)
            # get the methods of this function
            for meth in methods(f)
                # if this method has no args then try and evaluate it
                if meth.nargs == 1
                    # suppress STDOUT to prevent loads of stuff spamming the console
                    @test (@suppress_out begin
                        f()
                    end;
                    true)
                end
            end
        end
    end
end

@run_package_tests
