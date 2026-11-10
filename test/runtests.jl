# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Test
using TestItems
using TestItemRunner
using OpticSim

using Aqua

@info "Running Aqua tests"
@testset "Aqua" begin
    Aqua.test_all(OpticSim)
end

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
    using DataFrames
end

@testsnippet TEST_DATA begin
    include("TestData/TestData.jl")
end

@info "Running package tests"
@run_package_tests
