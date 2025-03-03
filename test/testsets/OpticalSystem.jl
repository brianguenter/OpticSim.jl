# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testitem "OpticalSystem" begin
    using DataFrames
    using Unitful
    using AGFFileReader

    function doubleconvex(
        ::Type{T}=Float64;
        temperature::Unitful.Temperature=AGFFileReader.TEMP_REF_UNITFUL,
        pressure::T=convert(T, PRESSURE_REF)
    ) where {T<:Real}
        AxisymmetricOpticalSystem{T}(
            DataFrame(
                SurfaceType=["Object", "Standard", "Standard", "Image"],
                Radius=[Inf64, 60, -60, Inf64],
                # OptimizeRadius = [false, true, true, false],
                Thickness=[Inf64, 10.0, 57.8, missing],
                # OptimizeThickness = [false, true, true, false],
                Material=[Air, Examples_N_BK7, Air, missing],
                SemiDiameter=[Inf64, 9.0, 9.0, 15.0]
            );
            temperature,
            pressure
        )
    end

    @testset "Single threaded trace makes sure function executes properly" begin
        conv = doubleconvex()
        rays = Emitters.Sources.CompositeSource(OpticSim.translation(-OpticSim.unitZ3()), repeat([Emitters.Sources.Source()], 10000))
        trace(conv, rays)
        @test true #just want to verify that the trace function executed properly
    end
end
