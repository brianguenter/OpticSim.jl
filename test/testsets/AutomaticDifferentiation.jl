# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Optimization needs gradients of merit functions with respect to lens parameters. This checks that ForwardDiff can
# differentiate through building a system and tracing rays, and that the gradient is correct.
@testitem "ForwardDiff through trace" begin
    using ForwardDiff, StaticArrays, AGFFileReader
    using DataFrames: DataFrame
    using OpticSim.Emitters

    # double convex lens; free parameters: both radii, the lens thickness and the distance to the detector
    function lens(x::AbstractVector{T}) where {T}
        r1, r2, thickness, back = x
        AxisymmetricOpticalSystem{T}(
            DataFrame(
                SurfaceType=["Object", "Standard", "Standard", "Image"],
                Radius=[T(Inf), r1, r2, T(Inf)],
                Thickness=[T(Inf), thickness, back, missing],
                Material=[AGFFileReader.Air, AGFFileReader.Examples_N_BK7, AGFFileReader.Air, missing],
                SemiDiameter=[T(Inf), T(9.0), T(9.0), T(15.0)],
            ),
            1, 1,
        )
    end

    const rays = collect(Sources.Source(
        transform=OpticSim.Geometry.translation(0.0, 2.0, 10.0),
        origins=Origins.Hexapolar(3, 6.0, 6.0),
        directions=Directions.Constant(0.0, sind(3.0), -cosd(3.0)),
        spectrum=Spectrum.DeltaFunction(0.55),
    ))

    # squared RMS spot size on the detector
    function spotsize2(x::AbstractVector{T}) where {T}
        sys = lens(x)
        points = SVector{3,T}[]
        for r in rays
            rT = OpticalRay(SVector{3,T}(OpticSim.origin(r)), SVector{3,T}(OpticSim.direction(r)), T(OpticSim.power(r)), T(OpticSim.wavelength(r)))
            res = trace(sys, rT, test=true)
            res isa LensTrace && push!(points, point(res))
        end
        centroid = sum(points) / length(points)
        return sum(p -> sum(abs2, p - centroid), points) / length(points)
    end

    x0 = [60.0, -60.0, 10.0, 57.8]
    g = ForwardDiff.gradient(spotsize2, x0)
    fd = map(eachindex(x0)) do i
        h = 1e-6 * abs(x0[i])
        e = zeros(length(x0))
        e[i] = h
        (spotsize2(x0 .+ e) - spotsize2(x0 .- e)) / 2h
    end
    @test all(isfinite, g)
    @test isapprox(g, fd, rtol=1e-5)
end

# Aspheric surfaces are intersected with a triangle mesh followed by Newton iteration. Their gradients used to be NaN
# near the axis (sqrt(0) in the mesh sampling) and zero for aspheric terms equal to zero (zero terms were dropped).
@testitem "ForwardDiff through aspheric lens" begin
    using ForwardDiff, StaticArrays, AGFFileReader
    using DataFrames: DataFrame

    # plano-convex lens with an even aspheric front surface; free: conic, 4th and 6th order terms (sag at the 9mm edge)
    function lens(x::AbstractVector{T}) where {T}
        conic, a4, a6 = x
        AxisymmetricOpticalSystem{T}(
            DataFrame(
                SurfaceType=["Object", "Aspheric", "Standard", "Image"],
                Radius=[T(Inf), T(40.0), T(Inf), T(Inf)],
                Conic=[missing, conic, zero(T), missing],
                Parameters=[missing, ["4" => a4 / 9^4, "6" => a6 / 9^6], Pair{String,T}[], missing],
                Thickness=[T(Inf), T(6.0), T(75.0), missing],
                Material=[AGFFileReader.Air, AGFFileReader.Examples_N_BK7, AGFFileReader.Air, missing],
                SemiDiameter=[T(Inf), T(9.0), T(9.0), T(15.0)],
            ),
            1, 1,
        )
    end

    # rays from the axis out to 6mm, including the on-axis ray
    heights = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    function spotsize2(x::AbstractVector{T}) where {T}
        sys = lens(x)
        ys = T[]
        for h in heights
            res = trace(sys, OpticalRay(SVector{3,T}(0, h, 10), SVector{3,T}(0, 0, -1), one(T), T(0.55)), test=true)
            res isa LensTrace && push!(ys, point(res)[2])
        end
        return sum(abs2, ys) / length(ys)
    end

    for x0 in ([0.0, 0.0, 0.0], [-0.6, 0.01, 0.002])
        g = ForwardDiff.gradient(spotsize2, x0)
        fd = map(eachindex(x0)) do i
            h = 1e-6
            e = zeros(length(x0))
            e[i] = h
            (spotsize2(x0 .+ e) - spotsize2(x0 .- e)) / 2h
        end
        @test all(isfinite, g)
        @test all(!iszero, g) # zero-valued aspheric terms still have a derivative
        @test isapprox(g, fd, rtol=1e-4)
    end
end
