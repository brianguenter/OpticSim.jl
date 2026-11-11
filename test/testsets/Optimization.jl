# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testitem "Optimization: axial real ray trace" begin
    using DataFrames
    using OpticSim: AGFFileReader, axialrayheights, paraxialtrace
    glass = AGFFileReader.Examples_N_BK7
    air = AGFFileReader.Air
    λ = 0.55
    n = AGFFileReader.index(glass, λ)
    # a plano-convex lens, convex side first, then a flat stop 20 after the back
    lens = DataFrame(SurfaceType=["Object", "Standard", "Standard", "Stop", "Image"],
        Radius=[Inf, 50.0, Inf, Inf, Inf], Thickness=[Inf, 5.0, 20.0, 60.0, missing],
        Material=[air, glass, air, air, missing], SemiDiameter=[Inf, 15.0, 15.0, 15.0, 15.0])

    # by hand: a ray at height 10 meets the sphere at height 10, refracts by Snell's law there and again at the flat back
    h = 10.0
    z1 = 50 - sqrt(50^2 - h^2)       # where it meets the sphere, after the vertex
    i = asin(h / 50)                 # angle of incidence: the normal points at the center of curvature
    u = i - asin(sin(i) / n)         # the ray's angle to the axis inside the glass, toward the axis
    y2 = h - (5 - z1) * tan(u)       # at the flat back
    u2 = asin(n * sin(u))            # after refracting there
    y3 = y2 - 20 * tan(u2)           # at the stop
    @test axialrayheights(lens, h, λ) ≈ [h, y2, y3]

    # small heights agree with the paraxial trace
    y, _ = paraxialtrace(1e-4, 0.0, [1 / 50, 0.0], [5.0, 20.0], [n, 1.0], 1.0)
    @test axialrayheights(lens, 1e-4, λ)[3] ≈ y rtol = 1e-6

    # a ray beyond the sphere misses it
    @test all(isnan, axialrayheights(lens, 51.0, λ))
end

@testitem "Optimization: clear apertures that pass the axial beam" begin
    using DataFrames
    using OpticSim: AGFFileReader, axialrayheights, passaxialbeam
    glass = AGFFileReader.Examples_N_BK7
    air = AGFFileReader.Air
    lens(stop) = DataFrame(SurfaceType=["Object", "Standard", "Standard", "Stop", "Image"],
        Radius=[Inf, 50.0, Inf, Inf, Inf], Thickness=[Inf, 5.0, 20.0, 60.0, missing],
        Material=[air, glass, air, air, missing], SemiDiameter=[Inf, 5.0, 5.0, stop, 15.0])
    # a large stop: the surfaces grow to the beam's height at each, the stop keeps its size
    p = passaxialbeam(lens(15.0), 10.0, [0.55])
    heights = axialrayheights(lens(15.0), 10.0, 0.55)
    @test p.SemiDiameter[2:3] ≈ heights[1:2]
    @test p.SemiDiameter[4] == 15.0
    # a small stop limits the beam: the surfaces pass only the ray through the stop's edge
    p = passaxialbeam(lens(6.0), 10.0, [0.55])
    @test p.SemiDiameter[2] < 10.0
    @test axialrayheights(lens(6.0), p.SemiDiameter[2], 0.55)[3] ≈ 6.0 atol = 1e-6
    # semi-diameters never shrink
    @test passaxialbeam(lens(15.0), 1.0, [0.55]).SemiDiameter[2:3] == [5.0, 5.0]
end
