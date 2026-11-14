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

@testitem "Optimization: sequential spot merit gradient matches ForwardDiff" begin
    using DataFrames, ForwardDiff
    using OpticSim: AGFFileReader, pupilrays, spotmerit, LensOptimizationProblem, spotmerit_and_gradient
    glass = AGFFileReader.Examples_N_BK7
    air = AGFFileReader.Air
    # an aspheric (conic with 4th and 6th order terms) singlet, a stop in air, then a plano-convex lens with a flat back
    # x = front curvature, curvature of the second lens, front conic, front 4th order term, first thickness, back focus
    function presc(x::AbstractVector{T}) where {T}
        DataFrame(SurfaceType=["Object", "Aspheric", "Standard", "Stop", "Standard", "Standard", "Image"],
            Radius=T[Inf, 1/x[1], -80.0, Inf, 1/x[2], Inf, Inf],
            Conic=Union{Missing,T}[missing, x[3], zero(T), missing, zero(T), zero(T), missing],
            Parameters=Union{Missing,Vector{Pair{String,T}}}[missing, ["4" => x[4], "6" => T(1e-9)], Pair{String,T}[], missing, Pair{String,T}[], Pair{String,T}[], missing],
            Thickness=Union{Missing,T}[Inf, x[5], 5.0, 4.0, 3.0, x[6], missing],
            Material=[air, glass, air, air, glass, air, missing],
            SemiDiameter=T[Inf, 10.0, 10.0, 6.0, 8.0, 8.0, 15.0])
    end
    x0 = [1 / 40, 1 / 30, -0.5, -2e-6, 4.0, 40.0]
    bounds = [(0.0, 0.05), (0.0, 0.05), (-2.0, 1.0), (-1e-5, 1e-5), (2.0, 8.0), (20.0, 60.0)]
    build(x::AbstractVector{T}) where {T} = AxisymmetricOpticalSystem{T}(presc(x), 1, 1)
    # two fields, two wavelengths; the stop vignettes part of each beam
    rays = [pupilrays(θ, 14.0, 0.0, [0.5, 0.6], [1.0, 0.8]; rings=5) for θ in (0.0, 3.0)]
    merit(x::AbstractVector{T}) where {T} = spotmerit(presc(x), rays, T)
    lens = LensOptimizationProblem("test", "", string.(1:6), x0, bounds, build, merit, [0.0, 3.0], rays; prescription=presc)
    for x in (x0, x0 .+ [0.002, -0.003, 0.2, 1e-6, 0.5, -3.0])
        m, g = spotmerit_and_gradient(lens, x)
        @test m ≈ merit(x) rtol = 1e-10
        @test g ≈ ForwardDiff.gradient(merit, x) rtol = 1e-8
        @test first(spotmerit_and_gradient(lens, x; withgradient=false)) ≈ m rtol = 1e-12
    end

    # rays aimed at a pupil that moves with the variables: the gradient includes the rays' dependence on x
    aimed(x) = [pupilrays(θ, 14.0, x[5] / 2, [0.5, 0.6], [1.0, 0.8]; rings=5) for θ in (0.0, 3.0)]
    aimedmerit(x::AbstractVector{T}) where {T} = spotmerit(presc(x), aimed(x), T)
    aimedlens = LensOptimizationProblem("test", "", string.(1:6), x0, bounds, build, aimedmerit, [0.0, 3.0], aimed(x0); aimedrays=aimed, prescription=presc)
    m, g = spotmerit_and_gradient(aimedlens, x0)
    @test m ≈ aimedmerit(x0) rtol = 1e-10
    @test g ≈ ForwardDiff.gradient(aimedmerit, x0) rtol = 1e-8

    # a prescription whose spot merit isn't `merit` is rejected
    @test_throws AssertionError LensOptimizationProblem("test", "", string.(1:6), x0, bounds, build, x -> 2merit(x), [0.0, 3.0], rays; prescription=presc)
end
