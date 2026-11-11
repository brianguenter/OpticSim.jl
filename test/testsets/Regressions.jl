# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# One test per bug found by the BugFinder review, so the bugs can't come back unnoticed.

@testsnippet RegressionSystems begin
    using DataFrames, StaticArrays, LinearAlgebra
    using OpticSim: AGFFileReader
    using OpticSim.Geometry
    using OpticSim.Emitters: Sources, Origins, Directions, Spectrum
    singlet() = AxisymmetricOpticalSystem(DataFrame(SurfaceType=["Object", "Standard", "Standard", "Image"],
        Radius=[Inf, -35.571, 35.571, Inf], Thickness=[Inf, 4.0, 44.748, missing],
        Material=[AGFFileReader.Air, AGFFileReader.Examples_N_SK16, AGFFileReader.Air, missing],
        SemiDiameter=[Inf, 7.054, 6.033, 15.0]))
    beam(n) = Sources.Source(transform=translation(0.0, 0.0, 10.0), origins=Origins.RectGrid(10.0, 10.0, n, n),
        directions=Directions.Constant(0.0, 0.0, -1.0), spectrum=Spectrum.DeltaFunction(0.55)) # one wavelength, so tracing is deterministic
end

@testitem "Regression: normal of an Interval" begin
    using StaticArrays
    using OpticSim: AGFFileReader
    glass = FresnelInterface{Float64}(AGFFileReader.Examples_N_BK7, AGFFileReader.Air)
    iv = surfaceintersection(Sphere(2.0, interface=glass), Ray([0.3, 0.2, -5.0], [0.0, 0.0, 1.0]))
    # normal(::Interval) returned the point instead of the normal
    @test OpticSim.normal(iv) == OpticSim.normal(closestintersection(iv))
    @test !(OpticSim.normal(iv) ≈ OpticSim.point(iv))
end

@testitem "Regression: CSG intersection of objects with disjoint bounding boxes" begin
    using StaticArrays
    using OpticSim.Geometry
    # separated along z only, which hit a typo in the bounding box overlap test
    a = leaf(Sphere(1.0), translation(0.0, 0.0, 5.0))
    b = leaf(Sphere(1.0), translation(0.0, 0.0, -5.0))
    obj = (a ∩ b)()
    @test surfaceintersection(obj, Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval
    @test !OpticSim.inside(obj, SVector(0.0, 0.0, 5.0))
    # separated along each axis
    for offset in ([3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0])
        x, y, z = offset
        @test (leaf(Sphere(1.0)) ∩ leaf(Sphere(1.0), translation(x, y, z)))() isa OpticSim.CSGTree
    end
end

@testitem "Regression: Spectrum.Measured" begin
    using DataFrames
    using OpticSim.Emitters.Spectrum: Measured, spectrumpower
    # unsorted wavelengths, and a first sample that isn't a multiple of the spacing
    s = Measured(DataFrame(Wavelength=[505, 405, 605], Power=[1.0, 0.5, 0.0]))
    @test s.low_wave_length == 405 && s.high_wave_length == 605
    @test spectrumpower(s, 455.0) ≈ 0.75
    @test spectrumpower(s, 555.0) ≈ 0.5
    @test spectrumpower(s, 605.0) ≈ 0.0
    @test spectrumpower(s, 700.0) === nothing
end

@testitem "Regression: copy and outpath of an optical system" setup = [RegressionSystems] begin
    sys = singlet()
    c = copy(sys)
    @test c isa AxisymmetricOpticalSystem && c !== sys
    @test size(detectorimage(c)) == size(detectorimage(sys))
    # saving the image needs a Makie backend, which the tests don't load: the error says so
    err = try
        trace(sys, beam(3); printprog=false, test=true, outpath=tempname() * ".png")
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("Makie", err.msg)
end

@testitem "Regression: traceMT and tracehitsMT match trace" setup = [RegressionSystems] begin
    single = trace(singlet(), beam(40); printprog=false, test=true)
    @test sum(single) > 0
    @test traceMT(singlet(), beam(40); printprog=false, test=true) == single
    hits = tracehitsMT(singlet(), beam(40); printprog=false, test=true)
    @test length(hits) == sum(single)
end

@testitem "Regression: eye models" begin
    using OpticSim.Geometry
    eye = OpticSim.Data.ArizonaEye()
    @test "SurfaceType" in names(eye)
    @test OpticSim.Data.ModelEye(OpticSim.LensAssembly(leaf(Sphere(1.0), translation(0.0, 0.0, 20.0))())) isa CSGOpticalSystem
end

@testitem "Regression: BSpline insertknots" begin
    using StaticArrays
    c = OpticSim.BSplineCurve{OpticSim.Euclidean,Float64,2,3}(OpticSim.KnotVector{Float64}([0.0, 0, 0, 0, 0.5, 1, 1, 1, 1]),
        [MVector(0.0, 0), MVector(1.0, 1), MVector(2.0, 0), MVector(3.0, 1), MVector(4.0, 0)])
    @test OpticSim.insertknots(c) isa OpticSim.BSplineCurve
end

@testitem "Regression: lens constructors" begin
    using StaticArrays
    using OpticSim: AGFFileReader
    glass = AGFFileReader.Examples_N_BK7
    # invalid input throws a clear error instead of logging and failing later
    @test_throws ArgumentError FresnelLens(glass, 0.0, 15.0, 1.0, 5.0, 0.8, aspherics=[4 => 0.001])
    # a hyperbolic surface can have a vertex radius smaller than the semidiameter
    @test ConicLens(glass, 0.0, 20.0, 0.0, -5.22, -1.35, 8.0, 5.3)() isa OpticSim.CSGTree
    # A decenter along one axis only was ignored when sizing a concave cap, so the cap didn't reach the far edge of the
    # aperture. The concave front's center of curvature is at (1.5, 0, 10).
    decentered = SphericalLens(glass, 0.0, -10.0, Inf, 6.0, 8.0; frontdecenter=(1.5, 0.0))()
    hit = closestintersection(surfaceintersection(decentered, Ray([-7.8, 0.0, 20.0], [0.0, 0.0, -1.0])))
    @test hit !== nothing && OpticSim.point(hit) ≈ [-7.8, 0.0, 10 - sqrt(100 - 9.3^2)]
    OpticSim.emptyintervalpool!()
end

@testitem "Regression: concave lens edge ray" begin
    using StaticArrays, LinearAlgebra
    using OpticSim: AGFFileReader
    # a steep ray entering the concave front bowl near the edge used to be reported as leaving the glass where it enters
    lens = SphericalLens(AGFFileReader.Examples_N_BK7, 0.0, -9.583029481514219, -12.115030518293151, 7.287463084135346, 7.287802713203115)()
    o = SVector(-8.172311986363153, -15.671583086896094, 8.374828628765762)
    r = Ray(o, normalize(SVector(4.817636397144307, 1.5617391170833512, -1.2596945078581279) - o))
    iv = surfaceintersection(lens, r)
    for t in range(0.0, 40.0, length=400)
        inintervals = iv isa EmptyInterval ? false : any(OpticSim.α(OpticSim.lower(i)) < t < (OpticSim.upper(i) isa Infinity ? Inf : OpticSim.α(OpticSim.upper(i))) for i in (iv isa Interval ? (iv,) : OpticSim.intervals(iv)))
        @test inintervals == OpticSim.inside(lens, OpticSim.point(r, t))
    end
    OpticSim.emptyintervalpool!()
end

@testitem "Regression: aspheric surfaces use their bounding cylinder" begin
    using OpticSim.Geometry
    s = AcceleratedParametricSurface(AsphericSurface(5.0, radius=20.0, conic=-0.5))
    m = which(surfaceintersection, (typeof(s), Ray{Float64,3}))
    @test occursin("AsphericSurface", string(m.file))
end

@testitem "Regression: Plane visualization axes" begin
    using StaticArrays
    p = Plane(SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 0.0), visvec=SVector(0.0, 2.0, 0.0), vishalfsizeu=1.0, vishalfsizev=1.0)
    @test !any(isnan, p.visuvec) && !any(isnan, p.visvvec)
end

@testitem "Regression: paraxial lens rotationvec" begin
    using StaticArrays
    a = ParaxialLensRect(10.0, 2.0, 1.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    b = ParaxialLensRect(10.0, 2.0, 1.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0); rotationvec=SVector(1.0, 0.0, 0.0))
    @test OpticSim.shape(a).uvec != OpticSim.shape(b).uvec
end

@testitem "Regression: thin grating" begin
    using StaticArrays
    using OpticSim: AGFFileReader, OpticalRay, processintersection, transmission, reflectance
    air = AGFFileReader.Air
    # the defaults couldn't be constructed, and transmission and reflectance were swapped
    g = ThinGratingInterface(SVector(1.0, 0.0, 0.0), 2.0, air, air)
    @test transmission(g, 0) ≈ 1 / 3 && reflectance(g, 0) == 0
    r = OpticalRay(SVector(0.0, 0.0, 5.0), SVector(0.0, 0.0, -1.0), 1.0, 0.55)
    n = SVector(0.0, 0.0, 1.0)
    p = SVector(0.0, 0.0, 0.0)
    # only the +1 order, transmitted: continues forward at the grating equation angle
    t = ThinGratingInterface(SVector(1.0, 0.0, 0.0), 2.0, air, air; transmission=[0.0, 0.0, 1.0])
    out = [processintersection(t, p, n, r, 20.0, 1.0, true) for _ in 1:50]
    dirs = [o[1] for o in out if o !== nothing]
    @test !isempty(dirs) && all(d -> d[3] < 0 && d[1] ≈ 0.55 / 2, dirs)
    # only the +1 order, reflected: comes back
    f = ThinGratingInterface(SVector(1.0, 0.0, 0.0), 2.0, air, air; reflectance=[0.0, 0.0, 1.0], transmission=[0.0, 0.0, 0.0])
    out = [processintersection(f, p, n, r, 20.0, 1.0, true) for _ in 1:50]
    dirs = [o[1] for o in out if o !== nothing]
    @test !isempty(dirs) && all(d -> d[3] > 0 && d[1] ≈ 0.55 / 2, dirs)
end

@testitem "Regression: ForwardDiff through the pole of a concave surface" begin
    using StaticArrays, ForwardDiff
    using OpticSim.Geometry
    # An on-axis ray hits a concave spherical surface (a SphericalCap) at its pole, where its polar uv coordinates are
    # singular. The normal used to be computed from them, and when the ray's derivative was a tiny rounding residue,
    # as in real traces, the normal's derivative came out NaN.
    D = ForwardDiff.Dual{:t,Float64,1}
    cap = SphericalCap(D(28.3), D(1.0), SVector{3,D}(0, 0, 1), SVector{3,D}(0, 0, 0))
    r = Ray(SVector{3,D}(D(0.0, ForwardDiff.Partials((1e-17,))), D(0.0), D(-10.0)), SVector{3,D}(0, 0, 1))
    hit = closestintersection(surfaceintersection(cap, r), false)
    @test OpticSim.normal(hit) ≈ [0.0, 0.0, 1.0]
    @test all(x -> all(isfinite, ForwardDiff.partials(x)), OpticSim.normal(hit))
    OpticSim.emptyintervalpool!(D)
end

@testitem "Regression: absorbing surfaces don't reflect" begin
    using StaticArrays
    using OpticSim: AGFFileReader, OpticalRay, processintersection
    # A ray in glass meets an absorbing surface (a blackened lens edge) beyond the critical angle. Total internal
    # reflection used to override the surface's zero reflectance, so stray rays bounced off lens edges.
    glass = AGFFileReader.Examples_N_BK7
    edge = FresnelInterface{Float64}(glass, AGFFileReader.Air, reflectance=0.0, transmission=0.0)
    n = SVector(1.0, 0.0, 0.0)                      # outward normal of the edge
    d = SVector(sind(80), 0.0, -cosd(80))            # 80° incidence from inside the glass: beyond the critical angle
    r = OpticalRay(SVector(-1.0, 0.0, 0.0), d, 1.0, 0.55)
    @test processintersection(edge, SVector(0.0, 0.0, 0.0), n, r, 20.0, 1.0, true) === nothing
    # the same ray at an ordinary glass-air surface is totally internally reflected
    surface = FresnelInterface{Float64}(glass, AGFFileReader.Air)
    @test processintersection(surface, SVector(0.0, 0.0, 0.0), n, r, 20.0, 1.0, true) !== nothing
end

@testitem "Regression: the ring filling a lens beyond its sphere absorbs" begin
    using StaticArrays
    using OpticSim: AGFFileReader
    # A spherical surface of radius 5 can't reach a semidiameter of 6, so a flat ring fills the gap to the rim. It's
    # outside the optical aperture and absorbs, like the barrel; it used to refract like the optical surface.
    lens = SphericalLens(AGFFileReader.Examples_N_BK7, 0.0, 5.0, Inf, 8.0, 6.0)()
    hit = closestintersection(surfaceintersection(lens, Ray([5.5, 0.0, 10.0], [0.0, 0.0, -1.0])), false)
    i = interface(hit)
    @test OpticSim.reflectance(i) == 0 && OpticSim.transmission(i) == 0
    OpticSim.emptyintervalpool!()
end

@testitem "Regression: a surface smaller than its lens ends in a flat step" begin
    using StaticArrays
    using OpticSim: AGFFileReader
    # The back surface's semidiameter, 4, is smaller than the lens's, 6. It ends at 4, with a flat absorbing ring from there
    # to the rim, as on a real lens; it used to continue its sphere out to the rim. Back radius 10 (concave): its sag at
    # height 4 is 10 - √84, so the ring is that far behind the back vertex at z = -8.
    lens = SphericalLens(AGFFileReader.Examples_N_BK7, 0.0, 30.0, 10.0, 8.0, 6.0; backsemidiameter=4.0)()
    ringz = -8.0 - (10.0 - sqrt(84.0))
    # a ray coming back from behind the lens at height 5 hits the ring, not the sphere (which would be at z ≈ -9.34)
    hit = closestintersection(surfaceintersection(lens, Ray([5.0, 0.0, -20.0], [0.0, 0.0, 1.0])), false)
    @test point(hit)[3] ≈ ringz
    i = interface(hit)
    @test OpticSim.reflectance(i) == 0 && OpticSim.transmission(i) == 0
    # inside the back surface's semidiameter the ray meets the sphere, an ordinary optical surface
    hit = closestintersection(surfaceintersection(lens, Ray([3.0, 0.0, -20.0], [0.0, 0.0, 1.0])), false)
    @test point(hit)[3] ≈ -8.0 - (10.0 - sqrt(91.0))
    @test OpticSim.transmission(interface(hit)) == 1
    OpticSim.emptyintervalpool!()
end

@testitem "Regression: optimization geometry check includes element rims and the stop" begin
    using DataFrames
    using OpticSim: AGFFileReader, surfacescross
    glass = AGFFileReader.Examples_N_BK7
    air = AGFFileReader.Air
    # one lens, then the stop: the lens's back surface (semidiameter 4) is smaller than its front (6). Its back surface is
    # concave, so the flat ring at its edge is 10 - √84 ≈ 0.83 behind the back vertex. A stop 2 behind the vertex clears
    # it; one 0.5 behind is cut by the rim of the back surface and the ring, outside the stop's hole.
    lens(gap) = DataFrame(SurfaceType=["Object", "Standard", "Standard", "Stop", "Image"],
        Radius=[Inf, 30.0, 10.0, Inf, Inf], Thickness=[Inf, 8.0, gap, 50.0, missing],
        Material=[air, glass, air, air, missing], SemiDiameter=[Inf, 6.0, 4.0, 3.0, 20.0])
    @test !surfacescross(lens(2.0))
    @test surfacescross(lens(0.5))
end
