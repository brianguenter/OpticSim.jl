# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Tests of the single-source sag kernels (src/Geometry/SagCurvature.jl): point, partials and normal of the sag surfaces
# come from one `sag` function per surface type, and surfacecurvature from its FastDifferentiation derivatives.

@testsnippet CurvatureHelpers begin
    using LinearAlgebra
    using StaticArrays
    using Random

    "Sag surfaces covering every kernel type, each with a sampler of normalized (u, v) parameters."
    function sagsurfaces()
        polar(rng) = (0.95 * rand(rng), 2π * rand(rng) - π)
        rect(rng) = (1.9 * rand(rng) - 0.95, 1.9 * rand(rng) - 0.95)
        return [
            ("even asphere", TestData.evenasphericsurface(), polar),
            ("odd asphere", TestData.oddasphericsurface(), polar),
            ("oddeven asphere", TestData.oddevenasphericsurface(), polar),
            ("conic", AsphericSurface(9.0; radius=40.0, conic=-1.2), polar),
            ("zernike 1", TestData.zernikesurface1(), polar),
            ("zernike 3", TestData.zernikesurface3(), polar),
            ("qtype", TestData.qtypesurface1(), polar),
            ("chebyshev 1", TestData.chebyshevsurface1(), rect),
            ("chebyshev 2", TestData.chebyshevsurface2(), rect),
            ("gridsag", TestData.gridsagsurfacebicubic(), polar),
        ]
    end
end

@testitem "surfacecurvature of sag surfaces matches finite differences of the normal" setup = [TestConstants, TEST_DATA, CurvatureHelpers] begin
    # Along a unit tangent dⱼ the normal changes by ∂n/∂sⱼ, and C[i,j] = -dᵢ⋅∂n/∂sⱼ (Weingarten). The normal comes from the
    # gradient kernel, independent of the Hessian kernel used by surfacecurvature.
    rng = Xoshiro(SEED)
    for (name, s, sampler) in sagsurfaces()
        for _ in 1:20
            u, v = sampler(rng)
            p = point(s, u, v)
            k = surfacecurvature(s, p)
            @test k.normal ≈ normal(s, u, v) atol = 1e-12
            @test abs(dot(k.d1, k.normal)) < 1e-12 && abs(dot(k.d2, k.normal)) < 1e-12
            @test cross(k.d1, k.d2) ≈ k.normal atol = 1e-12
            ε = 1e-6
            nat(q) = OpticSim.sagnormal(s, q[1], q[2])
            for (j, dj) in enumerate((k.d1, k.d2))
                dn = (nat(p + ε * dj) - nat(p - ε * dj)) / 2ε
                for (i, di) in enumerate((k.d1, k.d2))
                    @test k.C[i, j] ≈ -dot(di, dn) atol = 1e-6 * max(1, abs(k.C[i, j]))
                end
            end
        end
    end
end

@testitem "partials of sag surfaces match finite differences of point" setup = [TestConstants, TEST_DATA, CurvatureHelpers] begin
    rng = Xoshiro(SEED)
    h = 1e-6
    for (name, s, sampler) in sagsurfaces()
        for _ in 1:20
            u, v = sampler(rng)
            pu, pv = partials(s, u, v)
            @test pu ≈ (point(s, u + h, v) - point(s, u - h, v)) / 2h rtol = 1e-6 atol = 1e-7
            @test pv ≈ (point(s, u, v + h) - point(s, u, v - h)) / 2h rtol = 1e-6 atol = 1e-7
        end
    end
end

@testitem "regressions fixed by the sag kernels" setup = [TestConstants, TEST_DATA] begin
    using LinearAlgebra
    h = 1e-6
    # odd aspheres with an r¹ term: the radial derivative used to count the first term twice
    s = AsphericSurface(9.0; radius=-25.0, conic=0.3, aspherics=[(1, 1e-3), (3, 1e-4), (5, -3e-7)])
    @test partials(s, 0.4, 0.7)[1][3] ≈ (point(s, 0.4 + h, 0.7)[3] - point(s, 0.4 - h, 0.7)[3]) / 2h rtol = 1e-7
    @test partials(s, 0.0, 0.7)[1][3] ≈ 9.0 * 1e-3 # one sided radial slope at the axis: semidiameter × a₁
    # Zernike surfaces used to give NaN partials at ρ = 0 and ρ = 1, and the normal (0, 0, 1) at the center even with tilt terms
    z = TestData.zernikesurface1()
    for (ρ, ϕ) in ((0.0, 1.2), (1.0, 0.7), (1.0, -2.3))
        @test partials(z, ρ, ϕ)[1][3] ≈ (point(z, ρ + h, ϕ)[3] - point(z, ρ - h, ϕ)[3]) / 2h rtol = 1e-6
        @test all(isfinite, normal(z, ρ, ϕ))
    end
    rad = OpticSim.semidiameter(z)
    zx = (point(z, h, 0.0)[3] - point(z, h, Float64(π))[3]) / (2h * rad)
    zy = (point(z, h, π / 2)[3] - point(z, h, -π / 2)[3]) / (2h * rad)
    @test normal(z, 0.0, 0.0) ≈ normalize([-zx, -zy, 1.0]) atol = 1e-8
end

@testitem "bilinear GridSag interpolation" setup = [TestConstants] begin
    using StaticArrays
    using Random
    rng = Xoshiro(SEED)
    base = ChebyshevSurface(2.0, 3.0, [(1, 1, 0.01)]; radius=50.0)
    grid = zeros(5, 4, 4) # Nv × Nu × [z, ∂z/∂x, ∂z/∂y, ∂²z/∂x∂y]; bilinear interpolation uses only z
    grid[:, :, 1] = 0.01 * rand(rng, 5, 4)
    s = GridSagSurface(base, grid; interpolation=GridSagLinear)
    Nu, Nv = 4, 5
    for _ in 1:50
        u, v = 1.8 * rand(rng) - 0.9, 1.8 * rand(rng) - 0.9
        # bilinear interpolation of the grid, computed independently
        û, v̂ = (u + 1) / 2 * (Nu - 1), (v + 1) / 2 * (Nv - 1)
        i, j = min(floor(Int, û), Nu - 2), min(floor(Int, v̂), Nv - 2)
        κ, γ = û - i, v̂ - j
        g(a, b) = grid[b+1, a+1, 1]
        expected = point(base, u, v)[3] + κ * (γ * g(i + 1, j + 1) + (1 - γ) * g(i + 1, j)) + (1 - κ) * (γ * g(i, j + 1) + (1 - γ) * g(i, j))
        @test point(s, u, v)[3] ≈ expected rtol = 1e-12
    end
end

@testitem "surfacecurvature of implicit and planar surfaces" setup = [TestConstants] begin
    using LinearAlgebra
    using StaticArrays
    p = normalize(SVector(0.3, -0.5, 0.8)) * 7.0
    k = surfacecurvature(Sphere(7.0), p)
    @test k.C ≈ -I(2) / 7 atol = 1e-14 # convex, outward normal
    @test k.normal ≈ normalize(p)
    k = surfacecurvature(Cylinder(3.0), SVector(3.0, 0.0, 1.0))
    @test sort(eigvals(Matrix(k.C))) ≈ [-1 / 3, 0.0] atol = 1e-14
    k = surfacecurvature(Plane(0.0, 0.0, 1.0, 0.0, 0.0, 2.0), SVector(1.0, 1.0, 2.0))
    @test k.C ≈ zeros(2, 2) atol = 1e-15
    @test k.normal ≈ [0.0, 0.0, 1.0]
    cap = SphericalCap(5.0, π / 3)
    k = surfacecurvature(cap, point(cap, 0.4, 0.5))
    @test k.C ≈ I(2) / 5 atol = 1e-14 # the cap's normal points toward its center
    @test k.normal ≈ normal(cap, 0.4, 0.5)
    @test surfacecurvature(AsphericSurface(9.0; radius=-25.0), SVector(0.0, 0.0, 0.0)).C ≈ -I(2) / 25 atol = 1e-15
end

@testitem "surface ids identify the CSG leaf that was hit" setup = [TestConstants, TEST_DATA] begin
    using StaticArrays
    using LinearAlgebra
    sys = TestData.asphericsystem()
    r = OpticalRay(SVector(0.0, 0.0, 10.0), SVector(0.001, 0.001, -1.0), 1.0, 0.55)
    track = Vector{OpticSim.LensTrace{Float64,3}}(undef, 0)
    trace(sys, r, trackrays=track)
    allleaves = reduce(vcat, [leaves(e) for e in OpticSim.elements(sys.assembly) if e isa CSGTree])
    byid = Dict(l.id => l for l in allleaves)
    @test length(byid) == length(allleaves) # ids are unique
    lenshits = [OpticSim.intersection(t) for t in track if surfaceid(OpticSim.intersection(t)) != 0]
    @test length(lenshits) >= 2
    for hit in lenshits
        leaf = byid[surfaceid(hit)]
        k = surfacecurvature(leaf, point(hit))
        # the leaf's own normal at the hit agrees with the intersection normal up to orientation
        @test abs(dot(k.normal, normal(hit))) ≈ 1 atol = 1e-9
    end
end
