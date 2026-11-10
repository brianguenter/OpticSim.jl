# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Property based and randomized tests. The Supposition.jl properties shrink any failure to a minimal counterexample.
# The randomized geometry checks use fixed seeds, so a failure always reproduces: the message includes the seed.

@testmodule FuzzHelpers begin
    using OpticSim, OpticSim.Geometry, StaticArrays, LinearAlgebra, Random
    using OpticSim: lower, upper, intervals, α, inside, evalcsg, emptyintervalpool!

    randunit(rng) = normalize(randn(rng, SVector{3,Float64}))
    randtransform(rng) = Transform(rotmat(2π * rand(rng), 2π * rand(rng), 2π * rand(rng)), SVector{3,Float64}(2 .* randn(rng, 3)))

    function randprimitive(rng)
        k = rand(rng, 1:5)
        if k == 1
            return leaf(Sphere(0.5 + 2rand(rng)), randtransform(rng))
        elseif k == 2
            return leaf(Cylinder(0.5 + 2rand(rng)), randtransform(rng))
        elseif k == 3
            return leaf(Plane(randunit(rng), SVector{3,Float64}(randn(rng, 3))))
        elseif k == 4
            return OpticSim.transform(Cuboid(0.5 + rand(rng), 0.5 + rand(rng), 0.5 + rand(rng)), randtransform(rng))
        else
            return OpticSim.transform(HexagonalPrism(0.5 + rand(rng), 0.5 + rand(rng)), randtransform(rng))
        end
    end

    """A random CSG tree of spheres, cylinders, planes, cuboids and hexagonal prisms, combined with ∪, ∩ and -."""
    function randcsg(rng, depth)
        depth == 0 && return randprimitive(rng)
        a, b = randcsg(rng, rand(rng, 0:depth-1)), randcsg(rng, rand(rng, 0:depth-1))
        op = rand(rng, 1:3)
        return op == 1 ? a ∪ b : op == 2 ? a ∩ b : a - b
    end

    """The (α start, α end) of each interval in a ray intersection result."""
    function ranges(res)
        res isa EmptyInterval && return Tuple{Float64,Float64}[]
        ivs = res isa Interval ? [res] : collect(intervals(res))
        return map(ivs) do iv
            l, u = lower(iv), upper(iv)
            (l isa RayOrigin ? 0.0 : α(l), u isa Infinity ? Inf : α(u))
        end
    end

    """
        consistent(obj, ray, αmax; samples, tolerance) -> Union{Nothing,String}

    Checks that the intervals the ray intersection reports are sorted and disjoint, and that a point along the ray is
    inside `obj` exactly when it is within one of them. Returns nothing, or a description of the first disagreement.
    """
    function consistent(obj, ray, αmax; samples=300, tolerance=1e-5)
        rs = ranges(evalcsg(obj, ray))
        emptyintervalpool!()
        for i in 2:length(rs)
            rs[i-1][2] <= rs[i][1] || return "intervals overlap or are out of order: $rs"
        end
        bounds = [x for (a, b) in rs for x in (a, b)]
        for t in range(0.0, αmax, length=samples)
            any(x -> abs(t - x) < tolerance, bounds) && continue
            expected = any(a < t < b for (a, b) in rs)
            actual = inside(obj, point(ray, t))
            expected == actual || return "at α = $t the intervals say inside = $expected but inside() says $actual; intervals $rs"
        end
        return nothing
    end

    """Runs `check(seed)` for each seed and returns the first failures, as "seed => message" strings."""
    function failures(check, seeds; keep=5)
        found = String[]
        for seed in seeds
            message = try
                check(seed)
            catch e
                "threw " * sprint(showerror, e)
            end
            message === nothing || push!(found, "seed $seed: " * first(message, 400))
            length(found) >= keep && break
        end
        return found
    end
end

@testitem "Fuzz: CSG intervals agree with inside()" setup = [FuzzHelpers] begin
    using .FuzzHelpers: randcsg, randunit, consistent, failures
    using Random, StaticArrays
    # random CSG trees and rays: where the ray intersection says a point is inside the object, inside() must agree
    check(seed) = begin
        rng = Xoshiro(seed)
        obj = randcsg(rng, rand(rng, 1:3))()
        consistent(obj, Ray(SVector{3,Float64}(4 .* randn(rng, 3)), randunit(rng)), 20.0)
    end
    @test failures(check, 1:1500) == String[]
end

@testitem "Fuzz: lens intervals agree with inside()" setup = [FuzzHelpers] begin
    using .FuzzHelpers: consistent, failures
    using Random, StaticArrays, LinearAlgebra
    using OpticSim: AGFFileReader
    # random spherical, conic and decentered lenses, hit by rays from random directions, including through the edges
    function randlens(rng)
        sd = 2 + 8rand(rng)
        radius() = rand(rng) < 0.15 ? Inf : rand(rng, (-1, 1)) * sd * (1.05 + 4rand(rng))
        r1, r2 = radius(), radius()
        sag(r) = isinf(r) ? 0.0 : abs(r) - sqrt(r^2 - sd^2)
        t = sag(r1) + sag(r2) + 0.5 + 5rand(rng) # thick enough that the two surfaces can't cross
        glass = AGFFileReader.Examples_N_BK7
        k = rand(rng, 1:3)
        if k == 1
            return SphericalLens(glass, 0.0, r1, r2, t, sd)()
        elseif k == 2
            c1, c2 = isinf(r1) ? 0.0 : -rand(rng), isinf(r2) ? 0.0 : -rand(rng)
            return ConicLens(glass, 0.0, r1, c1, r2, c2, t, sd)()
        else
            d = (0.3randn(rng), 0.3randn(rng))
            return SphericalLens(glass, 0.0, r1, r2, t, sd; frontdecenter=d, backdecenter=d)()
        end
    end
    check(seed) = begin
        rng = Xoshiro(seed)
        lens = randlens(rng)
        o = SVector(8randn(rng), 8randn(rng), 20 + 5randn(rng))
        target = SVector(6randn(rng), 6randn(rng), -5rand(rng))
        consistent(lens, Ray(o, normalize(target - o)), 60.0)
    end
    @test failures(check, 1:600) == String[]
end

@testitem "Fuzz: accelerated surfaces agree with inside()" setup = [FuzzHelpers] begin
    using .FuzzHelpers: consistent, failures
    using Random, StaticArrays, LinearAlgebra
    # aspheric, Zernike, Q-type and Chebyshev surfaces, clipped to a finite solid the way a lens clips them
    check(seed) = begin
        rng = Xoshiro(seed)
        sd = 2 + 3rand(rng)
        R = rand(rng, (-1, 1)) * sd * (1.5 + 3rand(rng))
        k = rand(rng, 1:4)
        s = if k == 1
            AsphericSurface(sd, radius=R, conic=-rand(rng), aspherics=[(4, 1e-4randn(rng)), (6, 1e-6randn(rng))])
        elseif k == 2
            ZernikeSurface(sd, radius=R, conic=-rand(rng), zcoeff=[(1, 0.05randn(rng)), (4, 0.05randn(rng)), (7, 0.02randn(rng))])
        elseif k == 3
            QTypeSurface(sd, radius=R, conic=-rand(rng), αcoeffs=[(1, 0, 0.05randn(rng)), (2, 1, 0.02randn(rng))])
        else
            ChebyshevSurface(sd, sd, [(1, 0, 0.05randn(rng)), (0, 2, 0.05randn(rng))], radius=R)
        end
        side = k == 4 ? Cuboid(0.95sd, 0.95sd, 50.0) : leaf(Cylinder(0.95sd, 100.0))
        obj = (leaf(AcceleratedParametricSurface(s, 25)) ∩ side ∩ leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -8.0)))()
        o = SVector(6randn(rng), 6randn(rng), 10 + 3randn(rng))
        target = SVector(sd * randn(rng), sd * randn(rng), -4rand(rng))
        consistent(obj, Ray(o, normalize(target - o)), 40.0)
    end
    @test failures(check, 1:300) == String[]
end

@testitem "Fuzz: tracing random lens systems" begin
    using Random, StaticArrays, LinearAlgebra, DataFrames
    using OpticSim: AGFFileReader, OpticalRay, power, pathlength, intersection
    glasses = (AGFFileReader.Examples_N_BK7, AGFFileReader.Examples_N_SK16, AGFFileReader.Examples_N_SF2)
    air = AGFFileReader.Air
    function randsystem(rng)
        types = ["Object"]
        radii = [Inf]
        thicknesses = Union{Float64,Missing}[Inf]
        materials = Union{AGFFileReader.AbstractGlass,Missing}[air]
        semidiameters = [Inf]
        conics = Union{Float64,Missing}[missing]
        for _ in 1:rand(rng, 1:3)
            sd = 4 + 6rand(rng)
            r1 = rand(rng) < 0.2 ? Inf : rand(rng, (-1, 1)) * sd * (1.3 + 4rand(rng))
            r2 = rand(rng) < 0.2 ? Inf : rand(rng, (-1, 1)) * sd * (1.3 + 4rand(rng))
            sag(r) = isinf(r) ? 0.0 : abs(r) - sqrt(r^2 - sd^2)
            append!(types, ["Standard", "Standard"])
            append!(radii, [r1, r2])
            append!(thicknesses, [sag(r1) + sag(r2) + 1 + 4rand(rng), sag(r1) + sag(r2) + 1 + 10rand(rng)])
            append!(materials, [rand(rng, glasses), air])
            append!(semidiameters, [sd, sd])
            append!(conics, [rand(rng) < 0.3 ? -rand(rng) : 0.0, 0.0])
        end
        push!(types, "Image"); push!(radii, Inf); push!(thicknesses, missing); push!(materials, missing)
        push!(semidiameters, 30.0); push!(conics, missing)
        return AxisymmetricOpticalSystem(DataFrame(SurfaceType=types, Radius=radii, Thickness=thicknesses,
            Material=materials, SemiDiameter=semidiameters, Conic=conics))
    end

    function traceall()
        problems = String[]
        hits = 0
        for seed in 1:40
            rng = Xoshiro(seed)
            system = randsystem(rng)
            for _ in 1:150
                o = SVector(15randn(rng), 15randn(rng), 20 + 10rand(rng))
                d = normalize(SVector(0.6randn(rng), 0.6randn(rng), -1.0))
                lt = try
                    trace(system, OpticalRay(o, d, 1.0, 0.45 + 0.3rand(rng)); test=rand(rng, Bool))
                catch e
                    push!(problems, "seed $seed: trace threw " * first(sprint(showerror, e), 200))
                    continue
                end
                lt === nothing && continue
                hits += 1
                p = power(lt)
                isfinite(p) && 0 <= p <= 1 + 1e-12 || push!(problems, "seed $seed: power $p")
                all(isfinite, OpticSim.point(intersection(lt))) || push!(problems, "seed $seed: non-finite detector point")
                pathlength(lt) >= 0 || push!(problems, "seed $seed: negative path length")
            end
        end
        return problems, hits
    end
    problems, hits = traceall()
    @test problems == String[]
    @test hits > 100 # make sure the test really traces rays all the way to the detector
end

@testitem "Properties: transforms" begin
    using Supposition, StaticArrays, LinearAlgebra
    using OpticSim.Geometry
    using OpticSim.Geometry: matrix, rotmatbetween
    angle = Data.Floats{Float64}(; minimum=-10.0, maximum=10.0, nans=false, infs=false)
    coordinate = Data.Floats{Float64}(; minimum=-1.0, maximum=1.0, nans=false, infs=false)
    vec3 = map(v -> SVector{3,Float64}(v), Data.Vectors(coordinate; min_size=3, max_size=3))
    nonzero(v) = norm(v) > 1e-3
    # axis directions are the special cases most likely to go wrong
    axis = Data.SampledFrom([SVector(1.0, 0, 0), SVector(-1.0, 0, 0), SVector(0.0, 1, 0), SVector(0.0, -1, 0), SVector(0.0, 0, 1), SVector(0.0, 0, -1)])
    direction = Data.OneOf(filter(nonzero, vec3), axis)
    rot3(t) = matrix(t)[SOneTo(3), SOneTo(3)]
    db = Supposition.NoRecordDB()

    @check db = db max_examples = 300 function rotations_match_rotmat(a=angle)
        rot3(rotationX(a)) ≈ rotmat(a, 0.0, 0.0) && rot3(rotationY(a)) ≈ rotmat(0.0, a, 0.0) && rot3(rotationZ(a)) ≈ rotmat(0.0, 0.0, a)
    end
    @check db = db max_examples = 300 function rotmat_is_a_rotation(a=angle, b=angle, c=angle)
        R = rotmat(a, b, c)
        R' * R ≈ I && det(R) ≈ 1
    end
    @check db = db max_examples = 300 function frame_is_right_handed(f=direction)
        R = rot3(Transform(SVector(0.0, 0.0, 0.0), f))
        R' * R ≈ I && det(R) ≈ 1 && R[:, 3] ≈ normalize(f)
    end
    @check db = db max_examples = 300 function rotmatbetween_rotates_a_to_b(a=direction, b=direction, flip=Data.Booleans())
        b = flip ? -a : b # antiparallel vectors are a special case
        R = rotmatbetween(a, b)
        isapprox(R * normalize(a), normalize(b), atol=1e-9) && R' * R ≈ I && det(R) ≈ 1
    end
    @check db = db max_examples = 300 function inverse_undoes_transform(a=angle, b=angle, c=angle, t=vec3, p=vec3)
        tr = Transform(rotmat(a, b, c), 10t)
        isapprox(inv(tr) * (tr * p), p, atol=1e-9)
    end
    @test translation(Float32, 1.0, 2.0, 3.0) isa Transform{Float32}
end

@testitem "Properties: refraction and Fresnel" begin
    using Supposition, StaticArrays, LinearAlgebra
    using OpticSim: snell, fresnel, refractedray, reflectedray
    coordinate = Data.Floats{Float64}(; minimum=-1.0, maximum=1.0, nans=false, infs=false)
    unitvec = map(v -> normalize(SVector{3,Float64}(v)), filter(v -> norm(v) > 1e-3, Data.Vectors(coordinate; min_size=3, max_size=3)))
    index = Data.Floats{Float64}(; minimum=1.0, maximum=2.5, nans=false, infs=false)
    # an exponent for how close to grazing incidence the ray is
    grazing = Data.Integers(0, 15)
    db = Supposition.NoRecordDB()

    function incidence(n, t, g)
        tangent = cross(n, t)
        norm(tangent) < 1e-6 && return t
        return normalize(normalize(tangent) + (g == 0 ? 1.0 : 10.0^-g) * n)
    end

    @check db = db max_examples = 500 function fresnel_conserves_energy(n=unitvec, t=unitvec, g=grazing, n1=index, n2=index)
        d = incidence(n, t, g)
        si, st = snell(n, d, n1, n2)
        R, T = fresnel(n1, n2, si, st)
        0 <= si <= 1 && isfinite(R) && isfinite(T) && 0 <= R <= 1 && 0 <= T <= 1 && isapprox(R + T, 1, atol=1e-9)
    end
    @check db = db max_examples = 500 function refraction_obeys_snells_law(n=unitvec, t=unitvec, g=grazing, n1=index, n2=index)
        d = incidence(n, t, g)
        si, _ = snell(n, d, n1, n2)
        out = refractedray(n1, n2, n, d)
        if out === nothing
            n1 * si >= n2 - 1e-9 # total internal reflection only when it should happen
        else
            isapprox(norm(out), 1, atol=1e-9) &&
                isapprox(n1 * si, n2 * norm(cross(n, out)), atol=1e-7) && # Snell's law
                abs(dot(cross(n, d), out)) < 1e-7 &&                      # stays in the plane of incidence
                (dot(d, n) == 0 || sign(dot(out, n)) == sign(dot(d, n)))  # continues through the surface
        end
    end
    @check db = db max_examples = 500 function reflection_obeys_law(n=unitvec, d=unitvec)
        out = reflectedray(n, d)
        isapprox(norm(out), 1, atol=1e-9) && isapprox(dot(out, n), -dot(d, n), atol=1e-9)
    end
end
