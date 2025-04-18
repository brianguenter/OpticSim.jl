# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


@testitem "RayListSource" setup = [Dependencies] begin
    rays = Vector{OpticalRay{Float64,3}}(undef, 0)
    for i in 0:7
        λ = ((i / 7) * 200 + 450) / 1000
        r = OpticalRay(SVector(0.0, -3.0, 10.0), SVector(0.0, 0.5, -1.0), 1.0, λ)
        push!(rays, r)
    end
    raygen = Emitters.Sources.RayListSource(rays)

    for (i, ray) in enumerate(raygen)
        @assert ray == rays[i]
    end
end

#Angular power tests

@testitem "Lambertian" setup = [Dependencies] begin
    @test typeof(AngularPower.Lambertian()) === AngularPower.Lambertian{Float64}
    @test_throws MethodError AngularPower.Lambertian(String)

    @test apply(AngularPower.Lambertian(), Transform(), 1, Ray(zeros(3), ones(3))) === 1
end

@testitem "Cosine" setup = [Dependencies] begin
    @test AngularPower.Cosine().cosine_exp === one(Float64)
    @test_throws MethodError AngularPower.Cosine(String)

    @test apply(AngularPower.Cosine(), Transform(), 1, Ray(zeros(3), ones(3))) === 0.5773502691896258
end

@testitem "Gaussian" setup = [Dependencies] begin
    @test AngularPower.Gaussian(0, 1).gaussianu === 0
    @test AngularPower.Gaussian(0, 1).gaussianv === 1

    @test apply(AngularPower.Gaussian(0, 1), Transform(), 1, Ray(zeros(3), ones(3))) === 0.7165313105737892
end


#Directions tests

@testitem "Constant" setup = [Dependencies] begin
    @test Directions.Constant(0, 0, 0).direction === Geometry.Vec3(Int)
    @test Directions.Constant(Geometry.Vec3()).direction === Geometry.Vec3()
    @test Directions.Constant().direction === unitZ3()

    @test Base.length(Directions.Constant()) === 1
    @test Emitters.generate(Directions.Constant(), 0) === unitZ3()

    @test Base.iterate(Directions.Constant()) === (unitZ3(), 2)
    @test Base.iterate(Directions.Constant(), 2) === nothing
    @test Base.getindex(Directions.Constant(), 1) === unitZ3()
    @test Base.getindex(Directions.Constant(), 2) === unitZ3()
    @test Base.firstindex(Directions.Constant()) === 0
    @test Base.lastindex(Directions.Constant()) === 0
    @test Base.copy(Directions.Constant()) === Directions.Constant()
end

@testitem "RectGrid" setup = [Dependencies] begin
    @test Directions.RectGrid(unitX3(), 0.0, 0.0, 0, 0).uvec === unitY3()
    @test Directions.RectGrid(unitX3(), 0.0, 0.0, 0, 0).vvec === unitZ3()
    @test Directions.RectGrid(0.0, 0.0, 0, 0).uvec === unitX3()
    @test Directions.RectGrid(0.0, 0.0, 0, 0).vvec === unitY3()

    @test Base.length(Directions.RectGrid(0.0, 0.0, 2, 3)) === 6
    @test collect(Directions.RectGrid(ones(Geometry.Vec3), π / 4, π / 4, 2, 3)) == [
        [0.6014252485829169, 0.7571812496798511, 1.2608580392339483],
        [1.1448844430815608, 0.4854516524305293, 0.9891284419846265],
        [0.643629619770125, 1.0798265148205617, 1.0798265148205617],
        [1.225225479837374, 0.7890285847869373, 0.7890285847869373],
        [0.6014252485829169, 1.2608580392339483, 0.7571812496798511],
        [1.1448844430815608, 0.9891284419846265, 0.4854516524305293],
    ]
end

@testitem "UniformCone" setup = [Dependencies] begin
    @test Directions.UniformCone(unitX3(), 0.0, 0).uvec === unitY3()
    @test Directions.UniformCone(unitX3(), 0.0, 0).vvec === unitZ3()
    @test Directions.UniformCone(0.0, 0).uvec === unitX3()
    @test Directions.UniformCone(0.0, 0).vvec === unitY3()

    @test Base.length(Directions.UniformCone(0.0, 1)) === 1
end

@testitem "HexapolarCone" setup = [Dependencies] begin
    @test Directions.HexapolarCone(unitX3(), 0.0, 0).uvec === unitY3()
    @test Directions.HexapolarCone(unitX3(), 0.0, 0).vvec === unitZ3()
    @test Directions.HexapolarCone(0.0, 0).uvec === unitX3()
    @test Directions.HexapolarCone(0.0, 0).vvec === unitY3()

    @test Base.length(Directions.HexapolarCone(0.0, 1)) === 7
    @test Base.length(Directions.HexapolarCone(0.0, 2)) === 19
    @test collect(Directions.HexapolarCone(π / 4, 1)) == [
        [0.0, 0.0, 1.0],
        [0.7071067811865475, 0.0, 0.7071067811865476],
        [0.3535533905932738, 0.6123724356957945, 0.7071067811865476],
        [-0.3535533905932736, 0.6123724356957946, 0.7071067811865477],
        [-0.7071067811865475, 8.659560562354932e-17, 0.7071067811865476],
        [-0.35355339059327406, -0.6123724356957944, 0.7071067811865476],
        [0.3535533905932738, -0.6123724356957945, 0.7071067811865476],
    ]
end


#Origins tests
@testitem "Point" setup = [Dependencies] begin
    @test Origins.Point(Geometry.Vec3()).origin === Geometry.Vec3()
    @test Origins.Point(0, 0, 0).origin === Geometry.Vec3(Int)
    @test Origins.Point().origin === Geometry.Vec3()

    @test Base.length(Origins.Point()) === 1
    @test Emitters.visual_size(Origins.Point()) === 1
    @test Emitters.visual_size(Origins.Point()) === 1
    @test Emitters.generate(Origins.Point(), 1) === Geometry.Vec3()

    @test Base.iterate(Origins.Point(), 1) === (Geometry.Vec3(), 2)
    @test Base.iterate(Origins.Point(), 2) === nothing
    @test Base.getindex(Origins.Point(), 1) === Geometry.Vec3()
    @test Base.getindex(Origins.Point(), 2) === Geometry.Vec3()
    @test Base.firstindex(Origins.Point()) === 0
    @test Base.lastindex(Origins.Point()) === 0
    @test Base.copy(Origins.Point()) === Origins.Point()
end

@testitem "RectUniform" setup = [Dependencies] begin
    @test Origins.RectUniform(1, 2, 3).width === 1
    @test Origins.RectUniform(1, 2, 3).height === 2
    @test Origins.RectUniform(1, 2, 3).samples_count === 3

    @test Base.length(Origins.RectUniform(1, 2, 3)) === 3
    @test Emitters.visual_size(Origins.RectUniform(1, 2, 3)) === 2
end

@testitem "RectGrid" setup = [Dependencies] begin
    @test Origins.RectGrid(1.0, 2.0, 3, 4).ustep === 0.5
    @test Origins.RectGrid(1.0, 2.0, 3, 4).vstep === 2 / 3

    @test Base.length(Origins.RectGrid(1.0, 2.0, 3, 4)) === 12
    @test Emitters.visual_size(Origins.RectGrid(1.0, 2.0, 3, 4)) === 2.0

    @test collect(Origins.RectGrid(1.0, 2.0, 3, 4)) == [
        [-0.5, -1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.5, -1.0, 0.0],
        [-0.5, -0.33333333333333337, 0.0],
        [0.0, -0.33333333333333337, 0.0],
        [0.5, -0.33333333333333337, 0.0],
        [-0.5, 0.33333333333333326, 0.0],
        [0.0, 0.33333333333333326, 0.0],
        [0.5, 0.33333333333333326, 0.0],
        [-0.5, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 1.0, 0.0],
    ]
end

@testitem "RectJitterGrid" setup = [Dependencies] begin
    @test Origins.RectJitterGrid(1.0, 2.0, 3, 4, 1).width === 1.0
    @test Origins.RectJitterGrid(1.0, 2.0, 3, 4, 1).height === 2.0
    @test Origins.RectJitterGrid(1.0, 2.0, 2, 2, 3).ustep === 0.5
    @test Origins.RectJitterGrid(1.0, 2.0, 2, 2, 3).vstep === 1.0

    @test Base.length(Origins.RectJitterGrid(1.0, 2.0, 5, 6, 7)) === 210
    @test Emitters.visual_size(Origins.RectJitterGrid(1.0, 2.0, 5, 6, 7)) === 2.0
end

@testitem "Hexapolar" setup = [Dependencies] begin
    @test Origins.Hexapolar(1, 0, 0).nrings === 1
    @test_throws MethodError Origins.Hexapolar(1.0, 0, 0)

    @test Base.length(Origins.Hexapolar(1, 0, 0)) === 7
    @test Base.length(Origins.Hexapolar(2, 0, 0)) === 19
    @test Emitters.visual_size(Origins.Hexapolar(0, 1, 2)) === 4

    @test collect(Origins.Hexapolar(1, π / 4, π / 4)) == [
        [0.0, 0.0, 0.0],
        [0.7853981633974483, 0.0, 0.0],
        [0.39269908169872425, 0.6801747615878316, 0.0],
        [-0.392699081698724, 0.6801747615878317, 0.0],
        [-0.7853981633974483, 9.618353468608949e-17, 0.0],
        [-0.3926990816987245, -0.6801747615878315, 0.0],
        [0.39269908169872425, -0.6801747615878316, 0.0],
    ]
end

@testitem "Spectrum" setup = [Dependencies] begin
    @testset "Uniform" begin
        @test Spectrum.UNIFORMSHORT === 0.450
        @test Spectrum.UNIFORMLONG === 0.680
        @test Spectrum.Uniform(0, 1).low_end === 0
        @test Spectrum.Uniform(0, 1).high_end === 1
        @test Spectrum.Uniform().low_end === 0.450
        @test Spectrum.Uniform().high_end === 0.680
    end

    @testset "DeltaFunction" begin
        @test Emitters.generate(Spectrum.DeltaFunction(2)) === (1, 2)
        @test Emitters.generate(Spectrum.DeltaFunction(2.0)) === (1.0, 2.0)
        @test Emitters.generate(Spectrum.DeltaFunction(π)) === (true, π) # hopefully this is ok!
    end
end

#Sources tests
@testitem "Source" setup = [Dependencies] begin
    @test Sources.Source().transform === Transform()
    @test Sources.Source().spectrum === Spectrum.Uniform()
    @test Sources.Source().origins === Origins.Point()
    @test Sources.Source().directions === Directions.Constant()
    @test Sources.Source().power_distribution === AngularPower.Lambertian()
    @test Sources.Source().sourcenum === 0

    @test Sources.Source() === Sources.Source(
        Transform(),
        Spectrum.Uniform(),
        Origins.Point(),
        Directions.Constant(),
        AngularPower.Lambertian())

    @test Base.length(Sources.Source()) === 1
    @test Base.length(Sources.Source(; origins=Origins.Hexapolar(1, 0.0, 0.0))) === 7
    @test Base.length(Sources.Source(; directions=Directions.HexapolarCone(0.0, 1))) === 7
    @test Base.length(Sources.Source(;
        origins=Origins.Hexapolar(1, 0.0, 0.0),
        directions=Directions.HexapolarCone(0.0, 1)
    )) === 49

    @test Base.firstindex(Sources.Source()) === 0
    @test Base.lastindex(Sources.Source()) === 0
    @test Base.copy(Sources.Source()) === Sources.Source()
end

@testitem "CompositeSource" setup = [Dependencies] begin
    s() = Sources.Source(spectrum=Spectrum.Uniform(rng=Random.MersenneTwister(0)))
    tr = Transform()
    cs1 = Sources.CompositeSource(tr, [s()])
    cs2 = Sources.CompositeSource(tr, [s(), s()])
    cs3 = Sources.CompositeSource(tr, [s(), cs2])

    @test cs1.transform === tr

    @test cs1.uniform_length === 1
    @test cs2.uniform_length === 1
    @test cs3.uniform_length === -1

    @test cs1.total_length === 1
    @test cs2.total_length === 2
    @test cs3.total_length === 3

    @test cs1.start_indexes == []
    @test cs2.start_indexes == []
    @test cs3.start_indexes == [0, 1, 3]

    @test Base.length(cs1) === 1
    @test Base.length(cs2) === 2
    @test Base.length(cs3) === 3
end

