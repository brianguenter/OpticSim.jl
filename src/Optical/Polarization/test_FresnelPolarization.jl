using Test
include("FresnelPolarization.jl")

@testset "Fresnel Amplitude" begin
    # Normal incidence, nᵢ = nₜ
    nᵢ = 1.5
    nₜ = 1.5
    sinθᵢ = 0.0
    rₛ, tₛ, rₚ, tₚ = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    @test isapprox(rₛ, 0.0; atol=1e-12)
    @test isapprox(tₛ, 1.0; atol=1e-12)
    @test isapprox(rₚ, 0.0; atol=1e-12)
    @test isapprox(tₚ, 1.0; atol=1e-12)

    # Air to glass, normal incidence
    nᵢ = 1.0
    nₜ = 1.5
    sinθᵢ = 0.0
    rₛ, tₛ, rₚ, tₚ = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    expected_r = (nᵢ - nₜ) / (nᵢ + nₜ)
    expected_t = 2 * nᵢ / (nᵢ + nₜ)
    @test isapprox(rₛ, expected_r; atol=1e-12)
    @test isapprox(tₛ, expected_t; atol=1e-12)
    @test isapprox(rₚ, expected_r; atol=1e-12)
    @test isapprox(tₚ, expected_t; atol=1e-12)

    # Total internal reflection
    nᵢ = 1.5
    nₜ = 1.0
    sinθᵢ = nₜ / nᵢ + 1e-6
    rₛ, tₛ = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    @test isapprox(rₛ, 1.0; atol=1e-12)
    @test isapprox(tₛ, 0.0; atol=1e-12)
end

@testset "Fresnel Intensity" begin
    # Normal incidence, nᵢ = nₜ
    nᵢ = 1.5
    nₜ = 1.5
    cosθᵢ = 1.0
    cosθₜ = 1.0
    sinθᵢ = 0.0
    rₛ, tₛ, rₚ, tₚ = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    rₛ², tₛ², rₚ², tₚ² = fresnel_intensity(nᵢ, nₜ, cosθᵢ, cosθₜ)
    @test isapprox(rₛ², rₛ^2; atol=1e-12)
    @test isapprox(tₛ², tₛ^2; atol=1e-12)
    @test isapprox(rₚ², rₚ^2; atol=1e-12)
    @test isapprox(tₚ², tₚ^2; atol=1e-12)
end


@testset "Fresnel Amplitude with Complex Index" begin
    # Absorbing medium: nₜ complex
    nᵢ = 1.0
    nₜ = 1.5 + 0.1im
    sinθᵢ = 0.0
    rₛ, tₛ, rₚ, tₚ = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    # For normal incidence, r = (nᵢ - nₜ) / (nᵢ + nₜ)
    expected_r = (nᵢ - nₜ) / (nᵢ + nₜ)
    expected_t = 2 * nᵢ / (nᵢ + nₜ)
    @test isapprox(rₛ, expected_r; atol=1e-12)
    @test isapprox(tₛ, expected_t; atol=1e-12)
    @test isapprox(rₚ, expected_r; atol=1e-12)
    @test isapprox(tₚ, expected_t; atol=1e-12)

    # Metal: nₜ strongly complex
    nᵢ = 1.0
    nₜ = 0.2 + 3.0im
    sinθᵢ = 0.0
    rₛ, tₛ, rₚ, tₚ = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    expected_r = (nᵢ - nₜ) / (nᵢ + nₜ)
    expected_t = 2 * nᵢ / (nᵢ + nₜ)
    @test isapprox(rₛ, expected_r; atol=1e-12)
    @test isapprox(tₛ, expected_t; atol=1e-12)
    @test isapprox(rₚ, expected_r; atol=1e-12)
    @test isapprox(tₚ, expected_t; atol=1e-12)
end

@testset "Fresnel Intensity with Complex Index" begin
    nᵢ = 1.0
    nₜ = 1.5 + 0.1im
    cosθᵢ = 1.0
    # For normal incidence, cosθₜ can be computed as sqrt(1 - (nᵢ/nₜ * sinθᵢ)^2), sinθᵢ = 0
    cosθₜ = 1.0
    rₛ, tₛ, rₚ, tₚ = fresnel_amplitude(nᵢ, nₜ, 0.0)
    rₛ², tₛ², rₚ², tₚ² = fresnel_intensity(nᵢ, nₜ, cosθᵢ, cosθₜ)
    @test isapprox(rₛ², abs2(rₛ); atol=1e-12)
    @test isapprox(tₛ², abs2(tₛ); atol=1e-12)
    @test isapprox(rₚ², abs2(rₚ); atol=1e-12)
    @test isapprox(tₚ², abs2(tₚ); atol=1e-12)
end