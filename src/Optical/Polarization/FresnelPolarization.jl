"""
    fresnel(nᵢ, nₜ, sinθᵢ) -> (rₛ, tₛ, rₚ, tₚ)
where 
rₛ is amplitude of reflected s polarized light

tₛ is amplitude of transmitted s polarized light

rₚ is amplitude of reflected p polarized light

tₚ is amplitude of transmitted p polarized light

Fresnel equations from Chipman et al, "Polarized Light and Optical Systems", 2018, eqn 8.12-8.15 pg. 300
"""
function fresnel_amplitude(nᵢ::T, nₜ::T, sinθᵢ::T) where {T<:Real}
    if (sinθᵢ >= nₜ / nᵢ) # 100% reflectance, zero transmission
        return (one(T), zero(T))
    end

    cosθᵢ = sqrt(one(T) - sinθᵢ^2)
    #Fresnel equations from Chipman et al, "Polarized Light and Optical Systems", 2018, eqn 8.12-8.15 pg. 300

    n = nₜ / nᵢ
    #rₛ is amplitude of reflected s polarized light
    rₛ = (cosθᵢ - sqrt(n^2 - sinθᵢ^2)) / (cosθᵢ + sqrt(n^2 - sinθᵢ^2))
    #tₛ is amplitude of transmitted s polarized light
    tₛ = (2 * cosθᵢ) / (cosθᵢ + sqrt(n^2 - sinθᵢ^2))
    #rₚ is amplitude of reflected p polarized light
    rₚ = (n^2 * cosθᵢ - sqrt(n^2 - sinθᵢ^2)) / (n^2 * cosθᵢ + sqrt(n^2 - sinθᵢ^2))
    #tₚ is amplitude of transmitted p polarized light
    tₚ = (2 * n * cosθᵢ) / (n^2 * cosθᵢ + sqrt(n^2 - sinθᵢ^2))

    return rₛ, tₛ, rₚ, tₚ
end
export fresnel_amplitude

"""
    fresnel_intensity(nᵢ::T, nₜ::T, cosθᵢ::T, cosθₜ::T) -> Tuple{T,T,T,T}

    returns s and p polarization intensities with geometric correction for transmitted intensity.
"""
function fresnel_intensity(nᵢ::T, nₜ::T, cosθᵢ::T, cosθₜ::T) where {T<:Real}
    (rₛ, tₛ, rₚ, tₚ) = fresnel_amplitude(nᵢ, nₜ, sinθᵢ)
    rₛ² = rₛ^2
    rₚ² = rₚ^2
    tₛ² = tₛ^2 * (nₜ * cosθₜ) / (nᵢ * cosθᵢ) #transmitted intensity correction factor

    tₚ² = tₚ^2 * (nₜ * cosθₜ) / (nᵢ * cosθᵢ) #transmitted intensity correction factor

    return rₛ², tₛ², rₚ², tₚ²
end

export fresnel_intensity
