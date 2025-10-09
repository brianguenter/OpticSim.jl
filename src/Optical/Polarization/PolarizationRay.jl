# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    PolarizationRay{T,N} <: AbstractRay{T,N}

Ray incorporating polarization information.

**NOTE**: we use monte carlo integration to get accurate results on the detector, this means that all rays essentially hit the detector with power = 1 and some rays are thrown away at any interface to correctly match the reflection/transmission at that interface. For inspection purposes we also track the 'instantaneous' power of the ray in the `power` field of the `PolarizationRay`.

```julia
PolarizationRay(ray::Ray{T,N}, power::T, wavelength::T, opl=zero(T))
PolarizationRay(origin::SVector{N,T}, direction::SVector{N,T}, power::T, wavelength::T, opl=zero(T))
```

Has the following accessor methods:
```julia
ray(r::PolarizationRay{T,N}) -> Ray{T,N}
direction(r::PolarizationRay{T,N}) -> SVector{N,T}
origin(r::PolarizationRay{T,N}) -> SVector{N,T}
power(r::PolarizationRay{T,N}) -> T
wavelength(r::PolarizationRay{T,N}) -> T
pathlength(r::PolarizationRay{T,N}) -> T
sourcepower(r::PolarizationRay{T,N}) -> T
nhits(r::PolarizationRay{T,N}) -> Int
sourcenum(r::PolarizationRay{T,N}) -> Int
```
"""
struct PolarizationRay{T,N} <: AbstractRay{T,N}
    ray::Ray{T,N}
    power::T
    wavelength::T
    opl::T
    nhits::Int
    sourcepower::T
    sourcenum::Int

    function PolarizationRay(ray::Ray{T,N}, power::T, wavelength::T; opl::T=zero(T), nhits::Int=0, sourcenum::Int=0, sourcepower::T=power) where {T<:Real,N}
        return new{T,N}(ray, power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    function PolarizationRay(origin::SVector{N,T}, direction::SVector{N,T}, power::T, wavelength::T; opl::T=zero(T), nhits::Int=0, sourcenum::Int=0, sourcepower::T=power) where {T<:Real,N}
        return new{T,N}(Ray(origin, normalize(direction)), power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    function PolarizationRay(origin::AbstractArray{T,1}, direction::AbstractArray{T,1}, power::T, wavelength::T; opl::T=zero(T), nhits::Int=0, sourcenum::Int=0, sourcepower::T=power) where {T<:Real}
        @assert length(origin) == length(direction) "origin (dimension $(length(origin))) and direction (dimension $(length(direction))) vectors do not have the same dimension"
        N = length(origin)
        return new{T,N}(Ray(SVector{N,T}(origin), normalize(SVector{N,T}(direction))), power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    # Convenience constructor. Not as much typing
    PolarizationRay(ox::T, oy::T, oz::T, dx::T, dy::T, dz::T; wavelength=0.55) where {T<:Real} = PolarizationRay(SVector{3,T}(ox, oy, oz), SVector{3,T}(dx, dy, dz), one(T), T(wavelength)) #doesn't have to be inside struct definition but if it is then VSCode displays hover information. If it's outside the struct definition it doesn't.
end
export PolarizationRay

ray(r::PolarizationRay) = r.propagation_vector
direction(r::PolarizationRay{T,N}) where {T<:Real,N} = direction(ray(r))
origin(r::PolarizationRay{T,N}) where {T<:Real,N} = origin(ray(r))
power(r::PolarizationRay{T,N}) where {T<:Real,N} = r.power
wavelength(r::PolarizationRay{T,N}) where {T<:Real,N} = r.wavelength
pathlength(r::PolarizationRay{T,N}) where {T<:Real,N} = r.opl
nhits(r::PolarizationRay{T,N}) where {T<:Real,N} = r.nhits
sourcepower(r::PolarizationRay{T,N}) where {T<:Real,N} = r.sourcepower
sourcenum(r::PolarizationRay{T,N}) where {T<:Real,N} = r.sourcenum
export ray, power, wavelength, pathlength, nhits, sourcepower, sourcenum

function Base.print(io::IO, a::PolarizationRay{T,N}) where {T,N}
    println(io, "$(rpad("Origin:", 25)) $(origin(a))")
    println(io, "$(rpad("Direction:", 25)) $(direction(a))")
    println(io, "$(rpad("Power:", 25)) $(power(a))")
    println(io, "$(rpad("Source Power:", 25)) $(sourcepower(a))")
    println(io, "$(rpad("Wavelength (in air):", 25)) $(wavelength(a))")
    println(io, "$(rpad("Optical Path Length:", 25)) $(pathlength(a))")
    println(io, "$(rpad("Hits:", 25)) $(nhits(a))")
    if sourcenum(a) != 0
        println(io, "$(rpad("Source Number:", 25)) $(sourcenum(a))")
    end
end

function Base.:*(a::Transform{T}, r::PolarizationRay{T,N}) where {T,N}
    return PolarizationRay(a * ray(r), power(r), wavelength(r), opl=pathlength(r), nhits=nhits(r), sourcenum=sourcenum(r), sourcepower=sourcepower(r))
end
