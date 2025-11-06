# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


"""
Ray incorporating polarization information. Polarization rays are strictly 3D.

```julia
PolarizationRay(ray::Ray{T}, power::T, wavelength::T, opl=zero(T))
PolarizationRay(origin::SVector{N,T}, direction::SVector{N,T}, power::T, wavelength::T, opl=zero(T))
```

Has the following accessor methods:
```julia
ray(r::PolarizationRay{T}) -> Ray{T}
direction(r::PolarizationRay{T}) -> SVector{3,T}
origin(r::PolarizationRay{T}) -> SVector{3,T}
power(r::PolarizationRay{T}) -> T
wavelength(r::PolarizationRay{T}) -> T
pathlength(r::PolarizationRay{T}) -> T
sourcepower(r::PolarizationRay{T}) -> T
nhits(r::PolarizationRay{T}) -> Int
sourcenum(r::PolarizationRay{T}) -> Int
```
"""
struct PolarizationRay{T} <: AbstractRay{T,3}
    ray::Ray{T,3} 
    electric_field::SVector{3,Complex{T}}
    power::T
    wavelength::T
    opl::T
    nhits::Int
    sourcepower::T
    sourcenum::Int

    function PolarizationRay(ray::Ray{T,3}, electric_field::SVector{3,Complex{T}}, power::T, wavelength::T; opl::T=zero(T), nhits::Int=0, sourcenum::Int=0, sourcepower::T=power) where {T<:Real}
        return new{T}(ray, electric_field, power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    function PolarizationRay(origin::SVector{3,T}, electric_field::SVector{3,Complex{T}}, direction::SVector{3,T}, power::T, wavelength::T; opl::T=zero(T), nhits::Int=0, sourcenum::Int=0, sourcepower::T=power) where {T<:Real}
        return new{T}(Ray(origin, normalize(direction)), electric_field, power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    function PolarizationRay(origin::AbstractArray{T,1}, direction::AbstractArray{T,1}, electric_field::SVector{3,Complex{T}}, power::T, wavelength::T; opl::T=zero(T), nhits::Int=0, sourcenum::Int=0, sourcepower::T=power) where {T<:Real}
        @assert length(origin) == length(direction) "origin (dimension $(length(origin))) and direction (dimension $(length(direction))) vectors do not have the same dimension"
        N = length(origin)
        return new{T}(Ray(SVector{N,T}(origin), normalize(SVector{N,T}(direction))), electric_field, power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    # Convenience constructor. Not as much typing
    PolarizationRay(ox::T, oy::T, oz::T, dx::T, dy::T, dz::T, electric_field::SVector{3,Complex{T}}; wavelength=0.55) where {T<:Real} = PolarizationRay(SVector{3,T}(ox, oy, oz), SVector{3,T}(dx, dy, dz), electric_field, one(T), T(wavelength)) #doesn't have to be inside struct definition but if it is then VSCode displays hover information. If it's outside the struct definition it doesn't.
end
export PolarizationRay

electric_field(r::PolarizationRay{T}) where {T<:Real} = r.electric_field
ray(r::PolarizationRay) = r.ray
direction(r::PolarizationRay) = direction(ray(r))
origin(r::PolarizationRay) = origin(ray(r))
power(r::PolarizationRay) = r.power
wavelength(r::PolarizationRay) = r.wavelength
pathlength(r::PolarizationRay) = r.opl
nhits(r::PolarizationRay) = r.nhits
sourcepower(r::PolarizationRay) = r.sourcepower
sourcenum(r::PolarizationRay) = r.sourcenum
export ray, power, wavelength, pathlength, nhits, sourcepower, sourcenum

function Base.print(io::IO, a::PolarizationRay)
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

#not sure if this is ever used.
function Base.:*(a::Transform{T}, r::PolarizationRay)
    return PolarizationRay(a * ray(r), Geometry.rotate(a, r.electric_field), power(r), wavelength(r), opl=pathlength(r), nhits=nhits(r), sourcenum=sourcenum(r), sourcepower=sourcepower(r))
end
