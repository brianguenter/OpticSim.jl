# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    OpticalInterface{T<:Real}

Any subclass of OpticalInterface **must** implement the following:

```julia
processintersection(opticalinterface::OpticalInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, ::Bool, firstray::Bool = false) -> Tuple{SVector{N,T}, T, T}
insidematerial(i::OpticalInterface{T}) -> AGFFileReader.AbstractGlass
outsidematerial(i::OpticalInterface{T}) -> AGFFileReader.AbstractGlass
reflectance(i::OpticalInterface{T}) -> T
transmission(i::OpticalInterface{T}) -> T
```

See documentation for [`processintersection`](@ref) for details.
"""
abstract type OpticalInterface{T<:Real} end
export OpticalInterface

"""
    GlassID

A small handle for a glass, used by interfaces in place of the glass itself. A glass holds kilobytes of data that can't
be stored inline, so an interface that held one couldn't be plain data, and every intersection that stored its interface
had to be allocated on the heap. With glasses referred to by ID, interfaces and intersections are plain data.

`glassid(glass)` registers a glass, once, and returns its ID; `glass(id)` returns the glass. Air is always ID 0.
"""
struct GlassID
    id::Int32
end

# The registry of glasses. A lookup reads GLASSES[] without locking: registering a glass replaces the vector rather than
# resizing it, so a lookup on another thread always sees a complete vector.
const GLASSES = Ref(AGFFileReader.Glass[])
const GLASSIDS = Dict{AGFFileReader.Glass,Int32}()
const GLASSLOCK = ReentrantLock()

glassid(::AGFFileReader.AirType) = GlassID(0)
function glassid(g::AGFFileReader.Glass)
    lock(GLASSLOCK) do
        GlassID(get!(GLASSIDS, g) do
            glasses = copy(GLASSES[])
            push!(glasses, g)
            GLASSES[] = glasses
            Int32(length(glasses))
        end)
    end
end

"""Returns the glass with ID `id`: the air singleton or a Glass, a small union that compiles to a branch."""
glass(id::GlassID)::Union{AGFFileReader.AirType,AGFFileReader.Glass} = id.id == 0 ? AGFFileReader.Air : @inbounds GLASSES[][id.id]

"""
Valid modes for deterministic raytracing
"""
@enum InterfaceMode Reflect Transmit ReflectOrTransmit
export InterfaceMode, Reflect, Transmit, ReflectOrTransmit

######################################################################

"""
    NullInterface{T} <: OpticalInterface{T}

Interface which will be ignored totally by any rays, used only in construction of CSG objects.

```julia
NullInterface(T = Float64)
NullInterface{T}()
```
"""
struct NullInterface{T} <: OpticalInterface{T}
    NullInterface(::Type{T}=Float64) where {T<:Real} = new{T}()
    NullInterface{T}() where {T<:Real} = new{T}()
end

insidematerial(a::NullInterface) = AGFFileReader.Air
outsidematerial(a::NullInterface) = AGFFileReader.Air

reflectance(::NullInterface{T}) where {T<:Real} = zero(T)
transmission(::NullInterface{T}) where {T<:Real} = one(T)

export NullInterface

######################################################################

"""
    ParaxialInterface{T} <: OpticalInterface{T}

Interface describing an idealized planar lens, i.e. one that is thin and with no aberrations.

**In general this interface should not be constructed directly, the `ParaxialLensEllipse` and `ParaxialLensRect` functions should be used to create a [`ParaxialLens`](@ref) object directly.**

```julia
ParaxialInterface(focallength::T, centroid::SVector{3,T}, outsidematerial::Y)
```
"""
struct ParaxialInterface{T} <: OpticalInterface{T}
    focallength::T
    outsidematerial::GlassID
    centroid::SVector{3,T}
    function ParaxialInterface(focallength::T, centroid::SVector{3,T}, outsidematerial::AGFFileReader.AbstractGlass) where {T<:Real}
        return new{T}(focallength, glassid(outsidematerial), centroid)
    end
end
export ParaxialInterface

function Base.show(io::IO, a::ParaxialInterface{R}) where {R<:Real}
    print(io, "ParaxialInterface($(a.focallength), $(glassname(outsidematerial(a))))")
end

insidematerial(a::ParaxialInterface{T}) where {T<:Real} = AGFFileReader.Air
outsidematerial(a::ParaxialInterface{T}) where {T<:Real} = glass(a.outsidematerial)
reflectance(::ParaxialInterface{T}) where {T<:Real} = zero(T)
transmission(::ParaxialInterface{T}) where {T<:Real} = one(T)
opticalcenter(a::ParaxialInterface) = a.centroid
focallength(a::ParaxialInterface) = a.focallength

######################################################################

"""
    FresnelInterface{T} <: OpticalInterface{T}

Interface between two materials with behavior defined according to the [Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations), with a specified reflectance and transmission.
Assumes unpolarized light.

```julia
FresnelInterface{T}(insidematerial, outsidematerial; reflectance = 0, transmission = 1, interfacemode = ReflectOrTransmit)
```

The interfacemode can be used to trace rays deterministically. Valid values are defined in the InterfaceMode enum.
Reflect means that all values are reflected, Transmit means that all values are transmitted. ReflectOrTransmit will randomly
reflect and transmit rays with the distribution given by the reflection and transmission arguments. This is also the default.
In all cases the power recorded with the ray is correctly updated. This can be used to fake sequential raytracing. For
example a beamsplitter surface may be set to either Reflect or Transmit to switch between the two outgoing ray paths.

"""
struct FresnelInterface{T} <: OpticalInterface{T}
    insidematerial::GlassID
    outsidematerial::GlassID
    reflectance::T
    transmission::T
    interfacemode::InterfaceMode

    function FresnelInterface{T}(insidematerial::AGFFileReader.AbstractGlass, outsidematerial::AGFFileReader.AbstractGlass; reflectance::T=zero(T), transmission::T=one(T), interfacemode=ReflectOrTransmit) where {T<:Real}
        return new{T}(glassid(insidematerial), glassid(outsidematerial), reflectance, transmission, interfacemode)
    end

    # from glasses that are already registered, e.g. when an interface is built while tracing
    function FresnelInterface{T}(insidematerial::GlassID, outsidematerial::GlassID; reflectance::T=zero(T), transmission::T=one(T), interfacemode=ReflectOrTransmit) where {T<:Real}
        return new{T}(insidematerial, outsidematerial, reflectance, transmission, interfacemode)
    end

end
export FresnelInterface

function Base.show(io::IO, a::FresnelInterface{R}) where {R<:Real}
    print(io, "FresnelInterface($(glassname(insidematerial(a))), $(glassname(outsidematerial(a))), $(a.reflectance), $(a.transmission), $(a.interfacemode))")
end

insidematerial(a::FresnelInterface) = glass(a.insidematerial)
outsidematerial(a::FresnelInterface) = glass(a.outsidematerial)

reflectance(a::FresnelInterface{T}) where {T<:Real} = a.reflectance
transmission(a::FresnelInterface{T}) where {T<:Real} = a.transmission
interfacemode(a::FresnelInterface{T}) where {T<:Real} = a.interfacemode

transmissiveinterface(::Type{T}, insidematerial::X, outsidematerial::Y) where {T<:Real,X<:AGFFileReader.AbstractGlass,Y<:AGFFileReader.AbstractGlass} = FresnelInterface{T}(insidematerial, outsidematerial, reflectance=zero(T), transmission=one(T))
reflectiveinterface(::Type{T}, insidematerial::X, outsidematerial::Y) where {T<:Real,X<:AGFFileReader.AbstractGlass,Y<:AGFFileReader.AbstractGlass} = FresnelInterface{T}(insidematerial, outsidematerial, reflectance=one(T), transmission=zero(T))
opaqueinterface(::Type{T}=Float64) where {T<:Real} = FresnelInterface{T}(AGFFileReader.Air, AGFFileReader.Air, reflectance=zero(T), transmission=zero(T))
export opaqueinterface

######################################################################

"""
    ThinGratingInterface{T} <: OpticalInterface{T}

Interface representing an idealized thin grating. `period` is in microns, `vector` should lie in the plane of the surface.
Transmission and reflectance can be specified for an arbitrary number of orders up to 10, selected using the `maxorder` and `minorder` parameters.
If `nothing` then `reflectance` is assumed to be **0** and `transmission` is assumed to be **1**.

```julia
ThinGratingInterface(vector, period, insidematerial, outsidematerial; maxorder = 1, minorder = -1, reflectance = nothing, transmission = nothing)
```
"""
struct ThinGratingInterface{T} <: OpticalInterface{T}
    insidematerial::GlassID
    outsidematerial::GlassID
    vector::SVector{3,T}
    period::T
    maxorder::Int
    minorder::Int
    transmission::SVector{GRATING_MAX_ORDERS,T}
    reflectance::SVector{GRATING_MAX_ORDERS,T}

    function ThinGratingInterface(vector::SVector{3,T}, period::T, insidematerial::AGFFileReader.AbstractGlass, outsidematerial::AGFFileReader.AbstractGlass; maxorder::Int=1, minorder::Int=-1, reflectance::Union{Nothing,AbstractVector{T}}=nothing, transmission::Union{Nothing,AbstractVector{T}}=nothing) where {T<:Real}
        @assert maxorder >= minorder
        norders = maxorder - minorder + 1
        @assert norders <= GRATING_MAX_ORDERS "Thin grating is limited to $GRATING_MAX_ORDERS orders"
        @assert ((reflectance === nothing) || length(reflectance) == norders) && ((transmission === nothing) || length(transmission) == norders)
        if reflectance !== nothing
            sreflectance = vcat(SVector{length(reflectance),T}(reflectance), zeros(SVector{GRATING_MAX_ORDERS - length(reflectance),T}))
        else
            sreflectance = zeros(SVector{GRATING_MAX_ORDERS,T})
        end
        if transmission !== nothing
            stransmission = vcat(SVector{length(transmission),T}(transmission), zeros(SVector{GRATING_MAX_ORDERS - length(transmission),T}))
        else
            # spread the transmitted power evenly over the orders in use
            stransmission = SVector{GRATING_MAX_ORDERS,T}(ntuple(i -> i <= norders ? one(T) / norders : zero(T), GRATING_MAX_ORDERS))
        end
        @assert zero(T) <= sum(sreflectance) <= one(T)
        @assert zero(T) <= sum(stransmission) <= one(T)
        @assert zero(T) <= sum(sreflectance) + sum(stransmission) <= one(T) + sqrt(eps(T))
        new{T}(glassid(insidematerial), glassid(outsidematerial), normalize(vector), period, maxorder, minorder, stransmission, sreflectance)
    end
end
export ThinGratingInterface

insidematerial(a::ThinGratingInterface) = glass(a.insidematerial)
outsidematerial(a::ThinGratingInterface) = glass(a.outsidematerial)

reflectance(a::ThinGratingInterface{T}, order::Int) where {T<:Real} = a.reflectance[order-a.minorder+1]
transmission(a::ThinGratingInterface{T}, order::Int) where {T<:Real} = a.transmission[order-a.minorder+1]

function Base.show(io::IO, a::ThinGratingInterface{R}) where {R<:Real}
    print(io, "ThinGratingInterface($(glassname(insidematerial(a))), $(glassname(outsidematerial(a))), $(a.vector), $(a.period), $(a.transmission), $(a.reflectance))")
end

######################################################################

"""
`ConvergingBeam`, `DivergingBeam` or `CollimatedBeam`, defines the behavior of a beam in a [`HologramInterface`](@ref).
"""
@enum BeamState ConvergingBeam DivergingBeam CollimatedBeam
export BeamState, ConvergingBeam, DivergingBeam, CollimatedBeam

"""
    HologramInterface{T} <: OpticalInterface{T}

Interface representing a _thick_ hologram (though geometrically thin).
The efficiency, `η`, is calculated using Kogelnik's coupled wave theory so is only valid for the first order.
If the zero order is included then it has efficiency `1 - η`.
Also assumes that the HOE was recorded under similar conditions to the playback conditions, `thickness` is in microns.

`BeatState` arguments can be one of `ConvergingBeam`, `DivergingBeam` and `CollimatedBeam`. In the first two cases `signalpointordir` and `referencepointordir` are 3D point in global coordinate space. For `CollimatedBeam` they are normalized direction vectors.

For reference, see:
- _Coupled Wave Theory for Thick Hologram Gratings_ - H Kogelnik, 1995
- _Sequential and non-sequential simulation of volume holographic gratings_ - M Kick et al, 2018

```julia
HologramInterface(signalpointordir::SVector{3,T}, signalbeamstate::BeamState, referencepointordir::SVector{3,T}, referencebeamstate::BeamState, recordingλ::T, thickness::T, beforematerial, substratematerial, aftermaterial, signalrecordingmaterial, referencerecordingmaterial, RImodulation::T, include0order  = false)
```
"""
struct HologramInterface{T} <: OpticalInterface{T}
    beforematerial::GlassID
    substratematerial::GlassID
    aftermaterial::GlassID
    signalpointordir::SVector{3,T}
    signalbeamstate::BeamState
    referencepointordir::SVector{3,T}
    referencebeamstate::BeamState
    recordingλ::T
    signalrecordingmaterial::GlassID
    referencerecordingmaterial::GlassID
    thickness::T
    RImodulation::T
    include0order::Bool

    function HologramInterface(signalpointordir::SVector{3,T}, signalbeamstate::BeamState, referencepointordir::SVector{3,T}, referencebeamstate::BeamState, recordingλ::T, thickness::T, beforematerial::AGFFileReader.AbstractGlass, substratematerial::AGFFileReader.AbstractGlass, aftermaterial::AGFFileReader.AbstractGlass, signalrecordingmaterial::AGFFileReader.AbstractGlass, referencerecordingmaterial::AGFFileReader.AbstractGlass, RImodulation::T, include0order::Bool=false) where {T<:Real}
        @assert substratematerial !== AGFFileReader.Air "Substrate can't be air"
        if signalbeamstate === CollimatedBeam
            signalpointordir = normalize(signalpointordir)
        end
        if referencebeamstate === CollimatedBeam
            referencepointordir = normalize(referencepointordir)
        end
        new{T}(glassid(beforematerial), glassid(substratematerial), glassid(aftermaterial), signalpointordir, signalbeamstate, referencepointordir, referencebeamstate, recordingλ, glassid(signalrecordingmaterial), glassid(referencerecordingmaterial), thickness, RImodulation, include0order)
    end
end
export HologramInterface

insidematerial(a::HologramInterface{T}) where {T<:Real} = glass(a.beforematerial)
outsidematerial(a::HologramInterface{T}) where {T<:Real} = glass(a.aftermaterial)
substratematerial(a::HologramInterface{T}) where {T<:Real} = glass(a.substratematerial)
signalrecordingmaterial(a::HologramInterface{T}) where {T<:Real} = glass(a.signalrecordingmaterial)
referencerecordingmaterial(a::HologramInterface{T}) where {T<:Real} = glass(a.referencerecordingmaterial)

function Base.show(io::IO, a::HologramInterface{R}) where {R<:Real}
    print(io, "HologramInterface($(glassname(insidematerial(a))), $(glassname(substratematerial(a))), $(glassname(outsidematerial(a))), $(a.signalpointordir), $(a.signalbeamstate), $(a.referencepointordir), $(a.referencebeamstate), $(a.recordingλ), $(a.thickness), $(a.RImodulation))")
end

######################################################################

# a few things need a concrete type to prevent allocations resulting from the ambiguities introduce by an abstract type (i.e. OpticalInterface{T})
# if adding new subtypes of OpticalInterface they must be added to this definition as well (and only paramaterised by T so they are fully specified)
# They must also be plain data (refer to glasses by GlassID): then intersections and intervals are plain data, which
# Julia stores inline and returns without allocating.
AllOpticalInterfaces{T} = Union{NullInterface{T},ParaxialInterface{T},FresnelInterface{T},HologramInterface{T},ThinGratingInterface{T}}
NullOrFresnel{T} = Union{NullInterface{T},FresnelInterface{T}}
