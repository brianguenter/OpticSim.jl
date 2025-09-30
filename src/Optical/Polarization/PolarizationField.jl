#Functions for creating polarization fields which can be applied to emitters or optical elements
#The field is defined over a rectangular domain in u,v space. This space is mapped onto the surface of the emitter or optical element. The polarization at a point is defined relative to the local coordinate frame of the surface. For example, a linear polarization field defined as SVector(1,0,0) will be horizontal polarization on a horizontal plane and radial polarization on a spherical surface.

struct PolarizationField{T} where {T<:Complex}
    field::Function  # function of (u,v) returning SVector{3,T} polarization vector
    umin::T
    umax::T
    vmin::T
    vmax::T
end

polarization(field::PolarizationField, u::T, v::T) where {T<:Real} = field.field(u, v)

#define some typical polarization fields

struct PolarizationSurface{F,S} where {F<:PolarizationField,S<:Surface}
    field::F
    surface::S
end

"""Need to rotate polarization vector based on local surface coordinates. Not sure if just rotating polarization vector is enough or if need to establish orthogonal coordinate frame for all our surface types. Most of them do have orthogonal coordinate frames defined but this is not required by the Surface interface."""
polarization(surface::PolarizationSurface, u::T, v::T) where {T<:Real} = polarization(surface.field, u, v)