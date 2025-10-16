#code to compute s,p from surfacenormal and raydirection

struct PolarizationMatrix{T} <: AbstractMatrix{T}
    P::SMatrix{3,3,Complex{T},9}
end
export PolarizationMatrix

"""
    ElectricField(eₓ::Complex{T}, e_y::Complex{T}, e_z::Complex{T}) where {T<:Real}
end
"""
struct ElectricField{T<:Real}
    E::SVector{3,Complex{T}}
end
export ElectricField

Base .* (a::PolarizationMatrix, b::ElectricField) = ElectricField(a.P * b.E)
Base .* (a::PolarizationMatrix, b::PolarizationMatrix) = PolarizationMatrix(a.P * b.P)

"""
Computes the transformation that takes local s,p coordinates into global coordinates.

Assumes that normal and incident_vector are unit vectors."""
PolarizationTransform(normal::SVector{3,T}, propagation_vector::SVector{3,T}) where {T<:Real} = begin
    if normal ⋅ propagation_vector < 0
        propagation_vector = -propagation_vector #put propagation vector on same side of surface as normal to make sure transformation matrix is right handed.
    end
    s = normalize(cross(normal, propagation_vector))
    p = normalize(cross(propagation_vector, s))

    P = SMatrix{3,3,Complex{T},9}(vcat(propagation_vector, s, p))

    return PolarizationMatrix{T}(P)
end

Base.adjoint(a::PolarizationMatrix{T}) where {T<:Real} = PolarizationMatrix(adjoint(a.P))

"""
    PolarizationRay