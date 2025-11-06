#code to compute s,p from surfacenormal and raydirection


"""
Computes the transformation that takes local s,p coordinates into global coordinates.

Assumes that normal and incident_vector are unit vectors."""
PolarizationTransform(normal::SVector{3,T}, propagation_vector::SVector{3,T}) where {T<:Real} = begin
    if normal ⋅ propagation_vector < 0
        propagation_vector = -propagation_vector #put propagation vector on same side of surface as normal to make sure transformation matrix is right handed.
    end
    s = normalize(cross(normal, propagation_vector))
    p = normalize(cross(propagation_vector, s))

    P = SMatrix{3,3,T,9}(vcat(transpose(propagation_vector), transpose(s'), transpose(p'))) #put vectors into rows of 3x3 matrix. In general the polarization matrix can have complex entries but this one is always real.

    return PolarizationMatrix{T}(P)
end

Base.adjoint(a::PolarizationMatrix{T}) where {T<:Real} = PolarizationMatrix(adjoint(a.P))

"""
    PolarizationRay