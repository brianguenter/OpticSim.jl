include("FresnelPolarization.jl")
include("PolarizationField.jl")
include("PolarizationRay.jl")


struct PolarizationMatrix{T} <: AbstractMatrix{T}
    data::SMatrix{3,3,T,9}
end