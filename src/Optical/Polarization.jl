struct PolarizationMatrix{T} <: AbstractMatrix{T}
    data::SMatrix{3,3,T,9}
end