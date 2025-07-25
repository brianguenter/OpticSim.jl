# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module OpticSim

import Unitful
using LinearAlgebra: eigen, svd, I, qr, dot, cross, norm, det, normalize, inv
import LinearAlgebra

using StaticArrays
using DataFrames: DataFrame

using Base: @.

using StringEncodings
using AGFFileReader

# included here to allow a call to the activate! during the initialization


include("constants.jl")
include("utilities.jl")





include("Geometry/Geometry.jl")
include("Optical/Optical.jl")

include("Data/Data.jl")

#initialize these caches here so they will get the correct number of threads from the load time environment, rather than the precompile environment. The latter happens if the initialization happens in the const definition. If the precompile and load environments have different numbers of threads this will cause an error.
function __init__()
    for _ in 1:Threads.nthreads()
        push!(threadedtrianglepool, Dict{DataType,TrianglePool}((Float64 => TrianglePool{Float64}())))
        push!(threadedintervalpool, Dict{DataType,IntervalPool}((Float64 => IntervalPool{Float64}())))
    end
end


################################################

# This can be used to track NaN, particularly in ForwardDiff gradients, causing problems
# e.g. Diagnostics.testoptimization(lens = Examples.doubleconvex(NaNCheck{Float64}), samples = 1)

# struct NaNCheck{T<:Real} <: Real
#     val::T
#     function NaNCheck{T}(a::S) where {T<:Real, S<:Real}
#         @assert !(T <: NaNCheck)
#         new{T}(T(a))
#     end
# end
# export NaNCheck
# Base.isnan(a::NaNCheck{T}) where{T} = isnan(a.val)
# Base.isinf(a::NaNCheck{T}) where{T} = isinf(a.val)
# Base.typemin(::Type{NaNCheck{T}}) where{T} = NaNCheck{T}(typemin(T))
# Base.typemax(::Type{NaNCheck{T}}) where{T} = NaNCheck{T}(typemax(T))
# Base.eps(::Type{NaNCheck{T}}) where {T} = NaNCheck{T}(eps(T))
# Base.decompose(a::NaNCheck{T}) where {T} = Base.decompose(a.val)
# Base.round(a::NaNCheck{T}, m::RoundingMode) where {T} = NaNCheck{T}(round(a.val, m))

# struct NaNException <: Exception end

# # (::Type{Float64})(a::NaNCheck{S}) where {S<:Real} = NaNCheck{Float64}(Float64(a.val))
# (::Type{T})(a::NaNCheck{S}) where {T<:Integer,S<:Real} = T(a.val)
# (::Type{NaNCheck{T}})(a::NaNCheck{S}) where {T<:Real,S<:Real} = NaNCheck{T}(T(a.val))
# Base.promote_rule(::Type{NaNCheck{T}}, ::Type{T}) where {T<:Number} = NaNCheck{T}
# Base.promote_rule(::Type{T}, ::Type{NaNCheck{T}}) where {T<:Number} = NaNCheck{T}
# Base.promote_rule(::Type{S}, ::Type{NaNCheck{T}}) where {T<:Number, S<:Number} = NaNCheck{promote_type(T,S)}
# Base.promote_rule(::Type{NaNCheck{T}}, ::Type{S}) where {T<:Number, S<:Number} = NaNCheck{promote_type(T,S)}
# Base.promote_rule(::Type{NaNCheck{S}}, ::Type{NaNCheck{T}}) where {T<:Number, S<:Number} = NaNCheck{promote_type(T,S)}

# for op = (:sin, :cos, :tan, :log, :exp, :sqrt, :abs, :-, :atan, :acos, :asin, :log1p, :floor, :ceil, :float)
#     eval(quote
#         function Base.$op(a::NaNCheck{T}) where{T}
#             temp = NaNCheck{T}(Base.$op(a.val))
#             if isnan(temp)
#                 throw(NaNException())
#             end
#             return temp
#         end
#     end)
# end

# for op = (:+, :-, :/, :*, :^, :atan)
#     eval(quote
#         function Base.$op(a::NaNCheck{T}, b::NaNCheck{T}) where{T}
#             temp = NaNCheck{T}(Base.$op(a.val, b.val))
#             if isnan(temp)
#                 throw(NaNException())
#             end
#             return temp
#         end
#     end)
# end

# for op =  (:<, :>, :<=, :>=, :(==), :isequal)
#     eval(quote
#         function Base.$op(a::NaNCheck{T}, b::NaNCheck{T}) where{T}
#             temp = Base.$op(a.val, b.val)
#             return temp
#         end
#     end)
# end

################################################

end # module
