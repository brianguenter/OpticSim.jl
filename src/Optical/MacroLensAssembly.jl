# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# The original macro-generated lens assembly. Superseded by TupleLensAssembly (see LensAssembly.jl), which is what
# LensAssembly(...) now creates. Kept only so the two can be benchmarked against each other; delete this file once that
# is no longer needed.
#
# Create one with MacroLensAssembly(elements...; axis). It behaves like any other LensAssembly.

abstract type MacroLensAssembly{T<:Real} <: LensAssembly{T} end

# in order to prevent type ambiguities lens assembly must be paramaterised by its elements, but we want lens assembly to have an arbitrary number of elements
# as such we have to use macros to create lens assembly structs of a given size, and methods to handle them in a type unambiguous way
# we may often end up with lots of simple surfaces (Rectangle, Ellipse) in the assembly, and we can store these more efficiently in Vectors
macro lensassembly_constructor(N)
    # generate the type definition for a lens assembly of size N
    typesig = [:($(Symbol("T$(i)")) <: Union{Surface{T},CSGTree{T}}) for i in 1:N]
    fields = [:($(Symbol("E$(i)"))::$(Symbol("T$(i)"))) for i in 1:N]
    name = Symbol("MacroLensAssembly$N")
    return esc(quote
        struct $(name){T,$(typesig...)} <: MacroLensAssembly{T}
            axis::SVector{3,T}
            rectangles::Vector{Rectangle{T}}
            ellipses::Vector{Ellipse{T}}
            hexagons::Vector{Hexagon{T}}
            paraxials::Vector{ParaxialLens{T}}
            $(fields...)
        end
    end)
end

macro lensassembly_intersection(N)
    # generates the _closestintersection() function for a lens assembly of size N, this should never be called directly, instead call closestintersection() in all cases
    type = :($(Symbol("MacroLensAssembly$(N)")))
    fieldnames = [Symbol("E$(i)") for i in 1:N]
    checks = [
        quote
            intvl = surfaceintersection(obj.$(fieldnames[i]), r)
            if !(intvl isa EmptyInterval)
                intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                if intsct !== nothing
                    αcurr = α(intsct)
                    if αcurr < αmin
                        αmin = αcurr
                        closest = intsct
                    end
                end
            end
        end for i in 1:N
    ]
    return esc(quote
        function _closestintersection(obj::$(type){T}, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
            αmin = typemax(T)
            closest = nothing
            for rect in obj.rectangles
                intvl = surfaceintersection(rect, r)
                if !(intvl isa EmptyInterval)
                    intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                    if intsct !== nothing
                        αcurr = α(intsct)
                        if αcurr < αmin
                            αmin = αcurr
                            closest = intsct
                        end
                    end
                end
            end
            for ell in obj.ellipses
                intvl = surfaceintersection(ell, r)
                if !(intvl isa EmptyInterval)
                    intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                    if intsct !== nothing
                        αcurr = α(intsct)
                        if αcurr < αmin
                            αmin = αcurr
                            closest = intsct
                        end
                    end
                end
            end
            for hex in obj.hexagons
                intvl = surfaceintersection(hex, r)
                if !(intvl isa EmptyInterval)
                    intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                    if intsct !== nothing
                        αcurr = α(intsct)
                        if αcurr < αmin
                            αmin = αcurr
                            closest = intsct
                        end
                    end
                end
            end
            for par in obj.paraxials
                intvl = surfaceintersection(par, r)
                if !(intvl isa EmptyInterval)
                    intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                    if intsct !== nothing
                        αcurr = α(intsct)
                        if αcurr < αmin
                            αmin = αcurr
                            closest = intsct
                        end
                    end
                end
            end
            $(checks...)
            return closest
        end
    end)
end

macro lensassembly_elements(N)
    # generates the _typed_elements() function for a lens assembly of size N, this should never be called directly, instead call typed_elements() in all cases
    type = :($(Symbol("MacroLensAssembly$(N)")))
    fieldnames = [Symbol("E$(i)") for i in 1:N]
    tuple = [:(obj.$(fieldnames[i])) for i in 1:N]
    return esc(quote
        function _typed_elements(obj::$(type){T}) where {T<:Real}
            return tuple($(tuple...))
        end
    end)
end

# pregenerate a small number at compile time
for N in 1:PREGENERATED_LENS_ASSEMBLY_SIZE
    type = :($(Symbol("MacroLensAssembly$(N)")))
    if !isdefined(OpticSim, type)
        @eval @lensassembly_constructor($N)
        @eval @lensassembly_intersection($N)
        @eval @lensassembly_elements($N)
    end
end

# first element is a separate argument so T is always bound (at least one element is required)
function MacroLensAssembly(first::Union{Surface{T},CSGTree{T},LensAssembly{T}}, rest::Vararg{Union{Surface{T},CSGTree{T},LensAssembly{T}}}; axis::SVector{3,T}=SVector{3,T}(0.0, 0.0, 1.0)) where {T<:Real}
    elements = (first, rest...)
    # make the actual object
    actual_elements = []
    rectangles = Vector{Rectangle{T}}(undef, 0)
    ellipses = Vector{Ellipse{T}}(undef, 0)
    hexagons = Vector{Hexagon{T}}(undef, 0)
    paraxials = Vector{ParaxialLens{T}}(undef, 0)
    for e in elements
        if e isa Rectangle{T}
            push!(rectangles, e)
        elseif e isa Ellipse{T}
            push!(ellipses, e)
        elseif e isa Hexagon{T}
            push!(hexagons, e)
        elseif e isa ParaxialLens{T}
            push!(paraxials, e)
        elseif e isa Surface{T} || e isa CSGTree{T}
            push!(actual_elements, e)
        else
            # unpack the lens assembly
            rectangles = vcat(rectangles, e.rectangles)
            ellipses = vcat(ellipses, e.ellipses)
            hexagons = vcat(hexagons, e.hexagons)
            for lae in OpticSim.typed_elements(e)
                push!(actual_elements, lae)
            end
        end
    end
    N = length(actual_elements)
    # make the methods for this N
    type = :($(Symbol("MacroLensAssembly$(N)")))
    if !isdefined(OpticSim, type)
        @eval @lensassembly_constructor($N)
        @eval @lensassembly_intersection($N)
        @eval @lensassembly_elements($N)
    end
    eltypes = typeof.(actual_elements)
    naxis = normalize(axis)
    return eval(:($(type){$T,$(eltypes...)}($naxis, $rectangles, $ellipses, $hexagons, $paraxials, $(actual_elements...))))
end

function typed_elements(ass::MacroLensAssembly{T}) where {T<:Real}
    # under certain circumstances we run into the world age problem (see https://discourse.julialang.org/t/how-to-bypass-the-world-age-problem/7012)
    # to avoid this we need to use invokelatest(), but this is type ambiguous and ruins performance, and in most cases isn't necessary
    # so we can use try/catch to only call it when we need (and warn appropriately), retaining performance in almost all cases
    # NEVER call _typed_elements directly, always use this method
    try
        _typed_elements(ass)
    catch MethodError
        @warn "First time call to LensAssembly of previously unused size, performance will be worse than normal." maxlog = 1
        Base.invokelatest(_typed_elements, ass)
    end
end

function closestintersection(ass::MacroLensAssembly{T}, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
    # see typed_elements for why the try/catch is needed
    # NEVER call _closestintersection directly, always use this method
    try
        _closestintersection(ass, r)
    catch MethodError
        @warn "First time call to LensAssembly of previously unused size, performance will be worse than normal." maxlog = 1
        Base.invokelatest(_closestintersection, ass, r)
    end
end

elements(ass::MacroLensAssembly{T}) where {T<:Real} = (ass.rectangles..., ass.ellipses..., ass.hexagons..., ass.paraxials..., typed_elements(ass)...)
