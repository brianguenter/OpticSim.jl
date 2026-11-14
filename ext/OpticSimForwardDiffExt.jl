# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Support for differentiating OpticSim with ForwardDiff. Loaded automatically when ForwardDiff is loaded.
module OpticSimForwardDiffExt

import OpticSim
import ForwardDiff
import LinearAlgebra

# strip the derivative part, recursively for nested duals (e.g. Hessians)
OpticSim.primalvalue(x::ForwardDiff.Dual) = OpticSim.primalvalue(ForwardDiff.value(x))

# The spot merit and its gradient from compiled sequential ray Jacobians (src/sequentialgradient.jl). The merit's
# gradient with respect to the lens parameters θ and to each ray's start comes from the Jacobians; the chain rule through
# x -> (θ, ray starts) is one forward mode pass through the prescription (and the aimed rays, when they depend on x),
# which traces no rays.
function OpticSim.spotmerit_and_gradient(lens::OpticSim.LensOptimizationProblem, x::AbstractVector{<:Real}; withgradient::Bool=true)
    lens.prescription === nothing && throw(ArgumentError("spotmerit_and_gradient needs a LensOptimizationProblem created with a `prescription` function"))
    x = Vector{Float64}(x)
    fieldrays = lens.aimedrays(x)
    merit, gθ, gstart = OpticSim.sequentialspotmerit(lens.prescription(x), fieldrays; withgradient)
    withgradient || return merit, nothing
    isnan(merit) && return merit, fill(NaN, length(x))
    wavelengths = unique(Float64(OpticSim.wavelength(r)) for rays in fieldrays for (r, _) in rays)
    aimed = fieldrays !== lens.rays
    function contracted(xd)
        s = LinearAlgebra.dot(gθ, OpticSim.sequentiallens(lens.prescription(xd), wavelengths).θ)
        if aimed
            k = 0
            for rays in lens.aimedrays(xd), (r, _) in rays
                k += 1
                s += LinearAlgebra.dot(view(gstart, 1:3, k), OpticSim.origin(r)) + LinearAlgebra.dot(view(gstart, 4:6, k), OpticSim.direction(r))
            end
        end
        return s
    end
    config = ForwardDiff.GradientConfig(contracted, x, ForwardDiff.Chunk(length(x))) # all derivatives in one pass
    return merit, ForwardDiff.gradient(contracted, x, config)
end

end # module OpticSimForwardDiffExt
