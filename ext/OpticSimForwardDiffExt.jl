# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Support for differentiating OpticSim with ForwardDiff. Loaded automatically when ForwardDiff is loaded.
module OpticSimForwardDiffExt

import OpticSim
import ForwardDiff

# strip the derivative part, recursively for nested duals (e.g. Hessians)
OpticSim.primalvalue(x::ForwardDiff.Dual) = OpticSim.primalvalue(ForwardDiff.value(x))

end # module OpticSimForwardDiffExt
