# Optical Interfaces

Every [`Surface`](@ref) must have an `OpticalInterface` associated with it to defined the behavior of any ray when it intersects that surface.

```@docs; canonical = false
OpticSim.OpticalInterface
OpticSim.NullInterface
FresnelInterface
ParaxialInterface
ThinGratingInterface
HologramInterface
MultiHologramInterface
```

The critical behavior of each interface is defined in the `processintersection` function:

```@docs; canonical = false
OpticSim.processintersection
```
