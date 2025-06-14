# Lenses and Other Optical Components

## Lenses

A number of helper functions are provided to make constructing simple lenses easier.
Firstly ordinary thick lenses:

```@docs; canonical = false
SphericalLens
ConicLens
AsphericLens
FresnelLens
```

As well as idealized lenses:

```@docs; canonical = false
ParaxialLens
```

## Other Components

We also have some holographic elements implemented, note that these have not been extensively tested and should not be treated as wholely accurate at this stage.

It is relatively simple to extend the existing code to add these kinds of specialized surfaces providing a paired [`OpticalInterface`](@ref) subclass is also defined. In this case the `WrapperSurface` can often serve as a suitable base for extension.

```@docs; canonical = false
WrapperSurface
ThinGratingSurface
HologramSurface
MultiHologramSurface
```

## Eye Models

Eye models are often very useful in simulation of head mounted display systems. We have two models implemented currently.

```@docs; canonical = false
OpticSim.Data.ModelEye
OpticSim.Data.ArizonaEye
```
