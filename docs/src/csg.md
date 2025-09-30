# CSG

## CSG Operations

OpticSim represents objects as boolean combinations of halfspaces, a method of constructing solid objects known as Constructive Solid Geometry (CSG). A halfspace is a surface which divides 3D space into 3 regions: inside the halfspace,outside the halfspace, on the surface dividing inside from outside. 

A simple example of a halfspace is an infinite oriented plane, where oriented means we have defined a direction for the normal to the plane surface. Another simple example is a sphere, or a cylinder.

If the halfspaces have mathematically well difined inside and outside functions then except for edge cases, which are unlikely to occur in an optics setting, then any boolean combination of halfspaces will also have a well defined inside and outside.

You can define optical elements by combining the basic halfspace types defined for you in OpticSim:
* sphere
* infinite cylinder
* plane
* prism
* spherical cap

Many of the surfaces already defined for you in OpticSim, such as sphere, prism,plane, cylinder, describe well defined solid objects, i.e., there is a mathematically defined inside and outside. 

Unfortunately, many of the surfaces common in optics do not by themselves define solid objects. None of the parametric surfaces, which includes all the optical asphere types, define solid objects. 

Parametric surfaces do have a positive and negative side though determined by the direction of the surface normal. You can create a well defined solid object by taking the intersection of a cylinder with an asphere.

There are three binary csg operations which can construct extremely complex objects from very simple primitives: union (``\cup``), intersection (``\cap``) and subtraction (i.e. difference).

This diagram shows the basic idea:
![CSG Tree visualization](https://upload.wikimedia.org/wikipedia/commons/8/8b/Csg_tree.png)

The code for this in our system would look this this:

```@example
using OpticSim # hide
cyl = Cylinder(0.7)
cyl_cross = cyl ∪ leaf(cyl, Geometry.rotationd(90, 0, 0)) ∪ leaf(cyl, Geometry.rotationd(0, 90, 0))

cube = Cuboid(1.0, 1.0, 1.0)
sph = Sphere(1.3)
rounded_cube = cube ∩ sph

result = rounded_cube - cyl_cross
```



```@docs; canonical = false
leaf
∪
∩
-
```

## Pre-made CSG Shapes

There are also some shortcut methods available to create common CSG objects more easily:

```@docs; canonical = false
BoundedCylinder
Cuboid
HexagonalPrism
RectangularPrism
TriangularPrism
Spider
```

## CSG Types

These are the types of the primary CSG elements, i.e. the nodes in the CSG tree.

```@docs; canonical = false
OpticSim.CSGTree
OpticSim.CSGGenerator
OpticSim.ComplementNode
OpticSim.UnionNode
OpticSim.IntersectionNode
OpticSim.LeafNode
```

## Additional Functions and Types

These are the internal types and functions used for geometric/CSG operations.

### Functions

```@docs; canonical = false
surfaceintersection(::CSGTree{T}, ::AbstractRay{T,N}) where {T<:Real,N}
inside(a::CSGTree{T}, x::T, y::T, z::T) where {T<:Real}
onsurface(a::CSGTree{T}, x::T, y::T, z::T) where {T<:Real}
```

### Intervals

```@docs; canonical = false
Interval
EmptyInterval
DisjointUnion
OpticSim.isemptyinterval
OpticSim.ispositivehalfspace
OpticSim.israyorigininterval
OpticSim.halfspaceintersection
OpticSim.closestintersection
OpticSim.IntervalPool
```

### Intersections

```@docs; canonical = false
OpticSim.IntervalPoint
RayOrigin
Infinity
Intersection
OpticSim.isinfinity
OpticSim.israyorigin
```
