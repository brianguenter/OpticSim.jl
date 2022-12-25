# See LICENSE in the project root for full license information.

using Test
using OpticSim
using OpticSim.Emitters
using OpticSim.Geometry
using OpticSim.Vis
using OpticSim.GlassCat
using OpticSim.Repeat

using StaticArrays

Examples_N_BK7 = GlassCat.Glass("$(@__MODULE__).Examples_N_BK7", 2, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653, 0.0, 0.0, NaN, NaN, 0.3, 2.5, 1.86e-6, 1.31e-8, -1.37e-11, 4.34e-7, 6.27e-10, 0.17, 20.0, -0.0009, 2.3, 1.0, 7.1, 1.0, 1, 1.0, [(0.3, 0.05, 25.0), (0.31, 0.25, 25.0), (0.32, 0.52, 25.0), (0.334, 0.78, 25.0), (0.35, 0.92, 25.0), (0.365, 0.971, 25.0), (0.37, 0.977, 25.0), (0.38, 0.983, 25.0), (0.39, 0.989, 25.0), (0.4, 0.992, 25.0), (0.405, 0.993, 25.0), (0.42, 0.993, 25.0), (0.436, 0.992, 25.0), (0.46, 0.993, 25.0), (0.5, 0.994, 25.0), (0.546, 0.996, 25.0), (0.58, 0.995, 25.0), (0.62, 0.994, 25.0), (0.66, 0.994, 25.0), (0.7, 0.996, 25.0), (1.06, 0.997, 25.0), (1.53, 0.98, 25.0), (1.97, 0.84, 25.0), (2.325, 0.56, 25.0), (2.5, 0.36, 25.0)], 1.5168, 2.3, 0.0, 0, 64.17, 0, 2.51, 0.0)

@testset "Interfacea" begin
    @testset "Fresnel" begin
        @testset "Diffusion" begin

        a = Array[[1500.0, 1500.0, 1.0], [2.0, 500.0, 0.0], [500.0, 300.0, 0.5], [500.0, 1500.0, 0.3]]
        #a = Array[[1500.0, 1500.0, 1.0]]
    for params in a

        dethalfsize = params[1]
        detdist = params[2]
        diff = params[3]
        println(dethalfsize, ", ", detdist, ", ",diff)

        interface = FresnelInterface{Float64}(Examples_N_BK7, Air; reflectance=1.0, transmission=0.0, diffusereflection=diff)
		
		bs_1 = OpticSim.transform(
                Cuboid(50.0, 100.0, 2.0, interface=interface),
                translation(0.0, 0.0, -detdist))
		
        detector = Ellipse(dethalfsize, dethalfsize, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0); interface = opaqueinterface())
        system = CSGOpticalSystem(LensAssembly(bs_1()), detector)
		
		src = Sources.Source(transform = translation(0.0,0.0, -1.0), origins = Origins.Point(),directions = Directions.UniformCone(-unitZ3(),0.0, 2000))

        pwr_src = 0 
        pwr_det = 0

        resetdetector!(system)
        for (i, r) in enumerate(src)
            pwr_src += r.sourcepower
            res = trace(system, r, test = true)
            if res !== nothing
                pwr_det += res.ray.power
            end
        end

        # Vis.drawtracerays(system, raygenerator=src ,trackallrays=true, test=true, numdivisions=100)
        a = dethalfsize
        b = detdist
        r = a/sin(atan(a/b))
        sectorsphere_area = (2*pi*r*(r-b))
        hemisphere_area = 2*pi*r^2

        rtr = pwr_det/pwr_src
        ar = diff*(sectorsphere_area/hemisphere_area)+(1.0-diff)
        println(rtr, "?~=", ar)

        @test isapprox(rtr, ar, rtol = 0.03, atol = 0.03)

    end

        end
    end
end