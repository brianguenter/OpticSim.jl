# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using OpticSim: area, Rectangle, ParaxialLensRect, virtualpoint
using Unitful.DefaultSymbols

@testitem "ParaxialAnalysis" begin
    @testset "Lens" begin
        lens = ParaxialLensRect(10.0, 100.0, 100.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        displaypoint = [0.0, 0.0, -8.0]
        r1 = Ray(displaypoint, [1.0, 0.0, 1.0])
        intsct = OpticSim.surfaceintersection(lens, r1)
        intsctpt = point(intsct)
        # compute the refracted ray and project this backwards to its intersection with the optical axis. This should be the same position as returned by virtualpoint()
        refrac, _, _ = OpticSim.processintersection(OpticSim.interface(lens), OpticSim.point(intsct), OpticSim.normal(lens), OpticalRay(r1, 1.0, 0.55), OpticSim.TEMP_REF, OpticSim.PRESSURE_REF, false)
        # compute intersection of ray with optical axis. this should match the position of the virtual point
        slope = refrac[1] / refrac[3]
        virtptfromslope = [0.0, 0.0, -intsctpt[1] / slope]
        @test isapprox(virtptfromslope, point(OpticSim.virtualpoint(lens, displaypoint)))
    end



    @testset "SphericalPolygon" begin
        using StaticArrays

        """creates a circular polygon that subtends a half angle of θ"""
        function sphericalcircle(θ, nsides=10)
            temp = MMatrix{3,nsides,Float64}(undef)
            for i in 0:1:(nsides-1)
                ϕ = i * 2π / nsides
                temp[1, i+1] = sin(θ) * cos(ϕ)
                temp[2, i+1] = cos(θ)
                temp[3, i+1] = sin(θ) * sin(ϕ)
            end
            return SphericalPolygon(SMatrix{3,nsides,Float64}(temp), SVector(0.0, 0.0, 0.0), 1.0)
        end

        oneeigthsphere() = SphericalTriangle(SMatrix{3,3,Float64}(
                0.0, 1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0),
            SVector(0.0, 0.0, 0.0),
            1.0)

        onesixteenthphere() = SphericalTriangle(SMatrix{3,3,Float64}(
                0.0, 1.0, 0.0,
                1.0, 1.0, 0.0,
                0.0, 0.0, 1.0),
            SVector(0.0, 0.0, 0.0),
            1.0)

        onesixteenthspherepoly() = SphericalPolygon(SMatrix{3,3,Float64}(
                0.0, 1.0, 0.0,
                1.0, 1.0, 0.0,
                0.0, 0.0, 1.0),
            SVector(0.0, 0.0, 0.0),
            1.0)

        threesidedpoly() = SphericalPolygon(SMatrix{3,3,Float64}(
                0.0, 1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0),
            SVector(0.0, 0.0, 0.0),
            1.0)


        foursidedpoly() = SphericalPolygon(SMatrix{3,4,Float64}(
                0.0, 1.0, 0.0,
                1.0, 1.0, -0.9,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0),
            SVector(0.0, 0.0, 0.0),
            1.0)

        #test that spherical triangle and spherical polygon area algorithms return the same values for the same spherical areas
        @test isapprox(OpticSim.area(oneeigthsphere()), OpticSim.area(threesidedpoly()))
        @test isapprox(OpticSim.area(onesixteenthphere()), OpticSim.area(onesixteenthspherepoly()))

        #halve the spherical area and make sure area function computes 1/2 the area
        @test isapprox(OpticSim.area(onesixteenthphere()), OpticSim.area(oneeigthsphere()) / 2.0)
        @test isapprox(OpticSim.area(sphericalcircle(π / 8.0, 1000)) / 4.0, OpticSim.area(sphericalcircle(π / 16.0, 1000)), atol=1e-2)
    end
end

