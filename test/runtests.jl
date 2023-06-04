using JRT
using Test

using Unitful
using UnitfulEquivalences
import GeometryBasics

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M,
    Bq, mBq, μBq 


import JRT: 
A_BI214_316Ti, A_TL208_316Ti, A_BI214_CU_LIM, A_TL208_CU_LIM,
PhysicalMaterial, RadioactiveMaterial,  Cylinder, PhysicalCylinder,
Ray, ray, 
clength, perimeter, area_barrel, area_endcap, volume, 
mass, a_bi214, a_tl208,in_endcaps, in_barrel, inside_cylinder,
transmittance_at_qbb, ray_intersection_with_cylinder, cylinder_equation,
cylinder_intersection_roots, cylinder_intersection_poly,
get_cylinder



                           
pti316  = PhysicalMaterial(7.87 * g/cm^3, 0.039 * cm^2/g)
ti316 = RadioactiveMaterial(pti316, A_BI214_316Ti, A_TL208_316Ti)

pcu  = PhysicalMaterial(8.96 * g/cm^3, 0.039 * cm^2/g)
cu12 = RadioactiveMaterial(pcu, A_BI214_CU_LIM, A_TL208_CU_LIM)


rmdct = Dict("ti316" => ti316,
             "cu12" => cu12)


cyl =Cylinder(5.0, 0., 10.0)
gcyl2, cylm = get_cylinder(cyl)

gcyl = GeometryBasics.Cylinder(GeometryBasics.Point{3, Float64}(1,2,3), 
                               GeometryBasics.Point{3, Float64}(2,3,4), 1.0)

#Ray moving along z axis, with no x,y components. 
#It must intersect endcups (either zmax or zmin)
rzp = Ray([1.0, 0.0, 0.0], [0.0, 0.0, 1.0])
rxp = ray_intersection_with_cylinder(rzp, cyl)

rzn = Ray([1.0, 0.0, 0.0], [0.0, 0.0, -1.0])
rxn = ray_intersection_with_cylinder(rzn, cyl)

ry = Ray([1.0, 0.0, 0.0], [0.5, 0.5, 1.0])
#ff = cylinder_intersection_poly(ry, cyl)
rxy = ray_intersection_with_cylinder(ry, cyl)

rz = Ray([1.0, 0.0, 0.0], [0., 0.1, 1.0])
rxz = ray_intersection_with_cylinder(rz, cyl)

# Ray intersecting the cylinder surface much further than zmax


@testset "JRT.jl" begin
    @test gcyl.origin == [1.0, 2.0, 3.0]
    @test gcyl.extremity == [2.0, 3.0, 4.0]
    @test gcyl.r == 1.0
     
    @test cyl.r == 5.0
    @test cyl.zmin == 0.0
    @test cyl.zmax == 10.0
    @test cyl.p0 == [0.0, 0.0, 0.0]
    @test cyl.p1 == [0.0, 0.0, 10.0]

    @test gcyl2.origin == cyl.p0
    @test gcyl2.extremity == cyl.p1
    @test gcyl2.r == cyl.r

    @test in_endcaps(cyl, [0.0, 0.0, 5.0]) == false
    @test in_endcaps(cyl, [0.0, 0.0, 0.0] )
    @test in_endcaps(cyl, [0.0, 0.0, 10.0] )

    @test rzp.e == [1.0, 0.0, 0.0]
    @test rzp.d == [0.0, 0.0, 1.0]

    @test rzn.e == [1.0, 0.0, 0.0]
    @test rzn.d == [0.0, 0.0, -1.0]

    #cylinder_intersection_roots returns the smaller positive root c
    #orresponding to the intersection between the ray and the cylinder. 
    #Since in this case the ray does not intersect the cylinder 
    #(the endcaps are not the cylinder surface) the root must be zero 
    @test cylinder_intersection_roots(rzp, cyl) == 0.0
    @test cylinder_intersection_roots(rzn, cyl) == 0.0
    @test cylinder_intersection_roots(ry, cyl) == 6.0

    #ray_intersection_with_cylinder will propagate to
    #cylinder and return the recomputed ray. In this case it must
    #be the top endcup of cylinder
    @test rxp.e == [0.0, 0.0, 10.0]
    @test rxp.d == [0.0, 0.0, 1.0]

    @test rxn.e == [0.0, 0.0, 0.0]
    @test rxn.d == [0.0, 0.0, -1.0]

    @test in_endcaps(cyl, rxp.e)
    @test in_endcaps(cyl, rxn.e)

    @test in_endcaps(cyl, rxy.e) == false 
    @test cylinder_equation(cyl, rxy.e) == 0.0
    @test in_barrel(cyl, rxy.e)

    @test in_endcaps(cyl, rxz.e)

    for (name, mat) in rmdct
         @test mat.m.μ ≈ mat.m.μovrρ * mat.m.ρ
     end

    # @test volume(cyl) ≈ π*mm^3 
    # @test surface(cyl) ≈ 2π*mm^2
    # @test endcap_surface(cyl) ≈ π*mm^2 
    # @test inner_volume(cyls) ≈ π*mm^3
    # @test shell_volume(cyls) ≈ 3π * mm^3 
    # @test inner_surface(cyls) ≈ 2π*mm^2  
    # @test inner_endcap_surface(cyls) ≈ π*mm^2
    # @test outer_surface(cyls) ≈ 4π * mm^2     
    # @test outer_endcap_surface(cyls) ≈ 4π * mm^2
    # @test thickness(cyls) ≈ 1.0mm
end
