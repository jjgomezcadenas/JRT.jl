using Unitful
using UnitfulEquivalences
using Roots

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M,
    Bq, mBq, μBq 


A_BI214_316Ti   =    1.0    * mBq/kg
A_TL208_316Ti   =    0.4   * mBq/kg
A_BI214_CU_LIM  =   12      * μBq/kg
A_TL208_CU_LIM  =    1.4    * μBq/kg

"""
Represents a cylinder 

"""
struct Cylinder
    r   :: Float64
    zmin:: Float64
    zmax:: Float64
    p0::Vector{Float64}
    p1::Vector{Float64}

    function Cylinder(r:: Float64, zmin::Float64, zmax::Float64)
        p0 = [0.0, 0.0, zmin] # point in one endcup
        p1 = [0.0, 0.0, zmax] # point in the other
        new(r, zmin, zmax, p0, p1)
    end
end

"""
A cylinder with units

"""
struct PhysicalCylinder
    cyl::Cylinder
    units::Unitful.Length
end

"""
Equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0

"""
cylinder_equation(c::Cylinder, P::Vector{Float64}) = P[1]^2 + P[2]^2 - c.r^2


"""
Length of cylinder

"""
clength(c::Cylinder) = c.zmax - c.zmin
clength(c::PhysicalCylinder) = clength(c.cyl) * c.units

"""
Perimeter of cylinder

"""
perimeter(c::Cylinder) = 2π * c.r
perimeter(c::PhysicalCylinder) = perimeter(c.cyl) * c.units


"""
Cylinder area (barrel)

"""
area_barrel(c::Cylinder) = perimeter(c) * clength(c)
area_barrel(c::PhysicalCylinder) = perimeter(c) * clength(c)


"""
Cylinder area (endcap)

"""
area_endcap(c::Cylinder) = π * c.r^2
area_endcap(c::PhysicalCylinder) = area_endcap(c.cyl) * c.units^2

"""
Cylinder area (total)

"""
area(c::Cylinder) = area_barrel(c) + area_endcap(c)
area(c::PhysicalCylinder) = area_barrel(c) + area_endcap(c)


"""
Cylinder volume 

"""
volume(c::Cylinder) = area_endcap(c) * clength(c)
volume(c::PhysicalCylinder) = area_endcap(c) * clength(c)


"""
Returns True if point is in end-caps of cyinder
 This means that the third coordinate of p (p[3], 
 is close to zmin or to zmax)

 """
function in_endcaps(c::Cylinder, p::Vector{Float64}) 
    isapprox(p[3], c.zmin, atol=1e-7) || isapprox(p[3], c.zmax, atol=1e-7)
end


"""
Returns True if point is in barrel of cyinder
 This means that the point verifies the equation of cylinder

 """
function in_barrel(c::Cylinder, p::Vector{Float64}) 
    isapprox(cylinder_equation(c, p), 0.0, atol= 1.0e-7)
end

"""
Returns true if the point is inside the cylinder

"""
function inside_cylinder(c::Cylinder, p::Vector{Float64})
	s1 =sqrt(p[1]^2 + p[2]^2) <= c.r
	s2 = c.zmin <= p[3] <= c.zmax
	s1 && s2
end


"""
Simple representation of a physical material

# Fields
- `ρ::typeof(1.0g/cm^3)`          : density (ρ)
- `μovrρ::typeof(1.0cm^2/g)`      : μ/ρ
- `μ::Unitful.Length^-1`          : attenuation coefficient (μ) 
- `λ::Unitful.Length`             : attenuation length 

"""
struct PhysicalMaterial
    ρ::typeof(1.0g/cm^3) 
    μovrρ::typeof(1.0cm^2/g) 
	μ::typeof(1.0/cm)
	λ::Unitful.Length

    function PhysicalMaterial(ρ::typeof(1.0g/cm^3),  μovrρ::typeof(1.0cm^2/g))
        μ = μovrρ*ρ
        λ = 1.0/μ
		new(ρ, μovrρ, μ, λ)
	end
end

"""
Simple representation of a radioactive material

# Fields
- `m::PhysicalMaterial`          : a radioactive material is a physical material
- `a_bi214::typeof(1.0*Bq/kg)`   : activity of Bi-214
- `a_tl208::typeof(1.0*Bq/kg)`   : activity of Tl-208

"""
struct RadioactiveMaterial
    m::PhysicalMaterial
    a_bi214::typeof(1.0*Bq/kg)
    a_tl208::typeof(1.0*Bq/kg)
    
end

"""
Mass of a cylinder filled with radioactive material mat
"""
function mass(c::PhysicalCylinder, mat::RadioactiveMaterial) 
    uconvert(g, volume(c) * mat.m.ρ)
end


"""
Activity (Bi214) of a cylinder filled with radioactive material mat
"""
function a_bi214(c::PhysicalCylinder, mat::RadioactiveMaterial) 
    uconvert(Bq, mat.a_bi214 * mass(c, mat))
end


"""
Activity (Tl208) of a cylinder filled with radioactive material mat
"""
function a_tl208(c::PhysicalCylinder, mat::RadioactiveMaterial) 
    uconvert(Bq, mat.a_tl208 * mass(c, mat))
end

"""
Fraction of gammas (events) that make it through a distance d
"""
function transmittance_at_qbb(rmt::RadioactiveMaterial, d::Float64) 
	l = uconvert(mm, rmt.m.λ)/mm
	exp(-d/l)
end

"""
Defines a ray (e.g, a segment that propagates in straight line)
r = e + t * d, 
where e in the initial point and d is the direction vector
"""
struct Ray
    e::Vector{Float64}
    d::Vector{Float64}
	#u::Vector{Float64}
	#function Ray(e::Vector{Float64}, d::Vector{Float64})
	#	u = d .- e
	#	new(e,d,u)
	#end
end

"""
Propagates ray r along the distance defined by t
"""
ray(r::Ray, t::Float64) =r.e + t * r.d

"""
Computes  the polynomial defined by the intersection roots 
between a ray and a cylinder. 
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0

"""
function cylinder_intersection_poly(r::Ray, cl::Cylinder)
    function gf()
		a = r.d[1]^2 + r.d[2]^2
		b = 2 * (r.e[1] * r.d[1] + r.e[2] * r.d[2])
		c = r.e[1]^2 + r.e[2]^2 - cl.r^2

		f(x) =  a * x^2 + b * x + c
	end
	gf()
end


"""
Computes intersection roots between a ray and a cylinder
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0

"""
function cylinder_intersection_roots(r::Ray, cl::Cylinder, eps::Float64=1e-9)
	
	f = cylinder_intersection_poly(r, cl)
	
	#println("a = ", -2*clength(cyl) * cyl.r, " b = ", 2*clength(cyl) * cyl.r )
	
	rts = find_zeros(f, -100*clength(cl) * cl.r, 100*clength(cl) * cl.r)
	#println("roots: rts = ", rts)
	if length(rts) == 0
		return 0
	else
		prts = [x for x in rts if x >eps]
		if length(prts) == 0
			return 0
		else
			return minimum(prts)
		end
	end
end

"""
Intersection between a ray and the end-cups of a cylinder

"""
function ray_intersection_with_cylinder_end_caps(r::Ray, cl::Cylinder, 
                                                 t::Float64)
 
    p = ray(r, t)
    if p[3] > cl.zmax
        t = (cl.zmax - r.e[3])/r.d[3]
    else
        t = (cl.zmin - r.e[3])/r.d[3]
	end

    t, ray(r,t)
end

"""
 Intersection between a ray and a cylinder
 
 """
function ray_intersection_with_cylinder(r::Ray, cl::Cylinder)
   
    t = cylinder_intersection_roots(r, cl)

	if t == 0 && r.d[3] >= 0.0  # intersects upper endcap
		return Ray([0.0,0.0,  cl.zmax], [0.0,0.0,1.0])
	elseif t == 0 && r.d[3] < 0.0  # intersects lower endcap
		return Ray([0.0,0.0, cl.zmin], [0.0, 0.0, -1.0])
	end
    
	P = ray(r, t)
    z = P[3]
    if z < cl.zmin || z > cl.zmax
        t, P = ray_intersection_with_cylinder_end_caps(r, cl, t)
		return Ray(P, r.d)
	else
		return Ray(P, r.d)
	end
end