### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ a298b03b-c7c4-4114-9984-d37e7df8d4c8
using Pkg; Pkg.activate("/Users/jjgomezcadenas/Projects/JRT")

# ╔═╡ 5158389a-fca1-11ed-05a1-c7345c9bd354
begin
	#using Revise
	using PlutoUI
	using CSV
	using DataFrames
	using DataFramesMeta
	using Plots
	#using Printf
	using Markdown
	using InteractiveUtils
	using Statistics
	using Random
	using Distributions
	using PDMats
	#using Chain
	using LinearAlgebra
	using Images
	using GLMakie
	#using Makie
	import GeometryBasics
	using Colors
	using Roots
end

# ╔═╡ 7886beee-8947-4cb2-bcd9-89728829b7cd
using SparseArrays

# ╔═╡ c8cec5b9-abb9-443d-9e41-bc8b7070792d


# ╔═╡ e36751ff-a521-4529-95c8-6bbfc3314e66
plotlyjs()

# ╔═╡ df3dde7a-9437-409f-bd51-db975f74b78c
PlutoUI.TableOfContents(title="JRT draft", indent=true)


# ╔═╡ fc869138-1cf7-463b-aca5-729741f8b1ff
GLMakie.activate!()

# ╔═╡ f230deec-f4c9-429a-89f5-92d8c08e8472
md"""
### Drawing a cylinder with Makie
"""

# ╔═╡ 7371f833-a7a8-482c-925b-730d8378c22e
gcyl = GeometryBasics.Cylinder(GeometryBasics.Point{3, Float64}(1,2,3), GeometryBasics.Point{3, Float64}(2,3,4), 1.0)

# ╔═╡ 59ab61e2-00ff-4a06-b732-04f7d65f7325
#cyL = GeometryBasics.mesh(gcyl)

# ╔═╡ 2c803aef-9a52-4108-a252-e63a95929c67
cyL =GeometryBasics.mesh(GeometryBasics.Tesselation(gcyl, 100))

# ╔═╡ c59e3a3d-1821-49a1-8d93-29250146330e
with_theme(theme_light()) do
    fig = Figure(resolution = (1200,800))
	axs = [Axis3(fig[1, i]; aspect=:data, perspectiveness=0.1) for i=1:2]
    mesh!(axs[1], cyL, color = (:dodgerblue, 0.55), transparency = true, shading = false)
    GLMakie.wireframe!(axs[2], cyL; color = :blue, linewidth = 2.5)
    fig
end

# ╔═╡ b853da86-3242-4a9f-8163-1e5bd0d6a907


# ╔═╡ c7583c88-4cfc-49e1-b61e-c299801be450
md"""
## Drawing a cylinder
- A cylinder of radius r, and height h, having z-axis as symmetry axis, is a stack of circles of radius r:

$\begin{align}
x= r \times \cos(u) \\
y=r \times sin(u) \\
z=v 
\end{align}$

where:
$u \in [0,2\pi]$
$v \in [0,h]$

- This allows to plot the cylinder as any parameterized surface.
- For example:
	- Define u as a range beteen 0 and 2π
	- Define v as a range between 0 and 10.0 (the length of the cylinder)
"""

# ╔═╡ 75472827-851c-451f-9510-a7e33fee1a8e
md"""
## Defining and drawing a cylinder
"""

# ╔═╡ b6082524-84ae-4d98-b593-7250f554c0c3
function draw_cube(np=100)
	r = LinRange(-1, 1, np)
	cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
	#GLMakie.contour(cube, alpha=0.5)
	cube_with_holes = cube .* (cube .> 1.4)
	#GLMakie.volume(cube_with_holes, algorithm = :iso, isorange = 0.05, isovalue = 1.7)

	GLMakie.volume(cube, transparency=true)
end

# ╔═╡ 6795ad55-06b8-4a77-9a3a-cf3c15a3741e
draw_cube()

# ╔═╡ c48c5f53-2632-4cfc-b07b-7fde9f1d8eb6
md"""
## Testing functions
"""

# ╔═╡ de0fb8fd-c757-4769-a2da-d874add0a553
md"""
### Check function **in_endcaps(c::Cylinder, p::Vector{Float64})**
"""

# ╔═╡ c95dcf93-40ff-4db0-a323-29991f5135dd
md"""
### Demonstrate propagation functions
"""

# ╔═╡ c7200fa5-31eb-458c-80ca-70f968edd9b0
md"""
#### Ray moving along z axis, with no x,y components. It must intersect endcups (either zmax or zmin)
"""

# ╔═╡ a89c0c98-ee2c-4103-8cd3-9be03a2b75ae
md"""
- Function **cylinder\_intersection\_roots** returns the smaller positive root corresponding to the intersection between the ray and the cylinder. Since in this case the ray does not intersect the cylinder (the endcups are not the cylinder surface) the root must be zero 
"""

# ╔═╡ aa445ed2-dc2d-4802-ab12-18505487ceef
md"""
- Function **ray\_intersection\_with_cylinder** will return the recomputed ray
"""

# ╔═╡ 6c2b1465-02ae-4768-a4f6-e8f461339b97
md"""
- Since this is a ray going in the zp direction, t=0, and the recomputed e point in the ray is in zmax and belongs to the upper end-cap. 
"""

# ╔═╡ 26d72ea7-31e5-4c55-af55-7666f7c47f1c
md"""
- We can now check that the point **inendcaps** is indeed in the end-caps
"""

# ╔═╡ b9c450a6-446e-4a5c-93ee-28a45cf96356
md"""
- A ray pointing into negative z should end up in the bottom end-cap
"""

# ╔═╡ 59138131-46c1-4ac8-be0b-42f5ca2e30b8
md"""
#### Ray intersecting the cylinder surface
"""

# ╔═╡ 02a7fd08-440c-4e0c-93c8-d42792642d11
md"""
- The ray is moving with large director cosines so it will intersect the cylinder surface.
"""

# ╔═╡ bca3de2e-74a6-423e-a660-3845c3e9657b
md"""
- We can visualize where the roots will be
"""

# ╔═╡ c97578f2-1d39-497a-9e79-31e041e77b7b
md"""
#### Ray intersecting the cylinder surface much further than zmax
"""

# ╔═╡ 68e9ab79-21a4-48c1-82c3-45d8862e70cf
md"""
### Demonstrate generation of photons
"""

# ╔═╡ ff45920d-225a-48df-a9e3-f567b81bd579


# ╔═╡ b4110420-ef1c-4e1b-868a-68bd9d943fea
md"""
- Points must be normalized
"""

# ╔═╡ 5e866260-a798-4c5f-b407-38d091d0ed07
md"""
- Points must be in a unit sphere
"""

# ╔═╡ be8dd99e-60f5-4d75-959c-cffbd5bfa4b8
md"""
### Demonstrate generation around one point and propagation to cylinder
"""

# ╔═╡ 3b43b3d6-2621-4940-8f6e-cf5a5c6a2093
begin
	ax = 0.5:0.01:1.50 
	ay = 0.5:0.01:1.50 
	af(ax,ay) = abs2(ax-ay) 
	#ap = plot(ax,ay,af,st=:surface) 
	plot(ax,ay,af,st=:wireframe)
end

# ╔═╡ f675ebdf-e3f6-4f6b-bdae-e59db854f63c
indx=4

# ╔═╡ 3dce337c-3f5b-4733-870d-15cf8197525b
function draw_rays(e, gtdf)
    tt = collect(range(0.0, 1.0, 100))
	frst = true
	
	x1 = e[1] .+ tt * gtdf[1,"gdx"]
	y1 = e[2] .+ tt * gtdf[1,"gdy"]
	z1 = e[3] .+ tt * gtdf[1,"gdz"]

	p = plot(x1, y1, z1, legend=false)
    for r in eachrow(gtdf)
		xi = e[1] .+ tt * r.gdx
		yi = e[2] .+ tt * r.gdy
		zi = e[3] .+ tt * r.gdz
		if frst
			p = plot(xi, yi, zi)
			frst = false
		else
			p = plot!(p, xi, yi, zi,  legend=false)
		end
	end
	p
end

# ╔═╡ ce9da661-5f2f-4e5d-825e-5d029c274066
md"""
# Types
"""

# ╔═╡ efe04df4-cb56-4916-96f7-296dc37492ac
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
        #mag, v, n1, n2 = unit_vectors_()
        #P, P2, P3 = surfaces_(r, mag, v, n1, n2)
	end
end

# ╔═╡ 3c398d08-ef94-4202-9b87-62aa4fba07ed
cyl =Cylinder(5.0, 0., 10.0)

# ╔═╡ 50f7c64a-4fa6-4705-b2e9-896ef545029f
cyl

# ╔═╡ c7d4775a-f585-4852-a172-337655cb21ed
stdsk = Cylinder(1000.0, 0., 50.0)

# ╔═╡ 053952a1-0558-447c-a3f5-ae30c2d2f7be
"""
Defines a ray (e.g, a segment that propagates in straight line)
r = e + t * d, 
where e in the initial point and d is the direction vector
"""
struct Ray
    e::Vector{Float64}
    d::Vector{Float64}
	u::Vector{Float64}
	function Ray(e::Vector{Float64}, d::Vector{Float64})
		u = d .- e
		new(e,d,u)
	end
end


# ╔═╡ 9e193e81-0147-46d2-943f-1526dfb26e50
rzp = Ray([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]) # starts in (1,0,0), moves towards z+

# ╔═╡ 5de835fa-284f-4dd0-ba1d-033d4c50da5a
rzn = Ray([1.0, 0.0, 0.0], [0.0, 0.0, -1.0]) # starts in (1,0,0), moves towards z+

# ╔═╡ 40cff386-97de-4e6c-8d22-216302831893
ry = Ray([1.0, 0.0, 0.0], [0.5, 0.5, 1.0])

# ╔═╡ 4d17bf89-d01b-4f1b-9dc7-9dcd00ee1d8a
rz = Ray([1.0, 0.0, 0.0], [0., 0.1, 1.0])

# ╔═╡ 2facda59-dc5a-4951-8d6f-1d121690aad2
"""
Propagates ray r along the distance defined by t
"""
ray(r::Ray, t::Float64) =r.e + t * r.d


# ╔═╡ 991a026a-bfe2-4bc8-a200-31dd25a881d1
md"""
# Functions
"""

# ╔═╡ 539faa01-299b-4c5b-8113-fc21501dd84d
"""
Takes a Cylinder type and returns a GeometryBasics cylinder and mesh. 
"""
function get_cylinder(cl::Cylinder, np=100)
	p0 = GeometryBasics.Point{3, Float64}(cl.p0[1],cl.p0[2],cl.p0[3]) 
	p1 = GeometryBasics.Point{3, Float64}(cl.p1[1],cl.p1[2],cl.p1[3]) 
	gcyl = GeometryBasics.Cylinder(p0, p1, cl.r)
	cyL =GeometryBasics.mesh(GeometryBasics.Tesselation(gcyl, np))
	gcyl, cyL
end

# ╔═╡ ef077643-f5b1-4244-af7e-8111ede962a0
gcyl2, cylm = get_cylinder(cyl)

# ╔═╡ 63c8a7db-740a-4abf-981f-3c02be41493a
function dmcyl(cyl, np=100,wfrm=true)
	gcyl2, cylm = get_cylinder(cyl,np)
	with_theme(theme_light()) do
    fig = Figure(resolution = (1200,800))
	axs = [Axis3(fig[1, i]; aspect=:data, perspectiveness=0.1) for i=1:2]
    mesh!(axs[1], cylm, color = (:dodgerblue, 0.55), 
		 transparency = true, shading = false)
    GLMakie.wireframe!(axs[2], cylm; color = :blue, linewidth = 2.5)
    fig
	end
end

# ╔═╡ 911252fb-0b59-4d23-a90e-4edd293ea91f
function dmcyl1(cyl, np=100,wfm=true)
	gcyl2, cylm = get_cylinder(cyl,np)
	with_theme(theme_light()) do
	    fig = Figure(resolution = (1200,800))
		axs = [Axis3(fig[1, i]; aspect=:data, perspectiveness=0.1) for i=1:2]
	    if wfm
			mesh!(axs[1], cylm, color = (:dodgerblue, 0.55), 
			 	transparency = true, shading = false)
		else
	    	GLMakie.wireframe!(axs[1], cylm; color = :blue, linewidth = 2.5)
		end
	    fig
	end
end

# ╔═╡ 269e4b09-3505-40fe-8815-2bc847d02a99
"""
Takes a cylinder and returns a Makie drawing, either frame or wireframe
"""
function draw_cylinder(cl::Cylinder, np=100;  resolution = (1000,1000), 
                       col = :blue, clevel=0.5, linewidth = 2.5, wf=true)
	gcyl2, cylm = get_cylinder(cyl,np)

	with_theme(theme_light()) do
		fig = Figure(resolution = (1200,800))
		axs = Axis3(fig[1,1]; aspect=:data, perspectiveness=0.5) 
		if wf == false
    		GLMakie.mesh!(axs, cylm, color = (col, clevel), 
			  transparency = true, shading = false)
		else
    		GLMakie.wireframe!(axs, cylm; color = col, linewidth = linewidth)
		end
		fig
	end
end

# ╔═╡ 81d62300-1a88-4c9c-b050-8acd20857867
fig = draw_cylinder(cyl, wf=true)

# ╔═╡ c26c8511-afae-4447-9458-d1d21f565911
"""
Equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
cylinder_equation(c::Cylinder, P::Vector{Float64}) = P[1]^2 + P[2]^2 - c.r^2

# ╔═╡ deefebf9-07a3-4cce-b2d5-4b42644b5db7
"""
Length of cylinder
"""
clength(c::Cylinder) = c.zmax - c.zmin



# ╔═╡ e9d59365-ea79-4e97-b908-cbd528dc0804
"""
Perimeter of cylinder
"""
perimeter(c::Cylinder) = 2π * c.r



# ╔═╡ da14195e-414e-40ab-aca1-83d2d68b9172
"""
Cylinder area (barrel)
"""
area_barrel(c::Cylinder) = perimeter(c) * length(c)


# ╔═╡ f41fbfdf-4cee-4b82-b087-e224c70b7f6b
"""
Cylinder area (endcap)
"""
area_endcap(c::Cylinder) = π * c.r^2
        


# ╔═╡ 7f8c81be-0f2d-4214-b27d-a7470534d867
"""
Cylinder area (total)
"""
area(c::Cylinder) = area_barrel(c) + area_endcap(c)
       


# ╔═╡ 372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
"""
Cylinder volume 
"""
volume(c::Cylinder) = area_endcap(c) * length(c)


# ╔═╡ 0abeba90-ba75-4231-9dd8-d35400964c60
 """Returns True if point is in end-caps of cyinder
 This means that the third coordinate of p (p[3], is close to zmin or to zmax)
 """
in_endcaps(c::Cylinder, p::Vector{Float64}) = p[3] ≈ c.zmin || p[3] ≈ c.zmax

# ╔═╡ d06eaf38-f8bb-4732-9b8f-92f2a1053176
in_endcaps(cyl, [0.0, 0.0, 5.0]) # not in endcup

# ╔═╡ 2940f9d9-f2c3-4772-8208-2375370b7e61
in_endcaps(cyl, [0.0, 0.0, 0.0] )  # in bottom endcup

# ╔═╡ a7378f23-9b56-4a5f-8c4b-b51545e69633
in_endcaps(cyl, [0.0, 0.0, 10.0] )  # in top endcup

# ╔═╡ 464a3cf9-c030-4270-bae7-bf4339ed3fc8
 """Returns True if point is in barrel of cyinder
 This means that the point verifies the equation of cylinder
 """
in_barrel(c::Cylinder, p::Vector{Float64}) = cylinder_equation(c, p) ≈ 0.0

# ╔═╡ 58491933-49ec-4b70-ab1f-2301657ed9e8
"""
Returns true if the point is inside the cylinder
"""
function inside_cylinder(c::Cylinder, p::Vector{Float64})
	s1 =sqrt(p[1]^2 + p[2]^2) <= c.r
	s2 = c.zmin <= p[3] <= c.zmax
	s1 && s2
end

# ╔═╡ 4d53eaca-0cad-4f34-8205-216e9c0a16d9
md"""
To compute the intersection roots we write the ray as:

$r = \mathbf{e} + t \mathbf{d}$

and substitute coordinate by coordinate in the equation of the cylinder:

$F(x,y,z) = x^2 + y^2 - r^2 = 0$

thus:

$(e_1 + t d_1)^2 + (e_2 + t d_2)^2 -r^2 =0$

Then re-arrange

$(d_1^2 + d_2^2)t^2 + 2 t (e_1 d_1 + e_2 d_2) + e_1^2 + e_2^2 - r^2 =0$

So this is a polynomial of the form

$a t^2 + b t + c =0$

with:

$a = (d_1^2 + d_2^2)$

$b = 2 (e_1 d_1 + e_2 d_2)$

$c = e_1^2 + e_2^2 - r^2$

"""

# ╔═╡ 989cec01-548a-40e8-9a82-149dbfa7d367
"""
Computes  the polynomial defined by the intersection roots between a ray and a cylinder. Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
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

# ╔═╡ 7f4ab967-b5c2-46d9-90c4-ebbbdab0dbfe
ff = cylinder_intersection_poly(ry, cyl)

# ╔═╡ 3b4ea019-bba9-4ac1-b61a-fa251c9524d9
Plots.plot(-10:0.01:10, ff)

# ╔═╡ c2d3e54d-3ad5-46b6-9d28-be8acf1086b2
find_zeros(ff, -10.0, 10.0)  #finds two roots

# ╔═╡ c6f63bc1-6d4d-49b3-988f-ccaaaae562e1
"""
Computes intersection roots between a ray and a cylinder
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
function cylinder_intersection_roots(r::Ray, cl::Cylinder, eps::Float64=1e-9)
	
	f = cylinder_intersection_poly(r, cl)
	
	#println("a = ", -2*clength(cyl) * cyl.r, " b = ", 2*clength(cyl) * cyl.r )
	
	rts = find_zeros(f, -100*clength(cyl) * cyl.r, 100*clength(cyl) * cyl.r)
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

# ╔═╡ da50ac8a-8381-490b-a000-b6ba8b392888
cylinder_intersection_roots(rzp, cyl)

# ╔═╡ 71203991-d5c4-439f-8dd1-fe1b4cf3b875
cylinder_intersection_roots(rzn, cyl)

# ╔═╡ 96af216d-c0e3-4825-a11f-3aaa04c42d0b
cylinder_intersection_roots(ry, cyl) # we keep only the positive root

# ╔═╡ d137c26c-bf75-4899-ac0a-eb8e70c792bb
cylinder_intersection_roots(rz, cyl) 

# ╔═╡ c40ad59c-7687-4bba-b923-fa498295821c
"""
Intersection between a ray and the end-cups of a cylinder
"""
function ray_intersection_with_cylinder_end_caps(r::Ray, cl::Cylinder, t::Float64)
 
    p = ray(r, t)
    if p[3] > cl.zmax
        t = (cl.zmax - r.e[3])/r.d[3]
    else
        t = (cl.zmin - r.e[3])/r.d[3]
	end

    t, ray(r,t)
end

# ╔═╡ 8547b5f7-6b8a-4a9b-b0bf-5b1d8b08c072
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

# ╔═╡ 24cd55d7-6752-4e48-a3fa-7951cbff1b36
 rxp = ray_intersection_with_cylinder(rzp, cyl)

# ╔═╡ c3abda05-42ec-4759-bea4-cdf7ee92ed70
in_endcaps(cyl, rxp.e)

# ╔═╡ 366cf831-3558-44f2-bf09-6d9de74bf04d
 rxn = ray_intersection_with_cylinder(rzn, cyl)

# ╔═╡ 95a82d89-7b7f-4ce9-8f21-b527543c57f7
in_endcaps(cyl, rxn.e)

# ╔═╡ 21aeb736-8655-4747-bfa7-1ad9b5395fc3
 rxy = ray_intersection_with_cylinder(ry, cyl)

# ╔═╡ 5dcb797a-d7ef-463b-a629-336a8752c4ed
in_endcaps(cyl, rxy.e)

# ╔═╡ 92ef9d4b-ebf7-4247-a13e-aebe60cd79ad
cylinder_equation(cyl, rxy.e)

# ╔═╡ 80b17191-b750-4f31-9577-1f5214440429
in_barrel(cyl, rxy.e)

# ╔═╡ a143fe7e-e647-422e-8e0b-65099ff14d88
 rxz = ray_intersection_with_cylinder(rz, cyl)

# ╔═╡ 0dc2ca61-8ca9-43f2-ba71-87dcc2be08ec
in_endcaps(cyl, rxz.e)

# ╔═╡ 312154b9-aa06-4413-bcd7-7688fe168660
in_barrel(cyl, rxz.e)

# ╔═╡ 14583522-a4d2-4f37-9c8b-5e2f54ace834
"""
Generate three standard normally distributed numbers and normalize the vector.

We have to be careful in the case that the vector has a norm close to zero, in which we must worry about floating point precision by dividing by a very small number. This is the reason for the while loop.
"""
function vectors_spherical(npoints::Integer, eps=1e-4; seed=123)

	Random.seed!(seed)
	mean = zeros(3)
	C = ScalMat(3, 1.0) 
	d = MvNormal(mean, C)

	vx = zeros(npoints)
	vy = zeros(npoints)
	vz = zeros(npoints)
	
	for i = 1:npoints
		v = zeros(3)
	    while norm(v) < eps
	        v = rand(d, 1)
	    end
    
    	v = v / norm(v)
		vx[i] = v[1]
		vy[i] = v[2]
		vz[i] = v[3]
	end
	 DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end

# ╔═╡ 934e0612-9284-41f5-a58c-03e002094366
"""
Generate three standard normally distributed numbers and normalize the vector.

We have to be careful in the case that the vector has a norm close to zero, in which we must worry about floating point precision by dividing by a very small number. This is the reason for the while loop.
"""
function vectors_spherical2(npoints::Integer, eps=1e-4)

	vx = zeros(npoints)
	vy = zeros(npoints)
	vz = zeros(npoints)

	for i = 1:npoints
		v = zeros(3)
	    while norm(v) < eps
	        x = randn()  # random standard normal
	        y = randn()
	        z = randn()
	        v = [x, y, z]
	    end
    
    	v = v / norm(v)
		vx[i] = v[1]
		vy[i] = v[2]
		vz[i] = v[3]
	end
	 DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end

# ╔═╡ 1f6b02e3-ad55-44f3-a417-39285058aee5
gdf = vectors_spherical2(1000)

# ╔═╡ c80fa6be-29d5-4296-982a-f28a7a5bba90
@rtransform gdf begin
	:norm = norm([:vx, :vy, :vz])
end

# ╔═╡ cef79ccd-b9ac-4a63-9729-bfbc23906c6e
scatter(gdf.vx, gdf.vy, gdf.vz)

# ╔═╡ f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
"""
gtransport is a short name for
`generate photons in point p inside cylinder and propagate to cylinder surface'

"""
function gtransport(c::Cylinder, p::Vector{Float64}, nphotons::Integer=10)

    rdf = vectors_spherical2(nphotons) #
    gex = zeros(nphotons)
	gey = zeros(nphotons)
	gez = zeros(nphotons)
	gdx = zeros(nphotons)
	gdy = zeros(nphotons)
	gdz = zeros(nphotons)

    for i in 1:nphotons
    	d = [rdf[i,"vx"], rdf[i,"vy"], rdf[i,"vz"]]
        r = Ray(p,d)
        irr = ray_intersection_with_cylinder(r, c)
        gex[i] = irr.e[1]
		gey[i] = irr.e[2]
		gez[i] = irr.e[3]
		gdx[i] = irr.d[1]
		gdy[i] = irr.d[2]
		gdz[i] = irr.d[3]
	end
    gdf = DataFrame(id=1:nphotons, gex=gex, gey=gey, gez=gez, 
		                           gdx=gdx, gdy=gdy, gdz=gdz)
	leftjoin(rdf, gdf, on="ID"=>"id", matchmissing=:equal)
	#rdf,gdf
end

# ╔═╡ 47ff0204-6b05-4635-89fa-23dac3a68500
gtdf = gtransport(stdsk, [990.0,0.0,5.0], 20)

# ╔═╡ 23556eb1-2cad-47cf-b686-14d6f648ca9f
dgt = [gtdf[indx,"vx"], gtdf[indx,"vy"], gtdf[indx,"vz"]]

# ╔═╡ 35e8c307-1b9d-43cb-8de9-6a7776245b65
rt2 = Ray([990.0,0.0,5.0], dgt)

# ╔═╡ f3408829-adae-4c06-9fbd-ecb074795a1c
cylinder_intersection_roots(rt2, stdsk)

# ╔═╡ 557f9c90-6ca3-488d-84b5-f2fc285f3dc9
ray_intersection_with_cylinder(rt2, stdsk)

# ╔═╡ 09ec4f0d-bd08-4fee-83b4-40844a2b7236
begin
	scatter(gtdf.gex, gtdf.gey, gtdf.gez)
	plot!(gtdf.gex, gtdf.gey, gtdf.gez)
end

# ╔═╡ 1144066e-02ad-4d22-9d8d-e56abec9ce0f
draw_rays([10.0,10.0,100.0], gtdf)

# ╔═╡ f77c648e-275e-4103-905a-3927600d1aa5
size(gtdf)

# ╔═╡ f0fb5682-c081-44ab-84a5-7e418984f0df
md"""
## Drawing a cylinder
- A cylinder of radius r, and height h, having z-axis as symmetry axis, is a stack of circles of radius r:

$\begin{align}
x= r \times \cos(u) \\
y=r \times sin(u) \\
z=v 
\end{align}$

where:
$u \in [0,2\pi]$
$v \in [0,h]$

- This allows to plot the cylinder as any parameterized surface.
- For example:
	- Define u as a range beteen 0 and 2π
	- Define v as a range between 0 and 10.0 (the length of the cylinder)
"""

# ╔═╡ 5a74d977-a44e-4cba-b738-7bdd76056686
function plot_cylinder(c::Cylinder, ptheta=100, pz=100)
	r = c.r
	h = clength(c)
	m, n =ptheta, pz
	u = range(0, 2π , length=n)
	v = range(0, h, length=m)
	
	us = ones(m)*u'
	vs = v*ones(n)'
	#Surface parameterization
	X = r*cos.(us)
	Y = r*sin.(us)
	Z = vs
	GLMakie.surface(X, Y, Z)
end

# ╔═╡ 53a6f956-1957-4915-9900-3af1f44ed2fa
plot_cylinder(cyl)

# ╔═╡ 476bc4d8-7e22-4342-92ad-4896e47d437f
plot_cylinder(stdsk)

# ╔═╡ d9ef36e2-56f6-4655-98f2-37c2699363b5
"""
Given two vectors, jmgrid returns the julia equivalent to the meshgrid function.
If transp = true it returns the exact equivalent, which is the matrix transposed of
the usual matrix in julia. In Julia, contrary to Python, the first array index is a row and the second is a column (Julia has column-major arrays)
"""
function jmgrid(p0::Vector{<:Number}, p1::Vector{<:Number} ; transp=true)
	x = first.(Iterators.product(p0, p1))
	y = last.(Iterators.product(p0, p1))
	if transp
		collect(transpose(x)), collect(transpose(y))
	else
		x,y
	end
end

# ╔═╡ 32b53e87-034e-49db-842f-58f0c0002ed7
"""Normal to the cylinder barrel
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
then n = Grad(F)_P /Norm(Grad(F))_P
n = (2x, 2y, 0)/sqrt(4x^2 + 4y^2) = (x,y,0)/r (P)
"""
normal_to_barrel(c::Cylinder, P::Vector{Float64}) = [P[1], P[2], 0] ./ c.r



# ╔═╡ cb1d79b8-ada9-481e-8169-ab5976cb0759
"""
Give two points p0 and p1, thid function recturns;
mag: the magnitude of the vector defining the axis between p1 and p2
v: the unit vector in the direction of the axis defined by p0 and p1
n1, an n2: two unit vectors perpendicular to v and among themselves. 
"""
function unit_vectors_(p0::Vector{Float64}, p1::Vector{Float64})
	#vector in direction of axis
	v = p1 - p0

	#find magnitude of vector
	mag = norm(v,2)

	#unit vector in direction of axis
	v = v / mag

	# choose (1,0,0) as second axis unless is first axis
	not_v = [1.0, 0.0, 0.0]
	if v == not_v
		not_v = [0.0, 1.0, 0.0]
	end

	#make vector perpendicular to v and not v
	n1 = cross(v, not_v)

	#normalize n1
	n1 = n1 /norm(n1,2)

	#make unit vector perpendicular to v and n1
	n2 = cross(v, n1)

	mag, v,  n1, n2
end

# ╔═╡ 55a30f75-af85-4859-b1e9-b94b189e061f
function surfaces_(c::Cylinder, np=100)

	mag, v,  n1, n2 = unit_vectors_(c.p0, c.p1)
	
	
	#surface ranges over t from 0 to length of axis and 0 to 2*pi
	t = collect(range(0., mag, 2))
	theta = collect(range(0., 2π, np)) 
	rsample = collect(range(0., c.r, 2)) 

	#use meshgrid to make 2d arrays
	t, theta2 = jmgrid(t, theta)

	rsample,theta = jmgrid(rsample, theta)

	#generate coordinates for surface
	# "Tube"
	
	X,Y,Z = [c.p0[i] .+ v[i] * t .+ c.r * sin.(theta2) * n1[i] .+ c.r * cos.(theta2) *  n2[i] for i in 1:3]
	
	# "Bottom"
	
	X2, Y2, Z2 = [c.p0[i] .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	
	# "Top"
	X3, Y3, Z3 = [c.p0[i] .+ v[i]*mag .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	
	t, theta2, theta, rsample, (X,Y,Z), (X2, Y2, Z2), (X3, Y3, Z3)
end

# ╔═╡ 336ce641-c3e8-4fed-adc4-c8cdfec70ae7
cyl

# ╔═╡ 55ec1131-ef6c-4b99-9985-2b3a6d0f21e9
ts,theta2s,thetas, rs, P1, P2, P3 = surfaces_(cyl)

# ╔═╡ 252b84b2-aa7c-4d17-887b-241d0302e446
function draw_cylinderx(cl::Cylinder, np=100;  resolution = (1000,1000), 
                       col = :blue, clevel=0.5, linewidth = 2.5, wf=true)
	gcyl2, cylm = get_cylinder(cyl,np)
	ts,theta2s,thetas, rs, P1, P2, P3 = surfaces_(cyl)
	
	scene = Scene()
	cam3d!(scene)
	
	#with_theme(theme_light()) do
	#	fig = Figure(resolution = (1200,800))
	#	axs = Axis3(fig[1,1]; aspect=:data, perspectiveness=0.5) 
	#	if wf == false
    #cyld = GLMakie.wireframe!(scene, cylm; color = :grey90, linewidth = linewidth)
	
	sfs = GLMakie.surface!(scene, P1..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
	sfs2 = GLMakie.surface!(scene, P2..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
	sfs3 = GLMakie.surface!(scene, P3..., transparency = true, fxaa = true, overdraw=true, colormap=:greens)
	
	GLMakie.rotate!(sfs, Vec3f(1, 0, 0), 1.5)
	GLMakie.rotate!(sfs2, Vec3f(1, 0, 0), 1.5)
	GLMakie.rotate!(sfs3, Vec3f(1, 0, 0), 1.5)
	#GLMakie.rotate!(scene, Vec3f(1, 0, 0), 1.5)
	#GLMakie.rotate!(scene, Vec3f(0, 0, 1), 1.5)
	center!(scene)
    		
		#end
		#fig
	#end
	scene
end

# ╔═╡ 73d06337-6871-4e10-8c80-fac2346f8db4
draw_cylinderx(cyl, 200)

# ╔═╡ 2f10226d-093b-45eb-8894-e324a49a9663


# ╔═╡ 24a6ee8c-f761-49d2-96b8-320fbb35a4f3
function overd()
	scene = Scene()
	cam3d!(scene)
	sphere_plot = mesh!(scene, Sphere(Point3f(0), 0.5), color=:red)
	GLMakie.scale!(scene, 0.5, 0.5, 0.5)
	GLMakie.rotate!(scene, Vec3f(1, 0, 0), 0.5) # 0.5 rad around the y axis
	scene
end

# ╔═╡ cd62ccea-5adb-4a00-8f65-f72099bb36ee
overd()

# ╔═╡ 46ec8218-2ca1-45ba-8b55-cf79261c29fd
begin
GLMakie.surface(P1..., transparency = true, fxaa = true, overdraw=true, colormap=:greys)
GLMakie.surface(P2..., transparency = true, fxaa = true, overdraw=true, colormap=:greys)
end

# ╔═╡ 17e198a6-4f9a-40d6-95cb-270a4d3b14f4
GLMakie.surface(P2..., transparency = true, fxaa = true, colormap=:greys)

# ╔═╡ 023b42f9-c573-4055-9805-65adda40e0e0
function gridlaplacian(m, n)
    S = sparse(0.0I, n*m, n*m)
    linear = LinearIndices((1:m, 1:n))
    for i in 1:m
        for j in 1:n
            for (i2, j2) in ((i + 1, j), (i, j + 1))
                if i2 <= m && j2 <= n
                    S[linear[i, j], linear[i2, j2]] -= 1
                    S[linear[i2, j2], linear[i, j]] -= 1
                    S[linear[i, j], linear[i, j]] += 1
                    S[linear[i2, j2], linear[i2, j2]] += 1
                end
            end
        end
    end
    return S
end



# ╔═╡ 821b443e-5454-4aa9-a38b-c019c23c75d1
# d is used to denote the size of the data
d = 150

 # Sample centered Gaussian noise with the right correlation by the method
 # based on the Cholesky decomposition of the precision matrix


# ╔═╡ b124defe-cb13-45e3-b4fd-cf7c321023df
data = 0.1randn(d,d) + reshape(
        cholesky(gridlaplacian(d,d) + 0.003I) \ randn(d*d),
        d, d
)

# ╔═╡ a021735b-2e93-45bc-b635-00ebdc54f9df
begin
GLMakie.surface(data; shading=false, colormap = :deep)
#GLMakie.surface(data; shading=false, colormap = :deep)
end

# ╔═╡ Cell order:
# ╠═c8cec5b9-abb9-443d-9e41-bc8b7070792d
# ╠═a298b03b-c7c4-4114-9984-d37e7df8d4c8
# ╠═5158389a-fca1-11ed-05a1-c7345c9bd354
# ╠═e36751ff-a521-4529-95c8-6bbfc3314e66
# ╠═df3dde7a-9437-409f-bd51-db975f74b78c
# ╠═fc869138-1cf7-463b-aca5-729741f8b1ff
# ╠═f230deec-f4c9-429a-89f5-92d8c08e8472
# ╠═7371f833-a7a8-482c-925b-730d8378c22e
# ╠═59ab61e2-00ff-4a06-b732-04f7d65f7325
# ╠═2c803aef-9a52-4108-a252-e63a95929c67
# ╠═c59e3a3d-1821-49a1-8d93-29250146330e
# ╠═b853da86-3242-4a9f-8163-1e5bd0d6a907
# ╠═c7583c88-4cfc-49e1-b61e-c299801be450
# ╠═75472827-851c-451f-9510-a7e33fee1a8e
# ╠═3c398d08-ef94-4202-9b87-62aa4fba07ed
# ╠═53a6f956-1957-4915-9900-3af1f44ed2fa
# ╠═ef077643-f5b1-4244-af7e-8111ede962a0
# ╠═81d62300-1a88-4c9c-b050-8acd20857867
# ╠═6795ad55-06b8-4a77-9a3a-cf3c15a3741e
# ╠═b6082524-84ae-4d98-b593-7250f554c0c3
# ╠═c48c5f53-2632-4cfc-b07b-7fde9f1d8eb6
# ╠═de0fb8fd-c757-4769-a2da-d874add0a553
# ╠═50f7c64a-4fa6-4705-b2e9-896ef545029f
# ╠═d06eaf38-f8bb-4732-9b8f-92f2a1053176
# ╠═2940f9d9-f2c3-4772-8208-2375370b7e61
# ╠═a7378f23-9b56-4a5f-8c4b-b51545e69633
# ╠═c95dcf93-40ff-4db0-a323-29991f5135dd
# ╠═c7200fa5-31eb-458c-80ca-70f968edd9b0
# ╠═9e193e81-0147-46d2-943f-1526dfb26e50
# ╠═a89c0c98-ee2c-4103-8cd3-9be03a2b75ae
# ╠═da50ac8a-8381-490b-a000-b6ba8b392888
# ╠═aa445ed2-dc2d-4802-ab12-18505487ceef
# ╠═24cd55d7-6752-4e48-a3fa-7951cbff1b36
# ╠═6c2b1465-02ae-4768-a4f6-e8f461339b97
# ╠═c3abda05-42ec-4759-bea4-cdf7ee92ed70
# ╠═26d72ea7-31e5-4c55-af55-7666f7c47f1c
# ╠═b9c450a6-446e-4a5c-93ee-28a45cf96356
# ╠═5de835fa-284f-4dd0-ba1d-033d4c50da5a
# ╠═71203991-d5c4-439f-8dd1-fe1b4cf3b875
# ╠═366cf831-3558-44f2-bf09-6d9de74bf04d
# ╠═95a82d89-7b7f-4ce9-8f21-b527543c57f7
# ╠═59138131-46c1-4ac8-be0b-42f5ca2e30b8
# ╠═02a7fd08-440c-4e0c-93c8-d42792642d11
# ╠═40cff386-97de-4e6c-8d22-216302831893
# ╠═bca3de2e-74a6-423e-a660-3845c3e9657b
# ╠═7f4ab967-b5c2-46d9-90c4-ebbbdab0dbfe
# ╠═3b4ea019-bba9-4ac1-b61a-fa251c9524d9
# ╠═c2d3e54d-3ad5-46b6-9d28-be8acf1086b2
# ╠═96af216d-c0e3-4825-a11f-3aaa04c42d0b
# ╠═21aeb736-8655-4747-bfa7-1ad9b5395fc3
# ╠═5dcb797a-d7ef-463b-a629-336a8752c4ed
# ╠═92ef9d4b-ebf7-4247-a13e-aebe60cd79ad
# ╠═80b17191-b750-4f31-9577-1f5214440429
# ╠═c97578f2-1d39-497a-9e79-31e041e77b7b
# ╠═4d17bf89-d01b-4f1b-9dc7-9dcd00ee1d8a
# ╠═d137c26c-bf75-4899-ac0a-eb8e70c792bb
# ╠═a143fe7e-e647-422e-8e0b-65099ff14d88
# ╠═0dc2ca61-8ca9-43f2-ba71-87dcc2be08ec
# ╠═312154b9-aa06-4413-bcd7-7688fe168660
# ╠═68e9ab79-21a4-48c1-82c3-45d8862e70cf
# ╠═1f6b02e3-ad55-44f3-a417-39285058aee5
# ╠═ff45920d-225a-48df-a9e3-f567b81bd579
# ╠═b4110420-ef1c-4e1b-868a-68bd9d943fea
# ╠═c80fa6be-29d5-4296-982a-f28a7a5bba90
# ╠═5e866260-a798-4c5f-b407-38d091d0ed07
# ╠═cef79ccd-b9ac-4a63-9729-bfbc23906c6e
# ╠═be8dd99e-60f5-4d75-959c-cffbd5bfa4b8
# ╠═c7d4775a-f585-4852-a172-337655cb21ed
# ╠═476bc4d8-7e22-4342-92ad-4896e47d437f
# ╠═3b43b3d6-2621-4940-8f6e-cf5a5c6a2093
# ╠═47ff0204-6b05-4635-89fa-23dac3a68500
# ╠═f675ebdf-e3f6-4f6b-bdae-e59db854f63c
# ╠═23556eb1-2cad-47cf-b686-14d6f648ca9f
# ╠═35e8c307-1b9d-43cb-8de9-6a7776245b65
# ╠═f3408829-adae-4c06-9fbd-ecb074795a1c
# ╠═557f9c90-6ca3-488d-84b5-f2fc285f3dc9
# ╠═09ec4f0d-bd08-4fee-83b4-40844a2b7236
# ╠═1144066e-02ad-4d22-9d8d-e56abec9ce0f
# ╠═f77c648e-275e-4103-905a-3927600d1aa5
# ╠═3dce337c-3f5b-4733-870d-15cf8197525b
# ╠═ce9da661-5f2f-4e5d-825e-5d029c274066
# ╠═efe04df4-cb56-4916-96f7-296dc37492ac
# ╠═053952a1-0558-447c-a3f5-ae30c2d2f7be
# ╠═2facda59-dc5a-4951-8d6f-1d121690aad2
# ╠═991a026a-bfe2-4bc8-a200-31dd25a881d1
# ╠═539faa01-299b-4c5b-8113-fc21501dd84d
# ╠═63c8a7db-740a-4abf-981f-3c02be41493a
# ╠═911252fb-0b59-4d23-a90e-4edd293ea91f
# ╠═269e4b09-3505-40fe-8815-2bc847d02a99
# ╠═f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
# ╠═c26c8511-afae-4447-9458-d1d21f565911
# ╠═deefebf9-07a3-4cce-b2d5-4b42644b5db7
# ╠═e9d59365-ea79-4e97-b908-cbd528dc0804
# ╠═da14195e-414e-40ab-aca1-83d2d68b9172
# ╠═f41fbfdf-4cee-4b82-b087-e224c70b7f6b
# ╠═7f8c81be-0f2d-4214-b27d-a7470534d867
# ╠═372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
# ╠═0abeba90-ba75-4231-9dd8-d35400964c60
# ╠═464a3cf9-c030-4270-bae7-bf4339ed3fc8
# ╠═58491933-49ec-4b70-ab1f-2301657ed9e8
# ╠═4d53eaca-0cad-4f34-8205-216e9c0a16d9
# ╠═989cec01-548a-40e8-9a82-149dbfa7d367
# ╠═c6f63bc1-6d4d-49b3-988f-ccaaaae562e1
# ╠═c40ad59c-7687-4bba-b923-fa498295821c
# ╠═8547b5f7-6b8a-4a9b-b0bf-5b1d8b08c072
# ╠═14583522-a4d2-4f37-9c8b-5e2f54ace834
# ╠═934e0612-9284-41f5-a58c-03e002094366
# ╠═f0fb5682-c081-44ab-84a5-7e418984f0df
# ╠═5a74d977-a44e-4cba-b738-7bdd76056686
# ╠═d9ef36e2-56f6-4655-98f2-37c2699363b5
# ╠═32b53e87-034e-49db-842f-58f0c0002ed7
# ╠═cb1d79b8-ada9-481e-8169-ab5976cb0759
# ╠═55a30f75-af85-4859-b1e9-b94b189e061f
# ╠═336ce641-c3e8-4fed-adc4-c8cdfec70ae7
# ╠═55ec1131-ef6c-4b99-9985-2b3a6d0f21e9
# ╠═252b84b2-aa7c-4d17-887b-241d0302e446
# ╠═73d06337-6871-4e10-8c80-fac2346f8db4
# ╠═2f10226d-093b-45eb-8894-e324a49a9663
# ╠═24a6ee8c-f761-49d2-96b8-320fbb35a4f3
# ╠═cd62ccea-5adb-4a00-8f65-f72099bb36ee
# ╠═46ec8218-2ca1-45ba-8b55-cf79261c29fd
# ╠═17e198a6-4f9a-40d6-95cb-270a4d3b14f4
# ╠═7886beee-8947-4cb2-bcd9-89728829b7cd
# ╠═023b42f9-c573-4055-9805-65adda40e0e0
# ╠═821b443e-5454-4aa9-a38b-c019c23c75d1
# ╠═b124defe-cb13-45e3-b4fd-cf7c321023df
# ╠═a021735b-2e93-45bc-b635-00ebdc54f9df
