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

# ╔═╡ c8cec5b9-abb9-443d-9e41-bc8b7070792d


# ╔═╡ 110b6ec6-eb56-4281-940e-320ad49fa927
begin
	mean = zeros(3)
	C = ScalMat(3, 1.0) 
	d = MvNormal(mean, C)
	x = rand(d, 1)
	norm(x)
end

# ╔═╡ e36751ff-a521-4529-95c8-6bbfc3314e66
plotlyjs()

# ╔═╡ df3dde7a-9437-409f-bd51-db975f74b78c
PlutoUI.TableOfContents(title="JRT draft", indent=true)


# ╔═╡ a8a3afff-376e-4bb8-afdd-cb9d7803958b
md"""
## The meshgrid concept

- Consider the figure below. Suppose that we want to define a two-dimensional area formed by the coordinate points x = [0, 1, 2], y = [0, 1, 2]. The coordinates of the resulting 9 points are represented in the figure below. 
"""

# ╔═╡ 00b4f394-29af-4927-a47d-6956469cb2c8
mgrid = load("tutorial_numpy_function_meshgrid_02.jpeg")

# ╔═╡ b8954a27-a492-4084-9c85-6123e36a671e
md"""
- In julia is very simple to obtaind the matrix of tuples representing those point.
- First we define the function f(x,y) = (x,y)
- Then we apply it to the vector represeting x (a column vector of dimensions 3x1) times the adjoint of y (a raw vector of dimension 1x3) to obtain a 3x3 matrix
"""

# ╔═╡ 24b691e1-797a-44c0-a6c4-688fb05ce8e4
begin 
	xs = [0, 1, 2]
	ys = [0, 1, 2]	
	f(x, y) = (x, y)
	MC = f.(xs,ys')
end

# ╔═╡ 55a2c145-85c0-4fd0-90d6-90f55db07379
md"""
- Notice that in Julia (like in Fortran but unlike C and phython) c the first array index is a row and the second is a column (Julia has column-major arrays). Thus in the matrix of pairs above the first row describes the first column of the figure. 
"""

# ╔═╡ 8b8b88ed-e44a-450d-8999-395aa8c1bc0e
md"""
- We now want to obtain two arrays, each of size (3, 3) that include the x and y coordinates (separately) of the 9 points in question. This is also very easy.
"""

# ╔═╡ 61e07f7f-c1c7-4338-806c-17680d2a3b1d
ones(3) * xs'


# ╔═╡ dde96583-30e7-4e80-a927-9bbada1ce2d7
ys * ones(3)'

# ╔═╡ e8d11a5f-b2e7-49cf-801a-d5ecf8f14d0b
md"""
- Notice that the matrices read now by row.

- This is exactly what the numpy.meshgrid function does: it accepts as input the coordinates that define the hyperplane segment (it can be two dimensions or any other number of them) and returns arrays with the coordinates of those points

- Given two vectors, jmgrid returns the julia equivalent to the meshgrid function.

- If transp = true it returns the exact equivalent, which is the matrix transposed of the usual matrix in julia. 
"""

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

# ╔═╡ 653a5e64-cb0f-465e-83d1-5ce710ba17f3
u = collect(range(0, 2π , length=10))

# ╔═╡ 3c6c4864-7631-4b86-8496-785775c8f33d
v =collect(range(0, 10, length=10))


# ╔═╡ 728b4091-a2ab-462a-8880-d9c19604623c
md"""
- now obtain the matrix of first coordinates
"""

# ╔═╡ 668b8151-f69f-4bef-a5a9-c0835c18b2d6
	
	us = ones(10)*u'


# ╔═╡ dd08b21d-0d39-4cb7-b45e-f5cc3aa9572b
md"""
- And the matrix of second coordinates. 
"""

# ╔═╡ 2c8c91e7-003b-41c9-be46-7889e5b5692a
	vs = v*ones(10)'


# ╔═╡ 25a7b7f0-49d3-4f90-8245-fbfe656357b5
md"""
- Or directly using **jmgrid**
"""

# ╔═╡ c4681021-dfe4-4164-b272-87a7c8c4221a
md"""
- We can now do the surface parameterization
"""

# ╔═╡ 4d8a8ec6-7ac9-49f5-a9f5-429ba9ecabcd
begin
	r=5.0
	X = r*cos.(us)
	Y = r*sin.(us)
	Z = vs
end

# ╔═╡ 7a317c2f-aede-4802-897f-6ae907f01c91
X

# ╔═╡ f347eb0e-b197-437a-801f-bc3385ca7f7a
md"""
- Then one can plot the surfaces
"""
	
	

# ╔═╡ 194396b3-cdd0-4ddd-81e5-e9759184697f
Plots.surface(X, Y, Z, size=(600,600), cbar=:none, legend=false)

# ╔═╡ 75472827-851c-451f-9510-a7e33fee1a8e
md"""
## Defining and plotting a cylinder
"""

# ╔═╡ 64732133-48ad-4bc2-ac91-d5cf27f840e8
md"""
# Geometry Basics and Makie
"""

# ╔═╡ 0d23ed6d-8380-4de3-a81b-d28000f192ba
md"""
- Define a cylinder using GeometryBasics (it takes P0, P1 and radius)
"""

# ╔═╡ 1f3a6e11-8368-4181-838c-23dc2adc11db
cl = GeometryBasics.Cylinder(Point{3, Float64}(0,0,0), Point{3, Float64}(0,0,10), 5.0)

# ╔═╡ 7c37c89a-372c-47fc-9192-c55e2626f2bb
md"""
- Define a Sphere (takes a point for the center and a radius)
"""

# ╔═╡ 753bf5aa-daf6-4718-b059-0321347975a0
sphere =Sphere(Point3f(0), 1)

# ╔═╡ d88d0278-49cc-4044-8332-e35c0253fb9d
md"""
- Define the basic mesh and a finer mesh using teselation
"""

# ╔═╡ 26a44ca3-7622-4ca8-a621-49441df5ce70
msh1 = GeometryBasics.mesh(cl)

# ╔═╡ f94e907f-7ad2-4bba-8074-672abe609307
GeometryBasics.length(GeometryBasics.coordinates(msh1))

# ╔═╡ 17030aa1-c919-413d-8977-4493252a0182
msh2 =GeometryBasics.mesh(GeometryBasics.Tesselation(cl, 200))

# ╔═╡ 4693ff80-108e-45e6-8632-54473759960b
GeometryBasics.length(GeometryBasics.coordinates(msh2))

# ╔═╡ fd72c797-4dbc-4a55-ae4f-f5e7157f7cb9
md"""
- Draw using Makie
"""

# ╔═╡ 0929957e-f3b9-4cfa-bbd6-d93a2d1779fe
with_theme(theme_light()) do
    fig = Figure(resolution = (900,900))
    ax1 = Axis3(fig[1,1]; aspect = :data, perspectiveness = 0.5, azimuth=0.5)
	cmesh = GeometryBasics.mesh(GeometryBasics.Tesselation(cl, 64))
    mesh!(ax1, cmesh; color = (:dodgerblue, 0.75))
    fig
end

# ╔═╡ de0fb8fd-c757-4769-a2da-d874add0a553
md"""
- Check function **in_endcaps(c::Cylinder, p::Vector{Float64})**
"""

# ╔═╡ c95dcf93-40ff-4db0-a323-29991f5135dd
md"""
Check funding roots
"""

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

# ╔═╡ fb3ac3dd-4737-4543-912d-b34fa881850c
cyl

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


# ╔═╡ 40cff386-97de-4e6c-8d22-216302831893
ry = Ray([1.0, 0.0, 0.0], [0.5, 0.5, 1.0])

# ╔═╡ 2facda59-dc5a-4951-8d6f-1d121690aad2
ray(r::Ray, t::Float64) =r.e + t * r.d


# ╔═╡ 991a026a-bfe2-4bc8-a200-31dd25a881d1
md"""
# Functions
"""

# ╔═╡ c26c8511-afae-4447-9458-d1d21f565911
"""
Equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
cylinder_equation(c::Cylinder, P::Vector{Float64}) = P[1]^2 + P[2]^2 - c.r^2

# ╔═╡ 32b53e87-034e-49db-842f-58f0c0002ed7
"""Normal to the cylinder barrel
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
then n = Grad(F)_P /Norm(Grad(F))_P
n = (2x, 2y, 0)/sqrt(4x^2 + 4y^2) = (x,y,0)/r (P)
"""
normal_to_barrel(c::Cylinder, P::Vector{Float64}) = [P[1], P[2], 0] ./ c.r



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
function in_endcaps(c::Cylinder, p::Vector{Float64})
   p[3] ≈ c.zmin || p[3] ≈ c.zmax
end

# ╔═╡ d06eaf38-f8bb-4732-9b8f-92f2a1053176
in_endcaps(cyl, [0.0, 0.0, 5.0]) # not in end-caps

# ╔═╡ 2940f9d9-f2c3-4772-8208-2375370b7e61
in_endcaps(cyl, [0.0, 0.0, 0.0] )  # in bot

# ╔═╡ a7378f23-9b56-4a5f-8c4b-b51545e69633
in_endcaps(cyl, [0.0, 0.0, 10.0] )  # in top

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

# ╔═╡ c2d3e54d-3ad5-46b6-9d28-be8acf1086b2
find_zeros(ff, -10.0, 10.0)

# ╔═╡ 3b4ea019-bba9-4ac1-b61a-fa251c9524d9
Plots.plot(-100:0.01:100, ff)

# ╔═╡ c6f63bc1-6d4d-49b3-988f-ccaaaae562e1
"""
Computes intersection roots between a ray and a cylinder
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
function cylinder_intersection_roots(r::Ray, cl::Cylinder, eps::Float64=1e-9)
	
	f = cylinder_intersection_poly(r, cl)
	
	#println("a = ", -2*clength(cyl) * cyl.r, " b = ", 2*clength(cyl) * cyl.r )
	
	rts = find_zeros(f, -2*clength(cyl) * cyl.r, 2*clength(cyl) * cyl.r)
	if length(rts) == 0
		return 0
	else
		prts = [x for x in rts if x >eps]
		return minimum(prts)
	end
end

# ╔═╡ 96af216d-c0e3-4825-a11f-3aaa04c42d0b
cylinder_intersection_roots(ry, cyl)

# ╔═╡ 72625d7a-5146-4eaf-a2c8-d9e451f38eaf
   """
   Intersection between a ray and the end-cups of a cylinder
   """

# ╔═╡ c40ad59c-7687-4bba-b923-fa498295821c
function ray_intersection_with_cylinder_end_caps(r::Ray, cl::Cylinder, t::Float64)
 
    p = ray(r, t)
    if p[3] > cl.zmax
        t = (cl.zmax - r.e[3])/r.d[3]
    else
        t = (cl.zmin - r.e[3])/r.d[3]
	end

    t, ray(r,t)
end

# ╔═╡ 19ad1152-6ad4-4cff-8769-798d1be02364
ray_intersection_with_cylinder_end_caps(Ray([1.0, 0.0, 0.0], [0.0, 1.0, 1.0]), cyl, 10.0)

# ╔═╡ 8547b5f7-6b8a-4a9b-b0bf-5b1d8b08c072
 """
 Intersection between a ray and a cylinder
 """
function ray_intersection_with_cylinder(r::Ray, cl::Cylinder)
   
    t = cylinder_intersection_roots(r, cl)

	if t == 0 && r.d[3] >= 0  # cuts upper endcap
		return zmax, ray(r, zmax)
	elseif t == 0 && r.d[3] < 0  # cuts lower endcap
		return -zmax, ray(r, -zmax)
	end
    
	P = ray(r, t)
    z = P[3]
    if z < cl.zmin || z > cl.zmax
        t, P = ray_intersection_with_cylinder_end_caps(r, cl, t)
	end
    t, P
end

# ╔═╡ 21aeb736-8655-4747-bfa7-1ad9b5395fc3
ray_intersection_with_cylinder(ry, cyl)

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

# ╔═╡ f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
"""
gtransport is a short name for
`generate photons in point p inside cylinder and propagate to cylinder surface'

"""
function gtransport(c::Cylinder, p::Vector{Float64}, nphotons::Integer=10)

    rdf = vectors_spherical(nphotons) #
    gx = zeros(nphotons)
	gy = zeros(nphotons)
	gz = zeros(nphotons)

    for i in 1:n
    	d = [rdf[i,"vx"], rdf[i,"vy"], rdf[i,"vz"]]
        r = Ray(p,d)
        t, P = ray_intersection_with_cylinder(r, c)
        gx[i] = Px
		gy[i] = Py
		gz[i] = Pz
	end
    gdf = DataFrame(id=1:n, gx=gx, gy=gy, gz=gz)
	leftjoin(rdf, gdf, on="ID"=>"id", matchmissing=:equal)
end

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
	Plots.surface(X, Y, Z, size=(600,600), cbar=:none,  legend=false)
end

# ╔═╡ 53a6f956-1957-4915-9900-3af1f44ed2fa
plot_cylinder(cyl)

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

# ╔═╡ 50e60fc5-22d3-4a1a-b720-ea05b53971f3
jmgrid(xs,ys)

# ╔═╡ a8fba5d3-1d45-41b5-9a9c-0f27cbcf96d1
jmgrid(xs,ys, transp =false)

# ╔═╡ e009d524-96c1-4d80-b49d-244ad134a148
us2, vs2 = jmgrid(u,v)

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
	if v .== not_v
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
function surfaces_(r::Float64, mag::Float64, 
	               v::Vector{Float64}, n1::Vector{Float64}, n2::Vector{Float64},                     np=100)

	
	#surface ranges over t from 0 to length of axis and 0 to 2*pi
	t = LinRange(0., mag, 2) 
	theta = LinRange(0., 2π, np) 
	rsample = LinRange(0., c.r, 2) 

	#use meshgrid to make 2d arrays
	t, theta2 = jmgrid(t, theta)

	rsample,theta = jmgrid(rsample, theta)

	#generate coordinates for surface
	# "Tube"
	X, Y, Z = [p0[i] + v[i] * t + r * sin.(theta2) * n1[i] + r * cos.(theta2) *  n2[i] for i in 1:3]
	
	# "Bottom"
	X2, Y2, Z2 = [p0[i] + rsample[i] * sin.(theta) * n1[i] + rsample[i] * ncos.(theta) * n2[i] for i in 1:3]
	
	# "Top"
	X3, Y3, Z3 = [p0[i] + v[i]*mag + rsample[i] * sin.(theta) * n1[i] + rsample[i] * cos.(theta) * n2[i] for i in 1:3]
	(X,Y,Z), (X2, Y2, Z2), (X3, Y3, Z3)
end

# ╔═╡ Cell order:
# ╠═c8cec5b9-abb9-443d-9e41-bc8b7070792d
# ╠═a298b03b-c7c4-4114-9984-d37e7df8d4c8
# ╠═5158389a-fca1-11ed-05a1-c7345c9bd354
# ╠═110b6ec6-eb56-4281-940e-320ad49fa927
# ╠═e36751ff-a521-4529-95c8-6bbfc3314e66
# ╠═df3dde7a-9437-409f-bd51-db975f74b78c
# ╠═a8a3afff-376e-4bb8-afdd-cb9d7803958b
# ╠═00b4f394-29af-4927-a47d-6956469cb2c8
# ╠═b8954a27-a492-4084-9c85-6123e36a671e
# ╠═24b691e1-797a-44c0-a6c4-688fb05ce8e4
# ╠═55a2c145-85c0-4fd0-90d6-90f55db07379
# ╠═8b8b88ed-e44a-450d-8999-395aa8c1bc0e
# ╠═61e07f7f-c1c7-4338-806c-17680d2a3b1d
# ╠═dde96583-30e7-4e80-a927-9bbada1ce2d7
# ╠═e8d11a5f-b2e7-49cf-801a-d5ecf8f14d0b
# ╠═50e60fc5-22d3-4a1a-b720-ea05b53971f3
# ╠═a8fba5d3-1d45-41b5-9a9c-0f27cbcf96d1
# ╠═c7583c88-4cfc-49e1-b61e-c299801be450
# ╠═653a5e64-cb0f-465e-83d1-5ce710ba17f3
# ╠═3c6c4864-7631-4b86-8496-785775c8f33d
# ╠═728b4091-a2ab-462a-8880-d9c19604623c
# ╠═668b8151-f69f-4bef-a5a9-c0835c18b2d6
# ╠═dd08b21d-0d39-4cb7-b45e-f5cc3aa9572b
# ╠═2c8c91e7-003b-41c9-be46-7889e5b5692a
# ╠═25a7b7f0-49d3-4f90-8245-fbfe656357b5
# ╠═e009d524-96c1-4d80-b49d-244ad134a148
# ╠═c4681021-dfe4-4164-b272-87a7c8c4221a
# ╠═4d8a8ec6-7ac9-49f5-a9f5-429ba9ecabcd
# ╠═7a317c2f-aede-4802-897f-6ae907f01c91
# ╠═f347eb0e-b197-437a-801f-bc3385ca7f7a
# ╠═194396b3-cdd0-4ddd-81e5-e9759184697f
# ╠═75472827-851c-451f-9510-a7e33fee1a8e
# ╠═3c398d08-ef94-4202-9b87-62aa4fba07ed
# ╠═53a6f956-1957-4915-9900-3af1f44ed2fa
# ╠═64732133-48ad-4bc2-ac91-d5cf27f840e8
# ╠═0d23ed6d-8380-4de3-a81b-d28000f192ba
# ╠═1f3a6e11-8368-4181-838c-23dc2adc11db
# ╠═7c37c89a-372c-47fc-9192-c55e2626f2bb
# ╠═753bf5aa-daf6-4718-b059-0321347975a0
# ╠═d88d0278-49cc-4044-8332-e35c0253fb9d
# ╠═26a44ca3-7622-4ca8-a621-49441df5ce70
# ╠═f94e907f-7ad2-4bba-8074-672abe609307
# ╠═17030aa1-c919-413d-8977-4493252a0182
# ╠═4693ff80-108e-45e6-8632-54473759960b
# ╠═fd72c797-4dbc-4a55-ae4f-f5e7157f7cb9
# ╠═0929957e-f3b9-4cfa-bbd6-d93a2d1779fe
# ╠═de0fb8fd-c757-4769-a2da-d874add0a553
# ╠═50f7c64a-4fa6-4705-b2e9-896ef545029f
# ╠═d06eaf38-f8bb-4732-9b8f-92f2a1053176
# ╠═2940f9d9-f2c3-4772-8208-2375370b7e61
# ╠═a7378f23-9b56-4a5f-8c4b-b51545e69633
# ╠═c95dcf93-40ff-4db0-a323-29991f5135dd
# ╠═40cff386-97de-4e6c-8d22-216302831893
# ╠═fb3ac3dd-4737-4543-912d-b34fa881850c
# ╠═7f4ab967-b5c2-46d9-90c4-ebbbdab0dbfe
# ╠═c2d3e54d-3ad5-46b6-9d28-be8acf1086b2
# ╠═96af216d-c0e3-4825-a11f-3aaa04c42d0b
# ╠═3b4ea019-bba9-4ac1-b61a-fa251c9524d9
# ╠═19ad1152-6ad4-4cff-8769-798d1be02364
# ╠═21aeb736-8655-4747-bfa7-1ad9b5395fc3
# ╠═ce9da661-5f2f-4e5d-825e-5d029c274066
# ╠═efe04df4-cb56-4916-96f7-296dc37492ac
# ╠═053952a1-0558-447c-a3f5-ae30c2d2f7be
# ╠═2facda59-dc5a-4951-8d6f-1d121690aad2
# ╠═991a026a-bfe2-4bc8-a200-31dd25a881d1
# ╠═c26c8511-afae-4447-9458-d1d21f565911
# ╠═32b53e87-034e-49db-842f-58f0c0002ed7
# ╠═deefebf9-07a3-4cce-b2d5-4b42644b5db7
# ╠═e9d59365-ea79-4e97-b908-cbd528dc0804
# ╠═da14195e-414e-40ab-aca1-83d2d68b9172
# ╠═f41fbfdf-4cee-4b82-b087-e224c70b7f6b
# ╠═7f8c81be-0f2d-4214-b27d-a7470534d867
# ╠═372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
# ╠═0abeba90-ba75-4231-9dd8-d35400964c60
# ╠═4d53eaca-0cad-4f34-8205-216e9c0a16d9
# ╠═989cec01-548a-40e8-9a82-149dbfa7d367
# ╠═c6f63bc1-6d4d-49b3-988f-ccaaaae562e1
# ╠═72625d7a-5146-4eaf-a2c8-d9e451f38eaf
# ╠═c40ad59c-7687-4bba-b923-fa498295821c
# ╠═8547b5f7-6b8a-4a9b-b0bf-5b1d8b08c072
# ╠═14583522-a4d2-4f37-9c8b-5e2f54ace834
# ╠═934e0612-9284-41f5-a58c-03e002094366
# ╠═f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
# ╠═5a74d977-a44e-4cba-b738-7bdd76056686
# ╠═d9ef36e2-56f6-4655-98f2-37c2699363b5
# ╠═cb1d79b8-ada9-481e-8169-ab5976cb0759
# ╠═55a30f75-af85-4859-b1e9-b94b189e061f
