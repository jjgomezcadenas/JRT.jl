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
	#using Plots
	#using Printf
	using Markdown
	using InteractiveUtils
	using Statistics
	using Random
	using Distributions
	using PDMats
	using Chain
	using LinearAlgebra
	using Images
	using GLMakie
	#using Makie
	import GeometryBasics
	using Colors
	using Roots
	using Unitful
	using UnitfulEquivalences
end

# ╔═╡ e24e5710-04b9-49c3-b3f6-2b136a1744eb
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M,
    Bq, mBq, μBq 

# ╔═╡ e36751ff-a521-4529-95c8-6bbfc3314e66
begin
	A_BI214_316Ti   =    1.0    * mBq/kg
	A_TL208_316Ti   =    0.4   * mBq/kg
	A_BI214_CU_LIM  =   12      * μBq/kg
	A_TL208_CU_LIM  =    1.4    * μBq/kg
end

# ╔═╡ df3dde7a-9437-409f-bd51-db975f74b78c
PlutoUI.TableOfContents(title="JRT draft", indent=true)


# ╔═╡ fc869138-1cf7-463b-aca5-729741f8b1ff
GLMakie.activate!()

# ╔═╡ c7583c88-4cfc-49e1-b61e-c299801be450
md"""
## Drawing a cylinder
"""

# ╔═╡ c48c5f53-2632-4cfc-b07b-7fde9f1d8eb6
md"""
## Testing functions
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

# ╔═╡ 8a45b62d-cbfb-4f8c-9c61-fad0e3ca53da
md"""
- Generate photons around one poit and transport photons to cylinder surfaces
"""

# ╔═╡ 128d5ed6-3a09-427e-b50d-159d0288d1dd
np=2000

# ╔═╡ 5cca07ea-8740-4899-bca7-629bcdb8d46b


# ╔═╡ 775a2389-f62d-4341-b0e6-976d150207b1
md"""
### Generate points inside cylinder
"""

# ╔═╡ 2a795370-8b35-4db5-8f5c-c03b3efaa1de
md"""
## Attenuation calculation (steel) 
"""

# ╔═╡ ea0d7d4d-5dcf-4d7f-b8bc-dd4856957df5
md"""
1. Define the cylinder describing the disk
2. throw random gammas in cylinder
3. Propagate gammas. Consider only gammas that hit the upper end-cup
4. Compute the weight of the gammas (in terms of attenuation length)
"""

# ╔═╡ b2a2c523-7461-4637-87c5-7ecd321dc7a5
md"""
1. Cylinder: r = 500 mm, h=50 mm
"""

# ╔═╡ e500854c-1714-44c6-8c76-107950d65060
zsteel = 50.0

# ╔═╡ 92e67cdb-f4c7-4ada-9939-332a1fa3546d
md"""
2. Generate random gammas in cylinder
"""

# ╔═╡ f1255b5c-d28a-4e48-9457-104a36d4596d
ngam = 2000

# ╔═╡ 7a70c468-156b-47c4-a2af-2a6c30034191
md"""
3. Propagate gammas. Consider only gammas that hit the upper end-cup
"""

# ╔═╡ f329a7b4-75db-46f7-b137-5ed2ebaea571
md"""
## Attenuation calculation (adding copper shield
"""

# ╔═╡ 9d72493b-1fbb-46f3-a670-081ea8484af9
md"""
1. Define the cylinder describing the copper shield (zmin(copper) = zmax(steel))
2. Take the gammas arriving to the shield (previous step) and propagate. Consider only gammas that hit the upper end-cup
3. Compute the weight of the gammas (in terms of attenuation length)
"""

# ╔═╡ 354cb5bc-c5a9-4f6e-94f0-5b12390b8922
zcu = 120.0

# ╔═╡ f211e72a-bc1d-4798-a2b8-6cb39663f4a5
md"""
## Putting it all together
"""

# ╔═╡ 326a7877-e625-468a-8a7a-e9ff36e2f6fb
md"""
## Activity of reinforcement plate
"""

# ╔═╡ fb431946-28d1-4042-9b18-aa4439e48d84
md"""
## Making holes in plates
"""

# ╔═╡ 9148483a-c144-4a10-8606-d73ab695abaa
1200/76

# ╔═╡ ce9da661-5f2f-4e5d-825e-5d029c274066
md"""
# Types
"""

# ╔═╡ 24b1b8ef-10ee-4d96-8121-827234233873
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

# ╔═╡ 6458b9a6-e0b1-40cf-87bd-3e3f01990525
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

# ╔═╡ 4559bbcc-4e78-4e29-b0b2-0063be1335f9
begin
	pti316  = PhysicalMaterial(7.87 * g/cm^3, 0.039 * cm^2/g)
	mti316 = RadioactiveMaterial(pti316, A_BI214_316Ti, A_TL208_316Ti)
	pcu  = PhysicalMaterial(8.96 * g/cm^3, 0.039 * cm^2/g)
	cu12 = RadioactiveMaterial(pcu, A_BI214_CU_LIM, A_TL208_CU_LIM)

end


# ╔═╡ 478e01bc-d764-4af7-a6f4-661cdba4820d
function transmittance_at_qbb(rmt::RadioactiveMaterial, d::Float64) 
	l = uconvert(mm, rmt.m.λ)/mm
	exp(-d/l)
end

# ╔═╡ 3a182a84-25f7-4bb2-bb54-6e81c53142e7
begin
	trcu1l = transmittance_at_qbb(cu12, 70.0)
	trcu2l = transmittance_at_qbb(cu12, 120.0)
	trfe1l = transmittance_at_qbb(mti316, 50.0)
end

# ╔═╡ 2b73f038-b92e-473c-8201-e2c12913cebf
md"""
- **The approximate calculation falls within a factor 2 of the exact calculation**
- Attenuation length of steel = $(mti316.m.λ)
- Attenuation length of cu = $(cu12.m.λ)
- Transmitance through Ti316 (50 mm) = $trfe1l
- Transmitance through Cu (70 mm) = $trcu1l
- Transmitance through both (70 mm + 50 mm) = $(trcu1l * trfe1l)
- Transmitance through  (120 mm cu) = $(trcu2l)
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
cyl =Cylinder(5.0, 0., 10.0) # draws only the barrel 

# ╔═╡ c7d4775a-f585-4852-a172-337655cb21ed
stdsk = Cylinder(100.0, 0., 20.0)

# ╔═╡ f99d7b04-2fc6-48c5-a846-601e2e270691
gp = [0.0,0.0,(stdsk.zmax + stdsk.zmin)/2]

# ╔═╡ ac011074-fe75-4296-a926-a59f236825c8
Ti316 = Cylinder(500.0, 0., zsteel)

# ╔═╡ 29425331-b2bb-4fe2-9909-9ae74dafb349
cucyl = Cylinder(500.0, zsteel, zcu)

# ╔═╡ 2f01a131-d184-4ea8-97c1-89309716383f
struct PhysicalCylinder
	cyl::Cylinder
	units::Unitful.Length
end

# ╔═╡ 5e8c2718-7168-4ab6-8080-039a8c4d8806
dskti316 = PhysicalCylinder(Ti316, 1.0mm)

# ╔═╡ 50fa5b52-ab70-4746-a11f-15b09c735807
dsktcu = PhysicalCylinder(cucyl, 1.0mm)

# ╔═╡ 053952a1-0558-447c-a3f5-ae30c2d2f7be
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


# ╔═╡ 40cff386-97de-4e6c-8d22-216302831893
ry = Ray([1.0, 0.0, 0.0], [0.5, 0.5, 1.0])

# ╔═╡ 4d17bf89-d01b-4f1b-9dc7-9dcd00ee1d8a
rz = Ray([1.0, 0.0, 0.0], [0., 0.1, 1.0])

# ╔═╡ 991a026a-bfe2-4bc8-a200-31dd25a881d1
md"""
# Functions
"""

# ╔═╡ 2facda59-dc5a-4951-8d6f-1d121690aad2
"""
Propagates ray r along the distance defined by t
"""
ray(r::Ray, t::Float64) =r.e + t * r.d


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



# ╔═╡ 6f6cd50f-d4fd-48d8-8b53-7edf9112c45b
"""
Length of physical cylinder
"""
clength(c::PhysicalCylinder) = clength(c.cyl) * c.units

# ╔═╡ 9c3755c2-d7da-4bfe-8a7a-3af76af4aec0
clength(dskti316)

# ╔═╡ e9d59365-ea79-4e97-b908-cbd528dc0804
"""
Perimeter of cylinder
"""
perimeter(c::Cylinder) = 2π * c.r



# ╔═╡ c66ca1b8-1ae8-4b53-b828-88cac5d8eae9
"""
Perimeter of physical cylinder
"""
perimeter(c::PhysicalCylinder) = perimeter(c.cyl) * c.units

# ╔═╡ da14195e-414e-40ab-aca1-83d2d68b9172
"""
Cylinder area (barrel)
"""
area_barrel(c::Cylinder) = perimeter(c) * clength(c)


# ╔═╡ 0f55b2c3-db7b-4aad-bd8d-c678a5db0b02
"""
Physical Cylinder area (barrel)
"""
area_barrel(c::PhysicalCylinder) = perimeter(c) * clength(c)

# ╔═╡ f41fbfdf-4cee-4b82-b087-e224c70b7f6b
"""
Cylinder area (endcap)
"""
area_endcap(c::Cylinder) = π * c.r^2
        


# ╔═╡ 8da08d4a-f0cf-4b36-b634-b9f8716a3184
"""
Physical Cylinder area (endcap)
"""
area_endcap(c::PhysicalCylinder) = area_endcap(c.cyl) * c.units^2

# ╔═╡ 7f8c81be-0f2d-4214-b27d-a7470534d867
"""
Cylinder area (total)
"""
area(c::Cylinder) = area_barrel(c) + area_endcap(c)
       


# ╔═╡ e1f26c48-0c0c-480c-9574-96b29b65ccf7
"""
Physical Cylinder area (total)
"""
area(c::PhysicalCylinder) = area_barrel(c) + area_endcap(c)

# ╔═╡ f6ceb86d-dfa5-429d-a076-92b836471942
area(dskti316)

# ╔═╡ 372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
"""
Cylinder volume 
"""
volume(c::Cylinder) = area_endcap(c) * clength(c)


# ╔═╡ 834027c9-0a39-450f-9996-a1ec56f080a9
"""
Physical Cylinder volume 
"""
volume(c::PhysicalCylinder) = area_endcap(c) * clength(c)

# ╔═╡ 38879bba-124c-4a06-898c-613e63642240
volume(dskti316)

# ╔═╡ 9022e8da-fe82-4631-8e83-549dad785f9a
"""
Mass of a cylinder filled with radioactive material mat
"""
mass(c::PhysicalCylinder, mat::RadioactiveMaterial) = uconvert(g, volume(c) * mat.m.ρ)

# ╔═╡ 0a42e845-6f53-4ba8-9826-4165030adb96
mti316kg = uconvert(kg, mass(dskti316, mti316))

# ╔═╡ 8569dc46-642f-4709-9a78-0292f5aa1613
"""
Activity (Bi214) of a cylinder filled with radioactive material mat
"""
a_bi214(c::PhysicalCylinder, mat::RadioactiveMaterial) = uconvert(Bq, mat.a_bi214 * mass(c, mat))


# ╔═╡ 557eb397-665a-4475-a204-05887ce0aa45
abi214ti316Bq = a_bi214(dskti316, mti316)

# ╔═╡ 1061c7ff-fbf4-4dd4-bfc3-fb2c38022d87
"""
Activity (Tl208) of a cylinder filled with radioactive material mat
"""
a_tl208(c::PhysicalCylinder, mat::RadioactiveMaterial) = uconvert(Bq, mat.a_tl208 * mass(c, mat))


# ╔═╡ 25b9e691-ef24-4e01-b9e9-91f31b2bf61b
atl208ti316Bq = a_tl208(dskti316, mti316)

# ╔═╡ 0abeba90-ba75-4231-9dd8-d35400964c60
 """Returns True if point is in end-caps of cyinder
 This means that the third coordinate of p (p[3], is close to zmin or to zmax)
 """
in_endcaps(c::Cylinder, p::Vector{Float64}) = isapprox(p[3], c.zmin, atol=1e-7) || isapprox(p[3], c.zmax, atol=1e-7) 

# ╔═╡ 464a3cf9-c030-4270-bae7-bf4339ed3fc8
 """Returns True if point is in barrel of cyinder
 This means that the point verifies the equation of cylinder
 """
in_barrel(c::Cylinder, p::Vector{Float64}) = isapprox(cylinder_equation(c, p), 0.0, atol= 1.0e-7)

# ╔═╡ f675ebdf-e3f6-4f6b-bdae-e59db854f63c
function test_extrap(gtdf)
	cextp = []
	CE = []
	ZMX = []
	for indx in 1:size(gtdf)[1]
		dgt = [gtdf[indx,"gex"], gtdf[indx,"gey"], gtdf[indx,"gez"]]
		append!(ZMX, min(abs(dgt[3]-stdsk.zmin), abs(dgt[3]-stdsk.zmax)))
		append!(CE, cylinder_equation(stdsk, dgt))
		append!(cextp,in_endcaps(stdsk, dgt) || in_barrel(stdsk, dgt))
		
	end
	#df1 = insertcols(gtdf, :cextp => cextp)
	#df2 = insertcols(df1, :ce => CE)
	#df3 = insertcols(df2, :zm => ZMX)
	all(y->y==true, cextp)
end

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

# ╔═╡ b78b7410-886f-47db-b8a3-b8a7feeec341
lines(-10:0.01:10, ff)

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

# ╔═╡ 24cc705d-07dd-49bf-9542-49ce3b7c4855
"""
Transport gammas through the copper shield

"""
function transport_gammas_cs(c::Cylinder, df::DataFrame)
	nphotons = size(df)[1] # these are the points in df
    
    vex = zeros(nphotons)
	vey = zeros(nphotons)
	vez = zeros(nphotons)
	gex = zeros(nphotons)
	gey = zeros(nphotons)
	gez = zeros(nphotons)
	gdx = zeros(nphotons)
	gdy = zeros(nphotons)
	gdz = zeros(nphotons)

    for i in 1:nphotons
    	d = [df[i,"gdx"], df[i,"gdy"], df[i,"gdz"]] # direction
		p = [df[i,"gex"], df[i,"gey"], df[i,"gez"]] # position
        r = Ray(p,d)
        irr = ray_intersection_with_cylinder(r, c)
		vex[i] = p[1]
		vey[i] = p[2]
		vez[i] = p[3]
        gex[i] = irr.e[1]
		gey[i] = irr.e[2]
		gez[i] = irr.e[3]
		gdx[i] = irr.d[1]
		gdy[i] = irr.d[2]
		gdz[i] = irr.d[3]
	end
    DataFrame(vex=vex, vey=vey, vez=vez,
		      gex=gex, gey=gey, gez=gez,
		      gdx=gdx, gdy=gdy, gdz=gdz)
end

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
function vectors_spherical2(npoints::Integer, eps=1e-5)

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
GLMakie.scatter(gdf.vx, gdf.vy, gdf.vz)

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
gtdf = gtransport(stdsk, gp, np)

# ╔═╡ 23556eb1-2cad-47cf-b686-14d6f648ca9f
tgdt = test_extrap(gtdf)

# ╔═╡ 610fdd44-266d-4d1f-8c20-d6ad6ddb6af3
"""
Transport gammas in dataframe df to the surface of the cylinder

"""
function transport_gammas(c::Cylinder, df::DataFrame)
	nphotons = size(df)[1] # these are the points generated in cylinder
    
    vex = zeros(nphotons)
	vey = zeros(nphotons)
	vez = zeros(nphotons)
	gex = zeros(nphotons)
	gey = zeros(nphotons)
	gez = zeros(nphotons)
	gdx = zeros(nphotons)
	gdy = zeros(nphotons)
	gdz = zeros(nphotons)

    for i in 1:nphotons
		rdf = vectors_spherical2(1) # throw just one random
    	d = [rdf[1,"vx"], rdf[1,"vy"], rdf[1,"vz"]] # direction
		p = [df[i,"gex"], df[i,"gey"], df[i,"gez"]]
        r = Ray(p,d)
        irr = ray_intersection_with_cylinder(r, c)
		vex[i] = p[1]
		vey[i] = p[2]
		vez[i] = p[3]
        gex[i] = irr.e[1]
		gey[i] = irr.e[2]
		gez[i] = irr.e[3]
		gdx[i] = irr.d[1]
		gdy[i] = irr.d[2]
		gdz[i] = irr.d[3]
	end
    DataFrame(vex=vex, vey=vey, vez=vez,
		      gex=gex, gey=gey, gez=gez,
		      gdx=gdx, gdy=gdy, gdz=gdz)
end

# ╔═╡ 96f53108-a469-4449-9563-0d4e193cebcc
"""
Generate points inside a cylinder.
"""
function points_in_cylinder(cyl::Cylinder, npoints::Integer; seed=12345)
	Random.seed!(seed)
	r = cyl.r * sqrt.(rand(Float64, npoints)) # random points along r
	theta = 2π * rand(Float64, npoints) # random points along theta
	tr = 2π * rand(Float64, npoints)
	vx =  r .* cos.(theta)
	vy =  r .* sin.(theta)
	#vx = cyl.r * rand(Uniform(-1,1), npoints)
	#vy = cyl.r * rand(Uniform(-1,1), npoints)
	vz = clength(cyl) * rand(Float64, npoints)
	DataFrame(ID=1:npoints, gex=vx, gey=vy, gez=vz)
end

# ╔═╡ de326cb0-7258-4ce1-8fa3-50473d72801b
pcyl = points_in_cylinder(stdsk, 2000)

# ╔═╡ c7ccfaf5-e02e-4061-a3f1-e51cde623ade
gti = points_in_cylinder(Ti316, ngam)

# ╔═╡ fc8e6016-81e3-4798-83e1-f74f77e25489
tgti = transport_gammas(Ti316, gti)

# ╔═╡ 8839a652-22bb-439a-8878-1c17b69293a8
tigz = @chain tgti begin
    @rsubset isapprox(Ti316.zmax, :gez, atol=1e-7)
    @rtransform :d = norm([:gex-:vex, :gey-:vey, :gez-:vez])
	@rtransform :w = transmittance_at_qbb(mti316, :d)
end

# ╔═╡ c0fd4dd9-af77-48c3-8fc7-d329ec130a3d
scatter(tigz.d, tigz.w)

# ╔═╡ 01656238-3a59-4e9f-9895-0b00a1261f22
satt = sum(tigz.w) / size(tigz)[1]

# ╔═╡ 318c611f-a5ed-4b02-82ba-4956c93787e7
md"""
- The self-attenuation in the steel disk = $satt
"""

# ╔═╡ 41733909-c220-492d-a13c-d518e95159e7
cug = transport_gammas_cs(cucyl, tigz)

# ╔═╡ 5389442b-c311-4f35-a5b4-c3c1ca08a014
cugz = @chain cug begin
    @rsubset isapprox(cucyl.zmax, :gez, atol=1e-7)
    @rtransform :d = norm([:gex-:vex, :gey-:vey, :gez-:vez])
	@rtransform :w = transmittance_at_qbb(cu12, :d)
end

# ╔═╡ 2e9d61d5-78c3-4b22-9bef-3c4e57ba0d3e
scatter(cugz.d, cugz.w)

# ╔═╡ eb773df6-1aa1-4b7f-ad2b-1ba1440c99c5
cuatt = sum(cugz.w) / size(cugz)[1]

# ╔═╡ 8fc15be3-58eb-4300-abea-59c22216ff76
md"""
- The attenuation in the cu disk = $cuatt
"""

# ╔═╡ aa3e59b4-7a64-4ea5-8190-530cc632ca9a
toatt = sum(cugz.w)/ngam

# ╔═╡ 03a0dc11-4011-41a7-84dd-b753dc81771f
md"""
- The total attenuation (including geometry) = $toatt
"""

# ╔═╡ 7ec540e6-35d1-44b8-a904-2aa049914d54
tabi214ti316mBq = uconvert(mBq, abi214ti316Bq * toatt)

# ╔═╡ f2149fdb-4abd-475a-b351-5cbcf9318c11
tatl208ti316mBq = uconvert(mBq, atl208ti316Bq * toatt)

# ╔═╡ 9feb2372-b1de-4894-9997-b95a5f8a8b42
md"""
- Mass of reinforcement plate (kg) = $(round(mti316kg/kg, sigdigits=2))
- Acitivity (Bi214) of reinforcement plate (Bq) = $(round(abi214ti316Bq/Bq, sigdigits=2))
- Acitivity (Tl208) of reinforcement plate (Bq) = $(round(atl208ti316Bq/Bq, sigdigits=2))
- Transmited Acitivity (Bi214) of reinforcement plate (mBq) = $(round(tabi214ti316mBq/mBq, sigdigits=2))
- Transmited Acitivity (Tl208) of reinforcement plate (mBq) = $(round(tatl208ti316mBq/mBq, sigdigits=2))
"""

# ╔═╡ 4b54ca12-fca4-4542-9011-34d62cdbcc24
begin
	f = Figure()
	stephist(f[1, 1], tigz.d, bins = 10)
	stephist(f[1, 2], cugz.d, bins = 10)
	f
end

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
function cyl_unit_vectors(p0::Vector{Float64}, p1::Vector{Float64})
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

# ╔═╡ fdcde3b3-f80d-4852-9734-aa3b232e6b1d
"""
Parameterize the three surfaces of a cylinder
"""
function cylinder_surfaces_(c::Cylinder, np, nz, nr)

	mag, v,  n1, n2 = cyl_unit_vectors(c.p0, c.p1)
	
	#surface ranges over t from 0 to length of axis and 0 to 2*pi
	t = collect(range(0., mag, nz))
	theta = collect(range(0., 2π, np)) 
	rsample = collect(range(0., c.r, 2)) 

	#use meshgrid to make 2d arrays
	t, theta2 = jmgrid(t, theta)

	rsample,theta = jmgrid(rsample, theta)

	#generate coordinates for surface
	# "Tube"
	#X,Y,Z = [c.p0[i] .+ v[i] * t .+ c.r * sin.(theta2) * n1[i] .+ c.r * cos.(theta2) *  n2[i] for i in 1:3]
	
	X,Y,Z = [v[i] * t .+ c.r * sin.(theta2) * n1[i] .+ c.r * cos.(theta2) *  n2[i] for i in 1:3]
	
	# "Bottom"
	
	#X2, Y2, Z2 = [c.p0[i] .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	X2, Y2, Z2 = [rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	
	# "Top"
	#X3, Y3, Z3 = [c.p0[i] .+ v[i]*mag .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	X3, Y3, Z3 = [v[i]*mag .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	
	t, theta2, theta, rsample, (X,Y,Z), (X2, Y2, Z2), (X3, Y3, Z3)
end


# ╔═╡ 5a74d977-a44e-4cba-b738-7bdd76056686
"""
Shows explicitly how to parameterize and draw the cylinder barrel
"""
function draw_cylinder_barrel(c::Cylinder, ptheta=100, pz=100)
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

# ╔═╡ 252b84b2-aa7c-4d17-887b-241d0302e446
"""
Draw the three surfaces of a cylinder (barrel, bottom, top) using GLMakie
"""
function draw_cylinderx(cl::Cylinder, np=100, nz=10, 
	                    nr=10; viewa=[1,0,1], viewr=π,
                        showbarrel=true, showendcups=true)


	# Get the surfaces 
	ts,theta2s,thetas, rs, P1, P2, P3 = cylinder_surfaces_(cl, np, nz, nr)
	
	scene = Scene()
	cam3d!(scene)
	
	if showbarrel 
		sfs = GLMakie.surface!(scene, P1..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
	end
	if showendcups
		sfs2 = GLMakie.surface!(scene, P2..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
		sfs3 = GLMakie.surface!(scene, P3..., transparency = true, fxaa = true, overdraw=true, colormap=:greens)
	end
	
	#GLMakie.rotate!(scene, Vec3f(1, 0, 1), 1.5)
	GLMakie.rotate!(scene, Vec3f(viewa...), viewr)
	center!(scene)
	scene
end

# ╔═╡ 911d4179-909c-48f5-9a40-f23efe044498
draw_cylinderx(cyl, 200) #draws the three surfaces (cylinder is transparent)

# ╔═╡ 646c967d-145a-4f2a-9330-c4dbca8d78fc
draw_cylinderx(cyl, 200, showbarrel=true, showendcups=false,viewa=[1,0,0], viewr=π)

# ╔═╡ 4a047576-ab00-4821-9383-1ba945f39285
draw_cylinderx(cyl, 200, showbarrel=false, showendcups=true, viewa=[1,0,0], viewr=π)

# ╔═╡ 476bc4d8-7e22-4342-92ad-4896e47d437f
draw_cylinderx(stdsk, 200)

# ╔═╡ e515fae1-2ea9-4a0d-8064-0588b33e48e2
draw_cylinderx(stdsk, 200,showbarrel=false)

# ╔═╡ 7a1877d8-8ed1-4706-8027-810bd8a2e2f3
draw_cylinderx(Ti316, 200)

# ╔═╡ 72a529b6-01b2-44c4-b24c-139877cbb5d3
"""
Draw the three surfaces of a cylinder (barrel, bottom, top) using GLMakie
"""
function draw_cylinder2x(cl1::Cylinder, cl2::Cylinder, np=100, nz=10, 
	                    nr=10; viewa=[1,0,1], viewr=π,
                        showbarrel=true, showendcups=true)


	# Get the surfaces 
	_,_,_, _, P11, P21, P31 = cylinder_surfaces_(cl1, np, nz, nr)
	_,_,_, _, P12, P22, P32 = cylinder_surfaces_(cl2, np, nz, nr)
	scene = Scene()
	cam3d!(scene)
	
	if showbarrel 
		sfs = GLMakie.surface!(scene, P11..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
		sfs = GLMakie.surface!(scene, P12..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
	end
	if showendcups
		sfs2 = GLMakie.surface!(scene, P21..., transparency = true, fxaa = true, overdraw=false, colormap=:rainbow)
		sfs2 = GLMakie.surface!(scene, P21..., transparency = true, fxaa = true, overdraw=false, colormap=:greens)
		sfs3 = GLMakie.surface!(scene, P31..., transparency = true, fxaa = true, overdraw=true, colormap=:grays)
		sfs3 = GLMakie.surface!(scene, P32..., transparency = true, fxaa = true, overdraw=true, colormap=:bluesreds)
	end
	
	#GLMakie.rotate!(scene, Vec3f(1, 0, 1), 1.5)
	GLMakie.rotate!(scene, Vec3f(viewa...), viewr)
	center!(scene)
	scene
end

# ╔═╡ faa976e3-cf96-4141-ad41-8ed50052cf59
draw_cylinder2x(Ti316, cucyl, 200, viewa=[1,0,0], viewr=π, showbarrel=true, showendcups=false)

# ╔═╡ a41e2bc7-533c-40f7-b69f-94171750a44a
"""
Draw the three surfaces of a cylinder (barrel, bottom, top) using GLMakie
"""
function draw_rays_cylinderx(cl::Cylinder, gtdf, np=100, nz=10, 
	                    nr=10; viewa=[1,0,1], viewr=π,
                        showbarrel=true, showendcups=true, showrays=true)


	# Get the surfaces 
	ts,theta2s,thetas, rs, P1, P2, P3 = cylinder_surfaces_(cl, np, nz, nr)
	
	scene = Scene()
	cam3d!(scene)

	if showbarrel 
		sfs = GLMakie.surface!(scene, P1..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
	end
	if showendcups
		sfs2 = GLMakie.surface!(scene, P2..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
		sfs3 = GLMakie.surface!(scene, P3..., transparency = true, fxaa = true, overdraw=true, colormap=:greens)
	end
	if showrays
		GLMakie.scatter!(scene, gtdf.gex, gtdf.gey, gtdf.gez, color=:black, overdraw=false)
	end
	
	GLMakie.rotate!(scene, Vec3f(viewa...), viewr)
	center!(scene)
	scene
end

# ╔═╡ 09ec4f0d-bd08-4fee-83b4-40844a2b7236
draw_rays_cylinderx(stdsk, gtdf)

# ╔═╡ 7dbeff2f-2b28-4adf-b79a-d6b9b8c000d5
draw_rays_cylinderx(stdsk, gtdf, viewa=[1,0,0], viewr=2π)

# ╔═╡ 09e9e70d-0236-4c3e-be44-dad558229316
draw_rays_cylinderx(stdsk, pcyl, viewa=[1,0,0], viewr=2π, showbarrel=false, showendcups=true)

# ╔═╡ 1c19a47f-82f5-462b-8782-b9d392cbd5be
draw_rays_cylinderx(Ti316, gti, viewa=[1,0,0], viewr=2π, showbarrel=false, showendcups=true)

# ╔═╡ 78059f1f-5056-4da5-94b3-f995413ebbe4


# ╔═╡ 2f45109c-70d4-4050-99d8-5c999fb94ec3


# ╔═╡ 73d06337-6871-4e10-8c80-fac2346f8db4
function lines_in_3D()
    #seed!(123)
    n = 10
    x, y, z = randn(n), randn(n), randn(n)
	println("x =",x)
    aspect=(1, 1, 1)
    perspectiveness=0.5
    # the figure
    fig = Figure(; resolution=(1200, 500))
    ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
    ax2 = Axis3(fig[1, 2]; aspect, perspectiveness)
    ax3 = Axis3(fig[1, 3]; aspect=:data, perspectiveness)
    lines!(ax1, x, y, z; color=1:n, linewidth=3)
    scatterlines!(ax2, x, y, z; markersize=15)
    hm = meshscatter!(ax3, x, y, z; markersize=0.2, color=1:n)
    lines!(ax3, x, y, z; color=1:n)
    #Colorbar(fig[2, 1], hm; label="values", height=15, vertical=false,
     #   flipaxis=false, ticksize=15, tickalign=1, width=Relative(3.55/4))
    fig
end

# ╔═╡ 56e32497-e335-4ca4-b1c8-3ec9ebcef24f
lines_in_3D()

# ╔═╡ Cell order:
# ╠═a298b03b-c7c4-4114-9984-d37e7df8d4c8
# ╠═5158389a-fca1-11ed-05a1-c7345c9bd354
# ╠═e24e5710-04b9-49c3-b3f6-2b136a1744eb
# ╠═e36751ff-a521-4529-95c8-6bbfc3314e66
# ╠═df3dde7a-9437-409f-bd51-db975f74b78c
# ╠═fc869138-1cf7-463b-aca5-729741f8b1ff
# ╠═c7583c88-4cfc-49e1-b61e-c299801be450
# ╠═3c398d08-ef94-4202-9b87-62aa4fba07ed
# ╠═911d4179-909c-48f5-9a40-f23efe044498
# ╠═646c967d-145a-4f2a-9330-c4dbca8d78fc
# ╠═4a047576-ab00-4821-9383-1ba945f39285
# ╠═c48c5f53-2632-4cfc-b07b-7fde9f1d8eb6
# ╠═59138131-46c1-4ac8-be0b-42f5ca2e30b8
# ╠═02a7fd08-440c-4e0c-93c8-d42792642d11
# ╠═40cff386-97de-4e6c-8d22-216302831893
# ╠═bca3de2e-74a6-423e-a660-3845c3e9657b
# ╠═7f4ab967-b5c2-46d9-90c4-ebbbdab0dbfe
# ╠═b78b7410-886f-47db-b8a3-b8a7feeec341
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
# ╠═b4110420-ef1c-4e1b-868a-68bd9d943fea
# ╠═c80fa6be-29d5-4296-982a-f28a7a5bba90
# ╠═5e866260-a798-4c5f-b407-38d091d0ed07
# ╠═cef79ccd-b9ac-4a63-9729-bfbc23906c6e
# ╠═be8dd99e-60f5-4d75-959c-cffbd5bfa4b8
# ╠═c7d4775a-f585-4852-a172-337655cb21ed
# ╠═476bc4d8-7e22-4342-92ad-4896e47d437f
# ╠═e515fae1-2ea9-4a0d-8064-0588b33e48e2
# ╠═8a45b62d-cbfb-4f8c-9c61-fad0e3ca53da
# ╠═f99d7b04-2fc6-48c5-a846-601e2e270691
# ╠═128d5ed6-3a09-427e-b50d-159d0288d1dd
# ╠═47ff0204-6b05-4635-89fa-23dac3a68500
# ╠═f675ebdf-e3f6-4f6b-bdae-e59db854f63c
# ╠═23556eb1-2cad-47cf-b686-14d6f648ca9f
# ╠═5cca07ea-8740-4899-bca7-629bcdb8d46b
# ╠═09ec4f0d-bd08-4fee-83b4-40844a2b7236
# ╠═7dbeff2f-2b28-4adf-b79a-d6b9b8c000d5
# ╠═775a2389-f62d-4341-b0e6-976d150207b1
# ╠═de326cb0-7258-4ce1-8fa3-50473d72801b
# ╠═09e9e70d-0236-4c3e-be44-dad558229316
# ╠═2a795370-8b35-4db5-8f5c-c03b3efaa1de
# ╠═ea0d7d4d-5dcf-4d7f-b8bc-dd4856957df5
# ╠═b2a2c523-7461-4637-87c5-7ecd321dc7a5
# ╠═e500854c-1714-44c6-8c76-107950d65060
# ╠═ac011074-fe75-4296-a926-a59f236825c8
# ╠═7a1877d8-8ed1-4706-8027-810bd8a2e2f3
# ╠═92e67cdb-f4c7-4ada-9939-332a1fa3546d
# ╠═f1255b5c-d28a-4e48-9457-104a36d4596d
# ╠═c7ccfaf5-e02e-4061-a3f1-e51cde623ade
# ╠═1c19a47f-82f5-462b-8782-b9d392cbd5be
# ╠═7a70c468-156b-47c4-a2af-2a6c30034191
# ╠═fc8e6016-81e3-4798-83e1-f74f77e25489
# ╠═4559bbcc-4e78-4e29-b0b2-0063be1335f9
# ╠═8839a652-22bb-439a-8878-1c17b69293a8
# ╠═c0fd4dd9-af77-48c3-8fc7-d329ec130a3d
# ╠═01656238-3a59-4e9f-9895-0b00a1261f22
# ╠═318c611f-a5ed-4b02-82ba-4956c93787e7
# ╠═f329a7b4-75db-46f7-b137-5ed2ebaea571
# ╠═9d72493b-1fbb-46f3-a670-081ea8484af9
# ╠═354cb5bc-c5a9-4f6e-94f0-5b12390b8922
# ╠═29425331-b2bb-4fe2-9909-9ae74dafb349
# ╠═faa976e3-cf96-4141-ad41-8ed50052cf59
# ╠═41733909-c220-492d-a13c-d518e95159e7
# ╠═5389442b-c311-4f35-a5b4-c3c1ca08a014
# ╠═2e9d61d5-78c3-4b22-9bef-3c4e57ba0d3e
# ╠═4b54ca12-fca4-4542-9011-34d62cdbcc24
# ╠═eb773df6-1aa1-4b7f-ad2b-1ba1440c99c5
# ╠═8fc15be3-58eb-4300-abea-59c22216ff76
# ╠═aa3e59b4-7a64-4ea5-8190-530cc632ca9a
# ╠═03a0dc11-4011-41a7-84dd-b753dc81771f
# ╠═f211e72a-bc1d-4798-a2b8-6cb39663f4a5
# ╠═3a182a84-25f7-4bb2-bb54-6e81c53142e7
# ╠═2b73f038-b92e-473c-8201-e2c12913cebf
# ╠═326a7877-e625-468a-8a7a-e9ff36e2f6fb
# ╠═5e8c2718-7168-4ab6-8080-039a8c4d8806
# ╠═50fa5b52-ab70-4746-a11f-15b09c735807
# ╠═9c3755c2-d7da-4bfe-8a7a-3af76af4aec0
# ╠═f6ceb86d-dfa5-429d-a076-92b836471942
# ╠═38879bba-124c-4a06-898c-613e63642240
# ╠═0a42e845-6f53-4ba8-9826-4165030adb96
# ╠═557eb397-665a-4475-a204-05887ce0aa45
# ╠═25b9e691-ef24-4e01-b9e9-91f31b2bf61b
# ╠═7ec540e6-35d1-44b8-a904-2aa049914d54
# ╠═f2149fdb-4abd-475a-b351-5cbcf9318c11
# ╠═9feb2372-b1de-4894-9997-b95a5f8a8b42
# ╠═fb431946-28d1-4042-9b18-aa4439e48d84
# ╠═9148483a-c144-4a10-8606-d73ab695abaa
# ╠═ce9da661-5f2f-4e5d-825e-5d029c274066
# ╠═24b1b8ef-10ee-4d96-8121-827234233873
# ╠═6458b9a6-e0b1-40cf-87bd-3e3f01990525
# ╠═478e01bc-d764-4af7-a6f4-661cdba4820d
# ╠═efe04df4-cb56-4916-96f7-296dc37492ac
# ╠═2f01a131-d184-4ea8-97c1-89309716383f
# ╠═053952a1-0558-447c-a3f5-ae30c2d2f7be
# ╠═991a026a-bfe2-4bc8-a200-31dd25a881d1
# ╠═2facda59-dc5a-4951-8d6f-1d121690aad2
# ╠═539faa01-299b-4c5b-8113-fc21501dd84d
# ╠═269e4b09-3505-40fe-8815-2bc847d02a99
# ╠═f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
# ╠═610fdd44-266d-4d1f-8c20-d6ad6ddb6af3
# ╠═24cc705d-07dd-49bf-9542-49ce3b7c4855
# ╠═c26c8511-afae-4447-9458-d1d21f565911
# ╠═deefebf9-07a3-4cce-b2d5-4b42644b5db7
# ╠═6f6cd50f-d4fd-48d8-8b53-7edf9112c45b
# ╠═e9d59365-ea79-4e97-b908-cbd528dc0804
# ╠═c66ca1b8-1ae8-4b53-b828-88cac5d8eae9
# ╠═da14195e-414e-40ab-aca1-83d2d68b9172
# ╠═0f55b2c3-db7b-4aad-bd8d-c678a5db0b02
# ╠═f41fbfdf-4cee-4b82-b087-e224c70b7f6b
# ╠═8da08d4a-f0cf-4b36-b634-b9f8716a3184
# ╠═7f8c81be-0f2d-4214-b27d-a7470534d867
# ╠═e1f26c48-0c0c-480c-9574-96b29b65ccf7
# ╠═372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
# ╠═834027c9-0a39-450f-9996-a1ec56f080a9
# ╠═9022e8da-fe82-4631-8e83-549dad785f9a
# ╠═8569dc46-642f-4709-9a78-0292f5aa1613
# ╠═1061c7ff-fbf4-4dd4-bfc3-fb2c38022d87
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
# ╠═96f53108-a469-4449-9563-0d4e193cebcc
# ╠═d9ef36e2-56f6-4655-98f2-37c2699363b5
# ╠═32b53e87-034e-49db-842f-58f0c0002ed7
# ╠═cb1d79b8-ada9-481e-8169-ab5976cb0759
# ╠═fdcde3b3-f80d-4852-9734-aa3b232e6b1d
# ╠═5a74d977-a44e-4cba-b738-7bdd76056686
# ╠═252b84b2-aa7c-4d17-887b-241d0302e446
# ╠═72a529b6-01b2-44c4-b24c-139877cbb5d3
# ╠═a41e2bc7-533c-40f7-b69f-94171750a44a
# ╠═78059f1f-5056-4da5-94b3-f995413ebbe4
# ╠═2f45109c-70d4-4050-99d8-5c999fb94ec3
# ╠═73d06337-6871-4e10-8c80-fac2346f8db4
# ╠═56e32497-e335-4ca4-b1c8-3ec9ebcef24f
