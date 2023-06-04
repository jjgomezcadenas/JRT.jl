using GLMakie
using LinearAlgebra
import GeometryBasics
using Colors


"""
Given two vectors, jmgrid returns the julia equivalent 
to the meshgrid function.
If transp = true it returns the exact equivalent, 
which is the matrix transposed of
the usual matrix in julia. In Julia, contrary to Python, 
the first array index is a row and the second 
is a column (Julia has column-major arrays)

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

"""Normal to the cylinder barrel
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
then n = Grad(F)_P /Norm(Grad(F))_P
n = (2x, 2y, 0)/sqrt(4x^2 + 4y^2) = (x,y,0)/r (P)
"""
normal_to_barrel(c::Cylinder, P::Vector{Float64}) = [P[1], P[2], 0] ./ c.r


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


"""
Draw the three surfaces of two cylinders (barrel, bottom, top) using GLMakie
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

"""
Draw the points defined by gtdf together with 
three surfaces of a cylinder (barrel, bottom, top) using GLMakie

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