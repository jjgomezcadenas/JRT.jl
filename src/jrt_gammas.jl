using DataFrames
using Random
using Distributions
using PDMats
using Statistics

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

"""
Generate three standard normally distributed numbers and normalize the vector.
We have to be careful in the case that the vector has a norm close to zero, 
in which we must worry about floating point precision by dividing by 
a very small number. This is the reason for the while loop.

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

"""
Generate three standard normally distributed numbers 
and normalize the vector.

We have to be careful in the case that the vector 
has a norm close to zero, in which we must worry about 
floating point precision by dividing by a very small number. 
This is the reason for the while loop.

"""
function vectors_spherical2(npoints::Integer, eps=1e-5, seed=12345)
    Random.seed!(seed)
    
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