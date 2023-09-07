using GadgetIO
using LinearAlgebra
using Formatting

using Quadmath # add support for Float128 datatype

import Base.mod
function mod(x::Float128,y::Real)
    while true
        if sign(y)*x < 0
            x += y
        elseif sign(y)*x >= sign(y)*y
            x -= y
        else
            return x
        end
    end
end

global ORDER = 1 # order of soundwave
function make_soundwave( ; res = 64, # resolution
                         grid_type::String="regular", # type of grid. So far only regular, but should include other grids as well!
                         oname::String="soundwave_64.ic", # name of file to create
                         longx=1, longy=0.75, longz=0.75, # size of box
                         amplitude = 1e-4, # amplitude of density perturbation
                         doubleprecision=false, # use doubleprecision datatypes
                         arepo=false) # use proper center of mass for Arepo

    fout = open(oname, create=true, write=true)

    if grid_type == "regular"
        pos = make_regular_grid(res,
                                longx=longx, longy=longy, longz=longz)
    elseif grid_type == "hcp"
        pos = make_hcp_grid(res,
                            longx=longx, longy=longy, longz=longz)
    end
    global dx_offset_collection = Dict{Float64, Float64}()
    rho,  part_mass = make_rho_and_perturb_grid!(pos, amplitude=amplitude, box_volume=longx*longy*longz, res=res)
    
    npart = length(pos[1,:])
    write_header(fout, SnapshotHeader(Int32.([npart, 0,0,0,0,0]), # npart
                                      Float64.([part_mass, 0,0,0,0,0]), # massarr
                                      Float64(0.0), # time
                                      Float64(0.0), # redshift
                                      Int32(0), # flag sfr
                                      Int32(0), # flag feedback
                                      UInt32.([npart, 0,0,0,0,0]), # npart total
                                      Int32(0), #flag_cooling
                                      Int32(1), #num_files
                                      Float64(1), #boxsize
                                      Float64(0.3), #omega0
                                      Float64(0.7), #omegalambda
                                      Float64(1.0), #hubbleparam
                                      Int32(0), # flag stellarage
                                      Int32(0), # flag metals
                                      UInt32.([0, 0,0,0,0,0]),
                                      Int32(0), # flag entropy instead of u
                                      Int32(doubleprecision) # flag doubleprecision
                                      ))
    if doubleprecision
        float_type=Float64
    else
        float_type=Float32
    end
    write_block(fout, float_type.(pos), "POS")

    # for Arepo velocity, internal energy, and density actually have to use particle centers (center of mass)
    pos = if arepo
        tmp = if grid_type == "hcp"
            ""
        else
            "regular_"
        end
        ids = read_block("../../test_runs/out_soundwave_arepo_"*tmp*sprintf1("%d",res)*"/snap_000","ID",parttype=0)
        read_block("../../test_runs/out_soundwave_arepo_"*tmp*sprintf1("%d",res)*"/snap_000","CMCE",parttype=0,info=InfoLine("CMCE",Float32,3,[1,0,0,0,0,0]))[:,sortperm(ids)]
    else
        pos
    end

    vel = make_vel(pos, amplitude=amplitude)
    write_block(fout, float_type.(vel), "VEL")

    ids = collect(1:length(pos[1,:]))
    write_block(fout, UInt32.(ids), "ID")

    u   = make_u(pos, amplitude=amplitude)
    write_block(fout, float_type.(u), "U")

    rho = if arepo
        make_rho(pos, amplitude=amplitude)
    else
        rho
    end
    write_block(fout, float_type.(rho), "RHO")

    hsml = Array{Float64}(undef, length(pos[1,:]))
    hsml .= 0.021784158 * 90/res # value MFM finds
    write_block(fout, float_type.(hsml), "HSML")
    
    close(fout)
    
end
    
function make_regular_grid(res::Integer=64; # desired resolution per boxlength
                           longx=1, longy=1, longz=1) # dimensions of box

    nx = round(Integer,res*longx)
    ny = round(Integer,res*longy)
    nz = round(Integer,res*longz)

    pos = Array{Float64, 2}(undef, 3, nx*ny*nz)
    
    for ix in 1:nx
        for iy in 1:ny
            for iz in 1:nz
                id = 1 + (ix-1) * ny*nz + (iy-1) * nz + (iz-1)
                pos[1,id] = longx * Float64(ix)/Float64(nx)
                pos[2,id] = longy * Float64(iy)/Float64(ny)
                pos[3,id] = longz * Float64(iz)/Float64(nz)
            end
        end
    end
    
    return pos
end
function make_hcp_grid(res::Integer=64; # desired resolution per boxlength
                       longx=1, longy=1, longz=1) # domensions of the box
    ntot = round(Integer, res*longx) * round(Integer, res*longy) * round(Integer, res*longz)
    if mod(ntot,2) != 0
        ntot -= 1
    end

    # find by combining spacings with Ntot
    r = (sqrt(2.0)*longx*longy*longz/8.0 /ntot)^(1.0/3.0)

    # spacings
    dx = 2.0*r
    dy = sqrt(3.0)*r
    dz = sqrt(6.0)*2.0/3.0*r
   
    # particle numbers
    np = zeros(Int, 3)
    #np[:] .= round(Int, ntot^(1/3))
    np[1] = round(Int, res*longx)
    np[2] = round(Int, res*longy)
    np[3] = round(Int, res*longz)
        
    # enforce periodicity
    dx += (longx-dx*np[1])/np[1]
    dy += (longy-dy*np[2])/np[2]
    dz += (longz-dz*np[3])/np[3]
       
    # diagnostics
    println("Particle numbers :")
    println("     Nx = ",np[1])
    println("     Ny = ",np[2])
    println("     Nz = ",np[3])
    println(" Total  = ",np[1]*np[2]*np[3])
    println(" Wanted = ",ntot)
    println(" Delta  = ",ntot-np[1]*np[2]*np[3])

    ntot = np[1]*np[2]*np[3]
    pos = Matrix{Float128}(undef, 3, ntot) # Float64???

    # particle positions
    x = Array{Float128,3}(undef,np[1],np[2],np[3]) # type???
    y = Array{Float128,3}(undef,np[1],np[2],np[3])
    z = Array{Float128,3}(undef,np[1],np[2],np[3])

    idxarr = collect(1:np[1])

    # A(0) plane
    # 0st row
    x[:,1,1] .= r .+ idxarr*dx
    y[:,1,1] .= r
    z[:,1,1] .= r
    
    # 1st row
    x[:,2,1] = idxarr .*dx
    y[:,2,1] .= r + dy
    z[:,2,1] .= r

    # A(0) plane
    # 0st row
    x[:,1,1] .= r .+ idxarr*dx
    y[:,1,1] .= r
    z[:,1,1] .= r
    
    # 1st row
    x[:,2,1] = idxarr .*dx
    y[:,2,1] .= r + dy
    z[:,2,1] .= r

    # A-plane
    # even rows
    for i in 2:2:np[2]-1
        x[:,i+1,1] = x[:,1,1]
        y[:,i+1,1] = y[:,1,1] .+ i*dy
        z[:,i+1,1] = z[:,1,1]
    end
    
    # odd rows
    for i in 3:2:np[2]-1
        x[:,i+1,1] = x[:,2,1]
        y[:,i+1,1] = y[:,2,1] .+ (i-1)*dy
        z[:,i+1,1] = z[:,2,1]
    end

    #B(1)-plane
    # 0rst row
    x[:,1,2] = r .+ idxarr*dx 
    y[:,1,2] .= 0 
    z[:,1,2] .= r + dz

    # 1st row
    x[:,2,2] = idxarr .*dx 
    y[:,2,2] .= dy
    z[:,2,2] .= r + dz

    # B-plane
    # even rows
    for i in 2:2:np[2]-1
        x[:,i+1,2] = x[:,1,2]
        y[:,i+1,2] = y[:,1,2] .+ i*dy
        z[:,i+1,2] = z[:,1,2]
    end
   
    # odd rows
    for i in 3:2:np[2]-1
        x[:,i+1,2] = x[:,2,2]
        y[:,i+1,2] = y[:,2,2] .+ (i-1)*dy
        z[:,i+1,2] = z[:,2,2]
    end
    
    # all planes
    # even planes
    for i in 2:2:np[3]-1
        x[:,:,i+1] = x[:,:,1]
        y[:,:,i+1] = y[:,:,1]
        z[:,:,i+1] = z[:,:,1] .+ i*dz
    end
    
    # odd planes
    for i in 3:2:np[3]-1
        x[:,:,i+1] = x[:,:,2]
        y[:,:,i+1] = y[:,:,2]
        z[:,:,i+1] = z[:,:,2] .+ (i-1)*dz
    end
    
    # make particle arrays
    #pos = Array{Float64}(undef,np[1],np[2],np[3]) # is an input variable here

    bin = 0
    for i in 1:np[1]
        for j in 1:np[2]
            for k in 1:np[3]
                pos[1,bin+1] = x[i,j,k] + dx/2
                pos[2,bin+1] = y[i,j,k] + dy/2
                pos[3,bin+1] = z[i,j,k] + dz/2
                bin+=1
            end
        end
    end
    println("d ",dy," ",dy," ",dz)
    
    return pos
end

function make_rho_and_perturb_grid!(pos::Matrix; # array of positions, will be changed here!
                                   amplitude = 0.01, # amplitude of the density perturbation
                                   box_volume = 1.0*0.75^2,
                                   res = 64,
                                   longx=1)

    mean_rho = 1
    
    rho = Array{Float64}(undef, length(pos[1,:]))
    rho .= mean_rho
    part_mass = Float64(box_volume * mean_rho / length(rho))

    delta_x = (box_volume / length(rho))^(1/3) # mean particle distance in a regular grid
    for i in 1:length(rho)
        pos[1,i] = mod(pos[1,i] - dx_offset(pos[1,i], res, amplitude=amplitude, delta_x=delta_x, longx=longx), longx)
        rho[i] = rho[i] + amplitude * sin( 2 * pi * ORDER * pos[1,i] )
    end
    
    return rho, part_mass
end
function dx_offset(x::AbstractFloat, res=64; longx=1, amplitude=0.01, difference_side::String="left", delta_x=1/64)
    mean_rho = 1
    drho = x -> amplitude * sin( 2 * pi * ORDER * x )

    if ! @isdefined dx_offset_collection
        global dx_offset_collection = Dict{typeof(x), typeof(x)}()
    elseif haskey(dx_offset_collection,x) # value already caluclated once
        return dx_offset_collection[x]
    end
    
    x_arr = sort(mod.( x .+ collect(0:(res-1)) .* longx/res, longx ))
    ddelta_x = drho.(x_arr) ./ mean_rho .* delta_x # desired change in associated particle x-length
    if difference_side == "both"
        dx_offset_matrix = diagm(1 => ones(res-1).*0.5, -1 => ones(res-1).*-0.5)
        dx_offset_matrix[1,res] = -0.5
        dx_offset_matrix[res, 1] = 0.5
    elseif difference_side == "left"
        dx_offset_matrix = diagm( 0 => -ones(res), -1 => ones(res-1))
        dx_offset_matrix[1,res] = 1
    end

    dx_offset_solution = Array{typeof(x)}(undef, length(x_arr))
    fix_point=1
    dx_offset_solution[fix_point] = 0
    for i in fix_point-1:-1:1
        dx_offset_solution[i] = dx_offset_solution[i+1] - ddelta_x[i+1]
    end
    for i in fix_point+1:length(x_arr)
        dx_offset_solution[i] = dx_offset_solution[i-1] + ddelta_x[i] # difference_side=="left"
    end
    # add to known solutions
    dx_offset_collection = merge(dx_offset_collection, Dict(x_arr .=> dx_offset_solution))

    return dx_offset_collection[mod(x,longx)] # return the desired value
    
end

function make_rho(pos::Matrix; # array of positions
                  amplitude = 0.01) # amplitude of the density perturbation

    mean_rho = 1.0
    
    rho = Array{Float64}(undef, length(pos[1,:]))
    rho .= mean_rho

    for i in 1:length(rho)
        rho[i] = rho[i] + amplitude * sin( 2 * pi * ORDER * pos[1,i] )
    end
    
    return rho

end
function make_vel(pos; # array of positions
                  amplitude = 0.01) # amplitude of the density perturbation

    vel = Array{Float64, 2}(undef, 3, length(pos[1,:]))
    for i in 1:length(vel[1,:])
        vel[1,i] = - amplitude * 0.6666666 * sin( 2 * pi * ORDER * pos[1,i] )
        vel[2,i] = 0.0
        vel[3,i] = 0.0
    end
    return vel
end

function make_u(pos; # array of positions
                amplitude = 0.01) # amplitude of the density perturbation

    u = Array{Float64}(undef, length(pos[1,:]))
    for i in 1:length(u)
        u[i] = 0.4 + amplitude * 0.26666666 * sin( 2 * pi * ORDER * pos[1,i] )
    end
    return u
end

for res in [32, 45, 64, 90, 128]
    make_soundwave(res=res, oname="soundwave_"*sprintf1("%d", res)*".ic", grid_type="hcp", doubleprecision=false)
    make_soundwave(res=res, oname="soundwave_"*sprintf1("%d", res)*"_arepo.ic", grid_type="hcp", doubleprecision=false, arepo=true)
end
