using GadgetIO

"""
    fill_background_with_particles(ic::String;
                                   res=nothing,box_size=nothing,ndim=3)
    
Adapt format 2 IC file to be used with AREPO, fill background with low resolution, low density and pressure particles to improve mesh reconstruction at the boundary.
"""
function fill_background_with_particles(ic::String; # name of ic file (format 2!)
                                        res=nothing, # default: get resolution from ic file
                                        box_size=nothing, # default: get from snapshot header
                                        ndim=3)
    
    new_ic = split(ic, ".")[1]*"_arepo"
    if length(split(ic,".")) > 1
        for i in 2:length(split(ic,"."))
            new_ic = new_ic*"."*split(ic,".")[i]
        end
    end
    println("read from "*ic*" and write to "*new_ic)
    
    pos = read_block(ic, "POS", parttype=0)
    # now we can test the resolution

    head = read_header(ic)
    if box_size == nothing
        box_size = head.boxsize
    end
    println("BoxSize ",box_size)
    
    n0_particles = length(pos[1,:])

    find_region_to_fill = find_region_to_fill_recursive
    
    if res == nothing # get average resolution of particles
        println("find resolution")
        tmp, n_part_to_add = find_region_to_fill(pos, res=128, box_size=box_size, ndim=ndim) # some average resolution to get an idea of fraction of space occupied
        res = round(Int, (n0_particles / (1 - n_part_to_add/128^ndim))^(1/ndim))
    end
    println("res = ",res," => ",2^round(Int,log2(res)))
    res = 2^round(Int,log2(res))
    
    println("find where to fill")
    overlay_particle_filter, n_part_to_add = find_region_to_fill(pos, res=res, box_size=box_size, ndim=ndim)
     
    n_tot = n0_particles + n_part_to_add
    new_pos = Array{Float64, 2}(undef, 3, n_tot)
    new_pos[:, 1:n0_particles] = pos

    println("fill pos, n_tot = ",n0_particles,"+",n_part_to_add," = ",n_tot)
    # now fill the new positions with regular grid
    id = 1
    if ndim == 3
        for i in 1:res
            for j in 1:res
                for k in 1:res
                    if overlay_particle_filter[i,j,k] == 0 
                        new_pos[1, n0_particles + id] = (i-0.5)/res * box_size
                        new_pos[2, n0_particles + id] = (j-0.5)/res * box_size
                        new_pos[3, n0_particles + id] = (k-0.5)/res * box_size
                        id += 1
                    end
                end
            end
        end
    elseif ndim == 2
        for i in 1:res
            for j in 1:res
                if overlay_particle_filter[i,j] == 0
                    new_pos[1, n0_particles + id] = (i-0.5)/res * box_size
                    new_pos[2, n0_particles + id] = (j-0.5)/res * box_size
                    new_pos[3, n0_particles + id] = 0
                    id += 1
                end
            end
        end
    else
        for i in 1:res
            if overlay_particle_filter[i] == 0
                new_pos[1, n0_particles + id] = (i-0.5)/res * box_size
                new_pos[2, n0_particles + id] = 0
                new_pos[3, n0_particles + id] = 0
                id += 1
            end
        end
    end
        
    # now we can add the other particle properties
    
    println("write HEAD")
    fout = open(new_ic, create=true, write=true)
    head.boxsize=Float64(box_size)
    head.npart[1] = Int32(n_tot)
    head.nall[1] = UInt32(n_tot)
    mass = head.massarr[1]
    head.massarr[1] = Float32(0)
    write_header(fout, head)

    println("write POS")
    for type in 1:5
        if head.npart[type+1] == 0
            continue
        end
        new_pos = hcat(new_pos, read_block(ic, "POS", parttype=type))
    end
    write_block(fout, Float32.(new_pos), "POS")

    println("write VEL")
    new_vel = Array{Float64, 2}(undef, 3, n_tot)
    new_vel[:,1:n0_particles] = read_block(ic, "VEL", parttype=0)
    new_vel[:,n0_particles+1:n_tot] .= 0 # we might want to adjust that
    for type in 1:5
        if head.npart[type+1] == 0
            continue
        end
        new_vel = hcat(new_vel, read_block(ic, "VEL", parttype=type))
    end
    write_block(fout, Float32.(new_vel), "VEL")

    println("write ID")
    new_ids = Array{UInt64}(undef, n_tot)
    new_ids[1:n0_particles] = read_block(ic, "ID", parttype=0)
    new_ids[n0_particles+1:n_tot] = collect(n0_particles+1:n_tot)
    for type in 1:5
        if head.npart[type+1] == 0
            continue
        end
        new_ids = vcat(new_ids, read_block(ic, "ID", parttype=type).+n_part_to_add)
    end
    
    write_block(fout, UInt32.(new_ids), "ID")

    println("write MASS")
    new_mass = Array{Float64}(undef, n_tot)
    new_mass[1:n0_particles] = read_block(ic, "MASS", parttype=0, info=InfoLine("MASS", Float32, 1, [1,1,1,1,1,1]) )
    if any(x->x==0, new_mass[1:n0_particles])
        new_mass[1:n0_particles].= mass
    end
    new_mass[n0_particles+1:n_tot] .= 1e-10 * box_size^ndim / res^ndim # rho = 1e-10
    for type in 1:5
        if head.npart[type+1] == 0
            continue
        end
        new_mass = vcat(new_mass, read_block(ic, "MASS", parttype=type, info=InfoLine("MASS", Float32, 1, [1,1,1,1,1,1]) ))
    end
    write_block(fout, Float32.(new_mass), "MASS")

    println("write U")
    new_u = Array{Float64}(undef, n_tot)
    new_u[1:n0_particles] = read_block(ic, "U", parttype=0)
    new_u[n0_particles+1:n_tot] .= 1e-10
    # u only present for type 0
    write_block(fout, Float32.(new_u), "U")

    println("write RHO")
    try
        new_rho = Array{Float64}(undef, n_tot)
        new_rho[1:n0_particles] = read_block(ic, "RHO", parttype=0)
        new_rho[n0_particles+1:n_tot] .= 1e-10
        # rho only present for type 0
        write_block(fout, Float32.(new_rho), "RHO")
    catch
        println("block rho not present, skip it")
    end
    close(fout)

end
function find_region_to_fill_iterative(pos ; res=128, box_size=1)
    overlay_particle_filter = Array{Float64, 3}(undef, res, res, res) # filter for overlay grid
    
    for i in 1:res # now search for relevant region
        x0 = (i-1)/res * box_size
        x1 = i/res * box_size
        println(i)
        for j in 1:res
            y0 = (j-1)/res * box_size
            y1 = j/res * box_size
            for k in 1:res
                z0 = (k-1)/res * box_size
                z1 = k/res * box_size
                part_in_cell =
                    (x0 .<= pos[1,:] .< x1) .&
                    (y0 .<= pos[2,:] .< y1) .&
                    (z0 .<= pos[3,:] .< z1)
                if sum(part_in_cell) > 0
                    overlay_particle_filter[i,j,k] = 1
                else
                    overlay_particle_filter[i,j,k] = 0
                end
            end
        end
    end
    n_part_to_add = sum(overlay_particle_filter)
    return overlay_particle_filter, n_part_to_add
end
function find_region_to_fill_recursive(pos ; res::Integer=128, box_size=1,
                                       x0=0, y0=0, z0=0,
                                       ndim=3) # much faster than iterative version!

    x1 = x0 + box_size
    y1 = y0 + box_size
    z1 = y0 + box_size
    
    part_in_cell =
        (x0 .<= pos[1,:] .< x1) .&
        (y0 .<= pos[2,:] .< y1) .&
        (z0 .<= pos[3,:] .< z1)

    if ndim == 3
        overlay_particle_filter = Array{Int64, 3}(undef, res, res, res) # filter for overlay grid
    elseif ndim ==2
        overlay_particle_filter = Array{Int64, 2}(undef, res, res) # filter for overlay grid
    else
        overlay_particle_filter = Array{Int64}(undef, res) # filter for overlay grid
    end
    if any(x->x!=0, part_in_cell)
        if res == 1
            overlay_particle_filter .= 1
        else # do it recursively
            half_res = round(Int,res/2)
            if ndim == 3
                for i in 1:2, j in 1:2, k in 1:2
                    overlay_particle_filter[1+(i-1)*half_res:i*half_res,
                                            1+(j-1)*half_res:j*half_res,
                                            1+(k-1)*half_res:k*half_res], = find_region_to_fill_recursive(pos[:,part_in_cell], res=half_res, box_size=box_size/2, x0=x0+(i-1)*box_size/2,y0=y0+(j-1)*box_size/2,z0=z0+(k-1)*box_size/2)
                end
            elseif ndim == 2
                for i in 1:2, j in 1:2
                    overlay_particle_filter[1+(i-1)*half_res:i*half_res,
                                            1+(j-1)*half_res:j*half_res], = find_region_to_fill_recursive(pos[:,part_in_cell], res=half_res, box_size=box_size/2, x0=x0+(i-1)*box_size/2,y0=y0+(j-1)*box_size/2, ndim=2)
                end
            else
                for i in 1:2
                    overlay_particle_filter[1+(i-1)*half_res:i*half_res], = find_region_to_fill_recursive(pos[:,part_in_cell], res=half_res, box_size=box_size/2, x0=x0+(i-1)*box_size/2, ndim=1)
                end
            end
        end
    else
        #println(x0, " ", x1, ",",y0," ",y1,",",z0," ",z1)
        overlay_particle_filter[:] .= 0
    end
    
    n_part_to_add = res^ndim - sum(overlay_particle_filter)
    return overlay_particle_filter, n_part_to_add
end

#fill_background_with_particles("kepler_disk/snap_000",box_size=8,res=128,ndim=2) # effective resolution: 16
#fill_background_with_particles("grav_freefall/freefall.ic",box_size=2,res=16) # effective resolution: 8
#fill_background_with_particles("hydrostatic_sphere/hydro_sphere.ic",box_size=32000.0,res=64)
