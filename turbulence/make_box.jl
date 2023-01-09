# This file contains routines to produce a gadget snapshot 
# of a periodic box with a Kolmogorov velocity power spectrum
# from a glass file. The turbulent velocity contains the fraction
# X_turb of the thermal energy in the box. Temperature and density
# are set to values reasonable for a galaxy cluster atmosphere. 
# Note: This is not a real GADGET snapshot as hsml is set to a 
# reasonable constant. Upon running GADGET is going to recalculate it.
#
# The package GadgetIO is needed.
#
# Note: you have to compile this twice
# based on make_data.pro
# only used for very high resolution

# make initial conditions for a turbulent box

using GadgetIO
using FFTW
#using Random
using RandomNumbers.MersenneTwisters # use the same Random number generator as for IDL! Nevertheless, the first random number is still different, check that! For now, instead use a proper initialization, to improve random behavior

using Formatting

function make_box(npart::Integer=128^3; # number of particles in grid
                  nfiles=1, # number of IC subfiles
                  odir="./", # where to save the IC file
                  X_turb=0.3) # E_turb / E_therm

    println("Npart = ", npart)

# Gadget units & chemistry
    m_unit = 1.989e43       # [10^10 Msol]
    l_unit = 3.085678e21    # [kpc]
    v_unit = 100000         # [km/s]
    
    H_frac = 0.76
    umol = 4.0/(5.0*H_frac+3.0)
    mp = 1.6726231e-24      # proton mass cgs
    k_boltz = 1.3806580e-16 # [cgs]

# input values
    fout_name = odir*"turb.ic"      # output filename
    
    boxsize = 3000          # [kpc]
    
    T = 1e7                 # [K]
    rho = 1e-27/(m_unit/l_unit^3.)  # [GADGET]
    mass = rho * boxsize^3

# output
    println("Writing to : "*fout_name)

# create files
    fout = []
    npart_file = []
    if nfiles == 1
        append!(fout,[open(fout_name, create=true, write=true)])
        append!(npart_file,[round(Integer, npart/nfiles)])
    else
        for i_file in 1:nfiles
            append!(fout,[open(fout_name*"."*sprintf1("%d",i_file-1), create=true, write=true)])
            append!(npart_file,[round(Integer, npart/nfiles)])
        end
        npart_file[nfiles] = npart-sum(npart_file[1:(nfiles-1)])
    end
    filerange=[]
    if nfiles == 1
        append!(filerange,[1:npart])
    else
        append!(filerange,[1:npart_file[1]])
        for i in 2:nfiles
            append!(filerange,[sum(npart_file[1:(i-1)])+1:sum(npart_file[1:i])])
        end
    end
    
# make positions
    pos = Array{Float32,2}(undef,3, npart)
    println("make positions hcp")
    make_positions_hcp!(pos,boxsize, boxsize, boxsize, npart)
    println("done making positions")

    for i_file in 1:nfiles
# make data structures
        head = SnapshotHeader(Int32.([npart_file[i_file],0,0,0,0,0]), #npart = 
	                      Float64.([mass/npart,0,0,0,0,0]), #massarr = 
	                      Float64(0), #time = 
	                      Float64(0), #redshift = 
	                      Int32(0), #flag_sfr = 
	                      Int32(0), #flag_feedback = 
	                      UInt32.([npart,0,0,0,0,0]), #parttotal = 
	                      Int32(0), #flag_cooling = 
	                      Int32(nfiles), #num_files = 
	                      Float64(boxsize), #boxsize = 
	                      Float64(0.3), #omega0 = 
	                      Float64(0.7), #omegalambda = 
	                      Float64(0.7), #hubbleparam =
                              Int32(0), #flag_stellarage = 
                              Int32(0), #flag_metals = 
                              UInt32.([0,0,0,0,0,0]), #npartTotalHighWord =
                              Int32(0), #flag_entropy_instead_u = 
                              Int32(0), #flag_doubleprecision =
                              Int32(0), #flag_ic_info =
                              Float32(0.), #lpt_scalingfactor =
                              Int32.(zeros(12))) #fill =

        println("write header "*sprintf1("%d",i_file-1))
        write_header(fout[i_file], head)
        println("write pos "*sprintf1("%d",i_file-1)*sprintf1(" %d",length(pos[:,filerange[i_file]])/3))
        write_block(fout[i_file], pos[:,filerange[i_file]], "POS")
        
        head = nothing # don't need that anymore
    end
    
# thermal energy
    mass_cgs = mass * m_unit
    T = u2t(T,inv=true) # to be consistent with what IDL did so far.
    Etherm = mass_cgs/umol/mp * 3.0/2.0 *k_boltz * T

    println("Thermal Energy in Box [cgs] = ",Etherm)

# make velocities (the hard part) 
    Etherm /= m_unit *v_unit^2         # cgs to gadget
    v2 = npart*2.0*Etherm/mass*X_turb      # total v^2 of particles
    println("Total amplitude of (squared) velocitities = ",v2)
    
    ngrid = round(Int, npart^(1.0/3.0))

    vgrid = make_vel(ngrid, boxsize, v2)
    
    # Ncells = Npart 
    cellsize = boxsize / ngrid

    vel	 = Array{Float32, 2}(undef,3,npart)
    
    vel[1,:] = idlNGP(pos./cellsize, vgrid[1,:,:,:])
    vel[2,:] = idlNGP(pos./cellsize, vgrid[2,:,:,:])
    vel[3,:] = idlNGP(pos./cellsize, vgrid[3,:,:,:])

# write remaining blocks
    println("write vel")
    for i_file in 1:nfiles
        write_block(fout[i_file], vel[:,filerange[i_file]], "VEL")
    end
    vel = nothing # don't need that anymore
    vgrid = nothing # don't need that anymore
    
# IDs
    id	 = Array{UInt32}(undef,npart)
    id = UInt32.(collect(1:npart))
    println("write id")
    for i_file in 1:nfiles
        write_block(fout[i_file], id[filerange[i_file]], "ID")
    end
    id = nothing # don't need that anymore

# internal energy
    u	 = Array{Float32}(undef,npart)
    println("calculate temperature")
    u[:] .= T
    println("done calculating temperature")

    println("write u")
    for i_file in 1:nfiles
        write_block(fout[i_file], u[filerange[i_file]], "U")
    end
    u = nothing # don't need that anymore
    
# hsml
    hsml = Array{Float32}(undef,npart)
    hsml[:] .= 75.51 # what GADGET finds on average for this glass
    println("write hsml")
    for i_file in 1:nfiles
        write_block(fout[i_file], hsml[filerange[i_file]], "HSML")
    end
    hsml = nothing # don't need that anymore
    
# density
    dens = Array{Float32}(undef,npart)
    dens[:] .= rho
    println("write rho")
    for i_file in 1:nfiles
        write_block(fout[i_file], dens[filerange[i_file]], "RHO")
    end
    rho = nothing # don't need that anymore

    for i_file in 1:nfiles
        close(fout[i_file])
    end

end

# make velocity grid (this is where the magic happens :-)
function make_vel(ngrid, boxsize, amp)
 
    #Random.seed!(14041981)
    r = MT19937() #(14041981)
    
    kmin = 2*pi/(boxsize)       # box mode
    kmax = pi*ngrid/boxsize     # Nyquist mode

    vel = zeros(3,ngrid, ngrid, ngrid)	
    kmag = zeros(ngrid,ngrid,ngrid)
    cdata = zeros(ComplexF64, 3,ngrid,ngrid,ngrid)
    cdata_rl = zeros(ngrid,ngrid,ngrid)
    cdata_im = zeros(ngrid,ngrid,ngrid)
    iconj = 0
    jconj = 0
    kconj = 0

    for axes in 0:2
        for i in 0:ngrid-1
            for j in 0:ngrid-1
                for k in 0:floor(Int, ngrid/2)
                    # Generate k value	first (thank you Volker)
                    # Define conjugated indizes of the grid
		    if i != 0
                        iconj = ngrid - i
		    else
                        iconj = 0
                    end
		    if j != 0
                        jconj = ngrid - j
                    else
                        jconj = 0
                    end
		    if k != 0
                        kconj = ngrid - k
		    else
                        kconj = 0
                    end
                    
                    # Define grid
		    if i <= ngrid/2.
                        kx = i * kmin
		    else
                        kx = -iconj * kmin
                    end
	            
		    if j <= ngrid/2.
                        ky = j * kmin
		    else
                        ky = -jconj * kmin
                    end
	            
		    if k <= ngrid/2.
                        kz = k * kmin
		    else
                        kz = -kconj * kmin
                    end
                    
		    kmag[i+1,j+1,k+1]  = sqrt(kx^2 + ky^2 + kz^2)
	            
                    
		    if kmag[i+1,j+1,k+1] > kmax
		        continue    # Only do a sphere in k space
	            end
                    if i+j+k == 0
                        continue    # no DC current
                    end
                    # if (i eq ngrid/2) or (j eq ngrid/2) or (k eq ngrid/2) then $
                    #     continue    ; no DC current
                    
                    if i > ngrid/2 
                        continue    # these are done via symmetry
                    end
                    
                    # Power spectrum P(k)
                    Pk = kmin .* kolmog_3D(kmag[i+1,j+1,k+1],kmax,kmin)
                    
                    # Generate normal distributed random numbers with dispersion Pk
                    # using Box Mueller method
		    A =  sqrt( -log(rand(r,Float64)) * Pk ) 
		    phase = 2.0*pi*rand(r,Float64)
                    
                    # Cutting off all except the ~70 largest modes (Bauer&Springel2012)
                    if kmag[i+1,j+1,k+1] < (6.25/4)*kmin
                        A = 0
                    elseif kmag[i+1,j+1,k+1] > (12.57/4)*kmin
                        A = 0
                    end
	            
                    # Set power so we get a real vel after inverse FFT
                    if i > 0    # grid is hermitian in i>ngrid/2
                        cdata_rl[i+1,j+1,k+1] = A * cos(phase)
		        cdata_im[i+1,j+1,k+1] = A * sin(phase)
                        
                        cdata_rl[iconj+1,jconj+1,kconj+1] = cdata_rl[i+1,j+1,k+1]
    		        cdata_im[iconj+1,jconj+1,kconj+1] = -1*cdata_im[i+1,j+1,k+1]
                    else  # i = 0 needs special treatment
                        if j == 0 # first row
                            
                            if k > ngrid/2.
                                continue
                            end
                            
                            cdata_rl[i+1,j+1,k+1] = A * cos(phase)
		    	    cdata_im[i+1,j+1,k+1] = A * sin(phase)
                            
                            cdata_rl[i+1,j+1,kconj+1] = cdata_rl[i+1,j+1,k+1]
    			    cdata_im[i+1,j+1,kconj+1] = -1*cdata_im[i+1,j+1,k+1]
                        else  # j != 0 here
                            if j > ngrid/2.     # rest of the plane
                                continue
                            end
                            
                            cdata_rl[i+1,j+1,k+1] = A * cos(phase)
		    	    cdata_im[i+1,j+1,k+1] = A * sin(phase)
                            
                            cdata_rl[i+1,jconj+1,kconj+1] = cdata_rl[i+1,j+1,k+1]
    		    	    cdata_im[i+1,jconj+1,kconj+1] = -1*cdata_im[i+1,j+1,k+1]
                        end
                    end
	        end
            end
        end
        cdata[axes+1,:,:,:] = cdata_rl .+ im .* cdata_im
        data = ifft(cdata[1,:,:,:]) ./ ngrid^3
        
# check if we got the symmetries correct
        for i in 1:ngrid^3
            if abs(imag(data[i])) > 1e-7*abs(real(data[i]))
                println("symmetry failed!")
                return
            end
        end

# set velocities
	vel[axes+1,:,:,:] =  real(data)

    end
    
# norm to total amplitude because kolmog_3D is wrong
    norm = sqrt(sum(vel[1,:,:,:].^2 +vel[2,:,:,:].^2 +vel[3,:,:,:].^2 ))
    
    vel .*=  sqrt(amp) / norm
    cdata .*=  sqrt(amp) / norm  # Parseval's theorem -> ngrid^3

    return vel

end

function kolmog_3D(k, kmax, kmin) 
# here we require 1 = int^kmax_kmin dk P_0 * 4 pi k^2 k^(-11/3)
    norm = (6*pi * ( kmin^(-2.0/3.0) - kmax^(-2.0/3.0) ) )^(-1) 
    Pk = norm * k^(-11.0/3.0)
    return Pk
end

# sample Grid at Pos via NGP the IDL way
function idlNGP(pos, ingrid)

    ngrid = round(Integer,length(ingrid)^(1/3))
    
    grid = reshape(ingrid,(ngrid,ngrid,ngrid) )
    
    npos = round(Integer, length(pos)/3)

    u = reshape(pos[1,:],npos)
    v = reshape(pos[2,:],npos)
    w = reshape(pos[3,:],npos)

    i = floor.(Integer,u .- 1/2) #.+ 1
    j = floor.(Integer,v .- 1/2) .+ 1
    k = floor.(Integer,w .- 1/2) .+ 1

    ngp = Array{Float64, 1}(undef,npos)
    for i_part in 1:npos
        ngp[i_part] = grid[i[i_part],j[i_part],k[i_part]]
    end
    return ngp #reshape(ngp,(npos,3))
end


# convert internal energy to temperature for GADGET
function u2t(input; inv::Bool=false, xH::Float64=0.76, uvel::Float64=1e5, gamma::Float64=5/3) # rad=rad, ??
	
    bk=1.380658e-16   	    #k_boltzmann in cgs
    prtn=1.672623e-24		#m_proton in g 

    yhelium = ( 1. - xH ) / ( 4 * xH ) 

    mean_mol_weight = (1. + 4. * yhelium) / (1. + 3. * yhelium + 1) 

    if inv
	T = input
	u = T /( (gamma-1) * uvel^2 * prtn * mean_mol_weight / bk)
        return u
    else
        u = inpput
	T = u * (gamma-1) * uvel^2 * prtn * mean_mol_weight / bk
        return T
    end
end

# optimal hcp particle distribution in a periodic 
# box, with ntot being a rough upper bound of particles to
# distribute. return pos and ntot
# ntot values like X^3 recommended for cubic boxes
function make_positions_hcp!( pos, lx, ly, lz, ntot)

    # ntot is better divisible by two
    if mod(ntot, 2) != 0
        ntot-=1
    end

    # find by combining spacings with Ntot
    r = (sqrt(2.0)*lx*ly*lz/8.0 /ntot)^(1.0/3.0)

    # spacings
    dx = 2.0*r
    dy = sqrt(3.0)*r
    dz = sqrt(6.0)*2.0/3.0*r
   
    # particle numbers
    np = zeros(Int, 3)
    np[:] .= round(Int, ntot^(1/3))
	
    # enforce periodicity
    dx += (lx-dx*np[1])/np[1]
    dy += (ly-dy*np[2])/np[2]
    dz += (lz-dz*np[3])/np[3]
       
    # diagnostics
    println("Particle numbers :")
    println("     Nx = ",np[1])
    println("     Ny = ",np[2])
    println("     Nz = ",np[3])
    println(" Total  = ",np[1]*np[2]*np[3])
    println(" Wanted = ",ntot)
    println(" Delta  = ",ntot-np[1]*np[2]*np[3])

    ntot = np[1]*np[2]*np[3]

    # particle positions
    x = Array{Float64,3}(undef,np[1],np[2],np[3])
    y = Array{Float64,3}(undef,np[1],np[2],np[3])
    z = Array{Float64,3}(undef,np[1],np[2],np[3])

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

    # A-plane
    # even rows
    for i in 2:2:np[1]-1
        x[:,i+1,1] = x[:,1,1]
        y[:,i+1,1] = y[:,1,1] .+ i*dy
        z[:,i+1,1] = z[:,1,1]
    end
    
    # odd rows
    for i in 3:2:np[1]-1
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
    for i in 2:2:np[1]-1
        x[:,i+1,2] = x[:,1,2]
        y[:,i+1,2] = y[:,1,2] .+ i*dy
        z[:,i+1,2] = z[:,1,2]
    end
   
    # odd rows
    for i in 3:2:np[1]-1
        x[:,i+1,2] = x[:,2,2]
        y[:,i+1,2] = y[:,2,2] .+ (i-1)*dy
        z[:,i+1,2] = z[:,2,2]
    end
    
    # all planes
    # even planes
    for i in 2:2:np[2]-1
        x[:,:,i+1] = x[:,:,1]
        y[:,:,i+1] = y[:,:,1]
        z[:,:,i+1] = z[:,:,1] .+ i*dz
    end
    
    # odd planes
    for i in 3:2:np[2]-1
        x[:,:,i+1] = x[:,:,2]
        y[:,:,i+1] = y[:,:,2]
        z[:,:,i+1] = z[:,:,2] .+ (i-1)*dz
    end
    
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

    # randomize
    #Random.seed!(14041981)
    r = MT19937() #(14041981)
    bin = 0
    for ii in 1:np[1]
        for ij in 1:np[2]
            for ik in 1:np[3]
                i = rand(r,1:ntot) #round(Int, rand()*ntot + 1)
                j = rand(r,1:ntot) #round(Int, rand()*ntot + 1)
                
                pos[1,i], pos[1,j] = pos[1,j], pos[1,i]
                pos[2,i], pos[2,j] = pos[2,j], pos[2,i]
                pos[3,i], pos[3,j] = pos[3,j], pos[3,i]
                bin+=1
            end
        end
    end
end

for res in [64, 128, 256, 512, 1024] 
    if res==1024
        nfiles = 8
    else
        nfiles = 1
    end
    make_box(res^3, nfiles=nfiles, odir="./sph_"*sprintf1("%d",res)*"/")
end
