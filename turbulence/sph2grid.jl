using GadgetIO

#using Plots
using LaTeXStrings

ENV["GKSwsystype"]="nul"

using Distributed
addprocs(4)

using SPHKernels
using SPHtoGrid

using PyPlot

using Formatting

hydro_scheme="sph"
res=128
snap_num = 16
snap_dir = "../../test_runs/out_turb_"*hydro_scheme*"_"*sprintf1("%d",res)*"/"
snap = snap_dir*"snap_"*sprintf1("%03d",snap_num)

print(snap)

data = Dict( "POS"  => read_block( snap, "POS",  parttype=0 ),
       	     "RHO"  => read_block( snap, "RHO",  parttype=0 ),
	     "HSML" => read_block( snap, "HSML", parttype=0 ),
	     "MASS" => read_block( snap, "MASS", parttype=0 )
     	   )
size=6e5
par = mappingParameters(center=Float32[size,size,size],
			x_size=size,
			y_size=size,
			z_size=size,
			Npixels=4000)

# define the kernel
k = WendlandC6()

# density in pixel units
rho_2d = Float32.(density_2D.(data["RHO"], par.pixelSideLength))

# define weights
weights = ones(Float32, length(data["RHO"]))


# mapping loop
density_map = sphMapping( data["POS"], data["HSML"], data["MASS"], data["RHO"],
	      		  rho_2d, weights, show_progress=true,
			  param=par, kernel=k,
			  parallel=true, # run on more than 1 core
			  reduce_image=false, filter_particles=true
			)

# write to file
#cd(snap_dir)

image_file = hydro_scheme*"_"*sprintf1("%d",res)*"/rho.fits"
write_fits_image(image_file, density_map, par, snap=snap_num, units="g/cm^2")

imshow(density_map)
savefig(hydro_scheme*"_"*sprintf1("%d",res)*"/rho.png")
