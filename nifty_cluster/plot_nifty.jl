using GadgetIO

#using Plots
using LaTeXStrings

ENV["GKSwsystype"]="nul"

using Formatting
using PyPlot
using Statistics
using PyCall

using FITSIO

using PyPlotHelper # private helper to format plots


"""
    plot_surface_density_fits(; fitsname=nothing, ax=nothing, vmin=3e-4, vmax=2, name="", scaling=3)

Plot surface density from fits file
"""
function plot_surface_density_fits(; fitsname=nothing, ax=nothing, vmin=3e-4, vmax=2, name="", scaling=3)
    if !isfile(fitsname)
        println("file "*fitsname*" does not exist, skip this plot")
        ax.remove()
        return
    end

    density_map=read(FITS(fitsname)[2])

    conversion_factor = (1e15 * 1.989e33 / 3.085678e24^2)
    density_map = density_map .* conversion_factor
    
    im = ax.imshow(log10.(transpose(density_map)),cmap=get_colormap("density"),vmin=log10(vmin),vmax=log10(vmax))

    phi = 0:0.01:2*pi
    rvir = size(density_map)[1] / 2 / scaling
    x = @. cos(phi)*rvir + rvir*scaling
    y = @. sin(phi)*rvir + rvir*scaling
    ax.plot(x,y, color="red", linestyle="dashed")

    ax.set_xlim([0,2*rvir*scaling])
    ax.set_ylim([0,2*rvir*scaling])

    ax.set_xticks([rvir*(scaling-1),rvir*scaling,rvir*(scaling+1)])
    ax.set_yticks([rvir*(scaling-1),rvir*scaling,rvir*(scaling+1)])
    
    ax.set_xticklabels([L"-R_{vir}","0",L"R_{vir}"])
    ax.set_yticklabels([L"-R_{vir}","0",L"R_{vir}"])

    ax.text(10, 40, name, fontsize=PyPlotHelper.title_font_size, color="yellow")

    return im
    
end

"""
    find_half_mass_radius(pos,center,mass;n_points::Int64=1000,return_all::Bool=false,r_min::Real=-1,r_max::Real=-1)

Calcualte half-mass radius
"""
function find_half_mass_radius(pos,center,mass;n_points::Int64=1000,return_all::Bool=false,r_min::Real=-1,r_max::Real=-1)

    pos=deepcopy(pos)
    total_mass = sum(mass)

    pos1_mean = center[1]
    pos2_mean = center[2]
    pos3_mean = center[3]

    
    for i in 1:size(pos[1,:])[1]
        pos[1,i] = sqrt( (pos[1,i]-pos1_mean)^2 + (pos[2,i]-pos2_mean)^2 + (pos[3,i]-pos3_mean)^2)
    end

    if r_min == -1
        r_min = minimum(pos[1,:][pos[1,:] .!= 0])
    end
    if r_max == -1
        r_max = maximum(pos[1,:])
    end
    println(r_min," ",r_max)

    data = Dict( "r"  => Array{output_type}(undef,n_points),
                 "m"  => Array{output_type}(undef,n_points)
               )
    
    for i in 1:n_points
        r = exp( log(r_min) + (log(r_max)-log(r_min)) * i/n_points )
        function lessr(x)
            if x<r
                return true
            else
                return false
            end
        end
        enclosed = map(lessr,pos[1,:])
        enclosed_mass = sum(mass[enclosed])
        if !return_all
            if enclosed_mass >= total_mass/2e0
                r12_lower = r
                r12_higher = exp( log(r_min) + (log(r_max)-log(r_min)) * (i+1)/n_points )
                return (r12_lower + r12_higher)/2.0
            end
        end
        data["r"][i] = r
        data["m"][i] = enclosed_mass
    end
    return data
    
end

"""
    plot_radial_profile(snap::String, group::String; density::Bool=false,rho_ax=axes(1),temperature::Bool=false,T_ax=axes(1),plot_entropy::Bool=false,ent_ax=axes(1),nifty_comparison::Bool=false,color="red", use_sub_center::Bool=true,in_corner::Bool=false,hubble0 = 0.7,center=nothing, plot_nifty_in_bg::Bool=false,Tmax=15,rhomax=2e16, Smax=4e-7, r_min=-1, r_max=-1)

Create plots of radial density/temperature/entropy profile.
"""
function plot_radial_profile(snap::String, group::String; density::Bool=false,rho_ax=axes(1),temperature::Bool=false,T_ax=axes(1),plot_entropy::Bool=false,ent_ax=axes(1),nifty_comparison::Bool=false,color="red", use_sub_center::Bool=true,in_corner::Bool=false,hubble0 = 0.7,center=nothing, plot_nifty_in_bg::Bool=false,Tmax=15,rhomax=2e16, Smax=4e-7, r_min=-1, r_max=-1)

    if !isfile(snap) & !isfile(snap*".0")
        println("file "*snap*" does not exist, skip this plot")
        return
    end

    if nifty_comparison
        n_points = 50
    else
        n_points = 1000
    end
    
    pos = read_block( snap, "POS",  parttype=0, info=InfoLine("POS", output_type, 3, [1, 1, 1, 1, 0, 1]) ) #pos_info )
    if in_corner
        pos = @. mod(pos+5e5,1e6)
    end
    
    pos = pos/1e3 # kpc -> Mpc
    if use_sub_center
        pos_sub = read_subfind(group, "GPOS")
        if in_corner
            pos_sub = @. mod(pos_sub+5e5,1e6)
        end
        println("subfind central position: ",pos_sub[:,1]/1e3)
        center = pos_sub[:,1]/1e3
    elseif center==nothing
        center=[498.42328125, 496.4240625, 487.32759375] # good estimate at early times
        center = [500,500,500]
    end
    println("center (radial profile)")
    println(center)

    #mass = read_block( snap, "MASS",  parttype=0, info=mass_info)
    mass = ones(length(mass)) .* GadgetIO.read_header(snap).massarr[1]
    
    mass = mass*1e10 # 1e10 M_sun -> M_sun 

    if use_sub_center
        R_200 = read_subfind(group, "RCRI")[1] / 1e3
        R_vir = try
            read_subfind(group,"RVIR")[1] / 1e3
        catch
            read_subfind(group,"RTOP")[1] / 1e3
        end
        println("R200=",R_200," , R_vir=",R_vir)
    end
    
    if density
        println("plot density profile")
        
        data = find_half_mass_radius(pos,center,mass,return_all=true,r_min=r_min,r_max=r_max,n_points=n_points)
        rho = @. (data["m"][2:n_points]-data["m"][1:n_points-1])/(4*pi*data["r"][2:n_points]^2*(data["r"][2:n_points]-data["r"][1:n_points-1]))
        println(minimum(rho[@. !isnan(rho)]), " ",maximum(rho[@. !isnan(rho)]))
        rho_line, = rho_ax.loglog(data["r"][2:n_points],rho,"+",color=color)
        rho_ax.set_xlim([r_min,r_max])
        if rhomax != nothing
            ylim=[rhomax*1e-4,rhomax] #2e12,2e16
            rho_ax.set_ylim(ylim)
        else
            ylim = [minimum(rho), maximum(rho)]
        end
        
        if use_sub_center
            rho_ax.plot([R_200,R_200],ylim,linestyle="dashdot",color=color)
        end
        
        rho_ax.set_xlabel(L"R"*" ["*L"h^{-1}"*"Mpc]")
        rho_ax.set_ylabel(L"\rho_{gas}"*" ["*L"h^{-1}M_{\odot}/(h^{-3}"*"Mpc"*L"{}^3)]")
    end
    if temperature
        println("plot temperature profile")
        U = read_block( snap, "U",  parttype=0, info=u_info )
        data_U = find_half_mass_radius(pos,center,U.*mass,return_all=true,r_min=r_min,r_max=r_max,n_points=n_points) # mass weighted
        T = @. u2t((data_U["m"][2:n_points]-data_U["m"][1:n_points-1])/(data["m"][2:n_points]-data["m"][1:n_points-1])) # mass weighted
        println(minimum(T[@. !isnan(T)])*1.380658e-16*624150636.3094, " ", maximum(T[@. !isnan(T)])*1.380658e-16*624150636.3094)
        T_line, = T_ax.semilogx(data_U["r"][2:n_points],T.*1.380658e-16 .* 624150636.3094,"+",color=color) # kB*T in keV
        T_ax.set_xlabel(L"R"*" ["*L"h^{-1}"*"Mpc]")
        T_ax.set_ylabel(L"\rm{k}_{\rm{B}}T"*" [keV]")
        T_ax.set_xlim([r_min,r_max])
        if Tmax != nothing
            ylim=[0,Tmax]
            T_ax.set_ylim(ylim)
        else
            ylim = [minimum(T), maximum(T)].*1.380658e-16 .* 624150636.3094
        end
        
        if use_sub_center
            T_ax.plot([R_200,R_200],ylim,linestyle="dashdot",color=color)
        end
    end
    if plot_entropy
        println("plot entropy profile")
        
        rho = read_block( snap, "RHO",  parttype=0, info=rho_info )# ./ 3.0857e24^2 # in cgs
        mass = read_block( snap, "MASS", parttype=0, info=mass_info )
        #mass = ones(length(mass)) .* GadgetIO.read_header(snap).massarr[1]
        xH=0.76 # mass fraction
        proton_mass = 1.67262192369e-24 / 1.989e43 # g -> internal units 
        yHe = (1-xH)
        ne = rho ./ proton_mass * (xH + 2*yHe/4) # in kpc^-3 ( internal units) # assume full ionization
        println(minimum(ne), maximum(ne))

        U = read_block( snap, "U",  parttype=0, info=u_info )
        T = @. u2t(U)
        S = @. T/ne^(2/3)
        println(minimum(S),maximum(S))

        S = S .* 1e10 # increase numerical accuracy
        
        data_S = find_half_mass_radius(pos,center,S.*mass,return_all=true,r_min=r_min,r_max=r_max,n_points=n_points) # mass weighted
        S = @. (data_S["m"][2:n_points]-data_S["m"][1:n_points-1])/(data["m"][2:n_points]-data["m"][1:n_points-1])  # mass weighted
        S = S ./ 1e10 # convert back
        
        println(minimum(S[@. !isnan(S)])," ",maximum(S[@. !isnan(S)]))
        # convert kpc^2 ( internal units ) -> cm^2
        S = S .* 3.085678e21^2
        # T -> kB*T in keV
        S = S .*1.380658e-16 .* 624150636.3094
        S = S .* 1e10
        
        ent_line, = ent_ax.loglog(data_S["r"][2:n_points],S,"+",color=color)
        ent_ax.set_xlabel(L"R"*" ["*L"h^{-1}"*"Mpc]")
        ent_ax.set_ylabel(L"S=T/n_e^{2/3}"*" [keV "*L"h^{-2}"*"cm"*L"^2"*"]")
        ent_ax.set_xlim([r_min,r_max])

        if Smax != nothing
            ylim=[Smax*2e-3,Smax]
            ent_ax.set_ylim(ylim)
        else
            ylim = [minimum(S),maximum(S)]
        end
        
        if use_sub_center
            ent_ax.plot([R_200,R_200],ylim,linestyle="dashdot",color=color)
        end

    end
    lines = []
    if density
        push!(lines,rho_line)
    end
    if temperature
        push!(lines,T_line)
    end
    if plot_entropy
        push!(lines,ent_line)
    end
    return rho_line
end

"""
    u2t(U)

Convert internal energy to temperature
"""
function u2t(U)
    gamma=5/3
    xH=0.76
    uvel=1e5
    kB=1.380658e-16 # in cgs
    m_proton=1.672623e-24 # in g
    yHe=(1-xH)/(4*xH)
    mean_mol_weight=(1+4*yHe)/(1+3*yHe+1)

    T = U * (gamma-1)*uvel^2*m_proton*mean_mol_weight / kB
    return T
end


"""
    analyze_cluster(snap_num, hydros, sub_dirs;
                     outdir="./", cluster_name="nifty",
                     plot_density_fits=true,plot_profile::Bool=true,use_sub_center::Bool=false, plot_nifty_in_bg::Bool=false,
                     vmin=3e-4, vmax=2, fits=nothing,
                     Tmax=15, rhomax=2e16,Smax=4e-7, r_min=-1, r_max=-1,
                     hydro_names=nothing)

Create analysis plots for galaxy cluster, including radial profiles andsurface density.
"""
function analyze_cluster(snap_num, hydros, sub_dirs;
                         outdir="./", cluster_name="nifty",
                         plot_density_fits=true,plot_profile::Bool=true,use_sub_center::Bool=false, plot_nifty_in_bg::Bool=false,
                         vmin=3e-4, vmax=2, fits=nothing, # argumentes for plot_surface_density_fits
                         Tmax=15, rhomax=2e16,Smax=4e-7, r_min=-1, r_max=-1,
                         hydro_names=nothing)
    
    main_dirs = "../../test_runs/out_"*cluster_name*"_".*hydros.*sub_dirs #_switch5/"
    snap_dirs = @. main_dirs*"snapdir_"*sprintf1("%03d",snap_num)*"/"
    snaps = @. snap_dirs*"/snap_"*sprintf1("%03d",snap_num)
    group_dirs = @. main_dirs*"groups_"*sprintf1("%03d",snap_num)*"/"
    groups = @. group_dirs*"/sub_"*sprintf1("%03d",snap_num)

    if hydro_names==nothing
        hydro_names = hydros
    end

    if cluster_name == "nifty"
        in_corner = false
        hubble0 = 0.7
    else
        println("requested cluster name not supported yet. Assume cluster is in center of domain.")
        in_corner = false
        hubble0 = 0.72
    end
    if plot_profile
        pro_fig = figure(figsize=(12,4))
        style_plot(fig_width=12, print_columns=1)
        pro_gs = pro_fig.add_gridspec(1,3, left=0.08, right=0.98, bottom=0.15, top=0.98, hspace=0.05, wspace=0.25)
        pro_ax = pro_gs.subplots()
        rho_ax = pro_ax[1]
        T_ax = pro_ax[2]
        ent_ax = pro_ax[3]
    end

    if plot_density_fits
        n_width = minimum([2, length(main_dirs)])
        n_height = ceil(Int, length(main_dirs) / n_width)
        dens_fig = figure(figsize=(4*n_width,4*n_height))
        print_columns_density=1
        style_plot(fig_width=4*n_width, print_columns=print_columns_density)
        ratios = ones(n_width)
        append!(ratios,[0.05])
        if print_columns_density==2
            dens_top = 0.85
            dens_bottom = 0.1
            dens_right = 0.85
            dens_left = 0.1
            dens_space = 0.01
        elseif print_columns_density==1
            dens_top = 0.88
            dens_bottom = 0.07
            dens_right = 0.88
            dens_left = 0.07
            dens_space = 0.01
        end
        dens_gs = dens_fig.add_gridspec(n_height,n_width, hspace=dens_space, wspace=dens_space,
                                        left=dens_left,bottom=dens_bottom,top=dens_top,right=dens_right)
        dens_ax = dens_gs.subplots()1
        
    end

    #default values
    output_type=Float32
    pos_info = nothing
    rho_info = nothing
    hsml_info = nothing
    mass_info = nothing
    u_info = nothing
    ne_info = nothing

    all_lines = []
    for i in 1:length(hydros)

        if !use_sub_center # workaround if sub center not found corectly
            use_sub_center_here=false
            center = [499.49475, 493.14453125, 492.3893125] .* 1e3
            center_Mpc = center ./ 1e3
        else
            use_sub_center_here=true
            center=nothing
            center_Mpc=nothing
        end

        if length(snap_num) > 1
            i_snap = snap_num[i]
        else
            i_snap = snap_num
        end
            
        if plot_density_fits
            im = plot_surface_density_fits(fitsname=fits[i],ax=dens_ax[i],vmin=vmin,vmax=vmax,name=hydro_names[i])
            if i%n_width != 0
                dens_ax[i].set_xticklabels([])
            end
            if i > n_height
               dens_ax[i].set_yticklabels([])
            end 
            if i == 1
                cax = dens_fig.add_axes([dens_right+dens_space,dens_bottom,0.015,dens_top-dens_bottom])
                cbar = colorbar(im, cax=cax)
                cbar.set_label(L"\log~\rho"*" ["*L"10^{15}"*"h"*L"{}^{-1}M_{\odot}\cdot"*"Mpc"*L"{}^{-2}]")
                for j in length(main_dirs)+1:n_width*n_height
                    dens_ax[j].remove()
                end
            end
        end
        if plot_profile
            lines = plot_radial_profile(snaps[i],groups[i],density=true,rho_ax=rho_ax,temperature=true,T_ax=T_ax,plot_entropy=true,ent_ax=ent_ax,nifty_comparison=true,color=get_color(i),use_sub_center=use_sub_center_here,in_corner=in_corner,hubble0=hubble0,center=center_Mpc, plot_nifty_in_bg = plot_nifty_in_bg, Tmax=Tmax,rhomax=rhomax,Smax=Smax,r_min=r_min,r_max=r_max)
            if isfile(snaps[i]) || isfile(snaps[i]*".0")
                push!(all_lines,lines)
            end
        end
    end
    if plot_profile
        if plot_nifty_in_bg # data from Hopkins (2015)
            # density comparison
            arepo_data = [0.0244872401747  3.79725881255e+14
                          0.0538501175094  3.2726582962e+14
                          0.101891644342  2.75674979527e+14
                          0.164915429913  1.84738022891e+14
                          0.341636292891  9.96195030337e+13
                          0.646885322755  3.68319884332e+13
                          0.934350042247  1.49224818591e+13
                          1.55973922009  3.8703331019e+12]
            arepo_line, = rho_ax.plot(arepo_data[:,1], arepo_data[:,2], color="black", linestyle="dashed")
            g3music_data = [0.024737723069  3.18648044985e+15
                            0.0512539456361  1.44742485662e+15
                            0.0849908382119  7.98578991361e+14
                            0.151505184699  3.7113885523e+14
                            0.242366712858  1.68585877984e+14
                            0.387681612169  8.58568270324e+13
                            0.627720898703  3.68319884332e+13
                            1.07987628903  9.44420842329e+12
                            1.55998775196  3.22313239885e+12]
            g3music_line, = rho_ax.plot(g3music_data[:,1], g3music_data[:,2], color="red", linestyle="solid")
            # temperature comparison
            arepo_data = [0.0248362208632  11.503992951
                          0.0352055208846  11.0020941332
                          0.0434589935432  10.643330641
                          0.0634711056773  10.261704971
                          0.0905209756636  9.68780674058
                          0.118657331553  9.32945443865
                          0.165128661524  9.13960790073
                          0.267266206354  8.23038402232
                          0.393076865501  7.08036713415
                          0.577896279108  6.21851237242
                          0.708602236128  6.41201850356
                          0.905244414283  7.39825831559
                          1.23900947498  6.12768044644
                          1.61789980213  4.08834128791]
            T_ax.plot(arepo_data[:,1], arepo_data[:,2], color="black", linestyle="dashed")
            g3music_data = [0.0211935786243  4.10668037301
                            0.0365710806663  4.75878698877
                            0.0478895037926  5.19288053455
                            0.0607979140246  6.37118731184
                            0.0714452041834  6.9486217784
                            0.0876070628947  7.118114399
                            0.124118428548  7.0244452603
                            0.177937792339  7.07493942287
                            0.241652369205  7.19710404582
                            0.318553978606  7.12695498935
                            0.412559180462  6.81654747045]
            T_ax.plot(g3music_data[:,1], g3music_data[:,2], color="red", linestyle="solid")
            # entropy comparison
            arepo_data = [0.0236817477899  196.804201738
                          0.0284553795683  208.405217604
                          0.0371839193622  200.143989766
                          0.0487423436999  196.168216555
                          0.0580207058656  198.621808762
                          0.0646970963887  201.912101429
                          0.0772515397192  201.953728345
                          0.0919602812118  212.119777439
                          0.108788178717  218.301517234
                          0.19776127292  282.406903131
                          0.281117578609  317.969898996
                          0.435976365496  352.257495893
                          0.477164519102  366.950008236
                          0.593404397931  451.848147687
                          0.751991089598  652.261676811
                          0.983255481353  1121.97116931
                          1.17437887497  1462.60965861
                          1.33417854258  1474.79974968
                          1.4693006005  1499.2090285
                          1.59811415573  1542.74935333]
            ent_ax.plot(arepo_data[:,1], arepo_data[:,2], color="black", linestyle="dashed")
            g3music_data = [0.0252188883581  18.5854506233
                            0.0399892826919  30.9506249019
                            0.0570371032818  45.0504034391
                            0.0698393643679  57.5449788287
                            0.0898834950234  76.2557002781
                            0.129004152739  112.819757383
                            0.150274659163  131.742633259
                            0.200781152053  191.744369088
                            0.251265000329  248.956120324
                            0.313437660879  296.721553942
                            0.417400081183  358.029406367
                            0.528884875553  459.213700362
                            0.6849661314  644.26480974
                            0.833606040953  922.429378865
                            1.10697025384  1468.48226495
                            1.38096511256  1860.57275783
                            1.60838269216  1838.28687414]
            ent_ax.plot(g3music_data[:,1], g3music_data[:,2], color="red", linestyle="solid")
            my_legend = rho_ax.legend(all_lines,hydro_names, loc="lower left",bbox_to_anchor=[-0.02,0.01], borderpad=0, borderaxespad=0)
        end
        rho_ax.legend([arepo_line,g3music_line], ["AREPO","G3-MUSIC"], loc="upper right", title="Sembolini et al. (2016)")
        if plot_nifty_in_bg
            rho_ax.add_artist(my_legend)
        end
        
        pro_fig.show()
        pro_fig.savefig(outdir*"profile_comparison_"*sprintf1("%03d",maximum(snap_num))*".pdf")
        if close_all_figures
            close(pro_fig)
        end
    end
    if plot_density_fits
        dens_fig.show()
        dens_fig.savefig(outdir*"density_comparison_"*sprintf1("%03d",maximum(snap_num))*".png")
        if close_all_figures
            close(dens_fig)
        end
    end
end

cluster_name = "nifty"

plot_nifty_in_bg = false

test_dir = "../"

hydro_names = nothing
if cluster_name == "nifty" # nifty cluster
    i_snap = 81 # z=0
    hydros= ["mfm","sph","sph", "sph"]
    sub_dirs = ["/early/","/early/", "/older/early_msph/", "/older/early/"]
    outdir = test_dir*"/nifty_cluster/"
    fits = @. "../../test_runs/out_nifty_"*hydros*"/"*sub_dirs*"/fits/mymap_RHO."*sprintf1("%03d", i_snap)*".0.z.fits"
    plot_nifty_in_bg = true
    vmax = 3e-3
    vmin = 1e-4
    r_min=2e-2
    r_max=2e0
    Tmax = 15
    rhomax = 2e16
    Smax = 5e3
    hydro_names = ["MFM", "SPH"*L"(\alpha_{\rm{max}}^{\rm{cond}}=0.25)", "SPH"*L"(\kappa_{\rm{phys}})", "SPH"*L"(\alpha_{\rm{max}}^{\rm{cond}}=0)"]
end

final = true

if final
    use_sub_center=true
    plot_profile=true
else
    use_sub_center=false
    plot_profile=false
end

global output_type=Float32
global pos_info = nothing #InfoLine("POS", output_type, 3, [1, 1, 1, 1, 0, 1])
global rho_info = nothing #InfoLine("RHO", output_type, 1, [1, 0, 0, 0, 0, 0])
global hsml_info = nothing #InfoLine("HSML", output_type, 1, [1, 0, 0, 0, 0, 0])
global mass_info = nothing #InfoLine("MASS", output_type, 1, [1, 1, 1, 1, 1, 1])
global u_info = nothing #InfoLine("U", output_type, 1, [1, 0, 0, 0, 0, 0])
global ne_info = nothing #InfoLine("NE", output_type, 1, [1, 0, 0, 0, 0, 0])
        

global close_all_figures = false # if true, don't show figures longer than necesary, avoids to open a large number of figures at once
for snap_num in [i_snap]
    analyze_cluster(snap_num, hydros, sub_dirs, outdir=outdir, cluster_name=cluster_name,plot_density_fits=true,plot_profile=plot_profile,use_sub_center=use_sub_center, plot_nifty_in_bg=plot_nifty_in_bg, vmin=vmin, vmax=vmax, fits=fits,
                    Tmax=Tmax, rhomax=rhomax, Smax=Smax, r_min=r_min,r_max=r_max,
                    hydro_names=hydro_names)
end
