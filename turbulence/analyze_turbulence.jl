# plot velocity powerspectrum on a grid made by Sph2Grid
# and compare to IDL FFT. Note that differences at small k are due
# to the different data layout. IDL includes the redundant Hermitian
# part of the FFT data, which gives differences in the binning.

include("readgrid.jl")
using GadgetIO
using Formatting
using Statistics
using PyPlot
using LaTeXStrings
using LsqFit
using StatsBase
using LinearAlgebra
using CurveFit
using PyCall

using PyPlotHelper # private helper to format plots

global close_all_figures = false

global print_columns = 1.4 #2

global i_mfm = 1
global i_sph = 2
global i_arepo = 4
global i_64  = 6
global i_128 = 3
global i_256 = 1
global i_512 = 7

global t_sc = 10.59 # sound crossing time for turbulent box


"""
    method_comparison(test_dir::String,igrid1::Int64,igrid2::Int64=-1;normalized::Bool=false,split_legend::Bool=false, oname=nothing)

Create comparison plot of turbulent power spectrum for different methods.
"""
function method_comparison(test_dir::String,igrid1::Int64,igrid2::Int64=-1;normalized::Bool=false,split_legend::Bool=false, oname=nothing, plot_all_methods=true, plot_all_resolutions=true)

    if oname == nothing
        oname = "power_spectrum_method_comparison.pdf"
    end
    
    fig=figure(figsize=(12,8))
    style_plot(fig_width=12, print_columns=1)
    ax=fig.add_subplot()

    function method_color(i_method,i_res, plot_all_methods,plot_all_resolutions)
        if (plot_all_methods && plot_all_resolutions)
            return get_color(i_method)
        elseif plot_all_methods
            return get_color(i_method)
        elseif plot_all_resolutions
            if i_res == i_64
                return get_color(1)
            elseif i_res == i_128
                return get_color(2)
            elseif i_res == i_256
                return get_color(3)
            elseif i_res == i_512
                return get_color(5)
            end
            return get_color(i_res)
        end
    end
    function res_marker(i_method,i_res, plot_all_methods,plot_all_resolutions)
        if (plot_all_methods && plot_all_resolutions)
            return get_marker(i_res)
        elseif plot_all_methods
            return get_marker([6,3,0,1,7][i_method])
        elseif plot_all_resolutions
            return get_marker(i_res) # so far: v, *(default), P, p, 
        end
    end
    
    if !(plot_all_methods || plot_all_resolutions)
        println("Error: at least one of plot_all_methods and plot_all_resolutions has to be true!")
        return
    end
    
    if plot_all_methods
        # sph 128
        kol_line, sph_line = plot_powerspectrum(test_dir*"/out_turb_sph_128/",igrid1,label="sph_128",color=method_color(i_sph,i_128,plot_all_methods,plot_all_resolutions),marker=res_marker(i_sph,i_128,plot_all_methods,plot_all_resolutions),ax=ax,normalized=normalized)
        if plot_all_resolutions
            # sph 256
            sph256_line = plot_powerspectrum(test_dir*"/out_turb_sph_256/",igrid1,label="sph_256",color=method_color(i_sph,i_256,plot_all_methods,plot_all_resolutions),overplot=true,marker=res_marker(i_sph,i_256,plot_all_methods,plot_all_resolutions),ax=ax,normalized=normalized, scale=1e-2)
        end
    end
    if plot_all_resolutions
        # mfm 64 cubic
        mfm_64_line = plot_powerspectrum(test_dir*"/out_turb_mfm_64/",igrid1,label="mfm_64_cubic",color=method_color(i_mfm,i_64,plot_all_methods,plot_all_resolutions),overplot=plot_all_methods,marker=res_marker(i_mfm,i_64,plot_all_methods,plot_all_resolutions),ax=ax,normalized=normalized)
        if !plot_all_methods
            mfm_64_line = mfm_64_line[2] # ignore Kolmogorov line
        end
    end
    # mfm 128 cubic
    mfm_cubic_line = plot_powerspectrum(test_dir*"/out_turb_mfm_128_cubic/",igrid1,label="mfm_128_cubic",color=method_color(i_mfm,i_128,plot_all_methods,plot_all_resolutions),marker=res_marker(i_mfm,i_128,plot_all_methods,plot_all_resolutions),overplot=true,ax=ax,normalized=normalized)
    if plot_all_resolutions
        # mfm 256 cubic
        mfm_cubic256_line = plot_powerspectrum(test_dir*"/out_turb_mfm_256_cubic/",igrid1,label="mfm_256_cubic",color=method_color(i_mfm,i_256,plot_all_methods,plot_all_resolutions),overplot=true,marker=res_marker(i_mfm,i_256,plot_all_methods,plot_all_resolutions),ax=ax,normalized=normalized)
        # mfm 512 cubic
        mfm_cubic512_line = plot_powerspectrum(test_dir*"/out_turb_mfm_512/",igrid1,label="mfm_512_cubic",color=method_color(i_mfm,i_512,plot_all_methods,plot_all_resolutions),overplot=true,marker=res_marker(i_mfm,i_512,plot_all_methods,plot_all_resolutions),ax=ax,normalized=normalized)
    end
    if plot_all_methods
        # arepo_static_128
        arepo_static_line = plot_powerspectrum(test_dir*"/out_turb_arepo_static_128/",igrid1,label="arepo_static_128",color=method_color(5,i_128,plot_all_methods,plot_all_resolutions),marker=res_marker(5,i_128,plot_all_methods,plot_all_resolutions),overplot=true,ax=ax,normalized=normalized)
        # arepo_static_256
        #plot_powerspectrum(test_dir*"/out_turb_arepo_static_256/",igrid1,label="arepo_static_256",color="purple",overplot=true,marker="P",ax=ax,normalized=normalized) # ,scale=1
        # arepo_128
        arepo_line = plot_powerspectrum(test_dir*"/out_turb_arepo_128/",igrid1,label="arepo_128",color=method_color(i_arepo,i_128,plot_all_methods,plot_all_resolutions),marker=res_marker(i_arepo,i_128,plot_all_methods,plot_all_resolutions),overplot=true,ax=ax,normalized=normalized)
        gizmo_line = plot_powerspectrum(test_dir*"/out_turb_mfm_gizmo_128"*suffix*"/",igrid1,label="gizmo_128",color=method_color(i_gizmo,i_128,plot_all_methods,plot_all_resolutions),marker=res_marker(i_gizmo,i_128,plot_all_methods,plot_all_resolutions),overplot=true,ax=ax,normalized=normalized, scale=scale)
    end
    
    tight_layout()

    if plot_all_methods
        method_lines = [kol_line,mfm_cubic_line,gizmo_line,sph_line,arepo_line,arepo_static_line]
        method_descriptions = ["Kolmogorov", "MFM (OpenGadget3)", "MFM (gizmo)", "SPH (OpenGadget3)","Moving Mesh (AREPO)","Static Mesh (AREPO)"]
        method_legend = legend(method_lines,method_descriptions, bbox_to_anchor=[0.3,0.02], borderpad=0, borderaxespad=0, loc="lower left")
    end
    if plot_all_resolutions
        res_lines = [mfm_64_line,mfm_cubic_line, mfm_cubic256_line, mfm_cubic512_line]
        res_descriptions = [L"64^3",L"128^3",L"256^3", L"512^3"]
        res_legend = legend(res_lines,res_descriptions,title="resolution", bbox_to_anchor=[0.58,0.02], borderpad=0, borderaxespad=0, loc="lower left")
    end
    if (plot_all_methods && plot_all_resolutions)
        ax.add_artist(method_legend)
    end

    show()
    fig.savefig(oname)
    if close_all_figures
        close(fig)
    end
end
    

"""
    plot_powerspectrum(fdir::String="./", fnum::Int64=0; nodetails::Bool=false, overplot::Bool=false, color::String="red", label::String="powerspectrum", no_hsml::Bool=false, marker::String="*", scale::Number=1, ax=ax, normalized::Bool=false, pre_ylabel="",yrange=nothing)

Plot the turbulent power spectrum.
"""
function plot_powerspectrum(fdir::String="./", fnum::Int64=0; nodetails::Bool=false, overplot::Bool=false,color::String="red",label::String="powerspectrum",no_hsml::Bool=false,marker::String="*",scale::Number=1,ax=ax,normalized::Bool=false,pre_ylabel="",yrange=nothing)

    fname = fdir*"grid_"*sprintf1("%03d",fnum)

    head = readgrid(fname, "HEAD")

    gridsize = head.GridSize
    ngrid = head.NGrid
    cellsize = gridsize / ngrid

    kmin = 2*pi/(gridsize)      # box mode
    kmax = pi*ngrid/gridsize    # Nyquist mode

    #read file before plotting
    Pk = readgrid(fname, "VEL", Pk=true)
    #println(Pk)                    #for NAN control reasons
    bin_pos = readgrid(fname, "KPK")
    
    time=head.Time/t_sc    #divided by sound crossing time
        
    # plot
    #setcolors
    fname_short=sprintf1("%3d",fnum)
    if !overplot
        loglog([],[])
        xrange=[(0.9*kmin),(2*kmax)]
        xlim(xrange) #(1.2*kmax)])
        xlabel(L"k = 2π/L"*" [kpc"*L"{}^{-1}]")
        if normalized
            ylabel(pre_ylabel*L"k^{5/3}E(k)"*" [arb. units]")
        else
            ylabel(pre_ylabel*L"P(k)"*" [arb. units]")
        end
        if normalized && yrange==nothing 
            yrange=[1e-5,1e2]
        elseif !normalized && yrange==nothing
            yrange=[1e-10,1e1]
        end
        ylim(yrange)
    end

    #plot Kolmogorov_3D spectrum (-11/3)
    # Kolmogorov has a slope of -5/3, but here -11/3 due to 3d and integration over shells (additional k^2)
    spectrum=deepcopy(bin_pos)
    for i in 1:size(spectrum)[1]
        tmp = bin_pos[i]^(-11.0/3.0)/2.7e9
        spectrum[i] = tmp
    end
    if normalized
        normalization = spectrum
    else
        normalization = ones(Float64,size(spectrum))
    end
    if !overplot && normalized
        kol_line, = plot([minimum(bin_pos),xrange[2]], [1,1],"k-.",label="Kolmogorov")
    end
        
    #remove NANs and 0s through interpolation and plot corrected "Pk_good"
    Pk_good = deepcopy(Pk)
    for i in 2:size(Pk)[1]-1 # first and last element cannot be interpolated
        if isnan(Pk_good[i]) || Pk_good[i] == 0
            if Pk[i-1] != 0 && Pk[i+1] != 0
                Pk_good[i]=0.5*(Pk[i-1] + Pk[i+1])
            end
        end
    end
    #println(fname, Pk_good)
    line, = plot(bin_pos, scale.*Pk_good./normalization,marker=marker,color=color,label=label,markeredgewidth=0,linestyle=" ")

    if !overplot

        if normalized
            ypos_text = 3e1
        else
            ypos_text = 0.7e-3
        end
        
        if !no_hsml
            # plot average smoothing length
            sname=fdir*"/snap_"*sprintf1("%03d",fnum)
            
            hsml = read_block(sname, "HSML", parttype=0)
            av_hsml = mean(hsml)
            k_av_hsml=2.0*pi/av_hsml
            plot([k_av_hsml,k_av_hsml], yrange, "k-..")
            if !nodetails
                text(0.98*k_av_hsml, ypos_text, L"k_{\rm{SML}}^{128}", horizontalalignment="right",fontsize=PyPlotHelper.title_font_size) #e-9
            end
        end
        
        # plot Nyquist mode and box mode
        k_Nyquist=pi/cellsize
        k_box=2*pi/gridsize
        if !(0.95*k_av_hsml < k_Nyquist < 1.05*k_av_hsml)
            plot([k_Nyquist,k_Nyquist], yrange,"k--")
        end
        if !nodetails
            if 0.95*k_av_hsml < k_Nyquist < 1.05*k_av_hsml
                text(0.98*k_Nyquist, 0.4*ypos_text, L"=k_{\rm{Nyquist}}^{128}", horizontalalignment="right",fontsize=PyPlotHelper.title_font_size)
            else
                text(0.98*k_Nyquist, ypos_text, L"k_{\rm{Nyquist}}^{128}", horizontalalignment="right",fontsize=PyPlotHelper.title_font_size)
            end
        end
        plot([k_box,k_box], yrange, "k--")
        if !nodetails
            text(1.02*k_box, ypos_text, L"k_{\rm{box}}", horizontalalignment="left",fontsize=PyPlotHelper.title_font_size)
        end
        
        # plot k_seed_min and k_seed_max (see function make_vel @ make_data.pro for more information)
        k_seed_min=(6.25/4)*2*pi/gridsize
        k_seed_max=(12.57/4)*2*pi/gridsize
        plot([k_seed_min,k_seed_min], yrange, "k--")
        if !nodetails
            text(1.02*k_seed_min, ypos_text, L"k_{\rm{SEED,min}}", horizontalalignment="left",fontsize=PyPlotHelper.title_font_size)
        end
        plot([k_seed_max,k_seed_max], yrange, "k--")
        if !nodetails
            text(1.02*k_seed_max, ypos_text, L"k_{\rm{SEED,max}}", horizontalalignment="left",fontsize=PyPlotHelper.title_font_size)
        end
    end
    if !overplot
        return kol_line, line
    else
        return line
    end
end

"""
    decay_comparison(fdir::String="./"; secondary_legend=false, internal_energy=false, total_energy=false)

Create comparison plot of turbulent energy decay for different methods, initial turbulent energy fractions, and resolution.
"""
function decay_comparison(fdir::String="./";
                          secondary_legend=false, # legend describing the points
                          internal_energy=false, # plot internal energy
                          total_energy=false) # plot only total energy instead

    oname = "decay_comparison_".*["method","x_mfm","x_sph","x_arepo","resolution","x_combined"].*".pdf"
    
    comparison = [[["mfm_128_cubic_detailed","sph_128_no_visc","sph_128_detailed","sph_128_high_visc","arepo_128_detailed","arepo_static_128_detailed"]],
                  [["mfm_128_cubic_detailed","mfm_128_x1_cubic_detailed","mfm_128_x01_cubic_detailed","mfm_128_x003_cubic_detailed","mfm_128_x001_cubic_detailed","mfm_128_x0003_cubic","mfm_128_x0001_cubic"]],
                  [["sph_128_detailed","sph_128_x1_detailed","sph_128_x01_detailed","sph_128_x003_detailed","sph_128_x001_detailed","sph_128_x0003_detailed","sph_128_x0001"]],
                  [["arepo_128_detailed","arepo_128_x1_detailed","arepo_128_x01_detailed","arepo_128_x003_detailed","arepo_128_x001_detailed","arepo_128_x0003_detailed","arepo_128_x0001_detailed"]],
                  [["mfm_64","mfm_128_cubic_detailed","mfm_256_cubic","sph_128_detailed","sph_256"]],
                  [["sph_128_detailed","sph_128_x1_detailed","sph_128_x01_detailed","sph_128_x003_detailed","sph_128_x001_detailed","sph_128_x0003_detailed","sph_128_x0001"], ["mfm_128_cubic_detailed","mfm_128_x1_cubic_detailed","mfm_128_x01_cubic_detailed","mfm_128_x003_cubic_detailed","mfm_128_x001_cubic_detailed","mfm_128_x0003_cubic","mfm_128_x0001_cubic"]]
                  ]
    description = [[["MFM","SPH "*L"(\alpha_{\rm{max}}^{\rm{visc}}=0)","SPH "*L"(\alpha_{\rm{max}}^{\rm{visc}}=3)","SPH "*L"(\alpha_{\rm{max}}^{\rm{visc}}=10)","AREPO","AREPO(static)"]],
                   [string.([0.3, 0.1, 0.01, 0.003, 0.001, 0.0003, 0.0001])],
                   [string.([0.3, 0.1, 0.01, 0.003, 0.001, 0.0003, 0.0001])],
                   [string.([0.3, 0.1, 0.01, 0.003, 0.001, 0.0003, 0.0001])],
                   [[L"64^2"*" (MFM)","128^3",L"256^3", L"128^3"*" (SPH)",L"256^3"]], #,L"512^3"]]
                   [string.([0.3, 0.1, 0.01, 0.003, 0.001, 0.0003, 0.0001]), nothing]]
    titles = [["hydro method"], ["MFM: "*L"X_i = "], ["SPH: "*L"X_i = "], ["AREPO: "*L"X_i = "],["resolution = "], ["SPH: "*L"X_i=", "MFM"]]
    ylims = [[[4e-1, 1.1e0]],
             [[4e-1, 2e0]],
             [[4e-1, 2e0]],
             [[7e-1, 1.1e0]],
             [[4e-1,2e0]],
             [[4e-1, 1.6e0],[4e-1, 1.6e0]]]
    
    # read some general info that stays the same for all methods
    #times = []
    for igrid in 1:16
        fname = fdir*"/out_turb_mfm_128_cubic/snap_"*sprintf1("%03d",igrid)
        
        h = read_header(fname)
        #time = h.time/t_sc    #divided by sound crossing time
        #append!(times,[time])
    end
    masses = read_block(fdir*"/out_turb_mfm_128_cubic/snap_000", "MASS", parttype=0)
    
    for i_comp in [1,6] #1:length(oname)
        n_subplots = length(comparison[i_comp])
        fig = figure(figsize=(12,8*n_subplots))
        style_plot(fig_width=12, print_columns=print_columns)
        gs = fig.add_gridspec(n_subplots,1, left=0.15, bottom=0.1, top=0.99, right=0.99, wspace=0,hspace=0)
        ax = gs.subplots()
        if n_subplots==1
            ax=[ax]
        end

        for i_subplot in 1:n_subplots
            methods = comparison[i_comp][i_subplot]
            
            kin_lines = []
            if internal_energy
                u_lines = []
            end
            fit_lines = []
            
            for i in 1:length(methods)
                if i_comp == 5
                    masses = nothing
                end
                tmp_kin,tmp_u,tmp_fit = plot_energy_decay(fdir*"/out_turb_"*methods[i],ax[i_subplot],fitting=true,times=-1,masses=masses,color_kin=get_color(i),color_u=get_color(i),label=methods[i],total_energy=total_energy, internal_energy=internal_energy)
                append!(kin_lines,[tmp_kin])
                if internal_energy
                    append!(u_lines,[tmp_u])
                end
                append!(fit_lines,[tmp_fit])
            end
            ax[i_subplot].set_ylim(ylims[i_comp][i_subplot])
            if total_energy
                ax[i_subplot].set_ylim([0.9,1.1])
            end
            ax[i_subplot].set_ylabel(L"E/E_i^{\rm{fit}}")
            ax[i_subplot].set_xlabel(L"t/t_{\rm{sc}}")
            #tight_layout()

            if description[i_comp][i_subplot] != nothing
                method_legend = ax[i_subplot].legend(kin_lines,description[i_comp][i_subplot],title=titles[i_comp][i_subplot],loc="lower left",ncol=2, bbox_to_anchor=[-0.02,0.01], borderpad=0, borderaxespad=0)
                if secondary_legend
                    if total_energy
                        ax[i_subplot].legend([kin_lines[1],fit_lines[1]],[L"E_{tot}","fit"],loc="lower center")
                    elseif internal_energy
                        ax[i_subplot].legend([kin_lines[1],u_lines[1],fit_lines[1]],[L"E_{turb}", L"U","fit"],loc="lower center")
                    else
                        ax[i_subplot].legend([kin_lines[1],fit_lines[1]],[L"E_{turb}", "fit"],loc="lower center")
                    end
                    ax[i_subplot].add_artist(method_legend)
                end
            else
                method_legend = ax[i_subplot].legend([],[],title=titles[i_comp][i_subplot],loc="lower left")
            end
            if i_subplot < n_subplots
                ax[i_subplot].set_xlabel("")
                ax[i_subplot].set_xticklabels([])
            end
        end
        
        fig.savefig(oname[i_comp])
        if close_all_figures
            close(fig)
        end
    end
    
end
"""
    plot_energy_decay(fdir::String="./",ax=1; color_kin="red",marker_kin="*",linestyle_kin=" ", color_u="red",marker_u="+",linestyle_u=" ", fitting::Bool=false,small_range::Bool=false,linestyle_fit="dashed", times=nothing,masses=nothing, label::String="energy", plotting=true, internal_energy=false, total_energy=false)

Plot decay of turbulent energy.
"""
function plot_energy_decay(fdir::String="./",ax=1;
                           color_kin="red",marker_kin="*",linestyle_kin=" ",
                           color_u="red",marker_u="+",linestyle_u=" ",
                           fitting::Bool=false,small_range::Bool=false,linestyle_fit="dashed",
                           times=nothing,masses=nothing,
                           label::String="energy",
                           plotting=true,
                           internal_energy=false,
                           total_energy=false) # plot total energy instead of single components

    snaps=readdir(fdir)
    snaps=snaps[@. startswith(snaps,"snap_")]
    n_snaps = length(snaps)

    if times == -1 # estimate times assuming final time 16
        times = collect(1:n_snaps) * 16/n_snaps / t_sc # in units of sound crossing time
        need_times = false
    elseif times == nothing || length(times) != n_snaps
        times = []
        need_times = true
    else
        need_times = false
    end
    
    E_kin_tot = []
    u_tot = []
    
    for igrid in 1:n_snaps
        fname = fdir*"/"*snaps[igrid]
        vel_info = InfoLine("VEL", Float32, 3, [1, 1, 1, 1, 1, 1])
        v = read_block(fname, "VEL", parttype=0, info=vel_info)
        u_info = InfoLine("U", Float32, 1, [1, 0, 0, 0, 0, 0])
        u = read_block(fname, "U", parttype=0,info=u_info)
        
        if masses == nothing
            masses = read_block(fname, "MASS", parttype=0)
        end

        if need_times
            h = read_header(fname)
            time = h.time/t_sc    #divided by sound crossing time
        end
        
        E_kin = Array{Float64,1}(undef,size(u)[1])
        @. E_kin = 1/2*masses*(v[1,:]^2 + v[2,:]^2 + v[3,:]^2)

        append!(E_kin_tot,[sum(E_kin)])
        append!(u_tot,[sum(masses.*u)])
        if need_times
            append!(times,[time])
        end
    end
    println(minimum(E_kin_tot)," ",maximum(E_kin_tot))
    println(minimum(u_tot)," ",maximum(u_tot))
    if total_energy
        println(minimum(E_kin_tot .+ u_tot)," ",maximum(E_kin_tot .+ u_tot))
    end

    if total_energy
        total_energy_line, = ax.semilogy(times, (E_kin_tot.+u_tot) ./ (E_kin_tot[2]+u_tot[2]), marker=marker_kin, linestyle=linestyle_kin)
        return total_energy_line,total_energy_line,total_energy_line # return 3 times so the integration in the plotting is easier
    end
    if fitting
        if small_range
            fit_range = times .> 10
        else
            fit_range = 2:length(times)
        end
        a,b = exp_fit(times[fit_range], E_kin_tot[fit_range])
        time_err = sqrt(mean(@. (b - (log(E_kin_tot[fit_range])-log(a))/times[fit_range])^2))/b^2
        println(fdir)
        println("decay time = ",-1/b," +- ",time_err)
        f=open("decay_times.dat","a")
        write(f,label," ",string(-1/b),"\n")
        close(f)

        if plotting
            e_kin_line, = ax.semilogy(times, E_kin_tot./a,color=color_kin,marker=marker_kin,linestyle=linestyle_kin)
            if internal_energy
                u_line, = ax.semilogy(times, u_tot./u_tot[1],color=color_u,marker=marker_u,linestyle=linestyle_u)
            end
            
            fit_line, = ax.semilogy(times, a.*exp.(b.*times)./a,color=color_kin,linestyle=linestyle_fit)

            if internal_energy
                return e_kin_line, u_line, fit_line
            else
                return e_kin_line, e_kin_line, fit_line
            end
        else
            # calculate error for decay time from statistics
            return -1/b , time_err # decay time
        end
    else
        if plotting
            e_kin_line, = ax.semilogy(times, E_kin_tot./E_kin_tot[2],color=color_kin,marker=marker_kin,linestyle=linestyle_kin)
            if internal_energy
                u_line, = ax.semilogy(times, u_tot./u_tot[1],color=color_u,marker=marker_u,linestyle=linestyle_u)
            end

            if internal_energy
                return e_kin_line, u_line
            else
                return e_kin_line, e_kin_line
            end
        else
            return
        end
    end
end

#method_comparison("../../../test_runs/",16,normalized=true,split_legend=true)
method_comparison("../../../test_runs/",16,normalized=true,split_legend=true,plot_all_methods=false,oname="power_spectrum_res_comparison_pure.pdf")
method_comparison("../../../test_runs/",16,normalized=true,split_legend=true,plot_all_resolutions=false,oname="power_spectrum_method_comparison_pure.pdf")

decay_comparison("../../../test_runs/")
