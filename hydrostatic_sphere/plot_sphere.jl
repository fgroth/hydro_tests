using GadgetIO
using PyPlot
using PyPlotHelper # privat helper to format plots

using LaTeXStrings
using Formatting

using Statistics

"""
    plot_sphere(;oname="test", main_dir="../../test_runs/out_hydrostat_sphere_mfm/", snaps=[0,5,17,100,200], legend::Bool=false, nbins=100, ax=nothing, title=nothing, property="RHO", title_pos="upper")

Plot the radial profile of given `property`.

Plot all `snaps` inside `main_dir`. If `ax` is `nothing`,then create local figure and axis and save as `oname`.pdf
"""
function plot_sphere(;oname="test", main_dir="../../test_runs/out_hydrostat_sphere_mfm/", snaps=[0,5,17,100,200], legend::Bool=false, nbins=100, ax=nothing, title=nothing, property="RHO",title_pos="upper")

    if ax==nothing
        fig = figure(figsize=(4,4))
        ax = fig.add_subplot()
        style_plot(fig_width=4, print_columns=2)
        local_ax = true
    else
        local_ax = false
    end
    
    xrange = [14, 1e3]
    if property == "RHO"
        yrange = [1e-6,1e-4]
    elseif property == "U"
        yrange = [1e5,1.1e6]
    end

    colors = matplotlib.cm.plasma(0:1/length(snaps):1)
    for i in 1:length(snaps)
        color = colors[i, :]
        snap = main_dir*"snap_"*sprintf1("%03d",snaps[i])
        if !isfile(snap) && !isfile(snap*".0")
            continue # skip this
        end
        println("read "*snap)
        head = read_header(snap)
        println("time = ",head.time)

        mass = read_block(snap, "MASS", parttype=0)
        if minimum(mass) == 0 && maximum(mass) == 0
            if head.massarr[1] == 0
                head.massarr[1] = 0.0475000 # explicetely set the value I know from my IC
                println("use known mass")
            end
            mass = ones(Float64, length(mass)) .* head.massarr[1]
        end
        pos = read_block(snap, "POS", parttype=0)
        if property == "U"
            u = read_block(snap, "U", parttype=0)
        end

        pos0 = [0.0,0.0,0.0]
        for j in 1:3
            pos0[j] = mean(pos[j,:] .* mass)/mean(mass)
        end
        
        pos[1,:] = @. sqrt((pos[1,:]-pos0[1])^2 + (pos[2,:]-pos0[2])^2 + (pos[3,:]-pos0[3])^2)
        
        values = []
        radii = []
        for j in 1:nbins
            dx = (log10(xrange[2]) - log10(xrange[1]))/nbins
            x_lower = 10^(log10(xrange[1]) + (j-1)*dx)
            x_higher = 10^(log10(xrange[1]) + j*dx)
            v = 4/3*pi*x_higher^3
            
            bin = pos[1,:] .< x_higher #x_lower .< for slice

            if property == "RHO"
                rho_bin = sum(mass[bin]/v)
                append!(values, [rho_bin])
            elseif property == "U"
                u_bin = mean(u[bin])
                append!(values, [u_bin])
            end
            append!(radii, [(x_lower + x_higher)/2])
        end
        ax.loglog(radii, values, color=color)
    end
    
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    
    ax.set_xlabel(L"\log(r)")
    if property == "RHO"
        ax.set_ylabel(L"\rho")
    elseif property == "U"
        ax.set_ylabel(L"u")
    end

    if title != nothing
        ytitle_pos = 0.6*yrange[2]
        if title_pos=="lower"
            ytitle_pos = 1.3*yrange[1]
        end
        ax.text(1.12*xrange[1], ytitle_pos, title, fontsize=PyPlotHelper.title_font_size)
    end
    
    if local_ax
        fig.show()
        fig.savefig(oname*"_"*property*".pdf")
        if close_all_figures
            close(fig)
        end
    end
    
end

close_all_figures = false

"""
    create_comparison_plot( ; methods = ["mfm","sph","mfm_gizmo","arepo"], names = ["MFM", "SPH", "GIZMO", "AREPO"], snap_nums = 0:10:200, oname = "sphere.png", title_pos = "upper")

Create a plot comparing evolution of radial density and internal energy profiles for all `methods`, each at all `snap_nums`.

Use `names` as labels in legend, save file at `oname` 
"""
function create_comparison_plot( ; methods = ["mfm","sph","mfm_gizmo","arepo"],
                                 names = ["MFM", "SPH", "GIZMO", "AREPO"],
                                 snap_nums = 0:10:200,
                                 oname = "sphere.png",
                                 title_pos = "upper")
    
    fig = figure(figsize=(16,8))
    style_plot(fig_width=16, print_columns=1)
    gs = fig.add_gridspec(2,4,hspace=0,wspace=0, top=0.98, bottom=0.1, right=0.94, left=0.07)
    ax = gs.subplots()
    
    for i in 1:length(methods)
        if i == 1
            legend=true
        else
            legend=false
        end
        plot_sphere(main_dir="../../../test_runs/out_hydrostat_sphere_"*methods[i]*"/",ax=ax[1,i], title = names[i], property="RHO", title_pos=title_pos, snaps=collect(snap_nums))
        plot_sphere(main_dir="../../../test_runs/out_hydrostat_sphere_"*methods[i]*"/",ax=ax[2,i], property="U", snaps=collect(snap_nums))
        
        ax[1,i].set_ylim([0.9e-6,1e-4])
        ax[2,i].set_ylim([1e5,1.1e6])
        if i != 1
            ax[1,i].set_yticklabels([])
            ax[1,i].set_ylabel("")
            ax[2,i].set_yticklabels([])
            ax[2,i].set_ylabel("")
        end
        ax[1,i].set_xticklabels([])
        ax[1,i].set_xlabel("")
        ax[2,i].set_xlabel(L"r")
        
    end
    
    # add colormap to explain evolution. until t=10
    ax = fig.add_axes([0.94,0.1,0.015,0.98-0.1])
    cb = colorbar(matplotlib.cm.ScalarMappable(cmap="plasma"), ax, orientation="vertical")
    cb.set_ticks(collect(0:0.2:1))
    cb.set_ticklabels(collect(0:2:10))
    cb.set_label("t")
    
    fig.show()
    fig.savefig(oname)
end

#create_comparison_plot()
create_comparison_plot(methods=["mfm","mfm_testhll","mfm_testroe","mfm_testhllc"],
                       names=["MFM(Rs=exact)","MFM(Rs=HLL)","MFM(Rs=Roe)","MFM(Rs=HLLC)"],
                       oname="sphere_riemann.png",
                       title_pos="lower")
