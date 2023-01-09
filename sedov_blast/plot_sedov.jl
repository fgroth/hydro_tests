using GadgetIO

using PyPlot
using PyPlotHelper # private helper to format plots

using Formatting
using AnalyticMHDTestSolutions
using LaTeXStrings

"""
    plot_sedov( ; main_dir::String="../../../test_runs/out_sedov_", snap_num::Int64=6, method=["mfm", "sph","mfm_gizmo", "arepo"], names=["MFM", "SPH", "GIZMO", "AREPO"])

Create comparison plot of radial density profile for sedov blast wave test.
"""
function plot_sedov( ; main_dir::String="../../../test_runs/out_sedov_",
                    snap_num::Int64=6, method=["mfm", "sph","mfm_gizmo", "arepo"], names=["MFM", "SPH", "GIZMO", "AREPO"])

    fig = figure(figsize=(8,8))
    style_plot(fig_width=8, print_columns=2)
    rc("lines", markersize=3)
    gs = fig.add_gridspec(2,2, hspace=0,wspace=0)
    ax = gs.subplots()

    lines = []
    for i_method in 1:length(method)

        snap = main_dir * method[i_method] * "/snap_" * sprintf1("%03d",snap_num)
        pos = read_block(snap, "POS", parttype=0)
        rho = read_block(snap, "RHO", parttype=0)
        pos0 = [0.5, 0.5, 0.5]
        pos[1,:] = @. sqrt((pos[1,:]-pos0[1])^2 + (pos[2,:]-pos0[2])^2 + (pos[3,:]-pos0[3])^2)

        head = read_header(snap)
        println("t = ", head.time)
        
        # analytical solution
        par = SedovParameters(head.time, 10,
                        1e-6/(5/3-1)/1, 1)
        x = collect(minimum(pos[1,:]):(maximum(pos[1,:])-minimum(pos[1,:]))/4e3:maximum(pos[1,:]))
        sedov_ideal = solve(x, par)
        
        ax[i_method].plot(sedov_ideal.r, sedov_ideal.rho, color="black")
        println(maximum(sedov_ideal.rho))
        #end
        
        line = ax[i_method].scatter(pos[1,:], rho, marker=get_marker(1), color=get_color(1),s=1e-1)
        append!(lines,[line])

        ax[i_method].text(0.05,4,names[i_method],fontsize=PyPlotHelper.title_font_size)
        ax[i_method].set_xlim([0,0.5])
        ax[i_method].set_ylim([0,4.5])
    end

    ax[1,1].set_xticklabels([])
    ax[1,2].set_xticklabels([])
    ax[2,1].set_xlabel(L"r")
    ax[2,2].set_xlabel(L"r")
    
    ax[1,2].set_yticklabels([])
    ax[2,2].set_yticklabels([])
    ax[1,1].set_ylabel(L"\rho")
    ax[2,1].set_ylabel(L"\rho")
        
    fig.savefig("sedov.png")
    if close_all_figures
        close(fig)
    end
    
end

global close_all_figures = false

plot_sedov(snap_num=16)
