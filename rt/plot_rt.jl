using GadgetIO
using PyPlot
using Formatting
using PyPlotHelper # private helper to format plots

"""
    plot_rt(snap::String; ax=nothing,s=1e-1, vmin=0.5, vmax=2.5, problem = "rt")

Plot density distribution for kh or rt instability, set `problem` accordingly.
"""
function plot_rt(snap::String; ax=nothing,s=1e-1, vmin=0.5, vmax=2.5, problem = "rt")
    local_ax = false
    if ax == nothing
        fig,ax = subplots(111,figsize=(4,4))
        local_ax=true
    end

    if !isfile(snap) && !isfile(snap*".0")
        return
    end

    pos = read_block(snap, "POS", parttype=0)
    rho = read_block(snap, "RHO", parttype=0)

    println(minimum(rho),maximum(rho))
    
    ax.scatter(pos[1,:],pos[2,:],c=rho,cmap=get_colormap("density"),vmin=vmin,vmax=vmax,s=s)

    if problem == "rt"
        ax.plot([0, 1], [0.5, 0.5], color="green", linestyle="dashed")
    elseif problem == "kh"
        ax.plot([0,1],[0.25,0.25], color="green", linestyle="dashed")
        ax.plot([0,1],[0.75,0.75], color="green", linestyle="dashed")
    end
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    
    if local_ax
        fig.show()
    end
        
end

"""
    create_rt_plot(snaps,name; oname="rt_general.png", problem="rt", s=1e-1)

Create comparison plot for kh / rt instability.
"""
function create_rt_plot(snaps,name;
                        oname="rt_general.png",
                        problem="rt",
                        s=1e-1,
                        n_cols=2,
                        ypos_label=0.82,
                        print_columns=2)
    n_rows = ceil(Int, length(snaps)/n_cols)
    
    fig = figure(figsize=(4*n_cols,4*n_rows))
    style_plot(fig_width=4*n_cols, print_columns=print_columns)
    gs = fig.add_gridspec(n_rows,n_cols, left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=0.01, wspace=0.01)
    ax = gs.subplots()

    for i in 1:length(snaps)
        println("plot "*snaps[i])
        plot_rt(snaps[i],ax=ax[i],s=s[i], problem=problem)
        ax[i].text(0.05,ypos_label, name[i], fontsize=PyPlotHelper.title_font_size)
    end
    for i in length(snaps)+1:n_rows*n_cols
        ax[i].remove()
    end
    fig.show()
    fig.savefig(oname)
    if close_all_figures
        close(fig)
    end
end

problem = "rt" # "kh" / "rt"
main_dir="../../../test_runs/out_"*problem*"_"
if problem == "kh"
    ypos_label = 0.9
elseif problem == "rt"
    ypos_label = 0.9
end

comparison = "limiter" # "general" / "limiter"
if comparison == "general"
    methods = ["mfm", "sph", "mfm_gizmo", "arepo"]
    name = ["MFM", "SPH", "GIZMO", "AREPO"]
    s = [3e-1, 3e-1, 3e-1, 3e-1]
    n_cols = 2
    print_columns=1
elseif comparison == "limiter"
    methods = ["mfm", "mfm_arepo", "mfm_tvd"]
    name = ["", "", ""]
    s = [3e-1, 3e-1, 3e-1]
    n_cols = 3
    print_columns=2
    ypos_label=0.82
end

if problem == "rt"
    snap_num = 36
elseif problem == "kh"
    snap_num = 5
end

global close_all_figures=false

create_rt_plot(main_dir.*methods.*"/snap_".*sprintf1("%03d",snap_num),name,oname=problem*"_"*comparison*".png",problem=problem,s=s, n_cols=n_cols, ypos_label=ypos_label, print_columns=print_columns)
