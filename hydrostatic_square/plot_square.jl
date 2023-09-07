using GadgetIO
using PyPlot
using Formatting
using PyPlotHelper # private helper to format plots

"""
    plot_square(snap::String; ax=nothing,s=1e-3)

Plot density of the hydrostatic square test for given `snap`.

Set the marker size to `s`.
Plot into `ax`, if it is `nothing`, create local figure and axis.
"""
function plot_square(snap::String; ax=nothing,s=1e-3)
    if ax == nothing
        fig,ax = subplots(111,figsize=(4,4))
    end

    pos = read_block(snap, "POS", parttype=0)
    rho = read_block(snap, "RHO", parttype=0)

    ax.scatter(pos[1,:],pos[2,:],c=rho,cmap="plasma",vmin=0,vmax=5,s=s, rasterized=true)

    ax.plot([0.25, 0.25, 0.75, 0.75, 0.25], [0.25, 0.75, 0.75, 0.25, 0.25], color="green", linestyle="dashed")
    
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    
    if ax == nothing
        fig.show()
    end
        
end

"""
    create_square_plot(snaps, name; oname="square_general.png", n_cols=2, n_rows=nothing, s=4)

Create plot for hydrostatic square test, comparing different `snaps`.

Use `n_cols` columns and `n_rows` rows for setting up the gridspec. If `n_rows` is `nothing`, it is calculated from `n_cols`.
Use `name` as labels.
Set the marker size to `s`.
Save plot at `oname`
"""
function create_square_plot(snaps,name;
                            oname="square_general.png",
                            n_cols=2,
                            n_rows=nothing, # choose from n_rows and snaps
                            s=4,
                            print_columns=2)
    
    if n_rows == nothing
        n_rows = ceil(Int,length(snaps)/n_cols)
    elseif n_rows*n_cols < length(snaps)
        println("n_rows * n_cols has to be >=length(snaps)")
        return
    end
    
    fig = figure(figsize=(4*n_cols,4*n_rows), dpi=300)
    style_plot(fig_width=4*n_cols, print_columns=print_columns)
    gs = fig.add_gridspec(n_rows,n_cols, left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=0.01, wspace=0.01)
    ax = gs.subplots()

    if n_rows*n_cols == 1
        ax = [ax]
    end
    if length(s) == 1
        s = ones(n_rows*n_cols) .* s
    end
    
    for i in 1:length(snaps)
        println("plot "*snaps[i])
        plot_square(snaps[i],ax=ax[i],s=s[i])
        ax[i].text(0.05,0.82, name[i], fontsize=PyPlotHelper.title_font_size)
    end
    for i in length(snaps)+1:n_rows*n_cols
        ax[i].remove()
    end

    fig.show()
    if close_all_figures
        close(fig)
    end
    fig.savefig(oname)

end

main_dir="../../../test_runs/out_hydrostat_square_"

comparison = "general" #"general" / "resolution" / "limiter"
if comparison == "general"
    methods = ["mfm", "sph", "mfm_gizmo", "arepo"]
    name = ["MFM", "SPH", "GIZMO", "AREPO"]
    s = [4, 4, 4, 4]
    n_cols = 2
    print_columns=2
elseif comparison == "limiter"
    methods=["mfm", "mfm_springel2009", "mfm_tvd"]
    name = ["MFM(Li=GIZ)", "MFM(Li=ARE)","MFM(Li=TVD)"]
    s = [4,4,4]
    n_cols = 3
    print_columns=2
elseif comparison == "resolution"
    methods = ["mfm", "mfm_highres", "sph","sph_highres"]
    name = ["MFM", "", "SPH",""]
    s=[4, 0.1, 4, 0.1]
    n_cols = 2
    print_columns=2
end

snap_num = 20

snaps = main_dir .* methods .* "/snap_".*sprintf1("%03d",snap_num)

global close_all_figures=false
create_square_plot(snaps, name,s=s, oname="square_"*comparison*".pdf", n_cols=n_cols)

