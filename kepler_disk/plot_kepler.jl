using GadgetIO
using PyPlot
using Formatting
using LaTeXStrings
using PyPlotHelper # private helper to format plots

"""
    plot_density_kepler(snap::String; ax=nothing, vmin=0,vmax=2, part="all", plot_initial=true)

Plot density for kepler disk.
"""
function plot_density_kepler(snap::String; ax=nothing, vmin=0,vmax=2, part="all", plot_initial=true)
    if !isfile(snap)
        println("snapshot "*snap*" not present, skip this plot.")
        ax.remove()
        return nothing
    end
    pos = read_block(snap, "POS", parttype=0)
    rho = read_block(snap, "RHO", parttype=0)

    if ax == nothing
        fig=figure(figsize=(4,4))
        style_plot(fig_width=4, print_columns=2)
        ax = fig.add_subplot()
    end

    if part == "all"
        range = 1:length(pos[1,:])
    elseif part == "upper"
        range = pos[2,:] .> 4
    elseif part == "upper right"
        range = pos[2,:] .> 8 .- pos[1,:]
    elseif part == "upper left"
        range = pos[2,:] .> pos[1,:]
    elseif part == "lower"
        range = pos[2,:] .< 4
    elseif part == "lower right"
        range = pos[2,:] .< pos[1,:]
    elseif part == "lower left"
        range = pos[2,:] .< 8 .- pos[1,:]
    else
        println("ERROR: part has to be one of ['all', 'upper( left/right)', 'lower( left/right)']")
    end
    
    im = ax.scatter(pos[1,range],pos[2,range], c=rho[range], marker=".",vmin=vmin,vmax=vmax, s=1e-1, cmap=get_colormap("density"), rasterized=true)

    if plot_initial
        # overplot IC
        phi = collect(0:0.01:2*pi)
        r = 0.5
        center = [4,4]
        ax.plot(r.*cos.(phi).+center[1],r.*sin.(phi).+center[2],"r--")

        r = 2
        center = [4,4]
        ax.plot(r.*cos.(phi).+center[1],r.*sin.(phi).+center[2],"r--")
    end
    
    ax.set_ylim([1.9,6.1])
    ax.set_xlim([1.9,6.1])
    
    if ax == nothing
        fig.show()
    end

    return im 
end

combined = true

plot_dir = "./"
methods = ["mfm", "sph", "mfm_gizmo", "arepo"]
labels = ["MFM", "SPH", "GIZMO", "AREPO"]
main_dir = "/e/ocean2/users/fgroth/test_runs/"
dirs = main_dir.*"out_kepler_disk_".*methods
if combined
    snap_nums = [[5,5,5,5,5],[48,48,48,48,48]]
else
    snap_nums = [5,21,39,48]
end
close_all_figures = false

vmin = 0
vmax = 2

if combined
    n_cols = 2
    n_rows = ceil(Integer, length(dirs) / n_cols)
else
    n_rows = length(dirs)
    n_cols = length(snap_nums)
end
colorbar_space = 1.3
upper_extension=0.02
fig=figure(figsize=(4*n_cols+colorbar_space,4*n_rows+upper_extension), dpi=150) #, dpi=300)
style_plot(fig_width=4*n_cols+colorbar_space, print_columns=2)
gs = fig.add_gridspec(n_rows,n_cols, hspace=0.03, wspace=0.03*(4*n_cols)/(4*n_cols+colorbar_space),
                      left=0.1*(4*n_cols)/(4*n_cols+colorbar_space),right=0.99-(colorbar_space)/(4*n_cols+colorbar_space),top=0.99/(1+upper_extension),bottom=0.1/(1+upper_extension))
ax = gs.subplots()

global image = nothing
i_dir_text  = ones(Integer, length(snap_nums))

for i_dir in 1:length(dirs)
    if combined
        # upper left
        snap1 = dirs[i_dir]*"/snap_"*sprintf1("%03d",snap_nums[1][i_dir])
        this_image = plot_density_kepler(snap1, ax=ax[i_dir], part="upper left", vmin=vmin, vmax=vmax)
        # lower right
        snap2 = dirs[i_dir]*"/snap_"*sprintf1("%03d",snap_nums[2][i_dir])
        if this_image != nothing
            ax[i_dir].plot([0,10],[0,10], color="grey", linestyle="dashdot")
            plot_density_kepler(snap2, ax=ax[i_dir], part="lower right", plot_initial=false, vmin=vmin, vmax=vmax)
            # labels
            if image == nothing
                global image = this_image
            end
            ax[i_dir].set_yticks([3,4,5])
            if i_dir <= n_rows
                ax[i_dir].set_yticklabels([-1,0,1])
                ax[i_dir].set_ylabel(L"y")
            else
                ax[i_dir].set_yticklabels([])
            end
            ax[i_dir].set_xticks([3,4,5])
            if i_dir % n_rows == 0
                ax[i_dir].set_xticklabels([-1,0,1])
                ax[i_dir].set_xlabel(L"x")
            else
                ax[i_dir].set_xticklabels([])
            end
            ax[i_dir].text(1.95,1.95,labels[i_dir], fontsize=PyPlotHelper.title_font_size)
            ax[i_dir].text(1.95,5.6,L"t="*sprintf1("%.1f", read_header(snap1).time), fontsize=PyPlotHelper.title_font_size, color="grey")
            ax[i_dir].text(6.05,1.95,sprintf1("%.1f", read_header(snap2).time), fontsize=PyPlotHelper.title_font_size, horizontalalignment="right", color="grey")
        end
        
    else
    ax[i_dir,1].set_ylabel(labels[i_dir])
        for i_snap in 1:length(snap_nums)
            snap = dirs[i_dir]*"/snap_"*sprintf1("%03d",snap_nums[i_snap])
            this_image = plot_density_kepler(snap,ax=ax[i_dir,i_snap], vmin=vmin, vmax=vmax)
            if this_image != nothing
                if image == nothing
                    global image = this_image
                end
                if (i_snap != 1)
                    ax[i_dir,i_snap].set_yticks([])
                else
                    ax[i_dir,i_snap].set_yticks([3,4,5])
                    ax[i_dir,i_snap].set_yticklabels([-1,0,1])
                end
                if (i_dir != length(dirs))
                    ax[i_dir,i_snap].set_xticks([])
                else
                    ax[i_dir,i_snap].set_xticks([3,4,5])
                    ax[i_dir,i_snap].set_xticklabels([-1,0,1])
                end
            end
            if i_dir == i_dir_text[i_snap]
                if this_image == nothing
                    i_dir_text[i_snap] = i_snap + 1
                else
                    ax[i_dir,i_snap].text(2.5,5.5,L"t\approx"*sprintf1("%03f", read_header(snap).time))
                end
            end
        end
    end
end

if combined
    if n_rows*n_cols > length(dirs)
        for i_ax in n_rows*n_cols:-1:length(dirs)+1
            ax[i_ax].remove()
        end
    end
end
cb_ax = fig.add_axes([(4*n_cols)/(4*n_cols+colorbar_space)+0.02*(4*n_cols)/(4*n_cols+colorbar_space), 0.1/(1+upper_extension), 0.018, 0.89/(1+upper_extension)]) # ???
cb = colorbar(image, cb_ax, label=L"\rho")
cb.set_ticks([0,0.5,1,1.5,2])

fig.savefig(plot_dir*"kepler.pdf")
if close_all_figures
    close(fig)
end
