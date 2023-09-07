using GadgetIO
using PyPlot
using Formatting
using PyPlotHelper # private helper to format plots
using DelimitedFiles # to read in the exact solutions
using Statistics # to calculate mean
using LaTeXStrings

"""
    plot_zeldovich(snap::String; ax_v=nothing,ax_rho=nothing,ax_T=nothing,plot_exact::Bool=false,color="red",marker="+")

Plot zeldovich velocity, density and temperature profile.
"""
function plot_zeldovich(snap::String; ax_v=nothing,ax_rho=nothing,ax_T=nothing,plot_exact::Bool=false,color="red",marker="+")
    if !isfile(snap) && !isfile(snap*".0")
        return
    end
    
    if ax_v == nothing
        fig_v,ax_v = subplots(figsize=(4,4))
    end
    if ax_rho == nothing
        fig_rho,ax_rho = subplots(figsize=(4,4))
    end
    if ax_T == nothing
        fig_T0,ax_T = subplots(figsize=(4,4))
    end

    pos = read_block(snap, "POS", parttype=0) ./ 1e3
    rho = read_block(snap, "RHO", parttype=0)
    vel = read_block(snap, "VEL", parttype=0)
    u = read_block(snap, "U", parttype=0)

    # units and parameters
    m_unit = 1.989e43 #10^10 solar masses
    l_unit = 3.085678e21 #1.0 kpc /h
    v_unit = 1e5 #1 km/sec
    t_unit = l_unit / v_unit
    e_unit = m_unit * l_unit^2 / t_unit^2
    mp     = 1.6726231e-24 #[ g ]
    kcgs   = 1.380658e-16  #[ erg/K ]
    g1 = 5.0/3.0 - 1.0
    xH = 0.76
    xNe = 1.0
    nHe_fak = 3.0
    yhelium = (1. - xH) / (4. * xH)
    mu = (1. + 4. * yhelium) / (1. + nHe_fak * yhelium + xNe)
    
    if plot_exact
        data_exact=readdlm("zeldovich_exact.txt",Float64) # file from Hopkins (2015)
        ax_v.plot(data_exact[:,1],data_exact[:,4]./1e3)
        ax_rho.plot(data_exact[:,1],data_exact[:,2])
        ax_T.plot(data_exact[:,1],data_exact[:,3])
    end

    pos0 = mean(pos[1,:])
    # fold around periodicity
    pos[1,:] = @. mod(pos[1,:] + pos0,2*pos0)-pos0
    
    line = ax_v.scatter(pos[1,:], vel[1,:]./1e3,c=color,marker=marker,s=1, rasterized=true)

    T = @. g1 / kcgs * mp * u * e_unit / m_unit * mu
    ax_T.scatter(pos[1,:], log10.(T), c=color, marker=marker, s=1, rasterized=true)

    ax_rho.scatter(pos[1,:], log10.(rho./2.9e-8), c=color, marker=marker, s=1, rasterized=true)

    ax_rho.set_xlim([-32,32])
    ax_v.set_xlim([-32,32])
    ax_T.set_xlim([-32,32])
    ax_v.set_ylim([-2,2])
    ax_rho.set_ylim([-0.7,3.2])
    ax_T.set_ylim([-3,8.5])    
    
    if ax == nothing
        fig.show()
    end

    return line
        
end



main_dir="../../../test_runs/out_zeldovich_"

comparison = "general" # "general" / "switch"
if comparison == "general"
    methods = [["mfm"], ["sph"], ["mfm_gizmo"], ["arepo"]]
    name = ["MFM", "SPH", "GIZMO", "AREPO"]
    top=0.98
    bottom=0.08
    left=0.06
    right=0.99
    legend=false
    print_columns=1
elseif comparison == "switch"
    methods = [["mfm_pot0003", "mfm", "mfm_pot002"], ["mfm_kin0001", "mfm_kin0003", "mfm_kin0006"]]
    name = [L"U<\alpha_1 E_{pot}", L"U<\alpha_2 E_{kin}"]
    top=0.98
    bottom=0.08
    left=0.13
    right=0.99
    legend=true
    legend_title = [L"\alpha_1 = ", L"\alpha_2 ="]
    specific_names = [[L"3\cdot 10^{-3}", L"1\cdot 10^{-2}", L"2\cdot 10^{-2}"], [L"1\cdot10^{-3}", L"3\cdot10^{-3}", L"6\cdot 10^{-3}"]]
    print_columns=2
end

n_rows = 3 # v, rho, T
n_cols = length(methods)

fig = figure(figsize=(4*n_cols,4*n_rows*0.85), dpi=300)
style_plot(fig_width=4*length(methods), print_columns=print_columns)
gs = fig.add_gridspec(n_rows,n_cols, hspace=0.0, wspace=0.0, top=top, bottom=bottom, right=right, left=left)
ax = gs.subplots()

snap_num = 94
for i_col in 1:n_cols
    lines = []
    for i_inplot in 1:length(methods[i_col])
        println("plot "*main_dir*methods[i_col][i_inplot]*"/snap_"*sprintf1("%03d",snap_num))
        line = plot_zeldovich(main_dir*methods[i_col][i_inplot]*"/snap_"*sprintf1("%03d",snap_num),ax_v=ax[1,i_col],ax_rho=ax[2,i_col],ax_T=ax[3,i_col],plot_exact=true, color=get_color(i_inplot),marker=get_marker(i_inplot))
        append!(lines, [line])
    end
    if legend
        ax[1,i_col].legend(lines, specific_names[i_col], title = legend_title[i_col], loc="lower left", markerscale=5, bbox_to_anchor=[-0.02,0.01], borderpad=0, borderaxespad=0)
    end
    ax[1,i_col].text(31,1.95, name[i_col], fontsize=PyPlotHelper.title_font_size, horizontalalignment="right",verticalalignment="top")
    if i_col == 1
        ax[1,i_col].set_ylabel(L"\frac{v_x}{1000}"*" [km/s]")
        ax[2,i_col].set_ylabel(L"\log(\rho/\rho_0)")
        ax[3,i_col].set_ylabel(L"\log(T)"*" [K]")
    else
        ax[1,i_col].set_yticklabels([])
        ax[2,i_col].set_yticklabels([])
        ax[3,i_col].set_yticklabels([])
    end
    ax[3,i_col].set_xlabel(L"x~[h^{-1}"*"Mpc]")
end
fig.show()
fig.savefig("zeldovich_"*comparison*".pdf")
