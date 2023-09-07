using AnalyticMHDTestSolutions
using GadgetIO
using PyPlot
using LaTeXStrings
using Formatting

using PyPlotHelper # private helper to format plots

using PyCall
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

using Statistics

"""
    plot_density( ; data=nothing, snapdir::String="", snap_num::Integer=-1, ax=nothing,marker="+",color="red", no_labels::Bool=false, zoom_ax1=nothing, zoom_ax2=nothing)

Plot density profile for Sod shock tube test.
"""
function plot_density( ; data=nothing, snapdir::String="", snap_num::Integer=-1,
                       ax=nothing,marker="+",color="red",
                       no_labels::Bool=false,
                       zoom_ax1=nothing, zoom_ax2=nothing)
    if data == nothing
        if (snapdir  == "") || (snap_num == -1)
            println("ERROR: either data or snapdir and snap_num have to be provided")
            return
        end
        pos = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"POS",parttype=0)
        rho = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"RHO",parttype=0)
    else
        pos = data["pos"]
        rho = data["rho"]
    end
    line = ax.scatter(pos[1,:], rho, marker=marker, c=color, rasterized=true)
    if zoom_ax1 != nothing
        zoom_ax1.scatter(pos[1,:], rho, marker=marker, c=color, rasterized=true)
    end
    if zoom_ax2 != nothing
        zoom_ax2.scatter(pos[1,:], rho, marker=marker, c=color, rasterized=true)
    end
    
    if !no_labels
        ax.set_ylabel(L"\rho")
    end
    
    return line
end

"""
    plot_pressure( ; data=nothing, snapdir::String="", snap_num::Number=-1, ax=nothing,marker="+",color="red", gamma=1.4, no_labels=false, zoom_ax1=nothing, zoom_ax2=nothing)

Plot pressure profile for Sod shock tube test.
"""
function plot_pressure( ; data=nothing, snapdir::String="", snap_num::Number=-1,
                        ax=nothing,marker="+",color="red",
                        gamma=1.4,
                        no_labels=false,
                        zoom_ax1=nothing, zoom_ax2=nothing)
    if data == nothing
        if (snapdir  == "") || (snap_num == -1)
            println("ERROR: either data or snapdir and snap_num have to be provided")
            return
        end
        pos = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"POS",parttype=0)
        u = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"U",parttype=0)
        rho = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"RHO",parttype=0)
    else
        pos = data["pos"]
        u = data["u"]
        rho = data["rho"]
    end
    p = u.*(gamma-1).*rho
    ax.scatter(pos[1,:], p, marker=marker, c=color, rasterized=true)
    if zoom_ax1 != nothing
        zoom_ax1.scatter(pos[1,:], p, marker=marker, c=color, rasterized=true)
    end
    if zoom_ax2 != nothing
        zoom_ax2.scatter(pos[1,:], p, marker=marker, c=color, rasterized=true)
    end
    if !no_labels
        ax.set_ylabel(L"P")
    end
end

"""
    plot_velocity( ; data=nothing, snapdir::String="", snap_num::Number=-1, ax=nothing,marker="+",color="red", no_labels=false, zoom_ax1=nothing, zoom_ax2=nothing)

Plot velocity profile for Sod shock tube test.
"""
function plot_velocity( ; data=nothing, snapdir::String="", snap_num::Number=-1,
                        ax=nothing,marker="+",color="red",
                        no_labels=false,
                        zoom_ax1=nothing, zoom_ax2=nothing)
    if data == nothing
        if (snapdir  == "") || (snap_num == -1)
            println("ERROR: either data or snapdir and snap_num have to be provided")
            return
        end
        pos = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"POS",parttype=0)
        vel = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"POS",parttype=0)
    else
        pos = data["pos"]
        vel = data["vel"]
    end
    ax.scatter(pos[1,:],vel[1,:], marker=marker, c=color, rasterized=true)
    if zoom_ax1 != nothing
        zoom_ax1.scatter(pos[1,:], vel[1,:], marker=marker, c=color, rasterized=true)
    end
    if zoom_ax2 != nothing
        zoom_ax2.scatter(pos[1,:], vel[1,:], marker=marker, c=color, rasterized=true)
    end
    if !no_labels
        ax.set_ylabel(L"v")
    end
end

"""
    plot_entropy( ; data=nothing, snapdir::String="", snap_num::Number=-1, ax=nothing,marker="+",color="red", gamma=1.4, no_labels=false,  zoom_ax1=nothing,zoom_ax2=nothing)

Plot entropy profile for Sod shock tube test.
"""
function plot_entropy( ; data=nothing, snapdir::String="", snap_num::Number=-1,
                       ax=nothing,marker="+",color="red",
                       gamma=1.4,
                       no_labels=false,
                       zoom_ax1=nothing,zoom_ax2=nothing)
    if data == nothing
        if (snapdir  == "") || (snap_num == -1)
            println("ERROR: either data or snapdir and snap_num have to be provided")
            return
        end
        pos = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"POS",parttype=0)
        u = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"U",parttype=0)
        rho = read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"RHO",parttype=0)
    else
        pos = data["pos"]
        u = data["u"]
        rho = data["rho"]
    end
    A = u.*(gamma-1)./rho.^(gamma-1)
    ax.scatter(pos[1,:], A, marker=marker, c=color)
    if zoom_ax1 != nothing
        zoom_ax1.scatter(pos[1,:], A, marker=marker, c=color)
    end
    if zoom_ax2 != nothing
        zoom_ax2.scatter(pos[1,:], A, marker=marker, c=color)
    end
    if !no_labels
        ax.set_ylabel(L"A=P/\rho^{\gamma}")
    end
end

"""
    prepare_zoom(property::String, ax, par, solution)

Create and return zoom axes for `property` inside given axis `ax`.
"""
function prepare_zoom(property::String,
                      ax,
                      par, solution)

    pos_contact = solution.v34 * par.t + par.x_contact
    pos_shock = solution.vs * par.t + par.x_contact
    pos_raref = -solution.vt * par.t + par.x_contact
    pos_center = (pos_contact + pos_shock) / 2

    if property == "density"
        density_zoom_ax1 = inset_locator.inset_axes(ax, width="40%", height="40%", loc="upper right")
        density_zoom_ax1.set_xlim([pos_center+1, pos_center+5])
        density_zoom_ax1.set_ylim([min(solution.rho3*0.9, solution.rho4*0.9), max(solution.rho3*1.1, solution.rho4*1.1)])
        
        density_zoom_ax2 = inset_locator.inset_axes(ax, width="20%", height="25%", loc="lower left", borderpad=1.5)
        density_zoom_ax2.set_xlim([pos_raref-1.0,pos_raref+3.0])
        density_zoom_ax2.set_ylim(solution.rho3 .* [0.95, 1.05])
        return density_zoom_ax1, density_zoom_ax2
    elseif property == "pressure"
        pressure_zoom_ax1 = inset_locator.inset_axes(ax, width="49%", height="49%", loc="upper right")
        pressure_zoom_ax1.set_xlim([pos_shock-1,pos_shock+1])
        pressure_zoom_ax1.set_ylim([par.Pr-1.0, solution.P34+2.0])

        pressure_zoom_ax2 = inset_locator.inset_axes(ax, width="20%", height="25%", loc="lower left", borderpad=1.0)
        pressure_zoom_ax2.set_xlim([pos_contact-1.5,pos_contact+1.5])
        pressure_zoom_ax2.set_ylim([solution.P34-1.5, solution.P34+1.5])
        pressure_zoom_ax2.set_yticks([7,8,9])
        pressure_zoom_ax2.set_yticklabels(["","8","9"])
        pressure_zoom_ax2.set_xticks([88,89,90])
        pressure_zoom_ax2.set_xticklabels(["88","","90"])
        return pressure_zoom_ax1, pressure_zoom_ax2
    elseif property == "velocity"
        velocity_zoom_ax1 = inset_locator.inset_axes(ax, width="30%", height="30%", loc="lower center", borderpad=1.0)
        velocity_zoom_ax1.set_xlim([pos_shock-2,pos_shock+1])
        velocity_zoom_ax1.set_ylim([solution.v34-1.0, solution.v34+0.5])
        return velocity_zoom_ax1
    elseif property == "entropy"
        entropy_zoom_ax = nothing
        return entropy_zoom_ax
    end

end

"""
    comparison_plot(; plot_analytical=true,density=true,pressure=true,velocity=true,entropy=true, gamma=1.4,plot_zoom=false,oname=nothing, mach="1_5", hydro_methods=["mfm","sph","mfm_gizmo","arepo"], snap_num = 5)

Create comparison plot for Sod shock tube test.
"""
function comparison_plot(; plot_analytical=true,density=true,pressure=true,velocity=true,entropy=true, gamma=1.4,plot_zoom=false,oname=nothing,
                         mach="1_5", hydro_methods=["mfm","sph","mfm_gizmo","arepo"],
                         snap_num = 5)

    if typeof(mach) == String
        mach = [mach]
    end
    if typeof(hydro_methods) == String
        hydro_methods = [hydro_methods]
    end

    titlepos = nothing
    n_subplots=0
    if density
        n_subplots += 1
        i_density = n_subplots
        if titlepos == nothing
            titlepos = [70, 0.85]
        end
    end
    if pressure
        n_subplots += 1
        i_pressure = n_subplots
        if titlepos == nothing
            titlepos = [35, 45]
            pressure_yrange=[0,50]
        else
            pressure_yrange=[0,45]
        end
    end
    if velocity
        n_subplots += 1
        i_velocity = n_subplots
        if titlepos == nothing
            titlepos = [70, 30]
        end
    end
    if entropy
        n_subplots += 1
        i_entropy = n_subplots
        if titlepos == nothing
            titlepos = [45, 70]
        end
    end
    n_blocks = length(mach) * length(hydro_methods)

    height_frac = if n_subplots == 1
        1
    else
        0.61
    end
    fig = figure(figsize=(4*n_blocks,4*n_subplots*height_frac), dpi=300)
    style_plot(fig_width=4*n_blocks, print_columns=1)
    rc("lines", markersize=3)
    rc("legend", markerscale=2)
    bottom = 0.03 + 0.17/n_subplots
    gs = fig.add_gridspec(n_subplots,n_blocks, hspace=0, wspace=0, left=0.06, right=0.99, bottom=bottom, top=0.99)
    ax=gs.subplots(sharex=true, sharey="row")

    if n_subplots == 1 # make sure, we always have a matrix, not a vector
        ax_matrix = Matrix(undef, n_subplots, n_blocks)
        ax_matrix[1,:] .= ax
        ax = ax_matrix
    end
    if n_blocks == 1
        ax_matrix = Matrix(undef, n_subplots, n_blocks)
        ax_matrix[:,1] .= ax
        ax = ax_matrix
    end
    
    main_dir = "../../../test_runs/out_shock_"
    snapdirs = Matrix(undef, length(hydro_methods), length(mach))
    for i in 1:length(hydro_methods)
        for j in 1:length(mach)
            if (typeof(hydro_methods[i]) == String) & (typeof(mach[j]) == String)
                snapdirs[i,j] = [main_dir*hydro_methods[i]*"_M-"*mach[j]] # making this a list allows to do multiple lines per plot
            else
                snapdirs[i,j] = @. main_dir*hydro_methods[i]*"_M-"*mach[j] # already a list
            end
                      
        end
    end
    supply_info = Matrix{Bool}(undef, length(hydro_methods), length(mach))
    for i in 1:length(hydro_methods)
        supply_info[i,:] .= occursin("gizmo", hydro_methods[i])
    end

    names = Matrix{String}(undef, length(hydro_methods), length(mach))
    for i in 1:length(hydro_methods)
        for j in 1:length(mach)
            names[i,j] = "" # initialize
            if (j == 1) & (length(hydro_methods) != 1) # if more than 1 hydro method, include it's name, but only once
                if occursin("mfm", hydro_methods[i])
                    if occursin("gizmo", hydro_methods[i])
                        names[i,j] = "GIZMO"
                    else
                        names[i,j] = "MFM"
                    end
                elseif occursin("sph", hydro_methods[i])
                    names[i,j] = "SPH"
                elseif occursin("arepo", hydro_methods[i])
                    names[i,j] = "AREPO"
                end
                if length(mach) != 1 # prepare to patch together
                    names[i,j] = names[i,j]*", "
                end
            end
            if (i == 1) & (length(mach) != 1) # if more than 1 mach number, include it's value, but only once
                names[i,j] = names[i,j]*L"\mathcal{M}="*replace(mach[j],"_"=>".")
            end
        end
    end

    xmin = 35
    xmax = 105

    # create four columns: MFM, SPH, gizmo, Arepo
    for i_block in 1:length(snapdirs)

        println("plot ",snapdirs[i_block][1])
        # initialize zoom axes to nothing
        density_zoom_ax1 = density_zoom_ax2 = nothing
        pressure_zoom_ax1 = pressure_zoom_ax2 = nothing
        velocity_zoom_ax1 = velocity_zoom_ax2 = nothing
        entropy_zoom_ax1 = entropy_zoom_ax2 = nothing
        if plot_analytical || plot_zoom # prepare solution
            x = collect(xmin:(xmax-xmin)/4e3:xmax)
            # obtain input parameters
            rho = read_block(snapdirs[i_block][1]*"/snap_000", "RHO", parttype=0, info=get_info("RHO"))
            pos = read_block(snapdirs[i_block][1]*"/snap_000", "POS", parttype=0, info=get_info("POS"))
            U = read_block(snapdirs[i_block][1]*"/snap_000", "U", parttype=0, info=get_info("U"))
            time = read_header(snapdirs[i_block][1]*"/snap_"*sprintf1("%03d",snap_num)).time
            par = AnalyticMHDTestSolutions.SodParameters(
                rhor = Float64(mean(rho[75 .< pos[1,:] .< 135])), rhol = Float64(mean(rho[5 .< pos[1,:] .< 65])),
                Ur = Float64(mean(U[75 .< pos[1,:] .< 135])), Ul = Float64(mean(U[5 .< pos[1,:] .< 65])),
                t = time, γ_th = gamma )
            solution = solve(x, par)
            println("rho ",par.rhol," ",par.rhor)
            println("P ",par.Pl," ",par.Pr)
            println("Mach ",par.M)

            # insert zoom plots
            if density
                if plot_zoom
                    density_zoom_ax1, density_zoom_ax2 = prepare_zoom("density", ax[i_density, i_block],
                                                                     par, solution)
                end
            end
            if pressure
                if plot_zoom
                    pressure_zoom_ax1, pressure_zoom_ax2 = prepare_zoom("pressure", ax[i_pressure, i_block],
                                                              par, solution)
                end
            end
            if velocity
                if plot_zoom
                    velocity_zoom_ax1 = prepare_zoom("velocity", ax[i_velocity, i_block],
                                                     par, solution)
                end
            end
            if entropy
                if plot_zoom
                    entropy_zoom_ax = prepare_zoom("entropy", ax[i_entropy, i_block],
                                                   par, solution)
                end
            end
        end

        lines = []

        for i in 1:length(snapdirs[i_block]) # loop of lines shown per axis
            snapdir = snapdirs[i_block][i]
            if supply_info[i_block][i]
                output_type = Float32
                pos_info = InfoLine("POS", output_type, 3, [1, 1, 1, 1, 0, 1])
                rho_info = InfoLine("RHO", output_type, 1, [1, 1, 0, 0, 0, 0])
                vel_info = InfoLine("VEL", output_type, 3, [1, 1, 1, 1, 0, 1])
                u_info = InfoLine("U", output_type, 1, [1, 0, 0, 0, 0, 0])
            else
                pos_info = nothing
                rho_info = nothing
                vel_info = nothing
                u_info = nothing
            end
            data = Dict(
                "pos" => read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"POS",parttype=0, info=pos_info),
                "rho" => read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"RHO",parttype=0, info=rho_info),
                "vel" => read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"VEL",parttype=0, info=vel_info),
                "u" => read_block(snapdir*"/snap_"*sprintf1("%03d",snap_num),"U",parttype=0, info=u_info)
            )
            # only use a narrow tube
            slice = (0.45 .< data["pos"][2,:] .< 0.65) .&& (0.45 .< data["pos"][3,:] .> 0.65)
            data["pos"] = data["pos"][:,slice]
            data["rho"] = data["rho"][slice]
            data["vel"] = data["vel"][:,slice]
            data["u"] = data["u"][slice]
            
            if density
                line = plot_density(data=data,ax=ax[i_density,i_block], color=get_color(i), marker=get_marker(i), no_labels=(i_block!=1 || i!=1),
                                    zoom_ax1=density_zoom_ax1, zoom_ax2=density_zoom_ax2)
                append!(lines,[line])
                ax[i_density,i_block].set_ylim([0,1.1])
                ax[i_density,i_block].set_yticks([0,0.5,1])
                ax[i_density,i_block].set_yticklabels(["0","0.5","1.0"])
            end
            if pressure
                plot_pressure(data=data,ax=ax[i_pressure,i_block], color=get_color(i), marker=get_marker(i), no_labels=(i_block!=1 || i!=1),
                              zoom_ax1=pressure_zoom_ax1, zoom_ax2=pressure_zoom_ax2)
                ax[i_pressure,i_block].set_ylim(pressure_yrange)
                ax[i_pressure,i_block].set_yticks([0,20,40])
                ax[i_pressure,i_block].set_yticklabels(["0","20","40"])
            end
            if velocity
                plot_velocity(data=data,ax=ax[i_velocity,i_block], color=get_color(i), marker=get_marker(i), no_labels=(i_block!=1 || i!=1),
                              zoom_ax1=velocity_zoom_ax1, zoom_ax2=velocity_zoom_ax2)
                ax[i_velocity,i_block].set_ylim([-1,11])
            end
            if entropy
                plot_entropy(data=data,ax=ax[i_entropy,i_block], color=get_color(i), marker=get_marker(i), no_labels=(i_block!=1 || i!=1),
                             zoom_ax1=entropy_zoom_ax1, zoom_ax2=entropy_zoom_ax2)
                ax[i_entropy,i_block].set_ylim([0,110])
            end
        end

        if plot_analytical # overplot the solution
            if density
                ax[i_density,i_block].plot(x, solution.rho, color="black", linestyle="-")
                if plot_zoom
                    density_zoom_ax1.plot(x, solution.rho, color="black", linestyle="-")
                    inset_locator.mark_inset(ax[i_density,i_block], density_zoom_ax1, 3,4, alpha=0.8, zorder=1000.0, linestyle=":")

                    density_zoom_ax2.plot(x, solution.rho, color="black", linestyle="-")
                    inset_locator.mark_inset(ax[i_density,i_block], density_zoom_ax2, 1,4, alpha=0.8, zorder=1000.0, linestyle=":")
                end
            end
            if pressure
                ax[i_pressure,i_block].plot(x, solution.P, color="black", linestyle="-")
                if plot_zoom
                    pressure_zoom_ax1.plot(x, solution.P, color="black", linestyle="-")
                    inset_locator.mark_inset(ax[i_pressure,i_block], pressure_zoom_ax1, 3,4, alpha=0.8, zorder=1000.0, linestyle=":")
                    pressure_zoom_ax2.plot(x, solution.P, color="black", linestyle="-")
                    inset_locator.mark_inset(ax[i_pressure,i_block], pressure_zoom_ax2, 1,4, alpha=0.8, zorder=1000.0, linestyle=":")
                end
            end
            if velocity
                ax[i_velocity,i_block].plot(x, solution.v, color="black", linestyle="-")
                if plot_zoom
                    velocity_zoom_ax1.plot(x, solution.v, color="black", linestyle="-")
                    inset_locator.mark_inset(ax[i_velocity,i_block], velocity_zoom_ax1, 2,4, alpha=0.8, zorder=1000.0, linestyle=":")
                end
            end
            if entropy
                ax[i_entropy,i_block].plot(x, solution.P ./ solution.rho .^gamma, color="black", linestyle="-")
                if plot_zoom
                    # ???
                end
            end
        end

        ax[1,i_block].text(titlepos[1], titlepos[2], names[i_block], fontsize=PyPlotHelper.title_font_size)
        
        ax[n_subplots,i_block].set_xlabel(L"x")
        ax[1,i_block].set_xlim([xmin,xmax])
        
    end
    
    if oname == nothing
        oname = "shock"
        for m in mach
            oname = oname*"_M-"*m
        end
        oname = oname*".pdf"
    end
    fig.savefig(oname)
end

# full comparison

#comparison_plot(mach="1_5", plot_zoom=true)
#comparison_plot(mach="3_0", plot_zoom=true)
#comparison_plot(mach="10_0", plot_zoom=true)
#comparison_plot(mach="100_0", plot_zoom=true)

# orthogonal comparison

comparison_plot(mach=["1_5","3_0","10_0","100_0"], plot_zoom=false,
                density=true,pressure=true,velocity=true,entropy=true,
                hydro_methods="mfm",
                oname="shock_mfm.pdf", snap_num=25)
comparison_plot(mach="10_0", plot_zoom=true,
                density=false,pressure=true,velocity=false,entropy=false,
                hydro_methods=["mfm","sph","mfm_gizmo","arepo"],
                oname="shock_methods.pdf", snap_num=25)
