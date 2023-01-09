using GadgetIO

using PyPlot
using PyPlotHelper # private helper to format plots

using Formatting
using LaTeXStrings

using Statistics

"""
    plot_freefall(; dir::String="../../../test_runs/out_grav_freefall_", hydro_schemes=["mfm"], names=["MFM"])

Plot half-mass radius depending on time.
"""
function plot_freefall(; dir::String="../../../test_runs/out_grav_freefall_", hydro_schemes=["mfm"], names=["MFM"])

    fig = figure(figsize=(6,8))
    style_plot(fig_width=6, print_columns=2)
    gs = fig.add_gridspec(2,1,hspace=0.0, height_ratios=[1,4])
    ax = subplot(gs[2])
    ax_res = subplot(gs[1])
    ax_res.set_xticklabels([])

    lines = []

    for i_scheme in 1:length(hydro_schemes)
        snaps = readdir(dir * hydro_schemes[i_scheme])
        n_snaps = length(snaps[@. startswith(snaps, "snap_")])
        
        for i_snap in 1:n_snaps
            
            snap = dir * hydro_schemes[i_scheme] *"/snap_"* sprintf1("%03d", i_snap-1)

            head = read_header(snap)

            pos = read_block(snap, "POS", parttype=0, info=InfoLine("POS", Float32, 3, [1,1,1,1,1,1]) )
            mass = read_block(snap, "MASS", parttype=0, info=InfoLine("MASS", Float32, 1, [1,1,1,1,1,1]) )

            r12 = find_half_mass_radius(pos, mass)

            # theory
            if (i_scheme == 1) && (i_snap == 1)
                total_mass = sum(mass)

                r0 = find_half_mass_radius(pos, mass, full=true)
                rho0 = 3*total_mass / (4*pi* r0*r0*r0)
                tff = sqrt(3*pi/32/rho0)

                rdivr0 = collect(1e-4:1e-4:1)
                t = @. acos(sqrt(rdivr0)) + sqrt(rdivr0)*sqrt(1-rdivr0)

                println(r0)
                
                line_theory, = ax.plot(t*(tff*2/pi), rdivr0*r12, linestyle="solid", color="black")
                append!(lines, [line_theory])
                
                global t_ref = t*(tff*2/pi)
                global r_ref = rdivr0*r12
            end

            line = ax.scatter([head.time], [r12], marker=get_marker(i_scheme), color=get_color(i_scheme))
            if i_snap == 1
                append!(lines, [line])
            end

            # plot the residuals
            ax_res.plot([0,1.2], [0,0]) # reference line
            this_r_ref = minimum(r_ref[t_ref .<= head.time])
            ax_res.plot([head.time], [r12 / this_r_ref - 1], marker=get_marker(i_scheme), color=get_color(i_scheme))
        end
    end

    ax.set_xlim([0,1.2])
    ax.set_ylim([0,0.8])

    ax_res.set_ylim([-0.21,0.21])
    ax_res.set_xlim([0,1.2])
    ax_res.set_yticklabels(["","","0.0","0.2"])

    ax.set_ylabel(L"R_{1/2}")
    ax.set_xlabel(L"t")

    ax_res.set_ylabel("residuals")

    names = union(["analytical"], names)
    ax.legend(lines, names)
    fig.show()
    fig.savefig("grav_freefall.png")
    if close_all_figures
        close(fig)
    end
end

"""
    find_half_mass_radius(pos, mass; full::Bool=false)

Return half-mass radius for given particle positions `pos` and masses `mass`

For `full=true` return the full (maximum) radial extension. 
"""
function find_half_mass_radius(pos, mass; full::Bool=false) # for sure can be done faster.

    this_pos = deepcopy(pos)
    center = [mean(pos[1,:] .* mass), mean(pos[2,:] .* mass), mean(pos[3,:] .* mass)]./mean(mass)
    #println(center)
    pos[1,:] = @. sqrt((pos[1,:]-center[1])^2 + (pos[2,:]-center[2])^2 + (pos[3,:]-center[3])^2)
    if full
        pos = this_pos
        return maximum(pos[1,:])
    end
    total_mass = sum(mass)
    r_lower = -1
    r_higher = 1e99
    for r in pos[1,:]
        enclosed_mass = sum(mass[pos[1,:] .< r])
        if (enclosed_mass < 0.5*total_mass) && (r > r_lower)
            r_lower = r
        end
        if (enclosed_mass > 0.5*total_mass) && (r < r_higher)
            r_higher = r
        end
    end
    pos = this_pos
    return 0.5 * (r_lower + r_higher)
end

global close_all_figures = false

plot_freefall(hydro_schemes=["mfm", "sph", "mfm_gizmo", "arepo"], names=["MFM", "SPH", "MFM (GIZMO)", "AREPO"])
