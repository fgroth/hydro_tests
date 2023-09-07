using GadgetIO
using PyPlot
using PyPlotHelper # private helper to format plots
using Formatting
using Statistics
using CurveFit

"""
    residual(snap::String; debug::Bool=false, res=32, separate_amplitude_offset::Bool=false, separate_scatter_amplitude::Bool=false)

Return errors / residuals, comapring to soundwave.
"""
function residual(snap::String;
                  debug::Bool=false,
                  res=32,
                  separate_amplitude_offset::Bool=false,
                  separate_scatter_amplitude::Bool=false,
                  name="test",
                  arepo::Bool=false)
    
    time = read_header(snap).time
    x = if !arepo
        read_block(snap, "POS", parttype=0)[1,:]
    else
        read_block(snap, "CMCE", parttype=0, info=InfoLine("CMCE",Float32,3,[1,0,0,0,0,0]))[1,:]
    end
    rho = read_block(snap, "RHO", parttype=0)
    mass = read_block(snap, "MASS", parttype=0)
    # bin data per radius, so avoid scatter for fit
    xrho_bin = Array{Float64,2}(undef, round(Integer,res),2)
    for i_bin in 1:round(Integer,res)
        bin = (i_bin-1)/(res) .< x .<= i_bin/(res) 
        xrho_bin[i_bin,1] = mean(Float64.(x[bin]))
        xrho_bin[i_bin,2] = mean(Float64.(rho[bin]))
    end
    #rho0 = 0.5 * (maximum(xrho_bin[:,2]) + minimum(xrho_bin[:,2]))
    #rho0 = 1
    #rho0 = sum(mass) / (1*0.75*0.75)
    rho0 = mean(rho) # in principle, this would be slightly biased, as more weight to denser region, but givs the best estimate compared to other methods, as bias small for small amplitude.
    println("rho0 = ",rho0)
    u = read_block(snap, "U", parttype=0)
    u0 = mean(u)
    println("u0 = ",u0)

    vs = 0.66666666

    if separate_amplitude_offset
        x_rho = [Float64.(x) Float64.(rho)]
        coef0 = -4e-3
        #delta_rho = 0.5*(maximum(rho) - minimum(rho))
        delta_rho = 0.5*(maximum(xrho_bin[:,2]) - minimum(xrho_bin[:,2]))
        println("delta_rho ",delta_rho," = 0.5*( ",maximum(xrho_bin[:,2])-rho0," - ",minimum(xrho_bin[:,2])-rho0," )")
        # manual fit
        coef = fit_soundwave_offset(x_rho, coef0, rho0=rho0,res=res,time=time,delta_rho=delta_rho, debug=debug)
        offset_residual = coef
        if isfinite( offset_residual / time )
            vs += offset_residual / time
        end
        println("vs ",vs)
    end
    if separate_scatter_amplitude
        #delta_rho = 0.5*(maximum(rho) - minimum(rho)) # not the perfect choice, but good for a first test
        amplitude_residual = abs(delta_rho - 1e-4)
    else
        delta_rho = 1e-4
    end
    amplitude_scatter_residual = mean(abs.( rho .- rho_solution.(x, time, rho0=rho0, vs=vs, res=res, delta_rho=delta_rho) )) # actual residual
    println(amplitude_scatter_residual)
    
    if debug
        println("t = ",time)
        # mean substracted to allow easier comparison
        fig,ax = subplots(gridspec_kw=Dict("left"=>0.15,"right"=>0.99,"top"=>0.99,"bottom"=>0.1))
        println("rho0 ",rho0)
        ax.plot(x,rho .- rho0,"bo")
        ax.plot(xrho_bin[:,1],xrho_bin[:,2].-rho0, "b*") # binned data
        ax.plot(xrho_bin[:,1],rho_solution.(xrho_bin[:,1],time, vs=vs,res=res,delta_rho=delta_rho, rho0=rho0).-rho0, "ro") # delta rho from fit
        ax.plot(xrho_bin[:,1],rho_solution.(xrho_bin[:,1],time, vs=vs,res=res, rho0=rho0).-rho0, "r*") # delta_rho = 1e-4 (as it should be)
        ax.plot([0,1],[0,0],linestyle="dashed",color="black")
        ax.set_xlabel("x")
        ax.set_ylabel(L"\rho-\rho_0")
        if separate_amplitude_offset
            println("offset residual ",offset_residual)
        end
        if separate_scatter_amplitude
            println("amplitude residual ",amplitude_residual)
        end
        println("amplitude scatter residual ",amplitude_scatter_residual)
        fig.savefig(name*".png")
    end

    if separate_amplitude_offset
        if separate_scatter_amplitude
            return offset_residual, amplitude_scatter_residual, amplitude_residual
        else
            return offset_residual, amplitude_scatter_residual
        end
    end
    return amplitude_scatter_residual
end
"""
    rho_solution(x,t; rho0 = 1, delta_rho = 1e-4, k=2*pi/1, vs=0.66666666, res=64)

Return density of sinusodial soundwave with given parameters.
"""
function rho_solution(x,t;
                      rho0 = 1, delta_rho = 1e-4, # density and density perturbation
                      k=2*pi/1, # wavenumber (here defined by wavelength = size of box = 1)
                      vs=0.66666666, # soundspeed, should include any offset!
                      res=64) # to calculate offset due to the way how we create IC
    dx_ic = 1/res/2
    return rho0  + delta_rho * sin(k*(x+dx_ic+vs*t))
end
"""
    fit_soundwave_offset(xy, coef0; dcoef = 1e-2, iterations = 5, n_steps_per_direction = 10, rho0=nothing, delta_rho=nothing, res=64, time=time, debug::Bool=false)

Return fitted offset comparing to sinusodial soundwave with given parameters.
"""
function fit_soundwave_offset(xy, coef0;
                              dcoef = 1e-2,
                              iterations = 5,
                              n_steps_per_direction = 10,
                              rho0=nothing, #1,
                              delta_rho=nothing, #1e-4,
                              res=64,
                              time=time,
                              debug::Bool=false)
    # bin data per radius, so avoid scatter for fit
    xy_new = Array{Float64,2}(undef, res,2)
    for i_bin in 1:res
        bin = (i_bin-1)/(res) .< xy[:,1] .<= i_bin/(res) 
        xy_new[i_bin,1] = mean(xy[bin,1])
        xy_new[i_bin,2] = mean(xy[bin,2])
    end
    xy=xy_new
    if delta_rho == nothing
        delta_rho = 0.5*(maximum(xy[:;2]) - minimum(xy[2,:]))
    end
    rho0 = nothing
    if rho0 == nothing
        rho0 = 0.5*(maximum(xy[:,2]) + minimum(xy[:,2]))
        if debug
            println(xy)
            println("rho0 ",rho0)
        end
    end
    
    minimum_residual = mean(abs.(soundwave_residual_function(xy, [rho0, delta_rho,coef0], time=time,res=res)))
    coef_minimum_residual = coef0
    asymmetry = 1e99

    
    while iterations > 0
        if debug
            println("iteration ",iterations)
        end
        for i in -n_steps_per_direction:n_steps_per_direction
            residual = soundwave_residual_function(xy, [rho0, delta_rho, coef0+i/n_steps_per_direction * dcoef], time=time,res=res)
            residual = mean(abs.(residual))
            # check for symmetry in addition
            phase = mod.(xy[:,1] .+ (coef0+i/n_steps_per_direction * dcoef) .+ (1/res/2) .+ (0.6666666666*time), 1) # between 0 and 1, not 00 and 2pi!
            left_of_minimum = 0.25 .< phase .< 0.75
            right_of_minimum = (phase .< 0.25) .| (0.75 .< phase)
            left_residual = soundwave_residual_function(xy[left_of_minimum,:], [rho0, delta_rho, coef0+i/n_steps_per_direction * dcoef], time=time,res=res)
            left_residual = mean(left_residual)
            right_residual = soundwave_residual_function(xy[right_of_minimum,:], [rho0, delta_rho, coef0+i/n_steps_per_direction * dcoef], time=time,res=res)
            right_residual = mean(right_residual)
            
            if residual <= minimum_residual
                if debug
                    println("new minimum ",coef_minimum_residual," ",minimum_residual," ",residual)
                end
                minimum_residual = residual
                coef_minimum_residual = coef0+i/n_steps_per_direction * dcoef
            else
                #println(i, " ", coef0+i/n_steps_per_direction * dcoef, " ",dcoef, " ",coef0," ",residual-minimum_residual)
            end
            
        end
        iterations -= 1
        dcoef /= 5
        coef0 = coef_minimum_residual
    end

    return coef0
    # also vary the soundwave amplitude, as otherwise, some non-symmetric distribution might be preferred? Of make symmetry some criterion (though probably hard to implement, would need to check left and right wing, if mainly positive / negative errors)
end
"""
    soundwave_residual_function(x_rho, rho0_deltaRho_x0; time=2, res=64)

Return residual comparing to sinudodial soundwave with given parameters.
"""
function soundwave_residual_function(x_rho, rho0_deltaRho_x0; time=2, res=64)
    #unpack parameters
    x = x_rho[:,1]
    rho = x_rho[:,2]
    rho0 = rho0_deltaRho_x0[1]
    delta_rho = rho0_deltaRho_x0[2]
    x0 = rho0_deltaRho_x0[3]
    
    k = 2*pi/1
    vs = 0.6666666666 #* sqrt(u0/0.4)
    return @. rho0 + delta_rho * sin(k*(x + x0 + 1/res/2 + vs*time)) - rho
end
"""
    plot_residuals( ; main_dir = "/home/moon/fgroth/phd/test_runs/test_collection/test_runs/", hydro_methods = ["mfm", "sph", "mfm_gizmo", "arepo"], hydro_method_titles = ["MFM", "SPH", "GIZMO", "AREPO"], res = [32, 64, 128], snap_num=20, debug::Bool=false, separate_amplitude_offset::Bool=false, separate_scatter_amplitude::Bool=false)

Create plot comparing soundwave errors for different `hydro_methods`.
"""
function plot_residuals( ; main_dir = "/home/moon/fgroth/phd/test_runs/test_collection/test_runs/",
                         hydro_methods = ["mfm", "sph", "mfm_gizmo", "arepo"],
                         hydro_method_titles = ["MFM", "SPH", "GIZMO", "AREPO"],
                         res = [32, 64, 128],
                         snap_num=20,
                         debug::Bool=false,
                         separate_amplitude_offset::Bool=false,
                         separate_scatter_amplitude::Bool=false,
                         all_combined=false,
                         oname::String="soundwave_errors.png",
                         fit_convergence::Bool=false
                         )        
    if separate_amplitude_offset
        if separate_scatter_amplitude
            n_plots = 3
            fig_height=2.5
            i_offset=1
            i_amplitude=2
            i_scatter=3
        else
            n_plots = 2
            fig_height=1.75
            i_offset=1
            i_amplitude=2
        end
        bottom=0.08
    else
        n_plots = 1
        fig_height=1
        bottom=0.12
    end

    if all_combined
        separate_amplitude_offset=true
        separate_scatter_amplitude=true
        n_plots=2
        n2_plots=2
        fig_height=1.8
        fig_width=2
        i_offset=1
        i_amplitude=2
        i_scatter=3
        i_all=4
        bottom=0.1
        left=0.08
        print_columns=1.2
    else
        fig_width=1
        n2_plots=1
        left=0.195
        print_columns=2
    end
    
    fig = figure(figsize=(4*fig_width,4*fig_height))
    style_plot(fig_width=4*fig_width, print_columns=print_columns)
    gs = fig.add_gridspec(n_plots,n2_plots,
                          left=left, right=0.95, bottom=bottom, top=0.98,hspace=0,wspace=0.3)
    ax = gs.subplots()
    if n_plots*n2_plots == 1
        ax = [ax]
    end

    lines = Matrix(undef, length(hydro_methods), length(res))

    n_errors = 1 + separate_amplitude_offset + separate_scatter_amplitude + all_combined
    l1_norm = Array{Float64,3}(undef, length(hydro_methods), n_errors, length(res))

    for i_method in 1:length(hydro_methods)
        arepo = contains(hydro_methods[i_method], "arepo")
        for i_res in 1:length(res)
            snap = main_dir * "/out_soundwave_" * hydro_methods[i_method] * "_" * sprintf1("%d", res[i_res])*"/snap_"*sprintf1("%03d",snap_num)
            println(snap)
            # test if snap exists
            if !isfile(snap) && !isfile(snap*".0")
                println(snap*" does not exist")
                continue
            end
            residuals = residual(snap, debug=debug, res=res[i_res], separate_amplitude_offset=separate_amplitude_offset,separate_scatter_amplitude=separate_scatter_amplitude, name=hydro_methods[i_method]*"_"*sprintf1("%d",res[i_res]), arepo=arepo)
            if separate_amplitude_offset
                if separate_scatter_amplitude
                    lines[i_method, i_res], = ax[i_offset].loglog([res[i_res]], abs(residuals[1]), marker=get_marker(i_method), color=get_color(i_method), linestyle="")
                    l1_norm[i_method,i_offset, i_res] = abs(residuals[1])
                    ax[i_scatter].loglog([res[i_res]], residuals[2], marker=get_marker(i_method), color=get_color(i_method), linestyle="")
                    l1_norm[i_method,i_scatter, i_res] = residuals[2]
                    ax[i_amplitude].loglog([res[i_res]], residuals[3], marker=get_marker(i_method), color=get_color(i_method), linestyle="")       
                    l1_norm[i_method,i_amplitude, i_res] = residuals[3]
                    if all_combined
                        l1_norm[i_method,i_all,i_res] = residual(snap, debug=debug, res=res[i_res], separate_amplitude_offset=false,separate_scatter_amplitude=false, name=hydro_methods[i_method]*"_"*sprintf1("%d",res[i_res]), arepo=arepo)
                        ax[i_all].loglog([res[i_res]], l1_norm[i_method,i_all,i_res], marker=get_marker(i_method), color=get_color(i_method), linestyle="")
                    end
                else
                    lines[i_method, i_res], = ax[i_offset].loglog([res[i_res]], abs(residuals[1]), marker=get_marker(i_method), color=get_color(i_method), linestyle="")
                    l1_norm[i_method,i_offset, i_res] = abs(residuals[1])
                    ax[i_amplitude].loglog([res[i_res]], residuals[2], marker=get_marker(i_method), color=get_color(i_method), linestyle="")
                    l1_norm[i_method,i_amplitude, i_res] = residuals[2]
                end
            else
                lines[i_method, i_res], = ax[1].loglog([res[i_res]], residuals, marker=get_marker(i_method), color=get_color(i_method), linestyle="")
                l1_norm[i_method, 1, i_res] = residuals
            end
        end
    end
    for i_plot in 1:n_plots*n2_plots
        ax[i_plot].set_xticks([32,64,128])
        ax[i_plot].set_xticklabels([])
        #ax[i_plot].tick_params(axis="x", which="minor", bottom="off") #minorticks_off()
        ax[i_plot].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax[i_plot].xaxis.set_minor_locator(matplotlib.ticker.LogLocator(numticks=1,base=1e10)) # trick ro get rid of unwanted ticks
    end
    for i in 1:n2_plots
        ax[n_plots*i].set_xticklabels(["32","64","128"])
    end
    
    if separate_amplitude_offset
        ax[i_offset].set_ylabel("offset error")
        if separate_scatter_amplitude
            ax[i_scatter].set_ylabel("scatter")
            ax[i_amplitude].set_ylabel("amplitude error")
            for i in 1:n2_plots
                ax[n_plots*i].set_xlabel("resolution")
                if all_combined
                    ax[i_all].set_ylabel("L1(total)")
                end
            end
        else
            ax[i_amplitude].set_ylabel("amplitude/scatter error")
            ax[2].set_xlabel("resolution")
        end
    else
        ax[1].set_ylabel("L1")
        ax[1].set_xlabel("resolution")
    end

    if !fit_convergence # add lines for 1st and second order convergence to guide the eye
        xrange = [minimum(res),maximum(res)]
        if separate_amplitude_offset
            ax[i_offset].plot(xrange, xrange.^(-1) .* (1.0e-2*xrange[1]), linestyle="dashed", color="blue") # first order
            ax[i_offset].plot(xrange, xrange.^(-2) .* (2e-2*xrange[1]^2), linestyle="dotted", color="red") # second order
            # second line for SPH and Arepo
            ax[i_offset].plot(xrange, xrange.^(-1) .* (1.0e-3*xrange[1]), linestyle="dashed", color="blue") # first order
            if separate_scatter_amplitude
                _1st_order_line, = ax[i_scatter].plot(xrange, xrange.^(-1) .* (4e-6*xrange[1]), linestyle="dashed", color="blue") # first order
                _2nd_order_line, = ax[i_scatter].plot(xrange, xrange.^(-2) .* (3e-6*xrange[1]^2), linestyle="dotted", color="red") # second order
                # additional line for SPH and Arepo
                
                #ax[i_amplitude].plot(xrange, xrange.^(-1) .* (1e-5*xrange[1]), linestyle="dashed", color="blue") # first order
                ax[i_amplitude].plot(xrange, xrange.^(-2) .* (2e-5*xrange[1]^2), linestyle="dotted", color="red") # second order
                # additional line for SPH and Arepo
                ax[i_amplitude].plot(xrange, xrange.^(-1) .* (1.1e-6*xrange[1]), linestyle="dashed", color="blue") # first order
            end
        else
            _1st_order_line, = ax[1].plot(xrange, xrange.^(-1) .* (8e-6*xrange[1]), linestyle="dashed", color="blue") # first order
            _2nd_order_line, = ax[1].plot(xrange, xrange.^(-2) .* (1e-5*xrange[1]^2), linestyle="dotted", color="red") # second order
        end

        if all_combined
            convergence_labels = Vector{String}(undef,length(hydro_methods))
            for i_method in 1:length(hydro_methods)
                a, b = linear_fit(log.(res), log.(l1_norm[i_method, i_all,:]))
                ax[i_all].plot(res, exp.(a.+b.*log.(res)), color=get_color(i_method), linestyle="dashed")
                convergence_labels[i_method] = sprintf1("%1.1f",-b)
            end
            ax[i_all].legend(lines[:,1], hydro_method_titles.*L"\,(b=".*convergence_labels.*")", loc="upper right")
            ax[1].legend([_1st_order_line, _2nd_order_line], ["1st order", "2nd order"], title="convergence", loc="upper right")
        else
            method_legend = ax[1].legend(lines[:,1], hydro_method_titles, loc="upper right")
            ax[1].legend([_1st_order_line, _2nd_order_line], ["1st order", "2nd order"], title="convergence", loc="upper center")
            ax[1].add_artist(method_legend)
        end
    else # labels with convergence order
        convergence_labels = Vector{String}(undef,length(hydro_methods))
        for i_method in 1:length(hydro_methods)
            for i_error in 1:n_errors
                a, b = linear_fit(log.(res), log.(l1_norm[i_method, i_error,:]))
                ax[i_error].plot(res, exp.(a.+b.*log.(res)), color=get_color(i_method), linestyle="dashed")
                convergence_labels[i_method] = sprintf1("%1.1f",-b)
            end
        end
        # fill method names with spaces, so the convergence label aligns
        #max_length_name = maximum(length.(hydro_method_titles))
        #for i_method in 1:length(hydro_methods)
        #    while length(hydro_method_titles[i_method]) < max_length_name
        #        hydro_method_titles[i_method] *= " "
        #    end
        #end
        method_legend = ax[1].legend(lines[:,1], hydro_method_titles.*L"\,(b=".*convergence_labels.*")", loc="upper right")
    end
    #ax[1].set_ylim([2.1e-3,1.9e-2])
    #ax[2].set_ylim([0.3e-6, 2e-5])
    
    fig.savefig(oname)
end


plot_residuals(snap_num=30, res=[32,45,64,90,128], separate_amplitude_offset=false, separate_scatter_amplitude=false,
               hydro_methods = ["mfm", "sph", "mfm_gizmo","arepo"], hydro_method_titles = ["MFM", "SPH", "GIZMO","AREPO"],
               oname="soundwave_l1.pdf", fit_convergence=true)
plot_residuals(snap_num=30, res=[32,45,64,90,128], separate_amplitude_offset=true, separate_scatter_amplitude=true,
               hydro_methods=["mfm","sph","mfm_gizmo","arepo"], hydro_method_titles = ["MFM", "SPH", "GIZMO","AREPO"],
               oname="soundwave_errors.pdf")


# areopo:
plot_residuals(snap_num=30, res=[32,45,64,90,128], separate_amplitude_offset=true, separate_scatter_amplitude=true,all_combined=true,
               hydro_methods = ["arepo","arepo_improved"], hydro_method_titles = ["AREPO","AREPO(improved)"],
               oname="soundwave_arepo.pdf", debug=false)
