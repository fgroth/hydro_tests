using GadgetIO
using AnalyticMHDTestSolutions
using Formatting
using PyPlot

main_dir = "../../../test_runs/out_sedov_mfm/"
gamma=5/3

rho_mean = 1.0
Pr = 0.0

for i_snap in 0:20
    u = read_block(main_dir*"snap_"*sprintf1("%03d",i_snap), "U", parttype=0)
    rho = read_block(main_dir*"snap_"*sprintf1("%03d",i_snap), "RHO", parttype=0)
    
    
    p = u.*(gamma-1).*rho
    if i_snap == 0
        global Pr = minimum(p) # at snap 0 to get unperturbed background
    end

    # calculate the Mach number of the shock
    mach = AnalyticMHDTestSolutions.solveMach(maximum(p), Pr,
                                              Float64(maximum(rho)), rho_mean,
                                              gamma)
    mach_alt = AnalyticMHDTestSolutions.solveMach(maximum(p), Pr,
                                                  rho_mean, rho_mean,
                                                  gamma)
    
    println("snap ",i_snap," ",mach," ",mach_alt)
    plot([i_snap],[mach],"bo")
    plot([i_snap],[mach_alt],"ro")
end
