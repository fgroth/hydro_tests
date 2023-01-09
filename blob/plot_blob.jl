using GadgetIO
using PyPlot
using Formatting
using LaTeXStrings
using Statistics
using PyPlotHelper # private helper to format plots

"""
    plot_blob(snap::String; ax=nothing, x=1e-1, s="invrho", xlim=[0,6000], ylim=[0,2000])

Plot a density slice for the blob test.
"""
function plot_blob(snap::String; ax=nothing, x=1e-1, s="invrho", xlim=[0,6000], ylim=[0,2000])
    if ax == nothing
        fig,ax = subplots(111,figsize=(4,4))
    end
    
    if !isfile(snap) && !isfile(snap*".0")
        ax.remove()
        return
    end

    pos = read_block(snap, "POS", parttype=0)
    rho = read_block(snap, "RHO", parttype=0)

    rho_amb = 3.0515682e-8
    rho_blob = 3.1302562e-7
    vmin = 0.5*rho_amb
    vmax = 1.2*rho_blob # 4.5*rho_amb
    
    println("rho ",minimum(rho)," ",maximum(rho))

    slice = 1000-10 .< pos[2,:] .< 1000+10

    if s == "invrho"
        s = x * (rho[slice]/mean(rho[slice])) .^ (-1/3)
    end
    ax.scatter(pos[3,slice],pos[1,slice],c=rho[slice],cmap=get_colormap("density"),vmin=vmin,vmax=vmax,s=s)
    
    # IC: circle marking initial position
    rad = 0:0.01:2*pi+0.01
    ax.plot(200 .*cos.(rad).+2000,200 .*sin.(rad).+1000, color="green", linestyle="dashed")

    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
end

"""
    blob_comparison_plot(; methods=["mfm","sph","mfm_gizmo","arepo"],names=["MFM","SPH","GIZMO","AREPO"], oname="blob.png", snap_num=20, snap_num_insertion=nothing)

Create density comparison plot for the blob test.
"""
function blob_comparison_plot(; methods=["mfm","sph","mfm_gizmo","arepo"],names=["MFM","SPH","GIZMO","AREPO"], oname="blob.png", snap_num=20, snap_num_insertion=nothing)
    n_cols = 2
    n_rows = ceil(Int, length(methods)/n_cols)
    
    fig = figure(figsize=(4*n_cols, 2*n_rows))
    style_plot(fig_width=4*n_cols, print_columns=1)
    gs = fig.add_gridspec(n_rows, n_cols, left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=0.01, wspace=0.01)
    ax = gs.subplots()

    for i in 1:length(methods)
        snap = "../../../test_runs/out_blob_"*methods[i]*"/snap_"*sprintf1("%03d",snap_num)
        plot_blob(snap, ax=ax[i], xlim=[1500,5500])
        if snap_num_insertion != nothing
            insertion_snap = "../../../test_runs/out_blob_"*methods[i]*"/snap_"*sprintf1("%03d",snap_num_insertion)
            i_row = mod(i,n_rows) + 1
            i_col = floor(Int64,(i-1)/n_cols) + 1
            insertion_ax = fig.add_axes([i_col/n_cols*0.99-0.3, i_row/n_rows*0.99-0.2, 0.3, 0.2])
            plot_blob(insertion_snap, ax=insertion_ax, x=1e-2)
        end
        ax[i].text(1520,1700, names[i], fontsize=PyPlotHelper.title_font_size)
    end
    for i in length(methods)+1:n_rows*n_cols
        ax[i].remove()
    end
    
    fig.show()
    fig.savefig(oname)
end

"""
    plot_blob_decay(; methods=["mfm","sph","mfm_gizmo","arepo"],names=["MFM","SPH","GIZMO","AREPO"], oname="blob_decay.png")

Plot the decay of the fraction of the cloud surviving.
"""
function plot_blob_decay(; methods=["mfm","sph","mfm_gizmo","arepo"],names=["MFM","SPH","GIZMO","AREPO"], oname="blob_decay.png")
    fig = figure(figsize=(8,8))
    style_plot(fig_width=8, print_columns=2)
    ax = fig.add_subplot()

    lines = []
    
    # comparison lines from Hopkins (2015)
    tsph = [ 0.0320358578214  0.975703324808
             0.136998277572  0.918158567775
             0.235638603267  0.873657289003
             0.287973015293  0.854475703325
             0.334102771543  0.840664961637
             0.386190563182  0.837595907928
             0.434970509943  0.850639386189
             0.489731979748  0.872890025575
             0.587397567723  0.892071611253
             0.682424709014  0.883631713555
             0.737937783809  0.856777493606
             0.78770421212  0.80537084399
             0.840449658124  0.759335038363
             0.889945978391  0.72557544757
             0.939043008508  0.717902813299
             1.03693172921  0.722506393862
             1.13479696226  0.728644501279
             1.23577039511  0.731713554987
             1.33980505246  0.734782608696
             1.48669685265  0.737851662404
             1.58817527011  0.707928388747
             1.6405801451  0.684143222506
             1.70223915653  0.655754475703
             1.73639412287  0.624296675192
             1.83502270473  0.580562659847
             1.9547745185  0.556777493606
             2.06225794666  0.53452685422
             2.14506367764  0.524552429668
             2.26164857247  0.507672634271
             2.37547758234  0.47084398977
             2.45564486664  0.433248081841
             2.56027845921  0.397186700767
             2.66169815752  0.371099744246
             2.75680750561  0.357289002558
             2.92231327314  0.344245524297
             3.06331750091  0.331969309463
             3.16433790908  0.331969309463
             3.28387833394  0.32199488491
             3.34500887311  0.328132992327
             3.39378881988  0.341176470588
             3.46095960123  0.352685421995
             3.58356125059  0.342710997442
             3.69682655671  0.342710997442
             3.88958974894  0.348849104859
             3.94471527742  0.347314578005]
    mfm = [ 0.0384636463281  0.955754475703
            0.143520016702  0.892071611253
            0.276033456861  0.83452685422
            0.380490892009  0.809974424552
            0.508710005741  0.832992327366
            0.594271621692  0.842966751918
            0.680032882718  0.83989769821
            0.744729630983  0.813043478261
            0.79134088418  0.767774936061
            0.82569549559  0.723273657289
            0.856671799154  0.699488491049
            0.896737825565  0.681841432225
            0.936615950728  0.676470588235
            1.07419489535  0.687979539642
            1.14482619135  0.673401534527
            1.18802390521  0.651150895141
            1.29027741531  0.570588235294
            1.39629677958  0.443989769821
            1.46800850775  0.358823529412
            1.53638890339  0.291304347826
            1.60781877969  0.224552429668
            1.66054073803  0.180051150895
            1.74710057936  0.124808184143
            1.82113367086  0.0879795396419
            1.89813403622  0.0572890025575
            1.98425935592  0.0304347826087
            2.05486716426  0.0173913043478
            2.12538102197  0.0104859335038
            2.20504332168  0.00588235294118
            2.29997651234  0.00358056265985
            2.46225664179  0.00127877237852]
    sph_hopkins, = ax.plot(tsph[:,1],tsph[:,2],color="orange", linestyle="dashed")
    mfm_hopkins, = ax.plot(mfm[:,1],mfm[:,2],color="black", linestyle="solid")

    for i_method in 1:length(methods)
        snaps = readdir("../../../test_runs/out_blob_"*methods[i_method], join=true)
        println("plot "*methods[i_method])
        snaps = snaps[contains.(snaps, "snap_")]
        n_cloud0=1
        for i_snap in 1:length(snaps)
            rho = read_block(snaps[i_snap],"RHO",parttype=0)
            u = read_block(snaps[i_snap],"U",parttype=0)

            rho_blob = 3.1302562e-7
            u_amb = 105301.484
            cloud = (rho .> 0.64*rho_blob) .& (u .< 0.9*u_amb)
            n_cloud = sum(cloud)
            t_kh = 2

            if i_snap==1
                n_cloud0=n_cloud
            end
            time = read_header(snaps[i_snap]).time
            line, = ax.plot([time/t_kh],[n_cloud/n_cloud0],linestyle="",color=get_color(i_method),marker=get_marker(i_method))
            if i_snap == 1
                append!(lines,[line])
            end
        end
    end

    legend1 = ax.legend(lines,names,bbox_to_anchor=[0.68,0.68], borderpad=0, borderaxespad=0, loc="lower left")
    ax.legend([sph_hopkins,mfm_hopkins], ["SPH","MFM"], title="Hopkins (2015)", bbox_to_anchor=[0.68,0.66], borderpad=0, borderaxespad=0, loc="upper left")
    ax.add_artist(legend1)
    
    ax.set_xlim([0,4])
    ax.set_ylim([0,1])

    ax.set_xlabel(L"t/\tau_{\rm{KH}}")
    ax.set_ylabel(L"m_{\rm{cloud}}/m_{\rm{cloud}}(t=0)")

    fig.show()
    fig.savefig(oname)
end

blob_comparison_plot(snap_num_insertion=80)
plot_blob_decay()
