using Makie
using CairoMakie

include("load_data.jl")

# global nextletter = 'a'

function condfixtime_diffusion_split(N, N0::Int64, r=1)
    return (1/r)*4*N0*(1-N0/N)
end


function condfixtime_diffusion_branch(N,N0::Int64, r=1)
    return 4N0/r /(1-N0/N) /(1+N0/(N-1))
end

function branchrate_effective_split(r, N, N0, b)
    return 1/(1/r + 1/(2b)*log(N^2/(N0*(N-N0))))
end

function branchrate_effective_branch(r, N, N_0, b)
    return 1/(1/r + 1/(2b)*log(N/N_0))
end

function condfixtime_moran(N, b)
    return (N-1)^2 / (N*b)
end

function plot_panel!(ax, df, r, N, N0, μ, update, b, tmax=Inf; xticks=0:2:4, yticks=0:500:1000, i=1)
    lightcolor = "#54585E"
    darkcolor=:black
    newdf = df[(
        (df.N .== N) .&
        (df.N0 .== N0) .&
        (df.μ .≈ μ) .&
        (df.r .≈ r) .&
        (df.update .== update) .&
        (df.t .< tmax)
    ), :]
    for gdf in groupby(newdf, :simid)
        lines!(ax, gdf.t, gdf.fixedmutations_mean, linewidth=0.3, color=(lightcolor, 0.5))
    end
    newdf_av = combine(
        groupby(newdf, :t),
        :fixedmutations_mean => mean
    )
    lines!(ax, newdf_av.t, newdf_av.fixedmutations_mean_mean, linewidth=1.4, color=darkcolor)
    tmax = isinf(tmax) ? maximum(newdf_av.t) : tmax
    # if update == :asymmetric
    #     lines!(ax, 0:0.1:tmax, t->μ*b*t, linestyle=:dash, color=:black, linewidth=1, label="μb")
    # elseif update == :Moran
    #     lines!(ax, 0:0.1:tmax, t->2*μ*b*t, linestyle=:dot, color=:black, linewidth=1.2, label="2μb")
    # end
    colors = ["#F792F8", "#7ADB9D"]
    tcond = condfixtime_diffusion_split(N, N0, branchrate_effective_split(r, N, N0, b))
    tcond < tmax && vlines!(ax, tcond, color=colors[1], linewidth=1.2)
    if update == :Moran
        vlines!(ax, condfixtime_moran(N, b), color=colors[2], linewidth=1.2)
    end
    ax.xticks = xticks
    ax.yticks = yticks
    # letter = nextletter
    # letterstr = rich("$(letter). ", font=:bold)
    # global nextletter += 1
    ax.title = rich("N = $N, N",subscript("0")," = $N0")
    ax
end

function plot(df, r, μ, gl, update, xticks, yticks; tmax=1, b=48, yvisible=true, i=1)
    axes = [Axis(gl[i, j]) for i in 1:2, j in 1:2]
    plot_panel!(axes[1, 1], df, r, 20, 1, μ, update, b, tmax; xticks, yticks, i)
    plot_panel!(axes[1, 2], df, r, 100, 1, μ, update, b, tmax; xticks, yticks, i)
    plot_panel!(axes[2, 1], df, r, 20, 10, μ, update, b, tmax; xticks, yticks, i)
    plot_panel!(axes[2, 2], df, r, 100, 50, μ, update, b, tmax; xticks, yticks, i)

    for i in 1:2, j in 2:2 hideydecorations!(axes[i,j], ticklabels=true, ticks=false) end
    for i in [1,2] hidexdecorations!(axes[1, i], ticklabels=true, ticks=false) end
    colgap!(gl, 12)
    rowgap!(gl, 15)
    Label(gl[1:2, 0], "Mean fixed mutations", 
        rotation = pi/2, tellwidth=true, tellheight=false, justification=:left, halign=:left, 
        valign=:center, padding=(3,0,0,0), visible=yvisible)
    Label(gl[3,1:2], "Time [yrs]", padding = (0, 0, 0, 10), height=20, valign=:top, tellheight=true)

    return axes
end

function fullplot(dfshort, b, μ, savename=nothing, tmax_moran=Inf, tmax_asym=Inf)
    set_theme!(mytheme)
    fig = Figure()
    gtl = fig[1,1] = GridLayout(halign=:left)
    gtr = fig[1,2] = GridLayout(halign=:left)

    axes_tl = plot(dfshort, 10, μ, gtl, :Moran, 0:1:2, 0:4:8; b, 
        tmax=tmax_moran, i=1)
    axes_tr = plot(dfshort, 10, μ, gtr, :asymmetric, 0:2:4, 0:4:8; b, 
        tmax=tmax_asym, yvisible=false, i=2)
    linkaxes!(axes_tl...)
    linkaxes!(axes_tr...)

    # Legend(
    #     fig[0,:], 
    #     [LineElement(linestyle=:dash, linewidth=0.7, linecolor=:black),
    #         LineElement(linestyle=:dot, linewidth=0.7, linecolor=:black)],
    #     ["Moran", "MWF-split"],
    #     "Conditional fixation time:",
    #     titleposition=:left,
    #     orientation=:horizontal,
    #     tellheight=:true,
    #     halign=:left,
    #     titlegap=16,
    #     colgap=16,
    #     margin=(0,0,10,0)
    # )
    isnothing(savename) || save(savename, fig, pt_per_unit=1)


    return fig
end


mytheme = Theme(
    fontsize = 20, 
    rowgap = 10, 
    colgap = 10, 
    resolution = 28.8 .* (24,14),
    Axis = (
        titlefont=:bold,
        titlealign=:center,
        subtitlefont=:italic, 
        titlesize = 16,
        # titlegap=10,
        # yticklabelspace=16.0,
        xticks=WilkinsonTicks(2), 
        yticks=WilkinsonTicks(2), 
        xgridvisible=false, 
        ygridvisible=false,
        xticksize=4,
        yticksize=4,
        rightspinevisible=false,
        topspinevisible=false,
        titlegap=12,
        ylabelpadding=8,

    ),
    Legend = (titleposition=:top, tellheight=false, tellwidth=true, framevisible=false, rowgap=0),
    Scatter = (;markersize=10)
)


b = 120
r = 40
# df_stats_by_sim10yr = get_dataframe("data/fulldist_data/data_tenyrs").df_stats_by_sim
# df_stats_by_sim10yr = change_timescale!(df_stats_by_sim10yr, 48, b)
fig = fullplot(df_stats_by_sim10yr, b, 0.01, "plots/poster_trajectories.pdf",
    2, 5)


# fig

