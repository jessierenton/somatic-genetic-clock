using Makie
using CairoMakie

include("plot_init.jl")
include("load_data.jl")

#region Theory
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
#endregion

#region Plot functions
function plot_panel!(ax, df, r, N, N0, μ, update, b, tmax=Inf; xticks=0:2:4, yticks=0:500:1000, i=1)
    lightcolor=palette3[i]
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
    tcond = condfixtime_diffusion_split(N, N0, branchrate_effective_split(r, N, N0, b))
    tcond < tmax && vlines!(ax, tcond, linestyle=:dot, color=:black, linewidth=0.7)
    if update == :Moran
        vlines!(ax, condfixtime_moran(N, b), linestyle=:dash, color=:black, linewidth=0.7)
    end
    ax.xticks = xticks
    ax.yticks = yticks
    ax.title = rich("N = $N, N",subscript("0")," = $N0")
    ax
end

function plot(df, r, μ, gl, update, xticks, yticks, label; tmax=1, b=48, yvisible=true, i=1)
    axes = [Axis(gl[i, j]) for i in 1:2, j in 1:2]
    plot_panel!(axes[1, 1], df, r, 20, 1, μ, update, b, tmax; xticks, yticks, i)
    plot_panel!(axes[1, 2], df, r, 100, 1, μ, update, b, tmax; xticks, yticks, i)
    plot_panel!(axes[2, 1], df, r, 20, 10, μ, update, b, tmax; xticks, yticks, i)
    plot_panel!(axes[2, 2], df, r, 100, 50, μ, update, b, tmax; xticks, yticks, i)

    for i in 1:2, j in 2:2 hideydecorations!(axes[i,j], ticklabels=true, ticks=false) end
    for i in [1,2] hidexdecorations!(axes[1, i], ticklabels=true, ticks=false) end
    colgap!(gl, 12)
    rowgap!(gl, 10)
    Label(gl[1:2, 0], "Fixed SoGV", 
        rotation = pi/2, tellwidth=true, tellheight=false, justification=:left, halign=:left, 
        valign=:center, padding=(3,0,0,0), visible=yvisible)
        Label(gl[0, 1:2], label, valign = :bottom, halign=:left,
        font = :bold, tellwidth=false, justification=:left,
        padding = (0, 0, 0, 0), lineheight=0.5, tellheight=true
    )
    Label(gl[3,1:2], "Time (years)", padding = (0, 0, 0, -5), height=10, valign=:top)

    return axes
end

function fullplot(dfshort, b, μ, savename=nothing, tmax_moran=Inf, tmax_asym=Inf)
    set_theme!(mytheme; resolution = (72 .* (6.2,7.8)), Axis=(;titlehalign=:left, xticklabelspace=6.0))
    fig = Figure()
    gtl = fig[1,1] = GridLayout(halign=:left)
    gtr = fig[1,2] = GridLayout(halign=:left)
    gbl = fig[2,1] = GridLayout(halign=:left)
    gbr = fig[2,2] = GridLayout(halign=:left)

    title1 = rich(rich("a ", fontsize=8), "Symmetric division\n", rich("r = 2.5/year", 
        font=:bold_italic))
    title2 = rich(" \n", rich("r = 20/year", font=:bold_italic))
    title3 = rich(rich("b ", fontsize=8), "Asymmetric division\n", rich("r = 2.5/year", 
        font=:bold_italic))
    title4 = rich(" \n", rich("r = 20/year", font=:bold_italic))

    axes_tl = plot(dfshort, 2.5, μ, gtl, :Moran, 0:1:2, 0:4:8, title1; b, 
        tmax=tmax_moran, i=1)
    axes_tr = plot(dfshort, 20, μ, gtr, :Moran, 0:1:2, 0:4:8, title2; b, 
        tmax=tmax_moran, yvisible=false, i=4)
    axes_bl = plot(dfshort, 2.5, μ, gbl, :asymmetric, 0:2:4, 0:4:8, title3; b, 
        tmax=tmax_asym, i=1)
    axes_br = plot(dfshort, 20, μ, gbr, :asymmetric, 0:2:4, 0:4:8, title4; b, 
        tmax=tmax_asym, yvisible=false, i=4)
    linkaxes!(axes_tl..., axes_tr...)
    linkaxes!(axes_bl..., axes_br...)
    Legend(
        fig[0,:], 
        [LineElement(linestyle=:dash, linewidth=0.7, linecolor=:black),
            LineElement(linestyle=:dot, linewidth=0.7, linecolor=:black)],
        ["Moran conditional fixation time", "MWF-split conditional fixation time"],
        # "Conditional fixation time:",
        titleposition=:left,
        orientation=:horizontal,
        tellheight=:true,
        halign=:right,
        valign=:bottom,
        titlegap=14,
        colgap=16,
        margin=(0,0,-15,0),
        framevisible=false
    )
    isnothing(savename) || save(savename, fig, pt_per_unit=1)


    return fig
end
#endregion

#region Load data
b = 120
r = 40
df_stats_by_sim10yr = get_dataframe("data/fulldist_data/data_tenyrs").df_stats_by_sim
df_stats_by_sim10yr = change_timescale!(df_stats_by_sim10yr, 48, b)
#endregion

#Supplementary figure 6

fig = fullplot(df_stats_by_sim10yr, b, 0.01, "plots/supp_modelling_trajectories.pdf",
    2, 5)

for (N, N0, r) in [(100, 50, 2.5), (100, 50, 20), (20, 10, 2.5)] 
    cft = condfixtime_diffusion_split(N, N0, branchrate_effective_split(r, N, N0, b))
    println("N = $N, N0 = $N0, r = $r: $cft")
end

fig

