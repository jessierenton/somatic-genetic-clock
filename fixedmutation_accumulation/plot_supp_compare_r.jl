using Makie
using CairoMakie

include("load_data.jl")
include("plot_init.jl")

#region Plot functions
function plot_panel!(ax, df, N, N0, μ, update, b, tmax=Inf; rvals=nothing, yticks=0:3000:6000, ylims=(nothing, 7000))
    newdf = df[(
        (df.N .== N) .&
        (df.N0 .== N0) .&
        (df.μ .≈ μ) .&
        (df.update .== update) .&
        (df.t .< tmax)
    ), :]
    isnothing(rvals) || filter!(:r => r -> r in rvals, newdf)
    newdf_av = combine(
        groupby(newdf, [:r, :t]),
        :fixedmutations_mean => mean
    )
    sort!(newdf_av, :r)
    lineplts = [
        lines!(ax, gdf.t, gdf.fixedmutations_mean_mean, linewidth=0.6, label="$(round(Int64, gdf[1, :r]))")
            for gdf in groupby(newdf_av, :r)
    ]
    if update == :asymmetric
        lines!(ax, 0:maximum(newdf.t), t->μ*b*t, linestyle=:dash, color=:black, linewidth=1, label="μb")
    elseif update == :Moran
        lines!(ax, 0:maximum(newdf.t), t->2*μ*b*t, linestyle=:dot, color=:black, linewidth=1.2, label="2μb")
    end
    ax.ylabel = ""
    ax.xticks = 0:100:200
    xlims!(ax, nothing, 205)
    ylims!(ax, ylims...)
    ax.yticks = yticks

    ax.title = rich("N = $N, N",subscript("0")," = $N0")
    ax
end


function plot(df, gl, update, yticks, ylims, label; tmax=192, b=48, μ=1)
    axes = [Axis(gl[i, j]) for i in 1:2, j in 1:3]
    plot_panel!(axes[1, 1], df, 5, 1, μ, update, b, tmax; yticks, ylims)
    plot_panel!(axes[2, 1], df, 5, 2, μ, update, b, tmax; yticks, ylims)
    plot_panel!(axes[1, 2], df, 20, 1, μ, update, b, tmax; yticks, ylims)
    plot_panel!(axes[2, 2], df, 20, 10, μ, update, b, tmax; yticks, ylims)
    plot_panel!(axes[1, 3], df, 100, 1, μ, update, b, tmax; yticks, ylims)
    plot_panel!(axes[2, 3], df, 100, 50, μ, update, b, tmax; yticks, ylims)

    for i in 1:2, j in 2:3 hideydecorations!(axes[i,j], ticklabels=true, ticks=false) end
    for i in [1,2,3] hidexdecorations!(axes[1, i], ticklabels=true, ticks=false) end
    colgap!(gl, 12)
    rowgap!(gl, 10)
    Label(gl[1:2, 0], "Mean fixed SoGV", 
        rotation = pi/2, tellwidth=false, tellheight=false, justification=:left, halign=:right, valign=:center, padding=(2,0,0,0))
    Label(gl[0, 0:3], label, valign = :bottom, halign=:left,
        font = :bold, tellwidth=false, justification=:left,
        padding = (0, 0, 0, 0), height=25
    )
    axes[2,2].xlabel = "Time (years)"

    return axes
end

function fullplot(df_stats_by_sim, b, μ, ylims1, ylims2, yticks1, yticks2, savename=nothing)
    set_theme!(mytheme; resolution = (72 .* (6.2,7.4)), Axis=(;width=72*1.12), palette=(;color=palette3))
    fig = Figure()
    gltop = fig[1,1] = GridLayout(halign=:left)
    glbottom = fig[2,1] = GridLayout(halign=:left)
    plot(df_stats_by_sim, gltop, :Moran, yticks1, ylims1, rich(rich("a ", fontsize=8), "Symmetric division"); b, μ)
    plot(df_stats_by_sim, glbottom, :asymmetric, yticks2, ylims2, rich(rich("b ", fontsize=8), "Asymmetric division"); b, μ)

    leg = Legend(
        fig[1:2, 2], 
        [[LineElement(linecolor=palette3[i], linewidth=1) for i in 5:-1:1],
        [LineElement(linecolor=:black, linestyle=:dot, linewidth=1.2),
        LineElement(linecolor=:black, linestyle=:dash, linewidth=1)]], 
        [["r = $(r%1 == 0 ? Int64(r) : r)" for r in sort!(unique(df_stats_by_sim.r), rev=true)],
        ["2μbt", "μbt"]],
        ["Simulation", "Theory"], 
        tellwidth=true, tellheight=false, framevisible=true, margin=(10,0,0,0), rowgap=-2,
        orientation=:vertical, halign=:left
        
    )
    isnothing(savename) || save(savename, fig, pt_per_unit=1)
    return fig
end
#endregion

#region Load data
b = 120
df_long = get_dataframe("data/fulldist_data/data_long").df_stats_by_sim
df_long = change_timescale!(df_long, 48, b)
#endregion

#region Make plots (Figure S2 and S3)
μ = 0.01
fig1 = fullplot(df_long, b, μ, 
    (nothing, 480), (nothing, 320), 
    0:200:400, 0:150:300, 
    "plots/supp_modelling_2panel_mu0.01.pdf"
)

μ = 1
fig2 = fullplot(df_long, b, μ, 
    (nothing, 48000), (nothing, 32000), 
    0:20000:40000, 0:15000:30000, 
    "plots/supp_modelling_2panel_mu1.pdf"
)
#endregion