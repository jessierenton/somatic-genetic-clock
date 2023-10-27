using DataFrames
using CSV
using Makie
using CairoMakie

include("../fixed_mutation_accumulation/plot_init.jl")

#region Theory
function condfixtime_diffusion_split(N, N0::Int64, r=1)
    return (1/r)*4*N0*(1-N0/N)
end

function condfixtime_diffusion_branch_old(N0::Int64, r=1)
    return 4N0/r
end

function condfixtime_diffusion_branch(N, N0, r=1)
    # return 4*N0 / (1 - (N0/N)^2) / r
    return 4N0 / r / (1-N0/N) / (1 + N0/(N-1))
end

function branchrate_effective_split(r, N, N0, b)
    return 1/(1/r + 1/(2b)*log(N^2/(N0*(N-N0))))
end

function branchrate_effective_branch(r, N, N_0, b)
    return 1/(1/r + 1/(2b)*log(N/N_0))
end
#endregion

#region Axis level plot functions
function plot_panel_asymmetric_splitting!(ax, df, df_theory, r; std=true, N=100)
    newdf = df[
        (df.N .== N) .& (df.r .== r) .& (df.simulation .== "full (asymmetric, split)"), :]
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split_eff; color=Cycled(2))
    std && errorbars!(ax, newdf.N0, newdf.condfixtime_mean, newdf.condfixtime_std, color=:grey)
    scat = scatter!(ax, newdf.N0, newdf.condfixtime_mean; color=:black)
    return [lin, scat]
end

function plot_panel_moran!(ax, df, df_theory, r; std=true, N=100)
    newdf = df[
        (df.N .== N) .& (df.r .== r) .& (df.simulation .== "full (moran, split)"), :]
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin1 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_moran; color=Cycled(1))
    lin2 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split_eff; color=Cycled(2))
    std && errorbars!(ax, newdf.N0, newdf.condfixtime_mean, newdf.condfixtime_std, color=:grey)
    scat = scatter!(ax, newdf.N0, newdf.condfixtime_mean; color=:black)
    return [lin1, lin2, scat]
end

function plot_panel_asymmetric_branching!(ax, df, df_theory, r; std=true, N=100)
    newdf = df[
        (df.N .== N) .& (df.r .== r) .& (df.simulation .== "full (asymmetric, branch)"), :]
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_branch_eff; color=Cycled(3))
    std && errorbars!(ax, newdf.N0, newdf.condfixtime_mean, newdf.condfixtime_std, color=:grey)
    scat = scatter!(ax, newdf.N0, newdf.condfixtime_mean; color=:black)
    return [lin, scat]
end

function plot_theory(ax, df_theory, r; showylabel=true, N=100, tmax_branch=nothing)
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin1 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_moran; color=Cycled(1))
    lin2 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split; color=Cycled(2))
    lin3 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split_eff; color=Cycled(2), linestyle=:dash)
    if !isnothing(tmax_branch)
        newdf_theory = newdf_theory[newdf_theory.condfixtime_branch .< tmax_branch, :]
    end
    lin4 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_branch; color=Cycled(3))
    if !isnothing(tmax_branch)
        newdf_theory = newdf_theory[newdf_theory.condfixtime_branch_eff .< tmax_branch, :]
    end
    lin5 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_branch_eff; color=Cycled(3), linestyle=:dash)
    if showylabel ax.ylabel = "Conditional fixation time (years)" end
    return [lin1, lin2, lin3, lin4, lin5]
end
#endregion

#region Main plot
function full_plot(;std=true)
    rvals_theory = [2.5, 20]
    rvals_sim = [8.0]
    set_theme!(mytheme; 
        resolution = 72 .* (6.2,7.8),
        palette = (;color = colorpalette),
        Legend = (orientation=:vertical, rowgap=0),
        Scatter = (;markersize=5)
    )
    fig = Figure(Axis=(;linkxaxis=true))
    gplottop = fig[1,1] = GridLayout()


    ax1 = Axis(gplottop[1,1])
    plot_theory(ax1, df_theory, rvals_theory[1]; tmax_branch=120)
    ax2 = Axis(gplottop[1,2])
    plot_theory(ax2, df_theory, rvals_theory[2]; showylabel=false, tmax_branch=14)
    Label(gplottop[2, :], rich("Number of founder cells, N",subscript("0")), valign=:top, padding=(0,0,0,-6), tellheight=true)

    # ax1.title = "a. Theoretical results"
    ax1.title = ""
    ax1.subtitle = " r = $(rvals_theory[1])/year"
    ax1.titlealign = :left 
    ax1.yticklabelspace = 12.0

    # ax1.subtitlefont = :regular

    ax2.title = ""
    ax2.subtitle = " r = $(Int64(rvals_theory[2]))/year"
    ax2.titlealign = :left 
    # ax2.subtitlefont = :regular

    # ylims!(ax1, -3.5, 125)
    # ylims!(ax2, -0.4, 15)

    leg_group1 = [[LineElement(linecolor=c)] for c in reverse(colorpalette)]
    leg_group2 = [[LineElement(linecolor=:black, linestyle=ls)] for ls in (:solid, :dash)]
    # leg_labels1 = ["MWF (branching)", "MWF (splitting)", "Moran"]
    leg_labels1 = [
        rich("T", subscript("MWF-branch")),
        rich("T", subscript("MWF-split")), 
        rich("T", subscript("Moran")), 
    ]
    leg_labels2 = [rich("T(r",subscript("eff"),")"), "T(r)"]
    Legend(gplottop[1,3],
        [leg_group1, leg_group2],
        [leg_labels1, leg_labels2],
        ["",""],
        gridshalign=:left, groupgap=0, margin=(15,0,0,0), framevisible=true, tellheight=true,
        padding=(10,10,8,0)

    )
    Label(fig[0, 1:end, TopLeft()], rich(rich("a ", fontsize=8, font=:bold), rich("Theoretical results:  ", font=:bold), "N = 100"), font=:italic, halign=:left, valign=:top, tellwidth=false, tellheight=true)
    Label(fig[2, 1:end, TopLeft()], rich(rich("b ", fontsize=8, font=:bold), rich("Simulation results:  ", font=:bold), "r = 8/year"), font=:italic, halign=:left, valign=:top, tellwidth=false, tellheight=true)


    gplotbottom = fig[3,:] = GridLayout()

    axes1 = [Axis(gplotbottom[1, i], titlealign=:center) for i in 1:3, j in 1:length(rvals_sim)]
    axes2 = [Axis(gplotbottom[2, i], titlealign=:center) for i in 1:3, j in 1:length(rvals_sim)]
    axes3 = [Axis(gplotbottom[3, i], titlealign=:center) for i in 1:3, j in 1:length(rvals_sim)]

    Label(gplotbottom[1,4], "Module branching\n& asymmetric\ndivision", font=:bold, justification=:left, halign=:left, tellheight=false, tellwidth=true)
    Label(gplotbottom[2,4], "Module splitting\n& asymmetric\ndivision", font=:bold, justification=:left, halign=:left, tellheight=false, tellwidth=true)
    Label(gplotbottom[3,4], "Module splitting\n& symmetric\ndivision",  font=:bold, justification=:left, halign=:left, tellheight=false, tellwidth=true)


    for (ax1, ax2, ax3, (N, r)) in zip(axes1, axes2, axes3, Base.product([100, 20, 5], rvals_sim))
        plot_panel_asymmetric_branching!(ax1, df, df_theory, r; std, N)
        plot_panel_asymmetric_splitting!(ax2, df, df_theory, r; std, N)
        plot_panel_moran!(ax3, df, df_theory, r; std, N)
        for ax in (ax1, ax2) hidexdecorations!(ax, ticks=false) end
        ax1.subtitle = " N = $N"
        ax1.titlealign = :left
        for ax in [ax1, ax2, ax3]
            if N == 5
                ax.xticks = 1:4
            else
                ax.xticks = 0:floor(Int64,N/2):N
                xlims!(ax, -N*0.03, N+N*0.03)
            end
        end
    end
    for i in 1:3 linkxaxes!([axes[i] for axes in [axes1, axes2, axes3]]...) end
    Label(gplotbottom[4, 2], rich("Number of founder cells, N",subscript("0")), valign=:top, padding=(0,0,0,-6), tellheight=true, tellwidth=false)
    axes2[1].ylabel = "Conditional fixation time (years)"
    axes2[1].yticklabelspace = 12.0
    ylims!(axes1[1], -20, 250)
    ylims!(axes1[2], -2, 35)
    ylims!(axes1[3], -0.2, 6.2)
    ylims!(axes2[1], -0.6, 14)
    ylims!(axes2[2], -0.1, 2.9)
    ylims!(axes2[3], -0.03, 0.74)
    ylims!(axes3[1], -0.6, 14)
    ylims!(axes3[2], -0.1, 2.9)
    ylims!(axes3[3], -0.03, 0.74)
    axes1[1].yticks = [0,100, 200]
    axes1[2].yticks = [0,15,30]
    axes1[3].yticks = [0,2.5,5]
    axes2[1].yticks = [0,5,10]
    axes2[2].yticks = [0,1,2]
    axes2[3].yticks = [0,0.3, 0.6]
    axes3[1].yticks = [0,5,10]
    axes3[2].yticks = [0,1,2]
    axes3[3].yticks = [0,0.3, 0.6]


    rowsize!(fig.layout, 0, Auto(0.1))
    rowsize!(fig.layout, 1, Auto(1.0))
    rowsize!(fig.layout, 2, Auto(0.1))
    rowsize!(fig.layout, 3, Auto(10.0))

    display(fig)
    save("supp_condfixtimes_new.pdf", fig, pt_per_unit=1)
end
#endregion

#region Load data
df = CSV.read("df_8.csv", DataFrame)


b = 120
N = 100
rvals = [2.5, 8, 20]
df_theory = DataFrame([
    (;
        r, 
        N,
        N0, 
        condfixtime_split=condfixtime_diffusion_split(N, N0, r),
        condfixtime_split_eff=condfixtime_diffusion_split(N, N0, branchrate_effective_split(r, 100, N0, b)),
        condfixtime_branch=condfixtime_diffusion_branch(N, N0, r),
        condfixtime_branch_eff=condfixtime_diffusion_branch(N, N0, branchrate_effective_branch(r, 100, N0, b)),
        condfixtime_moran=(N-1)^2/(N*b)

    )
    for r in rvals for N in [3, 5, 20, 100] for N0 in 1:(N-1)
])
#endregion

colorpalette = [lime, coral, turquoise]

full_plot(std=false)
