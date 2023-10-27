using DataFrames
using CSV
using Makie
using CairoMakie

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

function plot_panel_asymmetric_splitting!(ax, df, df_theory, r; std=true, N=100, color=:black, marker=:circle)
    newdf = df[
        (df.N .== N) .& (df.r .== r) .& (df.simulation .== "full (asymmetric, split)"), :]
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split_eff; color=color)
    std && errorbars!(ax, newdf.N0, newdf.condfixtime_mean, newdf.condfixtime_std, color=:grey)
    scat = scatter!(ax, newdf.N0, newdf.condfixtime_mean; color=:black, marker=marker)
    return [lin, scat]
end

function plot_panel_moran!(ax, df, df_theory, r; std=true, N=100, color=:black, marker=:circle)
    newdf = df[
        (df.N .== N) .& (df.r .== r) .& (df.simulation .== "full (moran, split)"), :]
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin1 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_moran; color=color)
    std && errorbars!(ax, newdf.N0, newdf.condfixtime_mean, newdf.condfixtime_std, color=:grey)
    scat = scatter!(ax, newdf.N0, newdf.condfixtime_mean; color=:black, marker=marker)
    return [lin1, scat]
end

function plot_panel_asymmetric_branching!(ax, df, df_theory, r; std=true, N=100, color=:black, marker=:circle)
    newdf = df[
        (df.N .== N) .& (df.r .== r) .& (df.simulation .== "full (asymmetric, branch)"), :]
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_branch_eff; color=color)
    std && errorbars!(ax, newdf.N0, newdf.condfixtime_mean, newdf.condfixtime_std, color=:grey)
    scat = scatter!(ax, newdf.N0, newdf.condfixtime_mean; color=:black, marker=marker)
    return [lin, scat]
end

function plot_theory(ax, df_theory, r; showylabel=true, N=100)
    newdf_theory = df_theory[(df_theory[:, :r] .== r) .& (df_theory[:, :N] .== N), :]
    lin1 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_moran; color=Cycled(1))
    lin2 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split; color=Cycled(2))
    lin3 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_split_eff; color=Cycled(2), linestyle=:dash)
    lin4 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_branch; color=Cycled(3))
    lin5 = lines!(ax, newdf_theory.N0, newdf_theory.condfixtime_branch_eff; color=Cycled(3), linestyle=:dash)
    if showylabel ax.ylabel = "Conditional fixation\ntime [years]" end
    return [lin1, lin2, lin3, lin4, lin5]
end

function full_plot(;std=true)
 
    set_theme!(mytheme)
    colors = ["#427674","#7ADB9D", "#F792F8", "#F4BA59"]
    fig = Figure(Axis=(;linkxaxis=true))


    # leg_group1 = [[LineElement(linecolor=c)] for c in [green, yellow, blue]]
    # leg_group2 = [[LineElement(linecolor=:black, linestyle=ls)] for ls in (:solid, :dash)]
    # # leg_labels1 = ["MWF (branching)", "MWF (splitting)", "Moran"]
    # leg_labels1 = [
    #     rich("T", subscript("MWF-branch")),
    #     rich("T", subscript("MWF-split")), 
    #     rich("T", subscript("Moran")), 
    # ]
    # leg_labels2 = [rich("T(r",subscript("eff"),")"), "T(r)"]
    # Legend(gplottop[1,3],
    #     [leg_group1, leg_group2],
    #     [leg_labels1, leg_labels2],
    #     ["",""],
    #     gridshalign=:left, groupgap=0, margin=(15,0,0,0), framevisible=true, tellheight=true,
    #     padding=(10,10,8,0)

    # )


    gplotbottom = fig[1,1] = GridLayout()

    axes = [Axis(gplotbottom[1, i]) for i in 1:2]

    r = 8
    for (ax, N) in zip(axes, [100, 5])
        plot_panel_asymmetric_branching!(ax, df, df_theory, r; std, N, color=colors[1], marker=:star4)
        plot_panel_asymmetric_splitting!(ax, df, df_theory, r; std, N, color=colors[2], marker=:circle)
        plot_panel_moran!(ax, df, df_theory, r; std, N, color=colors[3], marker=:diamond)
        ax.title = "N = $N"
        if N == 5
            ax.xticks = 1:4
        else
            ax.xticks = 0:floor(Int64,N/2):N
            xlims!(ax, -N*0.03, N+N*0.03)
        end
    end
    ylims!(axes[1], -10, 125)
    ylims!(axes[2], -0.5, 5)
    axes[2].yticks = [0,2,4]
    Label(gplotbottom[2, 1:2], rich("Number of founder cells, N",subscript("0")), valign=:top, padding=(0,0,0,0), tellheight=true, tellwidth=false)
    axes[1].ylabel = "Time [yrs]"

    # # end

    display(fig)
    save("poster.pdf", fig, pt_per_unit=1)
end


df = CSV.read("df_all3.csv", DataFrame)
df[:, :N] .= 100
df_small = CSV.read("df_all_small_module.csv", DataFrame)
df = vcat(df_small, df)
df = vcat(df_small, df)
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

mytheme = Theme(
    fontsize = 20, 
    rowgap = 10, 
    colgap = 30, 
    resolution = 28.8 .* (22,12),
    Axis = (
        titlefont=:bold,
        titlealign=:center,
        subtitlefont=:italic, 
        titlesize = 18,
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

full_plot(std=false)