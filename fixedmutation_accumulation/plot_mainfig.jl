using Makie
using CairoMakie

include("load_data.jl")
include("plot_init.jl")

#region Axis level functions
function plot_panel!(ax, df, N, N0, μ, r, update, tmax, darkcolor=:black, lightcolor=:darkgrey; 
    extra_r=nothing, show_individual=false, expected_grad=nothing, show_xlabel=true, 
    show_ylabel=true, yticks=nothing, kwargs...)

    rvals = isnothing(extra_r) ? [r] : pushfirst!(extra_r, r)
    newdf = df[(
        (df.N .== N) .&
        (df.N0 .== N0) .&
        (df.μ .≈ μ) .&
        in.(df.r, (rvals, )) .&
        (df.update .== update) .&
        (df.t .< tmax)
    ), :]
    newdf1 = newdf[newdf.r .≈ r, :]
    show_individual && for gdf in groupby(newdf1, :simid)
        lines!(ax, gdf.t, gdf.fixedmutations_mean, linewidth=thin, color=(darkcolor, 0.3))
    end
    newdf_av = combine(
        groupby(newdf, [:t, :r]),
        :fixedmutations_mean => mean
    )
    isnothing(expected_grad) || lines!(ax, sort!(unique(newdf1.t)), t->t*expected_grad; color=:black, linestyle=:dash, linewidth=thin*2)
    mainlines = []
    for (i, gdf) in enumerate(groupby(newdf_av, :r))
        l = lines!(ax, gdf.t, gdf.fixedmutations_mean_mean, color=isnothing(extra_r) ? darkcolor : Cycled(i), linewidth=thick)
        push!(mainlines, l)
    end
    if show_xlabel ax.xlabel = "Time (years)" end
    if show_ylabel ax.ylabel = "Fixed SoGV" end
    if !isnothing(yticks) ax.yticks = yticks end
    return mainlines
end

function plot_panel_grad_singlept!(ax, df, N, N0, μ, r, update, tmax; expected_grad=nothing, 
    extra_r=nothing, show_individual=false, show_xlabel=true, show_ylabel=true, yticks=nothing, df_long_fit=nothing)
    rvals = isnothing(extra_r) ? [r] : pushfirst!(extra_r, r)
    newdf = df[(
        (df.N .== N) .&
        (df.N0 .== N0) .&
        (df.μ .≈ μ) .&
        in.(df.r, (rvals, )) .&
        (df.update .== update) .&
        (df.t .< tmax)
    ), :]
    newdf[:, :grad_est] = newdf.fixedmutations_mean ./ newdf.t
    if !isnothing(df_long_fit)
        newdf_long_fit = subset(df_long_fit, :N => (x-> x .== N), :r => (x -> in.(x, (rvals, ))), :update => (x -> x.== update))
    end
    newdf_av = combine(
        groupby(newdf, [:t, :r]),
        :grad_est => mean
    )
    show_individual && scatter!(ax, newdf.t, newdf.grad_est, color=:darkgrey, markersize=1.5)
    mainlines = []
    for (i, (key, gdf)) in enumerate(pairs(groupby(newdf_av, :r)))
        r = key.r
        l = lines!(ax, gdf.t, gdf.grad_est_mean, color=isnothing(extra_r) ? :black : Cycled(i), linewidth=thick)
        l_true = if !isnothing(df_long_fit)
            true_grad = newdf_long_fit[findfirst(==(3.0), newdf_long_fit.r), :gradient]
            lines!(ax, gdf.t, x->true_grad, color=Cycled(i), linewidth=0.2)
        end
        push!(mainlines, l)
    end
    lines!(ax, [0, maximum(newdf_av.t)], t->μ*122, color=:black, linestyle=:dash, linewidth=thin*2)
    isnothing(expected_grad) || lines!(ax, [0, maximum(newdf_av.t)], t->expected_grad, color=:black, linestyle=:dash, linewidth=0.8)


    if show_xlabel ax.xlabel = "Time (years)" end
    if show_ylabel ax.ylabel = "Fixed SoGV\nper year" end
    if !isnothing(yticks) ax.yticks = yticks end
    return mainlines
end


function general_panel!(gl, df)
    params = [(20, 1), (20, 10), (100, 2), (100, 50)]
    n = length(params)
    axes = [Axis(gl[1, i]) for i in 1:n]
    linkaxes!(axes...)
    lns = []
    for (ax, (N, N0)) in zip(axes, params)
        push!(lns, plot_panel!(ax, df, N, N0, 0.01, 5, :Moran, 20, palette1[1]; show_individual=true))
        push!(lns, plot_panel!(ax, df, N, N0, 0.01, 5, :asymmetric, 20, palette1[2]; show_individual=true))
        ax.title = rich("N = $N, N",subscript("0")," = $N0")
    end

    for ax in axes
        ylims!(ax, nothing, 55)
        ax.yticks = 0:25:50
        ax.xticks = 0:10:20
    end
    for i in 2:n hideydecorations!(axes[i], grid=false, ticks=false) end

    colgap!(gl, 10)
    return gl, lns
end

function plot_experiment_grid!(gl, df_fit_stats, df_long_fit, μ, b, params)

    axes = [Axis(gl[j,i]) for i in 1:2 for j in 1:2]
    linkaxes!(axes...)
    for ax in axes[1:2:3] hidexdecorations!(ax, ticks=false) end
    for ax in axes[3:4] hideydecorations!(ax, ticks=false) end
    updates = [:Moran, :asymmetric, :Moran, :asymmetric]
    for (ax, (N, N0), update) in zip(axes, params, updates)
        condition1 = :N => (x -> x .== N)
        condition2 = :N0 => (x -> x .== N0)
        condition3 = :update => (x -> x.== update)
        sdf = subset(df_fit_stats, condition1, condition2, condition3)
        sdf_long = subset(df_long_fit, condition1, condition2, condition3)
        xs = fill(3, length(sdf.r))
        xs[sdf.r .≈ 3.0] .= 1
        xs[sdf.r .≈ 5.0] .= 2
        barplot!(ax, xs, sdf.gradient_mean; color=palette2[xs], fillto=0.6)
        hlines!(ax, [b*μ], color=:black, linestyle=:dot, linewidth=0.7)
        errorbars!(ax, xs, sdf.gradient_mean, sdf.gradient_std, color=:black)
        scatter!(ax, xs.+0.22, sdf_long.gradient, color=:black, markersize=14, marker='*', strokewidth=0.2)
        ax.xticks = ([1,2,3], ["3", "5", "8"])
        ax.yticks = 0:0.2:1.0
        if update == :Moran 
            ax.title = rich("N = $(N), N",subscript("0")," = $(N0)")
        end            
    end

    Label(gl[1:2, 1, Left()], "Accumulation rate of fixed SoGV (per year)", 
        rotation = pi/2, tellwidth=true, tellheight=false, justification=:left, 
        halign=:left, valign=:center, padding=(0,30,0,0))
    Label(gl[2, 1:2, Bottom()], "Module formation rate, r (per year)", 
        tellwidth=true, tellheight=false, justification=:center, 
        halign=:center, valign=:center, padding=(0,0,0, 30))
    for (j, label) in enumerate(["Symmetric division", "Asymmetric division"])
        Label(gl[j, 2, Right()], label, rotation=-pi/2, tellheight=false, 
            padding=(4,0,0,0))
    end
    return gl
end

function plot_eelgrass_grid!(gl, df, μ, b, params; plotfunc=plot_panel!, ylabel="Fixed SoGV", yticks=nothing, ylabpad=25, df_long_fit=nothing)
    tmax = 20
    extra_r = [3,8]
    n = length(params)
    axes = [Axis(gl[i, j], palette=(; color=palette2)) for i in 1:2, j in 1:length(params)]
    linkaxes!(axes...)
    local lns
    for (i, (N, N0, r)) in enumerate(params)
        lns = plotfunc(axes[1, i], df, N, N0, μ, r, :Moran, tmax; extra_r, expected_grad=μ*b, show_xlabel=false, show_ylabel=false, yticks, df_long_fit)
        lns = plotfunc(axes[2, i], df, N, N0, μ, r, :asymmetric, tmax; extra_r, expected_grad=μ*b, show_xlabel=false, show_ylabel=false, yticks, df_long_fit)
        axes[1, i].title = rich("N = $(N), N",subscript("0")," = $(N0)")
        if i != 1
            for j in 1:2 
                hideydecorations!(axes[j, i], grid=false, ticks=false)
            end
        end
        hidexdecorations!(axes[1, i], grid=false, ticks=false)
    end
    Label(gl[1:2, 1, Left()], ylabel, 
        rotation = pi/2, tellwidth=true, tellheight=false, justification=:left, 
        halign=:left, valign=:center, padding=(0,ylabpad,0,0))
    Label(gl[2, 1:2, Bottom()], "Time (years)", 
        tellwidth=true, tellheight=false, justification=:center, 
        halign=:center, valign=:center, padding=(0,0,0, 30))

    for (j, label) in enumerate(["Symmetric division", "Asymmetric division"])
        Label(gl[j, 2, Right()], label, rotation=-pi/2, tellheight=false, 
            padding=(4,0,0,0))
    end
    return gl, lns
end

function eelgrass_panel_symmetric_asymmetric!(gl, df, df_exp_fit_stats, df_long_fit, b, params, params_exp)
    gl_left = gl[1,1] = GridLayout()
    gl_right = gl[1,2] = GridLayout()
    gl_right = plot_experiment_grid!(gl_right, df_exp_fit_stats, df_long_fit, μ, b, params_exp)
    gl_left, lns = plot_eelgrass_grid!(gl_left, df, μ, b, params; yticks=0:5:15)
    return gl, lns
end
#endregion

#region Figure level functions
function eelgrass_supp_plot()
    fig = Figure(resolution=(72 .* (3.5,3.5)))
    gl = fig[1,1] = GridLayout()
    params = [(7, 3, 5), (12, 4, 5)]
    gl, eel_lns = plot_eelgrass_grid!(gl, df2, μ, b_eelgrass, params; plotfunc=plot_panel_grad_singlept!, 
        ylabel="Accumulation rate of fixed SoGV (per year)", yticks=0:0.4:1.2, ylabpad=30)
    Legend(fig[0, 1], eel_lns, ["r = 3/year", "r = 5/year", "r = 8/year"],
        titlegap=14, patchlabelgap=6, halign=:right, valign=:bottom, tellwidth=false,
        tellheight=true, padding=(0,0,0,0)
    )
    return fig
end

function eelgrass_extremes_supp_plot()
    (N_1, N0_1), (N_2, N0_2) = (7, 1), (12, 6) 
    params = [(N_1, N0_1, 5), (N_2, N0_2, 5)]
    params_exp = [(N_1, N0_1), (N_1, N0_1), (N_2, N0_2), (N_2, N0_2)]

    fig = Figure(resolution=(72 .* (7, 4)))
    gl = fig[1,1] = GridLayout()
    gl, eel_lns = eelgrass_panel_symmetric_asymmetric!(gl, df2, df_fit_stats, df_long_fit, b_eelgrass, params, params_exp)

    Legend(gl[0, end], eel_lns, ["r = 3/year", "r = 5/year", "r = 8/year"],
        orientation = :horizontal, titleposition=:left, framevisible=false,
        titlegap=14, patchlabelgap=6, halign=:right, valign=:top, margin=(0,0,0,0),
        tellwidth=false, padding=(0,0,15,0)
    )
    lineheight = 2.6

    Label(gl[0, 1, Left()], rich("Eelgrass\n",rich("a", fontsize=8)), font=:bold, valign=:top, halign=:left, 
        justification=:left, tellwidth=false, tellheight=false, lineheight=lineheight,
        padding=(0,0,0,0))
    Label(gl[0, 2, Left()], rich(" \n", rich("b", fontsize=8)), font=:bold, valign=:top, halign=:left, 
        justification=:left, lineheight=lineheight, padding=(0,0,0,0), fontsize=8)

    colgap!(gl, 15)
    colgap!(gl, 1, 20)
    rowgap!(gl, 1, 0)
    rowsize!(fig.layout, 1, Auto(0.6))
    return fig

end

function main_plot()
    (N_1, N0_1), (N_2, N0_2) = (7, 3), (12, 4) 
    params = [(N_1, N0_1, 5), (N_2, N0_2, 5)]
    params_exp = [(N_1, N0_1), (N_1, N0_1), (N_2, N0_2), (N_2, N0_2)]


    fig = Figure(resolution=(72 .* (7, 6)))
    gl1 = fig[1,1] = GridLayout()
    gl2 = fig[2,1] = GridLayout()
    gl1, lns = general_panel!(gl1, df1)
    gl2, eel_lns = eelgrass_panel_symmetric_asymmetric!(gl2, df2, df_fit_stats, df_long_fit, b_eelgrass, params, params_exp)
    Legend(gl1[0, end], lns[1:2], ["Symmetric division", "Asymmetric division"],
        orientation = :horizontal, titleposition=:left, framevisible=false,
        titlegap=14, patchlabelgap=6, halign=:right, valign=:top, margin=(0,0,0,0),
        tellwidth=false, padding=(0,0,0,0)
    )
    Legend(gl2[0, end], eel_lns, ["r = 3/year", "r = 5/year", "r = 8/year"],
        orientation = :horizontal, titleposition=:left, framevisible=false,
        titlegap=14, patchlabelgap=6, halign=:right, valign=:top, margin=(0,0,0,0),
        tellwidth=false, padding=(0,0,15,0)
    )
    lineheight = 2.6
    Label(gl1[0, 1, Left()], rich("Generic clonal species\n", rich("a", fontsize=8)), font=:bold, valign=:center, halign=:left, 
        justification=:left, tellwidth=false, tellheight=false, lineheight=lineheight,
        padding=(0,0,0,0))
    Label(gl2[0, 1, Left()], rich("Eelgrass\n",rich("b", fontsize=8)), font=:bold, valign=:top, halign=:left, 
        justification=:left, tellwidth=false, tellheight=false, lineheight=lineheight,
        padding=(0,0,0,0))
    Label(gl2[0, 2, Left()], rich(" \n", rich("c", fontsize=8)), font=:bold, valign=:top, halign=:left, 
        justification=:left, lineheight=lineheight, padding=(0,0,0,0), fontsize=8)

    colgap!(gl1, 10)
    colgap!(gl2, 15)
    colgap!(gl2, 1, 20)
    rowgap!(gl2, 1, 0)
    rowsize!(fig.layout, 1, Auto(0.6))
    return fig
end
#endregion

#region Setup data
const b_general=120
const b_eelgrass=122
μ=0.0069
df1 = get_dataframe("data/fulldist_data/data_branching").df_stats_by_sim
df1 = change_timescale!(df1, 48, b_general)
df2 = get_stats_dataframe("data/averaged_data/data_branch_eelgrass_new")
df_samples = get_dataframe_samples("data/experimental_sampling_data/data")
df_stats = get_fixedmutation_stats(df_samples)
df_fit = get_gradient_estimates(df_stats)
df_fit_stats = get_gradient_estimates_stats(df_fit)
df_long = get_stats_dataframe(["data/averaged_data/data_branch_eelgrass_long/$i/" for i in 0:4])
df_long_fit = combine(groupby(df_long, [:N, :N0, :μ, :r, :update])) do gdf
    model = lm(@formula(fixedmutations_mean ~ t), gdf)
    return DataFrame(Dict(
        :gradient => coef(model)[2],
        :stderror_gradient => stderror(model)[2]
    ))
end
for df0 in (df_stats, df_fit, df_fit_stats, df_long_fit)
    sort!(df0, [:N, :N0, :r, :update], rev=[false, false, false, true])
end

set_theme!(mytheme, Axis=(xlabelpadding=3,))
#endregion

#region Make main plot (Figure 2)
fig = main_plot()
save("plots/main_modelling_v2.pdf", fig, pt_per_unit=1)
#endregion

#region Make extended data plot (Supplementary Figure 12)
fig2 = eelgrass_supp_plot()
save("plots/supp_modelling_eelgrass_grad_estimates.pdf", fig2, pt_per_unit=1)
#endregion

#region Make extended data plot (Supplementary Figure 11)
fig3 = eelgrass_extremes_supp_plot()
save("plots/supp_modelling_eelgrass_extremes.pdf", fig3, pt_per_unit=1)
#endregion