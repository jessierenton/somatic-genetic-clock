using Makie
using CairoMakie
using AlgebraOfGraphics

include("load_data.jl")
include("plot_init.jl")

function fullplot_split_branch(df_fit_split, df_fit_branch, μ, savename=nothing; b=120)
    set_theme!(
        mytheme; 
        resolution = (72 .* (6.5,5)),
        Legend = (position=:top, nbanks=2, framevisible=true, groupgap=25, colgap=10, halign=:center,
            patchsize=(12,15), markersize=7,padding=(10,10,5,0), titlevisible=false,
            gridshalign=:left, titlefont=:regular, framecolor=(:black, 0.5), framewidth=0.4
        ),
        Scatter = (;markersize=7,),
        palette = (;marker=[:circle, :rect], color=[coral, turquoise, lime, purple]),
        Lines = (;linewidth=1)

    )
    fig = Figure()
    gl_fig = fig[1,1] = GridLayout()
    gl_left = gl_fig[1,1] = GridLayout()
    gl_right = gl_fig[1,2] = GridLayout()
    leg_elements, leg_labels, leg_titles = [], [], []
    axes = []
    for (df_fit, gl) in ((df_fit_split, gl_left), (df_fit_branch, gl_right))
        sort!(df_fit, :)
        pltdata = AlgebraOfGraphics.data(df_fit[df_fit.μ .== μ, :]) 
        layer1 = pltdata *  mapping(:r => "Module formation rate, r (per year)", :gradient => "Accumulation rate of fixed SoGV (per year)")
        layer1 *= (visual(Lines) + visual(Scatter, alpha=1))
        layers = layer1
        plt = layers * mapping(        
            layout=:N => x -> "N = $x" |> sorter(getlabels(:N, df_fit)),
            marker=:update => (x->"$x" |> sorter(["Symmetric", "Asymmetric"])) => "Division",
            color=:N0 => nonnumeric => "N0", 
            )
        facet = (linkxaxes = :all, linkyaxes = :all)
        grid = draw!(gl, plt; facet, axis=(;xticks=0:20:40, yticks=[120*μ,240*μ], yautolimitmargin=(0.08, 0.12)))
        for g in grid 
            l = lines!(g.axis, 1:maximum(df_fit.r), x->2*μ*b, color=:black, linestyle=:dot, linewidth=1.2, label="2μb") 
            translate!(l, 0, 0, -5)
            l = lines!(g.axis, 1:maximum(df_fit.r), x->1*μ*b, color=:black, linestyle=:dash, linewidth=1, label="μb") 
            translate!(l, 0, 0, -5)
            push!(axes, g.axis)
        end
        leg_elements, leg_labels, leg_titles = AlgebraOfGraphics.compute_legend(grid)
    end
    newleg_elements = [[LineElement(linestyle=:dot, linewidth=1.2)], [LineElement(linestyle=:dash, linewidth=1)]]
    push!(leg_elements, newleg_elements)
    push!(leg_labels, ["2μb", "μb"])
    leg_titles = ["Division:", rich("N", subscript("0"),":"), "Theory:"]
    leg_titles = ["", "", ""]
    leg_labels[1][1] = "Symmetric"
    leg_labels2 = [rich("N", subscript("0"), " = $N") for N in leg_labels[2]]
    leg_labels = [leg_labels[1], leg_labels2, leg_labels[3]]
    leg_labels[1] = ["Symmetric division", "Asymmetric division"]
    leg_labels[1:2] = leg_labels[2:-1:1]
    leg_elements[1:2] = leg_elements[2:-1:1]
    Legend(gl_fig[0, 1:2], leg_elements, leg_labels,leg_titles)
    Label(gl_left[0,1], rich(rich("a ", fontsize=8), "Module splitting"), valign = :bottom, halign=:left,
    font = :bold, tellwidth=false, justification=:left,
    padding = (0, 0, 0, 0), height=25)
    Label(gl_right[0,1], rich(rich("b ", fontsize=8), "Module branching"), valign = :bottom, halign=:left,
    font = :bold, tellwidth=false, justification=:left,
    padding = (0, 0, 0, 0), height=25)
    colgap!(gl_fig, 1, 20)
    linkaxes!(axes...)
    isnothing(savename) || save(savename, fig, pt_per_unit=1)
    return fig
end

function fullplot(df_fit, μ, savename=nothing; b=120)
    sort!(df_fit, :)
    set_theme!(
        mytheme; 
        resolution = (72 .* (4.5,5)),
        Legend = (position=:top, nbanks=2, framevisible=true, groupgap=25, colgap=10, halign=:center,
            patchsize=(12,15), markersize=7,padding=(10,10,5,0), titlevisible=false,
            gridshalign=:left, titlefont=:regular, framecolor=(:black, 0.5), framewidth=0.4
        ),
        Scatter = (;markersize=7,),
        palette = (;marker=[:circle, :rect], color=[coral, turquoise, lime, purple]),
        Lines = (;linewidth=1)

    )
    pltdata = AlgebraOfGraphics.data(df_fit[df_fit.μ .== μ, :]) 
    layer1 = pltdata *  mapping(:r => "Module formation rate, r (per year)", :gradient => "Accumulation rate of fixed SoGV (per year)")
    layer1 *= (visual(Lines) + visual(Scatter, alpha=1))
    # layer2 = pltdata *  mapping(:n => "b/r", :gradient => "Gradient of regression line", :confint)
    # layer2 *= visual(Errorbars)
    layers = layer1
    plt = layers * mapping(        
        layout=:N => x -> "N = $x" |> sorter(getlabels(:N, df_fit)),
        marker=:update => (x->"$x" |> sorter(["Symmetric", "Asymmetric"])) => "Division",
        color=:N0 => nonnumeric => "N0", 
        )
    facet = (linkxaxes = :all, linkyaxes = :all)
    fig = Figure()
    gl = fig[1,1] = GridLayout()
    grid = draw!(gl, plt; facet, axis=(;xticks=0:20:40, yticks=[120*μ,240*μ], yautolimitmargin=(0.08, 0.12)))
    for g in grid 
        l = lines!(g.axis, 1:maximum(df_fit.r), x->2*μ*b, color=:black, linestyle=:dot, linewidth=1.2, label="2μb") 
        translate!(l, 0, 0, -5)
        l = lines!(g.axis, 1:maximum(df_fit.r), x->1*μ*b, color=:black, linestyle=:dash, linewidth=1, label="μb") 
        translate!(l, 0, 0, -5)
    end
    leg_elements, leg_labels, leg_titles = AlgebraOfGraphics.compute_legend(grid)
    newleg_elements = [[LineElement(linestyle=:dot, linewidth=1.2)], [LineElement(linestyle=:dash, linewidth=1)]]
    push!(leg_elements, newleg_elements)
    push!(leg_labels, ["2μb", "μb"])
    leg_titles = ["Division:", rich("N", subscript("0"),":"), "Theory:"]
    leg_titles = ["", "", ""]
    leg_labels[1][1] = "Symmetric"
    leg_labels2 = [rich("N", subscript("0"), " = $N") for N in leg_labels[2]]
    leg_labels = [leg_labels[1], leg_labels2, leg_labels[3]]
    leg_labels[1] = ["Symmetric division", "Asymmetric division"]
    leg_labels[1:2] = leg_labels[2:-1:1]
    leg_elements[1:2] = leg_elements[2:-1:1]
    Legend(fig[0, 1], leg_elements, leg_labels,leg_titles)
    isnothing(savename) || save(savename, fig, pt_per_unit=1)
    return fig
end

#region Load data
df_long_split = get_dataframe("data/fulldist_data/data_long").df_stats_by_sim
df_long_split = change_timescale!(df_long_split, 48, 120)
df_long_branch = get_dataframe("data/fulldist_data/data_long_branching").df_stats_by_sim
df_long_branch = change_timescale!(df_long_branch, 48, 120)
df_fit_split = get_grad_fits(df_long_split, "data/averaged_data/linear_fit_fulldist_long_split.csv")
df_fit_branch = get_grad_fits(df_long_branch, "data/averaged_data/linear_fit_fulldist_long_branch.csv")

#endregion

#region Make plot (Figure S3)
fig5 = fullplot_split_branch(df_fit_split, df_fit_branch, 1.0, "plots/supp_modelling_fits_mu1.pdf"; b=120)
#endregion

#region Other plots
# fig1 = fullplot(df_fit_split, 1.0, "plots/supp_modelling_fits_mu1_split.pdf")
# fig2 = fullplot(df_fit_split, 0.01, "plots/supp_modelling_fits_mu0.01_split.pdf")
# fig3 = fullplot(df_fit_branch, 1.0, "plots/supp_modelling_fits_mu1_branch.pdf")
# fig4 = fullplot(df_fit_branch, 0.01, "plots/supp_modelling_fits_mu0.01_branch.pdf")
#endregion