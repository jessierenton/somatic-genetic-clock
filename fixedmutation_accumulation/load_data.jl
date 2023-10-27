using CSV
using DataFrames
using StatsBase
using Statistics
using GLM
using JSON

function loaddata(filename)
    data = Dict()
    open(filename) do io
        data = JSON.parse(io, dicttype=Dict{Symbol, Any})
    end
    return data
end

loaddata(filename, dir) = loaddata("$dir/$filename")

function loaddata(filename, dirs::Vector)
    data_list = [loaddata(filename, dir) for dir in dirs]
    data = data_list[1]
    data[:seed] = Int64[data[:seed]]
    data[:fixedmutations] = reduce(vcat, [d[:fixedmutations] for d in data_list])
    return data
end

function get_module_timeseries(data, timestep, tmax=Inf)
    maxmoduleid = if !isinf(tmax)
        data[floor(Int64, tmax / timestep)][1][end]
    else
        data[end][1][end]
    end
    module_tvecs = Vector{Float64}[Float64[] for i in 1:maxmoduleid]
    module_mutvecs = Vector{Int64}[Int64[] for i in 1:maxmoduleid]
    for (i, vec) in enumerate(data)
        t = (i-1)*timestep
        t <= tmax || break
        for (moduleid, fixedmuts) in zip(vec[1], vec[2])
            push!(module_tvecs[moduleid], t)
            push!(module_mutvecs[moduleid], fixedmuts)
        end
    end
    return module_tvecs, module_mutvecs
end

getlabels(s::Symbol, df; suffix="", prefix=nothing) = (isnothing(prefix) ? string(s) : prefix) * " = " .* string.(sort(unique(df[:, s]))) .* suffix

function get_dataframe(dir)
    Nvals = Int64[]
    N0vals = Int64[]
    μvals = Float64[]
    rvals = Float64[]
    updatevals = Symbol[]
    simidvals = Int64[]
    tvals = Float64[]
    moduleidvals = Int64[]
    fixedmutationvals = Int64[]

    for filename in readdir(dir)
        first(filename) === '.' && continue
        data = loaddata("$dir/$filename")
        N = data[:input][:modulesize]
        N0 = data[:input][:branchinitsize]
        μ = data[:input][:μ]
        r = data[:input][:branchrate]
        update = data[:input][:moranrate] > 0 ? :Moran : :asymmetric
        tstep = data[:timestep]
        for (simid, simresult) in enumerate(data[:fixedmutations])
            for (t, (moduleids, fixedmuts)) in zip(0:tstep:length(simresult)*tstep, simresult)
                for (moduleid, fixedmut) in zip(moduleids, fixedmuts)
                    push!(Nvals, N)
                    push!(N0vals, N0)
                    push!(μvals, μ)
                    push!(rvals, r)
                    push!(updatevals, update)
                    push!(simidvals, simid)
                    push!(tvals, t)
                    push!(moduleidvals, moduleid)
                    push!(fixedmutationvals, fixedmut)
                end
           end
        end
    end
    df = DataFrame(
        N = Nvals,
        N0 = N0vals,
        μ = μvals,
        r = rvals,
        update = updatevals,
        simid = simidvals,
        t = tvals,
        moduleid = moduleidvals,
        fixedmutations = fixedmutationvals
    )
    df_stats_by_sim = combine(
        groupby(df, [:N, :N0, :μ, :r, :update, :simid, :t]), 
        :fixedmutations => mean,
        :fixedmutations => std
    )
    df_stats_by_sim[:, :pred] = ((df_stats_by_sim.update .== :Moran) .+ 1) .* df_stats_by_sim.μ * 48 .* df_stats_by_sim.t;

    df_stats = combine(
        groupby(df_stats_by_sim, [:N, :N0, :μ, :r, :update, :t]), 
        :fixedmutations_mean => mean,
        :fixedmutations_mean => std,
        :fixedmutations_std => mean
    )
    df_stats[:, :pred] = ((df_stats.update .== :Moran) .+ 1) .* df_stats.μ * 48 .* df_stats.t;


    return (;df, df_stats_by_sim, df_stats)
end

get_filenames(dir) = filter(x-> x[1] != '.', readdir(dir))
get_filenames(dir::Vector) = get_filenames(dir[1])


function get_stats_dataframe(dir)
    Nvals = Int64[]
    N0vals = Int64[]
    μvals = Float64[]
    rvals = Float64[]
    updatevals = Symbol[]
    simidvals = Int64[]
    tvals = Float64[]
    fixedmutation_means = Float64[]
    fixedmutation_vars = Float64[]

    for filename in get_filenames(dir)
        data = loaddata(filename, dir)
        N = data[:input][:modulesize]
        N0 = data[:input][:branchinitsize]
        μ = data[:input][:μ]
        r = data[:input][:branchrate]
        update = data[:input][:moranrate] > 0 ? :Moran : :asymmetric
        tstep = data[:timestep]
        for (simid, simresult) in enumerate(data[:fixedmutations])
            for (t, (fixedmuts_mean, fixedmuts_var)) in zip(0:tstep:(length(simresult)-1)*tstep, simresult)
                push!(Nvals, N)
                push!(N0vals, N0)
                push!(μvals, μ)
                push!(rvals, r)
                push!(updatevals, update)
                push!(simidvals, simid)
                push!(tvals, t)
                push!(fixedmutation_means, fixedmuts_mean)
                push!(fixedmutation_vars, isnothing(fixedmuts_var) ? 0 : fixedmuts_var)
           end
        end
    end
    df = DataFrame(
        N = Nvals,
        N0 = N0vals,
        μ = μvals,
        r = rvals,
        update = updatevals,
        simid = simidvals,
        t = tvals,
        fixedmutations_mean = fixedmutation_means,
        fixedmutations_var = fixedmutation_vars
    )
    return df
end


function change_timescale!(df, oldbirthrate, newbirthrate)
    df.t = df.t .* (oldbirthrate / newbirthrate)
    df.r = df.r .* (newbirthrate / oldbirthrate)
    df
end
function add_dimensionless!(df, df_stats, df_stats_by_sim; b=48)
    df[:, :τ] = Float64.(df.t * b)
    df[:, :n] = Int64.(b ./ df.r)

    df_stats[:, :τ] = Float64.(df_stats.t * b)
    df_stats[:, :n] = Int64.(b ./ df_stats.r)

    df_stats_by_sim[:, :τ] = Float64.(df_stats_by_sim.t * b)
    df_stats_by_sim[:, :n] = Int64.(b ./ df_stats_by_sim.r)
    return 
end


function get_dataframe_samples(dir)
    Nvals = Int64[]
    N0vals = Int64[]
    μvals = Float64[]
    rvals = Float64[]
    updatevals = Symbol[]
    simidvals = Int64[]
    cloneids = Int64[]
    fixedmuts = Int64[]
    tvals = Float64[]

    for filename in get_filenames(dir)
        data = loaddata("$dir/$filename")
        simid = parse(Int64,split(filename, "_")[end][1:end-5])
        N = data[:input][:modulesize]
        N0 = data[:input][:branchinitsize]
        μ = data[:input][:μ]
        r = data[:input][:branchrate]
        update = data[:input][:moranrate] > 0 ? :Moran : :asymmetric
        for (time, fixedmut_clones) in zip(data[:times], data[:fixedmutations])
            for (cloneid, fixedmut_modules) in enumerate(fixedmut_clones)
                for fixedmut in fixedmut_modules
                    push!(Nvals, N)
                    push!(N0vals, N0)
                    push!(μvals, μ)
                    push!(rvals, r)
                    push!(updatevals, update)
                    push!(simidvals, simid)
                    push!(cloneids, cloneid)
                    push!(tvals, time)
                    push!(fixedmuts, fixedmut)
                end
           end
        end
    end
    df = DataFrame(
        N = Nvals,
        N0 = N0vals,
        μ = μvals,
        r = rvals,
        update = updatevals,
        simid = simidvals,
        t = tvals,
        cloneid = cloneids,
        fixedmutations = fixedmuts
    )
    sort!(df, [:N, :N0, :μ, :r, :update, :simid, :t])
    return df
end

function get_fixedmutation_stats(df)
    grouped = groupby(df, [:N, :N0, :μ, :r, :update, :simid, :cloneid, :t])
    return combine(grouped, :fixedmutations => mean, :fixedmutations => std)
end

function get_gradient_estimates(df_stats)
    grouped = groupby(df_stats, [:N, :N0, :μ, :r, :update, :simid])
    df_fit = combine(grouped) do gdf
        model = lm(@formula(fixedmutations_mean ~ t), gdf)
        return DataFrame(Dict(
            :gradient => coef(model)[2],
            :intercept => coef(model)[1],
            :stderror_gradient => stderror(model)[2],
            :stderror_intercept => stderror(model)[1]
        ))
    end
    return df_fit
end

function get_gradient_estimates_stats(df_fit)
    grouped = groupby(df_fit, [:N, :N0, :μ, :r, :update])
    return combine(
        grouped, 
        :gradient => mean, 
        :gradient => std, 
        :intercept => mean,
        :intercept => std

    )
end

function get_grad_fits(df, filename=nothing)
    df_fit = combine(groupby(df, [:N, :N0, :μ, :r, :update])) do gdf
        model = lm(@formula(fixedmutations_mean ~ t), gdf)
        return DataFrame(Dict(
            :gradient => coef(model)[2],
            :confint => (confint(model)[2] - confint(model)[1])/2,
        ))
    end
    df_fit = sort!(df_fit, [:N, :N0, :μ, :r])
    isnothing(filename) || CSV.write(filename, df_fit)
    df_fit[df_fit.update .== :Moran , :update] .= :Symmetric
    df_fit[df_fit.update .== :asymmetric , :update] .= :Asymmetric
    df_fit.update = String.(df_fit.update)
    return df_fit
end
