using SomaticEvolution
using Random
using JSON
using StatsBase

function fixed_mutations_over_time(::Type{T}, input, time_clones_samples, seed, outputfile) where T<:SomaticEvolution.AbstractCell
    rng = Random.seed!(seed)
    clonalmuts_times = Vector{Vector{Int64}}[]
    for (time, clone_samples) in time_clones_samples
        input_t = newinput(input; tmax = time)
        clonalmuts_clones = Vector{Int64}[]
        for nsamples in clone_samples
            population = runsimulation(
                T,
                input_t,
                rng
            )
            sampled_modules = 
                if length(population.output) > nsamples
                    sample(rng, population.output, nsamples; replace=false)
                else
                    x
                end
                clonal_mutations.(sampled_modules)
            push!(clonalmuts_clones, clonal_mutations.(sampled_modules))
        end
        push!(clonalmuts_times, clonalmuts_clones)
    end
    inputdict = SomaticEvolution.inputdict(input)
    open("$outputfile.json", "w") do io
        JSON.print(
            io, 
            Dict(
                :input => inputdict, 
                :seed => seed, 
                :times => [ts[1] for ts in time_clones_samples], 
                :fixedmutations => clonalmuts_times
            ), 
            4
        )
    end
end

function main()
    birthrate = parse(Float64, ARGS[1])
    moranrate = parse(Float64, ARGS[2])
    asymmetricrate = parse(Float64, ARGS[3])
    modulesize = parse(Int64, ARGS[4])
    branchinitsize = parse(Int64, ARGS[5])
    μ = parse(Float64, ARGS[6])
    branchrate = parse(Float64, ARGS[7])
    modulebranching = Symbol(ARGS[8])
    maxmodules = parse(Int64, ARGS[9])
    id = parse(Int64, ARGS[10])
    tmax =  parse(Float64, ARGS[11])
    jobid = parse(Int64, ARGS[12])
    outdir = ARGS[13]
    seed = jobid + id

    input = MultilevelBranchingMoranInput(;
        maxmodules,
        birthrate,
        deathrate=0,
        moranrate,
        asymmetricrate,
        mutationdist=:poisson,
        μ,
        branchrate,
        branchinitsize,
        modulesize,
        tmax=tmax,
        modulebranching,
        moranincludeself=false
    )
    mkpath(outdir)
    outdir *= "/N$(modulesize)_Ninit$(branchinitsize)"
    outdir *= "_μ$(μ)_r$(branchrate)"
    outdir *= "_moranrate$(moranrate)_asymrate$(asymmetricrate)_$id"

    time_clones_samples = [[4, [2, 2, 2]], [17, [5, 6]]]

    fixed_mutations_over_time(SimpleTreeCell, input, time_clones_samples, seed, outdir)
end

main()