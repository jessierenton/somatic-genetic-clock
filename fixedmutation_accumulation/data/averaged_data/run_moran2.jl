using SomaticEvolution
using Random
using JSON
using StatsBase

function fixed_mutations_over_time(::Type{T}, input, timesteps, seed, outputfile, repeats; nsample=10) where T<:SomaticEvolution.AbstractCell

    data = Vector{Vector{Float64}}[]
    complete = 0
    Threads.@threads for i in 1:repeats
    # for i in 1:repeats
        rng = Random.seed!(seed + i)
        
        function getdata(population)
            N = length(population)
            sampled_modules = 
                if N > nsample
                    sample(rng, population, nsample; replace=false)
                else
                    population
                end
            clonalmuts = clonal_mutations.(sampled_modules)
            return [mean(clonalmuts), var(clonalmuts), N]
        end

        newdata = runsimulation_timeseries(
            T,
            input, 
            timesteps, 
            getdata,            
            rng
        )
        push!(data, newdata)
        complete += 1
        println("completed = $complete\n")
        flush(stdout)
    end
    inputdict = SomaticEvolution.inputdict(input)
    open("$outputfile.json", "w") do io
        JSON.print(io, Dict(:input => inputdict, :seed => seed, :timestep => step(timesteps), :fixedmutations => data), 4)
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
    tmax =  parse(Float64, ARGS[10])
    tstep =  parse(Float64, ARGS[11])
    repeats = parse(Int64, ARGS[12])
    outdir = ARGS[13]
    seed = rand(1:1000000)



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
    outdir *= "_moranrate$(moranrate)_asymrate$(asymmetricrate)"

    fixed_mutations_over_time(SimpleTreeCell, input, 0:tstep:tmax, seed, outdir, repeats)
end

main()