using SomaticEvolution
using Random
using JSON

function fixed_mutations_over_time(::Type{T}, input, timesteps, seed, outputfile, repeats) where T<:SomaticEvolution.AbstractCell
    data = Vector{Vector{Vector{Int64}}}[]
    complete = 0
    # Threads.@threads for i in 1:repeats
    for i in 1:repeats
        rng = Random.seed!(seed + i)
        newdata = runsimulation_timeseries(
            T,
            input, 
            timesteps, 
            x -> [[mod.id for mod in x], clonal_mutations.(x)],
            rng
        )
        push!(data, newdata)
        complete += 1
        println("completed = $complete")
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
    tmax =  parse(Float64, ARGS[8])
    tstep =  parse(Float64, ARGS[9])
    repeats = parse(Int64, ARGS[10])
    modulebranching = Symbol(ARGS[11])
    outdir = ARGS[12]
    seed = rand(1:1000000)



    input = MultilevelBranchingMoranInput(;
        maxmodules=100,
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
        modulebranching=modulebranching,
	moranincludeself=false
    )
    mkpath(outdir)
    outdir *= "/N$(modulesize)_Ninit$(branchinitsize)"
    outdir *= "_μ$(μ)_r$(branchrate)"
    outdir *= "_moranrate$(moranrate)_asymrate$(asymmetricrate)"

    fixed_mutations_over_time(SimpleTreeCell, input, 0:tstep:tmax, seed, outdir, repeats)
end

main()
