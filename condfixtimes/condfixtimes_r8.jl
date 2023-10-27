using SomaticEvolution
using Statistics
using Random
using Distributions
using StatsBase
using CSV
using DataFrames

function condfixtime_diffusion(N, N0::Int64, r=1; p=0)
    if p==0
        return (1/r)*4*N0*(1-N0/N)
    else
        return -1/(r*p) * 4*N0*(1-N0/N) * (1-p) * log(1-p)
    end
end

function condfixtime_diffusion(N, x0::Float64, r=1)
    return (1/r)*4*N*x0*(1-x0)
end


function wrightfisher_adjusted_split(N, N0, r=1, Nmutants0=1; discretetime=true, rng=Random.GLOBAL_RNG)
    t = 0
    Nmuts = Nmutants0
    while true
        if Nmuts == 0 || Nmuts == N
            return (fixed = Nmuts == N, t=t)
        end
        m = rand(rng, [N0, N-N0])
        sampled_muts = rand(rng, Hypergeometric(Nmuts, N-Nmuts, m))
        # sampled_muts = rand(rng, Binomial(m, Nmuts/N))
        Nmuts = rand(rng, Binomial(N, sampled_muts/m))
        if discretetime
            t += 1/r
        else
            t += rand(rng, Exponential(1/r))
        end
    end
end

function wrightfisher_adjusted_branch(N, N0, r=1, Nmutants0=1; discretetime=true, rng=Random.GLOBAL_RNG)
    t = 0
    Nmuts = Nmutants0
    while true
        if Nmuts == 0 || Nmuts == N
            return (fixed = Nmuts == N, t=t)
        end
        if rand(rng) < 0.5
            m = N0
            sampled_muts = rand(rng, Hypergeometric(Nmuts, N-Nmuts, m))
            Nmuts = rand(rng, Binomial(N, sampled_muts/m))
        end
        if discretetime
            t += 1/r
        else
            t += rand(rng, Exponential(1/r))
        end
    end
end

function wrightfisher_adjusted_with_growth(N, N0, r=1, b=nothing, Nmutants0=1; discretetime=true, rng=Random.GLOBAL_RNG)
    t = 0
    Nmuts = Nmutants0
    while true
        if Nmuts == 0 || Nmuts == N
            return (fixed = Nmuts == N, t=t)
        end
        Ncells = rand(rng, [N0, N-N0])
        Nmuts = rand(rng, Hypergeometric(Nmuts, N-Nmuts, Ncells))
        # growth phase
        while Ncells < N
            if rand(rng) < Nmuts/Ncells
                Nmuts += 1
            end
            if !isnothing(b)
                t += rand(rng, Exponential(1/(b*Ncells)))
            end
            Ncells += 1
        end
        if discretetime
            t += 1/r
        else
            t += rand(rng, Exponential(1/r))
        end
    end
end

function moran(N, b=1, Nmutants0=1; discretetime=true, rng=Random.GLOBAL_RNG)
    r = N*b
    t = 0
    cells = zeros(Int64, N)
    cells[1:Nmutants] .= 1
    while true
        if Nmuts == 0 || Nmuts == N
            return (fixed = Nmuts == N, t=t)
        end
        #kill cell
        popat!(cells, rand(rng, 1:N))
        #divide cell
        push!(cells, rand(rng, cells))
        if discretetime
            t += 1/r
        else
            t += rand(rng, Exponential(1/r))
        end
    end
end

function condfixtime(populationmodel, N, N0, r=1, Nmutants0=1; n=1000, discretetime=true, rng=Random.GLOBAL_RNG)
    i = 1
    tvec = Float64[]
    while i <= n
        fixed, t = populationmodel(N, N0, r, Nmutants0; discretetime, rng)
        if fixed
            push!(tvec, t)
            i += 1
        end
    end
    return (
        condfixtime_mean=mean(tvec), 
        condfixtime_std=std(tvec), 
        nsamples=length(tvec)
    )
end

function condfixtime_SE(N, N0, r, b, tmax; nmodules=1, repeats=3, update=:asymmetric, rng=Random.GLOBAL_RNG, modulebranching=:split)
    allfixtimes = reduce(vcat, [
        begin
            input = MultilevelBranchingMoranInput(;
                modulesize=N,
                fixedmu=false,
                birthrate=b,
                deathrate=0,
                moranrate=(update==:moran ? b : 0),
                asymmetricrate=(update==:asymmetric ? b : 0),
                clonalmutations=0,
                tmax=tmax,
                maxmodules=nmodules,
                branchrate=r,
                branchinitsize=N0,
                μ=0.01,
                modulebranching
            )
            condfixtimes = reduce(vcat, runsimulation_condfixtime(Cell, input, rng)[1])
        end
        for i in 1:repeats
    ])
    return (
        condfixtime_mean=mean(allfixtimes), 
        condfixtime_std=std(allfixtimes), 
        nsamples=length(allfixtimes)
    )
end

N_N0vals = vcat(
    [(3,1), (3,2)],
    [(5,i) for i in 1:4],
    [(20,i) for i in 1:2:19],
    [(100,i) for i in vcat([1], 10:10:90)]
)
rvals = [8]
repeats=20
b = 120
tmax_asym = 1000
tmax_moran = 200
df = DataFrame(
    r=Float64[], 
    N=Int64[],
    N0=Int64[], 
    simulation=String[], 
    condfixtime_mean=Float64[],
    condfixtime_std=Float64[],
    nsamples=Int64[]
)
Threads.@threads  for (r, (N,N0)) in collect(Base.product(rvals, (N_N0vals)))
    cft_SE_asym_branch = condfixtime_SE(N, N0, r, b, tmax_asym; nmodules=1, repeats=repeats, update=:asymmetric, modulebranching=:withoutreplacement_nomutations)
    cft_SE_asym_split = condfixtime_SE(N, N0, r, b, tmax_asym; nmodules=1, repeats=repeats, update=:asymmetric, modulebranching=:split)
    cft_SE_moran_split = condfixtime_SE(N, N0, r, b, tmax_moran; nmodules=1, repeats=repeats, update=:moran, modulebranching=:split)
    push!(df, (r, N, N0, "full (asymmetric, branch)", cft_SE_asym_branch...))
    push!(df, (r, N, N0, "full (asymmetric, split)", cft_SE_asym_split...))
    push!(df, (r, N, N0, "full (moran, split)", cft_SE_moran_split...))

    println("r = $r, N0 = $N0 complete")
end
sort!(df, [:r, :N, :N0, :simulation])
CSV.write("df_8.csv", df)
