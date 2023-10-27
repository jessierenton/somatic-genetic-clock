using Printf

N_N0vals = [(3,1), (5,1), (5,2), (20,1), (20,2), (20,10), (100,1), (100,2), (100,50)]
μvals = [0.01, 1]
rvals = [1, 2, 4, 8, 16]
#birthrate, moranrate, asymmetricrate
rates = [(48, 48, 0), (48, 0, 48)]

open("input.txt", "w") do io
    write(io, "#birthrate, moranrate, asymmetricrate, N, N0, μ, r\n")
    for ((N, N0), μ, r, (birthrate, moranrate, asymmetricrate)) in Base.product(N_N0vals, μvals, rvals, rates)
            write(io, @sprintf(
                "%2d  %2d  %2d  %3d  %3d  %.2f  %.1f\n",
                birthrate,
                moranrate,
                asymmetricrate,
                N,
                N0,
                μ,
                r
            )
        )
    end
end

open("input_smalltime.txt", "w") do io
    write(io, "#birthrate, moranrate, asymmetricrate, N, N0, μ, r, tmax, tstep\n")
    for ((N, N0), μ, r, (birthrate, moranrate, asymmetricrate)) in Base.product(N_N0vals, μvals, rvals, rates)
            write(io, @sprintf(
                "%2d  %2d  %2d  %3d  %3d  %.2f  %.1f  %2.1f  %1.2f\n",
                birthrate,
                moranrate,
                asymmetricrate,
                N,
                N0,
                μ,
                r,
		8 ./ r,
		0.16 ./ r
            )
        )
    end
end

