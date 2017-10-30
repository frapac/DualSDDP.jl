

function run_jointapprox(w)
    @assert 0 <= w <= 1
    sddp= SDDPInterface(sddpprimal.spmodel, sddpprimal.params,
                        SDDP.IterLimit(MAX_ITER),
                        verbose_it=10)
    # Replace model with inner approx
    init_jointmodeler!(sddp, sddpdual.bellmanfunctions, sddpprimal.bellmanfunctions, w)
    srand(2713)
    return @time SDDP.simulate(sddp, NSIMU)[1]
end

println("RUN JOINT APPROX")
run_jointapprox(0.1)
