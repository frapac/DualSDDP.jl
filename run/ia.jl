
println("RUN INNER APPROX")
# Build a new SDDP interface
sddp= SDDPInterface(model, params,
                    SDDP.IterLimit(MAX_ITER),
                    verbose_it=10)
# Replace model with inner approx
init_innermodeler!(sddp, sddpdual.bellmanfunctions)
#= println(sddp.solverinterface[2]) =#
srand(2713)
ci = @time SDDP.simulate(sddp, NSIMU)[1]
