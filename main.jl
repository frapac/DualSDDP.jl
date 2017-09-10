# Run dual SDDP on dual of dam1

include("sddp_optim.jl")

# params
MAXIT = 100
NSIMU = 10000

### Build primal problem
model, params = build_model()
sddpprimal = SDDPInterface(model, params,
                    SDDP.IterLimit(MAX_ITER),
                    verbose_it=10)
SDDP.init!(sddpprimal)

### Build dual problem
modeldual = buildual_dam1(build_noiselaws())
sddpdual = SDDPInterface(modeldual, params,
                    SDDP.IterLimit(MAX_ITER),
                    verbose_it=0)
initdual!(sddpdual)


### SDDP DUAL ####
# Run 1 combined iteration to init cuts in sddpdual
for iter in 1:20
    SDDP.iteration!(sddpprimal, sddpdual)
end

### RUN iterations in dual
lbdual = Float64[]
println("RUN DUAL SDDP")
lb = updateinitialstate!(sddpdual, [40. for i in 1:N_DAMS])
tic()
for iter in 1:MAXIT
    # Update initial costate
    lb = updateinitialstate!(sddpdual, [40. for i in 1:N_DAMS])
    # Run forward an backward pass
    SDDP.iteration!(sddpdual)

    # save current iterations
    push!(lbdual, lb)
    (iter % 10 == 0) && display(iter, lb)
end
texec = toq()
println("Dual exec time: ", texec)


### RUN iterations in primal
println("RUN PRIMAL SDDP")
tic()
for iter in 1:MAXIT
    SDDP.iteration!(sddpprimal)
end
texec = toq()
println("Primal exec time: ", texec)


### MONTE CARLO ESTIMATION

println("RUN MONTE CARLO")
c = SDDP.simulate(sddpprimal, NSIMU)[1]

lbprimal = sddpprimal.stats.lower_bounds[end]
ubdual = lbdual[end]

println("#"^70)
println("Results")
println("-------")
println("Primal LB:\t", lbprimal)
println("Dual UB:\t", -ubdual)
println("Monte Carlo:\t", mean(c))

