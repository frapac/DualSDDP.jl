################################################################################
# Run Monte-Carlo estimation
################################################################################

function tailcuts(Vs, ncuts)

    vtails = PolyhedralFunction[]

    for i in 1:endof(Vs)
        β, λ = Vs[i].betas, Vs[i].lambdas

        nc = min(ncuts, Vs[i].numCuts)

        v = PolyhedralFunction(β[1:nc], λ[1:nc, :])
        push!(vtails, v)
    end

    return vtails
end


function mc!(sddp::SDDPInterface, Vs, ncuts::Int, scen)
    # update Bellman value functions
    sddp.bellmanfunctions = tailcuts(Vs, ncuts)
    # reload JuMP Model
    SDDP.reload!(sddp)

    # run Monte Carlo
    cost = SDDP.simulate(sddp, scen)[1]

    return mean(cost), std(cost)
end


function computeprimalMC(sddpprimal, ti, to, dt; nscen=1000,seed=2713)
    V = copy(sddpprimal.bellmanfunctions)
    μmc = Float64[]
    σmc = Float64[]

    model, params = sddpprimal.spmodel, sddpprimal.params
    sddp = SDDPInterface(model, params, SDDP.IterLimit(to), verbose_it=10)

    srand(seed)
    scen = SDDP.simulate_scenarios(sddp.spmodel.noises, nscen)
    for it in ti:dt:to
        tic()
        #c, s = mc!(sddpprimal, v, it, scen)
        SDDP.reload!(sddp)
        c, s = mc!(sddp,V,it,scen)
        tmc = toq()
        @printf("%s: %.3e", it, c)
        @printf("\t%.3e", s)
        @printf("\t%.0fs\n", tmc)
        push!(μmc, c)
        push!(σmc, s)
    end
    return μmc, σmc
end


function computedualMC(sddpprimal,sddpdual, ti, to, dt; nscen=1000,seed=2713)
    Vdual = copy(sddpdual.bellmanfunctions)

    model, params = sddpprimal.spmodel, sddpprimal.params
    sddp = SDDPInterface(model, params, SDDP.IterLimit(to), verbose_it=10)

    v = copy(sddp.bellmanfunctions)
    μmc = Float64[]
    σmc = Float64[]

    srand(seed)
    scen = SDDP.simulate_scenarios(sddpprimal.spmodel.noises, nscen)

    for it in ti:dt:to
        SDDP.reload!(sddp)
        vtail = tailcuts(Vdual, 2*it)
        init_innermodeler!(sddp, vtail)
        ci = SDDP.simulate(sddp, scen)[1]

        @printf("%s: %.3e", it, mean(ci))
        @printf("\t%.3e\n", std(ci))
        push!(μmc, mean(ci))
        push!(σmc, std(ci))
    end

    return μmc, σmc
end
