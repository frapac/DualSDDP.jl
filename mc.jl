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
    cost = SDDP.simulate(sddpprimal, scen)[1]

    return mean(cost), std(cost)
end


function computeprimalMC(sddp, ti, to, dt, nscen=1000)
    v = copy(sddp.bellmanfunctions)
    μmc = Float64[]
    σmc = Float64[]

    srand(2713)
    scen = SDDP.simulate_scenarios(sddpprimal.spmodel.noises, nscen)
    for it in ti:dt:to
        c, s = mc!(sddp, v, it, scen)
        @printf("%s: %.3e", it, c)
        @printf("\t%.3e\n", s)
        push!(μmc, c)
        push!(σmc, s)
    end
    return μmc, σmc
end


function computedualMC(sddpdual, ti, to, dt, nscen=1000)
    Vdual = copy(sddpdual.bellmanfunctions)

    model, params = sddpprimal.spmodel, sddpprimal.params
    sddp = SDDPInterface(model, params, SDDP.IterLimit(MAX_ITER), verbose_it=10)

    v = copy(sddp.bellmanfunctions)
    μmc = Float64[]
    σmc = Float64[]

    srand(2713)
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
