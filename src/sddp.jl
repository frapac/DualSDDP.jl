################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# SDDP dual code.
################################################################################


"""Init primal SDDP interface"""
function initprimal(mpts; maxit=100)
    model = build_model(mpts)
    params = getparams()
    # init SDDP interface
    sddpprimal = SDDPInterface(model, params, SDDP.IterLimit(maxit),
                               verbose_it=0)
    SDDP.init!(sddpprimal)
    return sddpprimal
end


"""Init dual SDDP interface"""
function initdual(mpts::MPTS, sddp; maxit=100)
    modeldual = buildemptydual(mpts)
    sddpdual = SDDPInterface(modeldual, sddp.params,
                            SDDP.IterLimit(maxit),
                            verbose_it=0)
    initdual!(mpts, sddpdual)
    return sddpdual
end


"""Run SDDP primal alone."""
function runprimal!(sddpprimal; nbsimu=100, maxiterations=100,
                   Δsimu=100)
    println("RUN PRIMAL SDDP")
    ubp = []
    stdp = []
    scen = SDDP.simulate_scenarios(sddpprimal.spmodel.noises, nbsimu)

    tic()
    for iter in 1:maxiterations
        # run a single SDDP iteration with StochDynamicProgramming
        SDDP.iteration!(sddpprimal)

        # compute statistical upper bound
        if iter % Δsimu == 0
            cost = SDDP.simulate(sddpprimal, scen)[1]
            push!(ubp, mean(cost))
            push!(stdp, std(cost))
        end
        # reload JuMP model to avoid memory leak
        (iter % 10 == 0) && SDDP.reload!(sddpprimal)
    end

    texec = toq()
    println("Primal exec time: ", texec)

    return ubp, stdp
end


"""Run SDDP dual alone"""
function rundual!(sddpdual, sddpprimal; nbsimu=100, maxiterations=100,
                 Δsimu=100)
    ubd = Float64[]
    stdd = Float64[]
    X0 = sddpprimal.spmodel.initialState

    # Run 1 combined iteration to init cuts in sddpdual
    for iter in 1:1
        SDDP.iteration!(sddpprimal, sddpdual)
    end
    p0 = SDDP.get_subgradient(sddpprimal.bellmanfunctions[1], X0)


    ### RUN iterations in dual
    lbdual = Float64[]
    timedual = Float64[]
    println("RUN DUAL SDDP")
    lb = updateinitialstate!(sddpdual, X0)
    tic()
    for iter in 1:maxiterations
        tic()
        # Update initial costate
        lb = updateinitialstate!(sddpdual, X0)

        # Run forward an backward pass
        SDDP.iteration!(sddpdual)
        tdual = toq()

        # save current iterations
        push!(lbdual, lb)
        push!(timedual, tdual)
        (iter % 10 == 0) && displayit(iter, lb)
    end
    texec = toq()
    println("Dual exec time: ", texec)
    return lbdual, timedual
end


"""Run jointly primal and dual SDDPs."""
function runjoint!(sddpprimal, sddpdual; nbsimu=100, maxiterations=100,
                  Δsimu=100)
    # evolution of UB and LB
    lbdual = Float64[]
    timedual = Float64[]
    ubp = []
    stdp = []

    X0 = sddpprimal.spmodel.initialState
    # generate scenarios to estimate statistical UB
    scen = SDDP.simulate_scenarios(sddpdual.spmodel.noises, nbsimu)

    println("RUN DUAL SDDP")

    # Time counter
    tic()

    # Run SDDP iterations
    for iter in 1:maxiterations
        tic()
        # perform a mixed iteration between primal and dual SDDP
        td = SDDP.iteration!(sddpprimal, sddpdual)

        # update initial co-state
        lb = updateinitialstate!(sddpdual, X0)
        # run a CUPPS forward pass in dual
        SDDP.fwdcuts(sddpdual)

        # reupdate initial co-state
        lb = updateinitialstate!(sddpdual, X0)
        tdual = toq() + td

        # save current iterations
        push!(lbdual, lb)
        push!(timedual, tdual)

        # if specified, compute statistical UB
        if iter % Δsimu == 0
            cost = SDDP.simulate(sddpprimal, scen)[1]
            push!(ubp, mean(cost))
            push!(stdp, std(cost))
        end

        if (iter % 10 == 0)
            print("Pass n\° ", iter)
            @printf("\tLB: %.4e ", sddpprimal.stats.lowerbound)
            @printf("\tUB: %.4e \n", lb)
        end
        # reload JuMP model to avoid memory leak
        (iter % 10 == 0) && SDDP.reload!(sddpprimal)
    end
    texec = toq()
    println("Total exec time: ", texec)

    return lbdual, timedual, ubp, stdp
end
