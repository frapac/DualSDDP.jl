################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# SDDP dual code.
################################################################################


struct Trajectory
    x::Array{Float64, 2}
    Trajectory(x::Array{Float64, 2}) = new(x)
    Trajectory(x::Array{Float64, 3}) = new(reshape(x, (size(x, 1), size(x, 3))))
end

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


"""Run SDDP primal alone.
nbsimu is the number of MonteCarlo scenarios used to evaluate the upper-bound
every Δsimu iterations. Evaluations returned in ubp, with standard deviation in stdp,
and trajectories in trajs."""
function runprimal!(sddpprimal; nbsimu=100, maxiterations=100,
                   Δsimu=100)
    println("RUN PRIMAL SDDP")
    ubp = []
    stdp = []
    scen = SDDP.simulate_scenarios(sddpprimal.spmodel.noises, nbsimu)

    trajs = Trajectory[]

    tic()
    for iter in 1:maxiterations
        # run a single SDDP iteration with StochDynamicProgramming
        SDDP.iteration!(sddpprimal)
        push!(trajs, Trajectory(SDDP.simulate(sddpprimal, 1)[2]))

        # compute statistical upper bound
        if iter % Δsimu == 0
            cost = SDDP.simulate(sddpprimal, scen)[1]
            push!(ubp, mean(cost))
            push!(stdp, std(cost))
        end
        # reload JuMP model to avoid memory leak
        (iter % 10 == 0) && SDDP.reload!(sddpprimal)

        if (iter % 10 == 0)
            print("Pass n\° ", iter)
            @printf("\tLB-P: %.4e ", sddpprimal.stats.lowerbound)
        end
    end

    texec = toq()
    println("Primal exec time: ", texec)

    return ubp, stdp, trajs
end


"""Run SDDP dual alone
nbsimu is the number of MonteCarlo scenarios used to evaluate the upper-bound
every Δsimu iterations. Evaluations returned in ubd, with standard deviation in stdd,
and trajectories in trajs."""
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


"""Run jointly primal and dual SDDPs.
nbsimu is the number of MonteCarlo scenarios used to evaluate the upper-bound
every Δsimu iterations. Evaluations returned in ubp, with standard deviation in stdp,
and trajectories in trajs."""
function runjoint!(sddpprimal, sddpdual; nbsimu=100, maxiterations=100,
                  Δsimu=100)
    # evolution of UB and LB
    ubd = Float64[]
    timedual = Float64[]
    ubp_mc = Float64[]
    c_ia_mc =Float64[]
    stdp = Float64[]
    stdia = Float64[]

    trajp = Trajectory[]

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
        td, traj = SDDP.iteration!(sddpprimal, sddpdual)
        # save primal trajectory
        push!(trajp, Trajectory(traj))

        # update initial co-state
        updateinitialstate!(sddpdual, X0)
        # run a CUPPS forward pass in dual
        SDDP.fwdcuts(sddpdual)

        # reupdate initial co-state
        ub = updateinitialstate!(sddpdual, X0)
        tdual = toq() + td

        # save current iterations
        push!(ubd, ub)
        push!(timedual, tdual)

        # if specified, compute statistical UB
        if iter % Δsimu == 0
            # OA MC
            c_oa = SDDP.simulate(sddpprimal, scen)[1]
            push!(ubp_mc, mean(c_oa))
            push!(stdp, std(c_oa))
            # IA MC - TODO NOT WORKING
            # sddp_IA  = deepcopy(sddpprimal)
            # V_dual = copy(sddp_IA.bellmanfunctions);
            # # we init the JuMP model inside the primal SDDP object (we simulate in the primal, not in the dual!)
            # init_innermodeler!(sddp_IA, V_dual)
            # c_ia = SDDP.simulate(sddp_IA, scen)[1]
            # push!(c_ia_mc, mean(c_ia))
            # push!(stdia, std(c_ia))
        end

        if (iter % 10 == 0)
            print("Pass n\° ", iter)
            @printf("\tLB-P: %.4e ", sddpprimal.stats.lowerbound)
            @printf("\tUB-D: %.4e \n", ub)
        end
        # reload JuMP model to avoid memory leak
        (iter % 51 == 0) && SDDP.reload!(sddpprimal)
    end
    texec = toq()
    println("Total exec time: ", texec)


    return ubd, timedual, ubp_mc, stdp, trajp #, c_ia_mc, sdtia
end
