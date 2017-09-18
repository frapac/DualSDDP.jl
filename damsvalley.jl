

X0 = [STOCK0 for i in 1:N_DAMS]


"""Build matrix corresponding to the damsvalley"""
function getmatrix(ndams, t)
    I = eye(ndams)
    O = zeros(Float64, ndams, ndams)
    i = ones(Float64, ndams)
    o = zeros(Float64, ndams)

    A = I
    B = - [I I]
    bt = [o' o'; I I][1:ndams, :]
    B += bt

    c = vcat(-COST[t]*i, 0*i)

    D = vcat(I, -I, O, O, O, I, -I)
    E = [O O; O O; I O; -I O; O -I; B; -B]


    return A, B, D, E, c
end


"""Build dual problem of `dam1`."""
function buildemptydual(laws)
    # Damsvalley configuration:
    x0 = [-3156.06 for i in 1:N_DAMS]
    x_bounds = [(-1e4, 1e4) for i in 1:N_DAMS];
    u_bounds = [(-Inf, 0) for i in 1:7*N_DAMS]

    # write dual dynamic
    dynamic(t, x, u, w) = 0.
    constdual(t, x, u, w) = 0.
    cost_t(t, x, u, w) = 0.

    model = SDDP.LinearSPModel(TF,       # number of timestep
                               u_bounds, # control bounds
                               x0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                               info=:DH,
                               eqconstr=constdual
                              )
    SDDP.set_state_bounds(model, x_bounds)
    return model
end


"""Build model in Decision-Hazard."""
function build_model_dual(model, param, t)
    ndams = model.dimStates
    A, B, D, E, c = getmatrix(ndams, t)

    m = Model(solver=param.SOLVER)
    law = model.noises

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    ns = law[t].supportSize
    w = collect(law[t].support[:, :])
    πp = law[t].proba

    @variable(m,  x[i=1:nx] )
    @variable(m, model.ulim[i,t][1] <= u[i=1:nu, j=1:ns] <=  model.ulim[i,t][2])
    @variable(m, xf[i=1:nx, j=1:ns])
    @variable(m, alpha[1:ns] >= 0)

    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    # linking all samples to previous state
    @constraint(m, sum(πp[j]*(xf[:, j] - D'*u[:, j]) for j in 1:ns) .== x)

    # add admissible constraint in dual:
    for j=1:ns
        @constraint(m, c + B'*xf[:, j] - E'*u[:, j] .== 0)
    end

    #= # add objective as minimization of expectancy: =#
    i = ones(Float64, ndams)
    @objective(m, Min, sum(πp[j]*(-(dot(w[:, j], xf[:, j])+
                                    dot(vcat(80*i, 0*i, 40*i, 0*i, 0*i, 80*i - w[:, j], w[:, j]), u[:, j]))
                                    + alpha[j]) for j in 1:ns))

    # take care of final cost: final co-state must equal 0
    if t == model.stageNumber - 1
        for j=1:ns
            @constraint(m, xf[:, j] .== 0)
            @constraint(m, alpha[j] == 0)
        end
    end

    # store number of cuts for StochDynamicProgramming.jl
    m.ext[:ncuts] = 0

    return m
end


# Read costs along time:
const α = EFFICIENCY*COST_CONV
const COST = collect(Float64, readdlm("$DATA_SOURCE/generated_prices.txt") * α)

# Read intakes:
INTAKES = [readdlm("$DATA_SOURCE/generated_intakes_$i.txt") for i in 0:N_DAMS-1]


# Build noiselaws with input data:
function build_noiselaws()
    laws = Vector{SDDP.NoiseLaw}(TF-1)
    n_scenarios = size(INTAKES[1], 1)
    # Uniform probabilities:
    proba = 1./n_scenarios*ones(n_scenarios)
    for t in 1:TF-1
        alea = hcat([INTAKES[i][:, t] for i in 1:N_DAMS]...)'
        laws[t] = SDDP.NoiseLaw(alea, proba)
    end
    return laws
end


# Build dynamic model and SDDP solver:
function build_model()
    srand(1111)
    laws = build_noiselaws();
    i = ones(Float64, N_DAMS)
    c = vcat(-i, 0*i)

    A, B, D, E, _ = getmatrix(N_DAMS, 1)
    dynamic(t, x, u, w) = x + B*u + w
    cost_t(t, x, u, w) = COST[t] * dot(c, u)

    # Build bounds:
    x_bounds = [(STOCK_MIN, STOCK_MAX) for i in 1:N_DAMS];
    u_bounds = vcat([(CONTROL_MIN, CONTROL_MAX) for i in 1:N_DAMS], [(0., Inf) for i in 1:N_DAMS]);

    model = SDDP.LinearSPModel(TF,       # number of timestep
                               u_bounds, # control bounds
                               X0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                              )
    SDDP.set_state_bounds(model, x_bounds)

    params = SDDP.SDDPparameters(SOLVER, passnumber=FORWARD_PASS,
                                 gap=EPSILON, max_iterations=MAX_ITER,
                                 confidence=CONFIDENCE_LEVEL,
                                 montecarlo_in_iter=MONTE_CARLO_SIZE,
                                 montecarlo_final=FINAL_MONTE_CARLO_SIZE,
                                 prune_cuts=PRUNE_CUTS,
                                 pruning_algo="exact+",
                                 rho0=0, alpha=.95)
    params.compute_ub = UPPER_BOUND
    return model, params
end

