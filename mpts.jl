
using MPTS


"""Build matrix corresponding to the damsvalley"""
function MPTSmatrix(nzones, narcs, R)
    I = eye(nzones)
    O = zeros(Float64, nzones, nzones)
    Oq = zeros(Float64, nzones, narcs)
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)

    # x+ = Ax + Bu + Cw
    # dimA: nx x nx
    A = I
    # u = [uturb; uspill, utherm, urec1, urec2, q]
    # (we get rid of F because we do not use decomposition)
    # dimB: nx x nu
    B = - [I I O O O Oq]
    # dimC: nx x nw
    C = [I O]


    c = [o; o; CTHERM*i; CPENAL*i; CPENAL*i ; zeros(Float64, narcs)]

    # Dx + Eu <= Gw
    # we note nc the number of constraints
    # dimD: nc x nx
    D = [O; O; O; O; O; O;O; O; O; O; A; -A]

    # dimE: nc x nu
    E = [I O O O O Oq; # uturb max
        -I O O O O Oq; # uturb min
         O -I O O O Oq; # uspill min
         O O I O O Oq; # utherm max
         O O -I O O Oq; # utherm min
         O O O -I O Oq; # penal1 min
         O O O O -I Oq; # penal2 min
         O O O O O I; # q max
         O O O O O -I; # q min
         I O I I -I R; # u + Aq = w
         B; # xf max
         -B] # xf min
    # lots of constraints :(
    # nc = 12 * nzones

    balance = [I O I I -I R]
    Gt = [O I]
    G = [uturb*i; o; o; utherm*i; o; o; o; o; -qmax, qmax; o; xmax; o]


    return A, B, D, E, G, Gt, C, c, balance
end

"""Build dual problem of MPTS."""
function buildual_dam1(laws)
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
    A, B, D, E, G, Gt, C, c, balance = MPTSmatrix(NZONES, NARCS, R)
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)
    g = [UTURB;o;o;UTHERM;o;o;o;QMAX;-QMAX; Gt*w[: ,j];XMAX - C*w[:, j];C*w[:, j]]

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
    @objective(m, Min, sum(πp[j]*(-(dot(C*w[:, j], xf[:, j])+ dot(g, u[:, j])
                                    + alpha[j]) for j in 1:ns)))

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

function buildprimal_MPTS()
    laws = build_noiselaws();
    i = ones(Float64, NZONES)
    c = vcat(-i, 0*i)

    # R = buildincidence
    A, B, D, E, G, Gt, C, c, balance = MPTSmatrix(NZONES, NARCS, R)

    dynamic(t, x, u, w) = A*x+ B*u + C*w
    constr(t, x, u, w) = balance * u - Gt * w
    cost_t(t, x, u, w) = dot(c, u)

    # Build bounds:
    x_bounds = XBOUNDS
    u_bounds = UBOUNDS

    model = SDDP.LinearSPModel(TF,       # number of timestep
                               u_bounds, # control bounds
                               X0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                               eqconstr=constr
                              )
    SDDP.set_state_bounds(model, x_bounds)
end
