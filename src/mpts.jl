################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# Build primal and dual SP problems.
################################################################################

################################################################################
# PROBLEM DEFINITION

#= if NODES == 2 =#
#=     NAMES = [:GER, :FRA] =#
#= elseif NODES == 4 =#
#=     NAMES = [:FRA, :GER, :ESP, :UK] =#
#= elseif NODES == 7 =#
#=     NAMES = [:FRA, :GER, :ESP, :UK, :PT, :ITA, :SUI] =#
#= elseif NODES == 8 =#
#=     NAMES = [:FRA, :GER, :ESP, :UK, :PT, :ITA, :SUI, :BEL] =#
#= end =#

"Store MPTS problem inside a dedicated structure."
struct MPTS
    # zones' names
    names
    # state upper bound
    XMAX
    # turbinate upper bound
    UTURB
    # thermal prod. upper bound
    UTHERM
    # initial position
    X0
    # node-arc incidence matrix
    R
    # thermal cost along time
    CTHERM
    # maximum transport through edges
    QMAX
    # number of stagges
    nstages::Int
    # number of nodes
    nzones::Int
    # number of arcs
    narcs::Int
    # quantization size to discretize scenarios
    nbins::Int
    # noise laws
    laws
end
function MPTS(names, nstages, nbins)
    nzones = length(names)
    XMAX, UTURB, UTHERM, X0, R, CTHERM_RAW, QMAX = getglobalparams(names)
    narcs  = size(R, 2)
    CTHERM = CTHERM_RAW .+ 15*rand(nzones, nstages)
    laws = fitgloballaw(names, nstages, nbins, nscen=50)

    return MPTS(names, XMAX, UTURB, UTHERM, X0, R, CTHERM, QMAX, nstages,
                nzones, narcs, nbins, laws)
end

"Return primal cost at time `t` as a vector."
getcost(t::Int, mpts::MPTS) = [zeros(mpts.nzones);
                               zeros(mpts.nzones);
                               mpts.CTHERM[:, t];
                               DATA["CPENAL"]*ones(mpts.nzones);
                               DATA["CTRANS"]*ones(Float64, mpts.narcs)]


################################################################################
# MATRIX
"""Build matrix corresponding to the damsvalley"""
function MPTSmatrix(mpts::MPTS)
    nzones, narcs = size(mpts.R)

    I = eye(nzones)
    Iq = eye(narcs)
    O = zeros(Float64, nzones, nzones)
    Oq = zeros(Float64, nzones, narcs)
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)

    # x+ = Ax + Bu + Cw
    # dimA: nx x nx
    A = I
    # u = [uturb, uspill, utherm, urec1, urec2, q]
    # (we get rid of F because we do not use decomposition)
    # dimB: nx x nu
    B = - [I I O O Oq]
    # dimC: nx x nw
    C = [I O]

    # Dx + Eu <= Gw
    # we note nc the number of constraints
    # dimD: nc x nx
    Oq2 = zeros(Float64, narcs, nzones)
    D = [O; O; O; O; O;O; Oq2; Oq2; O; A; -A]

    # dimE: nc x nu
    E = [I O O O Oq; # uturb max
        -I O O O Oq; # uturb min
         O -I O O Oq; # uspill min
         O O I O Oq; # utherm max
         O O -I O Oq; # utherm min
         O O O -I Oq; # utherm min
         Oq2 Oq2 Oq2 Oq2 Iq; # q max
         Oq2 Oq2 Oq2 Oq2 -Iq; # q min
         I O I I mpts.R; # u + Rq = w
         B; # xf max
         -B] # xf min
    # lots of constraints :(
    # nc = 10 * nzones + 2 * narcs

    balance = [I O I I mpts.R]
    Gt = [O I]

    return A, B, D, E, Gt, C,  balance
end



################################################################################
# PRIMAL PROBLEM
"Build primal SP model."
function build_model(mpts::MPTS)
    nzones, narcs = size(mpts.R)

    A, B, D, E, Gt, C, balance = MPTSmatrix(mpts)

    dynamic(t, x, u, w) = A*x+ B*u + C*w
    constr(t, x, u, w) = balance * u - Gt * w
    cost_t(t, x, u, w) = dot(getcost(t, mpts), u)

    # Build bounds:
    x_bounds = [(0, xub) for xub in mpts.XMAX]
    u_bounds = vcat(
                [(0, uub) for uub in mpts.UTURB],
                [(0, Inf) for uub in mpts.UTURB],
                [(0, tub) for tub in mpts.UTHERM],
                [(0, Inf) for uub in mpts.UTURB],
                [(-f, f) for f in mpts.QMAX])

    function finalcost(model, m)
        alpha = m[:alpha]
        xf = m[:xf]
        z = @JuMP.variable(m, [1:nzones], lowerbound=0)
        @JuMP.constraint(m, z[i=1:nzones] .>= mpts.X0 - xf)
        @JuMP.constraint(m, alpha == DATA["COST_HF"]*sum(z))
    end

    model = SDDP.LinearSPModel(mpts.nstages,       # number of timestep
                               u_bounds, # control bounds
                               mpts.X0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               mpts.laws,
                               Vfinal=finalcost,
                               eqconstr=constr
                              )
    SDDP.set_state_bounds(model, x_bounds)
    return model
end


################################################################################
# DUAL PROBLEM
"""Build dual problem of MPTS."""
function buildemptydual(mpts::MPTS)
    nzones, narcs = size(mpts.R)
    # Damsvalley configuration:
    x0 = [-3156.06 for i in 1:nzones]
    x_bounds = [(-1e4, 1e4) for i in 1:nzones];
    u_bounds = vcat([(0, Inf) for i in 1:(6*nzones+2*narcs)],
                    [(-Inf, Inf) for i in 1:nzones],
                    [(0, Inf) for i in 1:(2*nzones)])

    # write dual dynamic
    dynamic(t, x, u, w) = 0.
    constdual(t, x, u, w) = 0.
    cost_t(t, x, u, w) = 0.

    model = SDDP.LinearSPModel(mpts.nstages,       # number of timestep
                               u_bounds, # control bounds
                               x0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               mpts.laws,
                               info=:DH,
                               eqconstr=constdual
                              )
    SDDP.set_state_bounds(model, x_bounds)
    return model
end


"""Build model in Decision-Hazard."""
function build_model_dual(mpts::MPTS, model, param, t)
    nzones, narcs = size(mpts.R)

    ndams = model.dimStates
    A, B, D, E, Gt, C, balance = MPTSmatrix(mpts)
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)

    getrhs(w, j) = -[mpts.UTURB;
                     o;
                     o;
                     mpts.UTHERM;
                     o;
                     o;
                     mpts.QMAX;
                     mpts.QMAX;
                     Gt*w[: ,j];
                     mpts.XMAX - C*w[:, j];
                     C*w[:, j]]

    m = Model(solver=param.SOLVER)
    law = model.noises

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    ns = law[t].supportSize
    w = collect(law[t].support[:, :])
    πp = law[t].proba

    @variable(m, x[i=1:nx])
    @variable(m, model.ulim[i,t][1] <= u[i=1:nu, j=1:ns] <=  model.ulim[i,t][2])
    @variable(m, -LIPSCHITZ <= xf[i=1:nx, j=1:ns] <= 0)
    @variable(m, alpha[1:ns])

    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    # linking all samples to previous state
    @constraint(m, sum(πp[j]*(xf[:, j] + D'*u[:, j]) for j in 1:ns) .== x)

    # add admissible constraint in dual:
    c = getcost(t, mpts)
    for j=1:ns
        @constraint(m, c + B'*xf[:, j] + E'*u[:, j] .== 0)
    end

    #= # add objective as minimization of expectancy: =#
    i = ones(Float64, ndams)
    @objective(m, Min, sum(πp[j]*(-(dot(C*w[:, j], xf[:, j]) + dot(getrhs(w, j), u[:, j]))
                                  + alpha[j]) for j in 1:ns))

    # take care of final cost:
    if t == model.stageNumber - 1
        for j=1:ns
            @constraint(m, xf[:, j] .<= 0)
            @constraint(m, xf[:, j] .>= -DATA["COST_HF"])
            @constraint(m, alpha[j] == dot(mpts.X0, xf[:, j]))
        end
    end

    # store number of cuts for StochDynamicProgramming.jl
    m.ext[:ncuts] = 0

    return m
end
