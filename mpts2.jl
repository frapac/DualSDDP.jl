################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# Build primal and dual SP problems.
################################################################################

################################################################################
# PROBLEM DEFINITION
# import data from MPTS
### DEFINE HERE TOTAL NUMBER OF STAGES
NSTAGES = 3
α = 1 #/30 #13 / NSTAGES
### DEFINE HERE NUMBER OF NODES TO CONSIDER
NODES = 2

if NODES == 2
    NAMES = [:GER, :FRA]
elseif NODES == 4
    NAMES = [:FRA, :GER, :ESP, :UK]
elseif NODES == 7
    NAMES = [:FRA, :GER, :ESP, :UK, :PT, :ITA, :SUI]
elseif NODES == 8
    NAMES = [:FRA, :GER, :ESP, :UK, :PT, :ITA, :SUI, :BEL]
end

begin const
    XMAX, UTURB, UTHERM, X0, R, CTHERM_RAW, QMAX = getglobalparams(NAMES)
    UTURB  *= α
    UTHERM *= α
    QMAX   *= α
    NZONES = length(CTHERM_RAW)
    NARCS  = size(R, 2)
    CTHERM = CTHERM_RAW .+ 15*rand(NZONES, NSTAGES)
    NBINS  = 10
    SCEN   = "4"

    COST_HF = 3000
    CPENAL = 3000
    CTRANS = 1
end

"Return primal cost at time `t` as a vector."
getcost(t::Int) = [zeros(NZONES); zeros(NZONES); CTHERM[:, t]; CPENAL*ones(NZONES); CTRANS*ones(Float64, NARCS)]


################################################################################
# SCENARIO TREE DEFINITION
"""Quantize scenarios with KMeans quantization.  """
function optquantiz(scenarios, nbins; scale=true)
    # fix seed for reproductability
    srand(2713)
    ntime, nscenarios, nperturb = size(scenarios)

    output = Vector{NoiseLaw}(ntime)

    for t in 1:ntime
        data = copy(reshape(scenarios[t, :, :], nscenarios, nperturb)')
        # Launch kmeans:
        if scale
            μ    = mean(data, 2)
            σ    = std(data, 2)
            data = data .- μ
            pos  = vec(σ .> 0)
            data[pos, :] = data[pos, :] ./ σ[pos]
        end
        quantiz = kmeans(data, nbins)

        center =  quantiz.centers
        if scale
            center = σ .* quantiz.centers .+ μ
        end
        output[t] = NoiseLaw(center, 1/nscenarios*quantiz.counts)
    end

    return output
end

"Fit global noise laws corresponding to countries specified in `names`."
function fitgloballaw(names, nstages, nbins, Δt; nscen=10)

    nzones = length(names)
    scenarios = zeros(Float64, nstages, nscen, 2*nzones)

    ntime = 50

    for (idx, name) in enumerate(names)
        # inflow (in GWh)
        inflow = 24*readcsv("data/scen/$name/inflow.txt")[1:nscen, 1:ntime]
        # demand (in GWh)
        demand = 24*readcsv("data/scen/$name/demands.txt")[1:nscen, 1:ntime]

        nscen = size(inflow, 1)
        _inflow = zeros(nscen, nstages)
        _demand = zeros(nscen, nstages)

        # get monthly average
        for nn in 1:nstages
            s1 = nn
            s2 = nn + 1
            _inflow[:, nn] = sum(inflow[:, s1:s2], 2)
            _demand[:, nn] = sum(demand[:, s1:s2], 2)
        end

        scenarios[:, :, idx] = _inflow'
        scenarios[:, :, idx+nzones] = _demand'
    end

    return optquantiz(scenarios, nbins)
end

################################################################################
# MATRIX
"""Build matrix corresponding to the damsvalley"""
function MPTSmatrix()
    nzones, narcs = size(R)

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
         I O I I R; # u + Rq = w
         B; # xf max
         -B] # xf min
    # lots of constraints :(
    # nc = 10 * nzones + 2 * narcs

    balance = [I O I I R]
    Gt = [O I]

    return A, B, D, E, Gt, C,  balance
end



################################################################################
# PRIMAL PROBLEM
"Build primal SP model."
function build_model()
    dt = 30
    laws = fitgloballaw(NAMES, NSTAGES, NBINS, dt, nscen=50)
    nzones, narcs = size(R)

    A, B, D, E, Gt, C, balance = MPTSmatrix()

    dynamic(t, x, u, w) = A*x+ B*u + C*w
    constr(t, x, u, w) = balance * u - Gt * w
    cost_t(t, x, u, w) = dot(getcost(t), u)

    # Build bounds:
    x_bounds = [(0, xub) for xub in XMAX]
    u_bounds = vcat(
                [(0, uub) for uub in UTURB],
                [(0, Inf) for uub in UTURB],
                [(0, tub) for tub in UTHERM],
                [(0, Inf) for uub in UTURB],
                [(-f, f) for f in QMAX])

    function finalcost(model, m)
        alpha = m[:alpha]
        xf = m[:xf]
        z = @JuMP.variable(m, [1:nzones], lowerbound=0)
        @JuMP.constraint(m, z[i=1:nzones] .>= X0 - xf)
        @JuMP.constraint(m, alpha == COST_HF*sum(z))
    end

    model = SDDP.LinearSPModel(NSTAGES,       # number of timestep
                               u_bounds, # control bounds
                               X0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                               Vfinal=finalcost,
                               eqconstr=constr
                              )
    SDDP.set_state_bounds(model, x_bounds)
    return model
end


################################################################################
# DUAL PROBLEM
"""Build dual problem of MPTS."""
function buildemptydual(laws)
    nzones, narcs = size(R)
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

    model = SDDP.LinearSPModel(NSTAGES,       # number of timestep
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
    nzones, narcs = size(R)

    ndams = model.dimStates
    A, B, D, E, Gt, C, balance = MPTSmatrix()
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)

    getrhs(w, j) = -[UTURB;o;o;UTHERM;o;o;QMAX;QMAX; Gt*w[: ,j];XMAX - C*w[:, j];C*w[:, j]]

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
    c = getcost(t)
    for j=1:ns
        @constraint(m, c + B'*xf[:, j] + E'*u[:, j] .== 0)
    end

    #= # add objective as minimization of expectancy: =#
    i = ones(Float64, ndams)
    @objective(m, Min, sum(πp[j]*(-(dot(C*w[:, j], xf[:, j]) + dot(getrhs(w, j), u[:, j]))
                                  + alpha[j]) for j in 1:ns))

    # take care of final cost: final co-state must equal 0
    if t == model.stageNumber - 1
        for j=1:ns
            @constraint(m, xf[:, j] .<= 0)
            @constraint(m, xf[:, j] .>= -COST_HF)
            @constraint(m, alpha[j] == dot(X0, xf[:, j]))
            #= @constraint(m, xf[:, j] .== 0) =#
            #= @constraint(m, alpha[j] == 0) =#
        end
    end

    # store number of cuts for StochDynamicProgramming.jl
    m.ext[:ncuts] = 0

    return m
end
