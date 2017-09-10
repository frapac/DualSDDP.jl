# Write dual version of problems in problem.jl

#------------------------------
# 1 dam valley
#------------------------------
"""Build dual problem of `dam1`."""
function buildual_dam1(laws)


    # Damsvalley configuration:
    x0 = [-3156.06 for i in 1:N_DAMS]

    x_bounds = [(-1e9, 1e9) for i in 1:N_DAMS];
    u_bounds = [(-Inf, 0) for i in 1:7*N_DAMS]

    B = - [1. 1.]
    # Discretization of controls for assessment:
    DU = [10, 8]

    D = [1; -1; -1; 1; 0; 0; 0]
    E = [0 0; 0 0; 1 1; -1 -1; 1 0; -1 0; 0 -1]

    # write dual dynamic
    function dynamic(t, x, u, w)
        return 0.
    end

    function constdual(t, x, u, w)
        return 0.
    end

    function cost_t(t, x, u, w)
        return 0.
    end

    model = SDDP.LinearSPModel(TF,       # number of timestep
                               u_bounds, # control bounds
                               STOCK_TARGET,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                               info=:DH,
                               eqconstr=constdual
                              )
    SDDP.set_state_bounds(model, x_bounds)
    return model
end


function getmatrix(t)
    B = - [1. 1.]
    D = [1; -1; -1; 1; 0; 0; 0]
    E = [0 0; 0 0; 1 1; -1 -1; 1 0; -1 0; 0 -1]
    c = [-COST[t]; 0]

    return B, D, E, c
end

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

    @variable(m, model.xlim[i,t][1] <= x[i=1:nx] <= model.xlim[i,t][2])
    @variable(m, model.ulim[i,t][1] <= u[i=1:nu, j=1:ns] <=  model.ulim[i,t][2])
    @variable(m, xf[i=1:nx, j=1:ns])
    @variable(m, alpha[1:ns])

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


"""Overwrite JuMP Model in `sddp` to consider dual version."""
function initdual!(sddp)
    ms = JuMP.Model[build_model_dual(sddp.spmodel, sddp.params, t) for t=1:sddp.spmodel.stageNumber-1]
    sddp.solverinterface = ms
end


"""Compute dual upper bound with `max_p <x0, p> - D_0(p)`."""
function dualbound(sddpdual, x)
    V = sddpdual.bellmanfunctions[1]

    m = Model(solver=sddpdual.params.SOLVER)

    @variable(m, θ)
    p = @variable(m, p[1:sddpdual.spmodel.dimStates])

    for i in 1:V.numCuts
        lambda = vec(V.lambdas[i, :])
        @constraint(m, V.betas[i] + dot(lambda, p) <= θ)
    end

    @objective(m, Min, θ - dot(x, p))
    solve(m)
    return getobjectivevalue(m), getvalue(p)
end


"""Change initial co-state for dual SDDP."""
function updateinitialstate!(sddpdual, x)
    lb, p0 = dualbound(sddpdual, x)
    sddpdual.spmodel.initialState = p0
    return lb
end

function display(it, lb)
    print("Pass n\° ", it)
    @printf("\tLower-bound: %.4e \n", lb)
end
