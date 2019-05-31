################################################################################
# DP solvers
################################################################################

import StochDynamicProgramming: StochDynamicProgramming.SPModel

"""
Vanilla DP solver

Solve SPModel defined inside `model`.

# Return
* `V::Array{Float64, 2}`, size = (size_grid, ntime)
    Value functions

"""
function solvedp!(model::SPModel, trajectories::Vector{Trajectory})
    TF = model.stageNumber
    law = model.noises

    nstates = length(trajectories)

    # store resulting flows in array for multiplier's update
    V = zeros(Float64, nstates, TF)

    #Compute final value functions
    for (ix, traj) in enumerate(trajectories)
        x = view(traj.x[TF, :], :)
        V[ix, TF] = 3000 * sum(max.(0, model.initialState - x))
    end

    # Loop over time:
    @inbounds for t = (TF-1):-1:1
        lpproblem = buildlp(model, V, t, trajectories)

        # loop over state
        @inbounds for (ix, traj) in enumerate(trajectories)
            state = traj.x[t, :]
            JuMP.fix.(lpproblem.ext[:state], state)
            cost = 0.

            @inbounds for iξ in 1:law[t].supportSize
                optcost = solvelp(lpproblem, law[t].support[:, iξ])
                cost += law[t].proba[iξ] * optcost
            end
            # update value function
            V[ix, t] = cost
        end
    end
    return V
end


################################################################################
# SOLVER INTERFACE
################################################################################
"""Solve DP LP problem with warmstart"""
function solvelp(lpproblem, ξ)::Float64
    w = lpproblem[:w]
    JuMP.fix.(w, ξ)
    status = solve(lpproblem)
    if status ∉ [:Optimal, :Suboptimal]
        println(lpproblem)
        error("OPT")
    end
    return getobjectivevalue(lpproblem)
end


"""Build LP Model at a given stage for a single node.

WARNING: this function is specific to the given model
# Arguments
* `model::SPModel`
    Node's SPModel
* `V::Vector{PolyhedralFunction}`, size=(ntime,)
    Vector of value functions
* `t::Int`
    Current timestep
* `xfgrid`
    Interpolation grid

# Return
* `JuMP.Model`

"""
function buildlp(model::SPModel, V, t::Int, trajectories)
    # Build build corresponding to current trajectories.
    # quick and dirty...
    xfgrid = Array[traj.x[t+1, :] for traj in trajectories]

    nd = length(xfgrid)
    m = Model(solver=SOLVER)

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    # define variables in JuMP:
    m.ext[:state] = @variable(m, x[i=1:nx] )
    @variable(m, y[i=1:nx] )
    @variable(m, v[1:nx])
    @variable(m, model.ulim[i][1] <= u[i=1:nu] <= model.ulim[i][2])
    @variable(m, w[1:nw])
    # weights to approximate value function, in simplex
    @variable(m, 0 <= α[1:nd] <= 1)
    @constraint(m, sum(α) == 1)

    # Supply == Demand.
    constr_t = get(model.equalityConstraints)
    @constraint(m,  constr_t(t, x, u, w) .== 0.)

    # code the dynamic in the model
    # ∑ α_i V_{t+1}_i - f_t(x, u, w) == 0
    # as x and w are given, we store them in the right hand side (here set to 0.)
    #= @constraint(m, y .== model.dynamics(t, x, u, w)) =#
    @constraint(m, -y + sum(α[i]*xfgrid[i] for i in 1:nd) .== 0)

    # Lipschitz regularization of DP equations.
    @constraint(m, v .>= y - model.dynamics(t, x, u, w))
    @constraint(m, v .>= -y + model.dynamics(t, x, u, w))
    # Define objective function
    @objective(m, Min, LIPSCHITZ*sum(v) + model.costFunctions(t, x, u, w) + sum(α[i]*V[i, t+1] for i in 1:nd))

    return m
end
