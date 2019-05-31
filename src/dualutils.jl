################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# Write dual version of problems in problem.jl
################################################################################

# Lipschitz constant
const LIPSCHITZ = 10e5

"Define Linear Bellman Operator model with JuMP."
function initdual!(mpts, sddp)
    ms = JuMP.Model[build_model_dual(mpts, sddp.spmodel, sddp.params, t) for t=1:sddp.spmodel.stageNumber-1]
    # set JuMP model inside SDDP object
    sddp.solverinterface = ms
end


"""Compute dual upper bound with `max_p <x0, p> - D_0(p)`."""
function dualbound(sddpdual, x)
    V = sddpdual.bellmanfunctions[1]

    m = Model(solver=sddpdual.params.SOLVER)

    @variable(m, θ)
    p = @variable(m, -LIPSCHITZ <= p[1:sddpdual.spmodel.dimStates] <= 0)

    for i in 1:V.numCuts
        lambda = vec(V.lambdas[i, :])
        @constraint(m, V.betas[i] + dot(lambda, p) <= θ)
    end

    @objective(m, Max, dot(x, p) - θ)
    status = solve(m)

    # if solution is not optimal, raise an error
    if status != :Optimal
        println(m)
        println(getvalue(p))
        error("Starting point not found in dual SDDP.")
    end
    return getobjectivevalue(m), getvalue(p)
end


"""Change initial co-state in interface of dual SDDP."""
function updateinitialstate!(sddpdual, x)
    ub, p0 = dualbound(sddpdual, x)
    sddpdual.spmodel.initialState = p0
    return ub
end


"""Get lowerbound if NaN"""
function getlowerbound(sddpdual, p)

    Vt = sddpdual.bellmanfunctions[1]

    m = Model(solver=sddpdual.params.SOLVER)
    @variable(m, alpha)

    for i in 1:Vt.numCuts
        lambda = vec(Vt.lambdas[i, :])
        @constraint(m, Vt.betas[i] + dot(lambda, p) <= alpha)
    end

    @objective(m, Min, alpha)
    solve(m)
    return getvalue(alpha)
end


function displayit(it, lb)
    print("Pass n\° ", it)
    @printf("\tLower-bound: %.4e \n", lb)
end
