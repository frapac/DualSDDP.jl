# Write dual version of problems in problem.jl

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
    p = @variable(m, -10e4 <= p[1:sddpdual.spmodel.dimStates] <= 0)

    for i in 1:V.numCuts
        lambda = vec(V.lambdas[i, :])
        @constraint(m, V.betas[i] + dot(lambda, p) <= θ)
    end

    @objective(m, Min, θ - dot(x, p))
    status = solve(m)
    if status != :Optimal
        println(m)
        println(getvalue(p))
        error("Bang")
    end
    return getobjectivevalue(m), getvalue(p)
end


"""Change initial co-state for dual SDDP."""
function updateinitialstate!(sddpdual, x)
    lb, p0 = dualbound(sddpdual, x)
    sddpdual.spmodel.initialState = p0
    return lb
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


function display(it, lb)
    print("Pass n\° ", it)
    @printf("\tLower-bound: %.4e \n", lb)
end
