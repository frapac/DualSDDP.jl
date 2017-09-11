
function build_inner_approx!(m, V)
    # m stores the initial model
    xf = m[:xf]
    alpha = m[:alpha]
    ncuts = V.numCuts

    # define simplex Λ
    @variable(m, eta[1:ncuts] >= 0.)
    @constraint(m, sum(eta) == 1.)

    # we build the inner approximation all in once
    @constraint(m, alpha == -sum(eta[i]*V.betas[i] for i in 1:ncuts))
    @constraint(m, xf .== sum(eta[i]*V.lambdas[i, :] for i in 1:ncuts))
end


function build_joint_approx!(m, Vdual, Vprimal, ω)
    # m stores the initial model
    xf = m[:xf]
    alpha = m[:alpha]

    # inner cost
    @variable(m, αi)
    # outer cost
    @variable(m, αo)

    # add inner approx
    ncuts = Vdual.numCuts

    # define simplex Λ
    @variable(m, eta[1:ncuts] >= 0.)
    @constraint(m, sum(eta) == 1.)

    # we build the inner approximation all in once
    @constraint(m, αi == -sum(eta[i]*Vdual.betas[i] for i in 1:ncuts))
    @constraint(m, xf .== sum(eta[i]*Vdual.lambdas[i, :] for i in 1:ncuts))

    # add outer approx
    for idx in 1:Vprimal.numCuts
        beta = Vprimal.betas[idx]
        lambda = Vprimal.lambdas[idx, :]
        @constraint(m, beta + dot(lambda, xf) <= αo)
    end

    # mix outer and inner approx.
    @constraint(m, alpha == ω*αi + (1 - ω)*αo)
end

function init_innermodeler!(sddp, Vs)
    for t in 1:sddp.spmodel.stageNumber - 2
        build_inner_approx!(sddp.solverinterface[t], Vs[t+1])
    end
end

function init_jointmodeler!(sddp, Vd, Vp, ω)
    for t in 1:sddp.spmodel.stageNumber - 2
        build_joint_approx!(sddp.solverinterface[t], Vd[t+1], Vp[t+1], ω)
    end
end
