
function build_inner_approx!(m, V)
    # m stores the initial model
    xf = m[:xf]
    alpha = m[:alpha]
    ncuts = V.numCuts

    # define simplex Λ
    @variable(m, eta[1:ncuts] >= 0.)
    @constraint(m, sum(eta) == 1.)

    # we build the inner approximation all in once
    @constraint(m, alpha == - sum(eta[i]*V.betas[i] for i in 1:ncuts))
    @constraint(m, xf .== sum(eta[i]*V.lambdas[i] for i in 1:ncuts))
end
