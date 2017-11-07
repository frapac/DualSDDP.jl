

"""Display control."""
function dispu(up, nscen, nzones)
    colors = ["k", "r", "b", "g", "c", "m", "y"]
    figure()
    for (ind, i) in enumerate(up*nzones+1:(up+1)*nzones)
        plot(u[:, 1:nscen, i], c=colors[ind])
    end
end


"""Display state."""
function dispx(nscen, nzones)
    colors = ["k", "r", "b", "g", "c", "m", "y"]
    figure()
    for (ind, i) in enumerate(1:nzones)
        plot(x[:, 1:nscen, i], c=colors[ind])
    end
end


"""Display convergence of SDDP in primal and in dual."""
function dispconv(; nscen=1000, ptol=.95, delta=UPPER_BOUND)
    ΔI= delta
    tol = √2 * erfinv(2*ptol - 1)
    tol2 = √2 * erfinv(2*.999 - 1)
    ubf = sddpprimal.stats.upper_bounds
    lprim = sddpprimal.stats.lower_bounds

    figure(figsize=(9,6))

    plot(lbdual, c="darkblue", lw=4, label="Dual UB")
    plot(lprim, c="darkred", lw=4, label="Primal LB")
    #= plot(ubf, c="k", lw=.1, label="Forward primal cost") =#


    plot(ΔI:ΔI:MAXIT, ubp, color="k", lw=1.5, label="Outer strat.", marker="s")
    plot(ΔI:ΔI:MAXIT, ubp + tol*stdp/√nscen, color="k", lw=1, linestyle="--")
    plot(ΔI:ΔI:MAXIT, ubp - tol*stdp/√nscen, color="k", lw=1, linestyle="--")
    #= plot(ΔI:ΔI:MAXIT, ubp + tol2*stdp/√nscen, color="k", lw=1, linestyle=":") =#
    #= plot(ΔI:ΔI:MAXIT, ubp - tol2*stdp/√nscen, color="k", lw=1, linestyle=":") =#
    fill_between(ΔI:ΔI:MAXIT, ubp - tol*stdp/√nscen, ubp + tol*stdp/√nscen,
                 color="grey", alpha=.4, label="Confidence ($(100*ptol)%)")
    #= fill_between(ΔI:ΔI:MAXIT, ubp - tol2*stdp/√nscen, ubp + tol2*stdp/√nscen, =#
    #=              color="grey", alpha=.1, label="Confidence ($(100*.999)%)") =#


    legend(loc=4)
    xlabel("Iterations")
    ylim(.99lbprimal, 1.01lbprimal)
    tight_layout()
end


"""Display evolution of inner strategy."""
function dispubd(;nscen=1000)
    tol = √2 * erfinv(2*.95 - 1)
    plot(ΔMC:ΔMC:MAXIT, ubd, color="olive", marker="o", linestyle="--",
        label="Inner strat.")
    plot(ΔMC:ΔMC:MAXIT, ubd + tol*stdd/√nscen, color="darkgreen", alpha=.5, lw=1, linestyle=":")
    plot(ΔMC:ΔMC:MAXIT, ubd - tol*stdd/√nscen, color="darkgreen", alpha=.5, lw=1, linestyle=":")
    #= fill_between(ΔMC:ΔMC:MAXIT, ubd - tol*stdd/√nscen, ubd + tol*stdd/√nscen, =#
    #=              color="green", alpha=.2, label="Confidence ($(100*.95)%)") =#
end


function dispres(res, di)
    nbit = size(res, 1)
    lbdual = res[:, 2]
    lprim = res[:, 1]
    lbprimal = lprim[end]
    ubp = res[di:di:nbit, 6]
    ubd = res[di:di:nbit, 7]

    figure(figsize=(9,6))

    plot(lbdual, c="darkblue", lw=4, label="Dual UB")
    plot(lprim, c="darkred", lw=4, label="Primal LB")

    plot(di:di:nbit, ubp, color="k", lw=1.5, label="Outer strat.", marker="s")
    plot(di:di:nbit, ubd, color="olive", marker="o", linestyle="--",
        label="Inner strat.")

    legend(loc=4)
    xlabel("Iterations")
    ylim(.99lbprimal, 1.01lbprimal)
    tight_layout()
end
