
res = zeros(Float64, MAXIT, 9)
gap = lbdual ./ sddpprimal.stats.lower_bounds[1:end] - 1

res[:, 1] = sddpprimal.stats.lower_bounds[1:end]
res[:, 2] = lbdual
res[:, 3] = gap
res[:, 4] = cumsum(sddpprimal.stats.exectime[1:end])
res[:, 5] = cumsum(timedual)
res[ΔMC:ΔMC:MAXIT, 6] = ubp
res[ΔMC:ΔMC:MAXIT, 7] = ubd
res[ΔMC:ΔMC:MAXIT, 8] = stdp
res[ΔMC:ΔMC:MAXIT, 9] = stdd

writecsv("res/conv_$(MAXIT)_$(NODES)_$(NSTAGES).csv", res)
