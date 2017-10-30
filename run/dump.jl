
res = zeros(Float64, MAXIT, 7)
gap = lbdual ./ sddpprimal.stats.lower_bounds[1:end] - 1

res[:, 1] = sddpprimal.stats.lower_bounds[1:end]
res[:, 2] = lbdual
res[:, 3] = gap
res[:, 4] = cumsum(sddpprimal.stats.exectime[1:end])
res[:, 5] = cumsum(timedual)
res[50:50:MAXIT, 6] = ubp
res[50:50:MAXIT, 7] = ubd

writecsv("res/conv_$(MAXIT)_$(NODES)_$(NSTAGES).csv", res)
