
res = zeros(Float64, MAXIT, 6)
gap = lbdual ./ sddpprimal.stats.lower_bounds[1:end] - 1

res[:, 1] = sddpprimal.stats.lower_bounds[1:end]
res[:, 2] = lbdual
res[:, 3] = gap
res[:, 4] = cumsum(sddpprimal.stats.exectime[1:end])
res[:, 5] = cumsum(sddpdual.stats.exectime)
res[:, 6] = sddpprimal.stats.upper_bounds[1:end]

writecsv("res/conv_$(MAXIT)_$(NODES)_$(NSTAGES).csv", res)
