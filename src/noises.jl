################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# SCENARIO TREE DEFINITION
################################################################################

"""Quantize scenarios with KMeans quantization.  """
function optquantiz(scenarios, nbins; scale=true)
    # fix seed for reproductability
    srand(2713)
    ntime, nscenarios, nperturb = size(scenarios)

    output = Vector{NoiseLaw}(ntime)

    for t in 1:ntime
        data = copy(reshape(scenarios[t, :, :], nscenarios, nperturb)')
        # Launch kmeans:
        if scale
            μ    = mean(data, 2)
            σ    = std(data, 2)
            data = data .- μ
            pos  = vec(σ .> 0)
            data[pos, :] = data[pos, :] ./ σ[pos]
        end
        quantiz = kmeans(data, nbins)

        center =  quantiz.centers
        if scale
            center = σ .* quantiz.centers .+ μ
        end
        output[t] = NoiseLaw(center, 1/nscenarios*quantiz.counts)
    end

    return output
end

"Fit global noise laws corresponding to countries specified in `names`."
function fitgloballaw(names, nstages, nbins; nscen=10)

    nzones = length(names)
    scenarios = zeros(Float64, nstages, nscen, 2*nzones)

    ntime = 50

    for (idx, name) in enumerate(names)
        # inflow (in GWh)
        inflow = 24*readcsv("data/scen/$name/inflow.txt")[1:nscen, 1:ntime]
        # demand (in GWh)
        demand = 24*readcsv("data/scen/$name/demands.txt")[1:nscen, 1:ntime]

        nscen = size(inflow, 1)
        _inflow = zeros(nscen, nstages)
        _demand = zeros(nscen, nstages)

        # get monthly average
        for nn in 1:nstages
            s1 = nn
            s2 = nn + 1
            _inflow[:, nn] = sum(inflow[:, s1:s2], 2)
            _demand[:, nn] = sum(demand[:, s1:s2], 2)
        end

        scenarios[:, :, idx] = _inflow'
        scenarios[:, :, idx+nzones] = _demand'
    end

    return optquantiz(scenarios, nbins)
end
