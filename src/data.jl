################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# Data for production-transport problem
################################################################################


const DATA = Dict(
            "ORDER"=> [:BEL, :ESP, :FRA, :GER, :ITA, :PT, :SUI, :UK],
            # states lower and upper bounds for each countries (in GWh)
            "XMIN"=> 1000*Float64[0, 0, 0,   0, 0, 0,   0,   0],
            "XMAX"=> 1000*Float64[0, 8, 10,  4, 4, 2,   10,  2],
            # initial position of dams in each countries
            "XINI"=> 1000*Float64[0, 6, 7.5, 3, 3, 1.5, 7.5, 1.5],
            # maximum turbinate (in GWh)
            "UMAX"=> 30*24*Float64[0, 9, 12,  3, 10.5, 1.5, 12,  1.5],
            # maximum thermal production (in GWh)
            "PMAX"=> 30*24*Float64[40, 40, 100, 100, 40, 20, 20, 60],
            # piecewise linear cost
            "czpl"=> Float64[10 40 70 90;
                             10 40 70 90;
                             5  15 30 45;
                             10 25 35 50;
                             10 40 70 90;
                             20 100 100 100;
                             20 100 100 100;
                             20 40 60 80],
            # graph's incidence matrix
            "INCIDENCE"=>30*24*Float64[0 0 2 1 0 0 0 0;
                                       0 0 1 0 0 1 0 0;
                                       2 1 0 3 1 0 2 1;
                                       1 0 3 0 0 0 1 0;
                                       0 0 1 0 0 0 1 0;
                                       0 1 0 0 0 0 0 0;
                                       0 0 2 1 1 0 0 0;
                                       0 0 1 0 0 0 0 0],
            "COST_HF" => 3000,
            "CPENAL" =>3000,
            "CTRANS" => 1
                )


getpos(name::Symbol)=findfirst(DATA["ORDER"] .== name)

"Load data corresponding to specified countries."
function getglobalparams(names)
    nnodes = length(names)

    pos = getpos.(names)

    xini = DATA["XINI"][pos]

    # get number of controls
    nc = 4*nnodes
    nv = 1

    # NETWORK TOPOLOGY
    A, fexch = buildincidence(DATA["INCIDENCE"][pos, pos])
    # bounds on the state (here: dam)
    s_bounds = [DATA["XMAX"][p] for p in pos]
    # bounds on controls
    uturb = [DATA["UMAX"][p] for p in pos]
    utherm = [DATA["PMAX"][p] for p in pos]
    # thermal cost
    ctherm = DATA["czpl"][pos, 1]

    return s_bounds, uturb, utherm, xini, A, ctherm, fexch
end
