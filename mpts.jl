
"""Build matrix corresponding to the damsvalley"""
function MPTSmatrix(nzones, narcs, R)
    I = eye(nzones)
    O = zeros(Float64, nzones, nzones)
    Oq = zeros(Float64, nzones, narcs)
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)

    # x+ = Ax + Bu + Cw
    # dimA: nx x nx
    A = I
    # u = [uturb; uspill, utherm, urec1, urec2, q]
    # (we get rid of F because we do not use decomposition)
    # dimB: nx x nu
    B = - [I I O O O Oq]
    # dimC: nx x nw
    C = [I O]


    c = [o; o; CTHERM*i; CPENAL*i; CPENAL*i ; zeros(Float64, narcs)]

    # Dx + Eu <= Gw
    # we note nc the number of constraints
    # dimD: nc x nx
    D = [O; O; O; O; O; O;O; O; O; O; A; -A]

    # dimE: nc x nu
    E = [I O O O O Oq; # uturb max
        -I O O O O Oq; # uturb min
         O -I O O O Oq; # uspill min
         O O I O O Oq; # utherm max
         O O -I O O Oq; # utherm min
         O O O -I O Oq; # penal1 min
         O O O O -I Oq; # penal2 min
         O O O O O I; # q max
         O O O O O -I; # q min
         I O I I -I R; # u + Aq = w
         B; # xf max
         -B] # xf min
    # lots of constraints :(
    # nc = 12 * nzones

    Gt = [O I]
    G = [uturb*i; o; o; utherm*i; o; o; o; o; -qmax, qmax; o; xmax; o]


    return A, B, D, E, G, Gt, C, c
end
