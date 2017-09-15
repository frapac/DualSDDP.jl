
"""Build matrix corresponding to the damsvalley"""
function getmatrix(nzones, narcs, R)
    I = eye(nzones)
    O = zeros(Float64, nzones, nzones)
    Oq = zeros(Float64, nzones, narcs)
    i = ones(Float64, nzones)
    o = zeros(Float64, nzones)

    A = I
    # u = [uturb; uspill, utherm, urec1, urec2, q]
    B = - [I I O O O Oq]
    C = [I O]


    c = [o; o; ctherm*i; cpenal*i; cpenal*i ; zeros(Float64, narcs)]

    D = [O; O; O; O; O; O; O; O; A; -A]
    E = [I O O O O Oq;
        -I O O O O Oq;
         O -I O O O Oq;
         O O I O O Oq;
         O O -I O O Oq;
         O O O -I O Oq;
         O O O O -I Oq;
         I O I I -I R;
         B;
         -B]

    Gt = [O I]
    G = [uturb*i; o; o; utherm*i; o; o; o; o; o; xmax; o]


    return A, B, D, E, c
end
