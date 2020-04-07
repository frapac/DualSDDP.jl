

using DualSDDP
using Test

@testset "SDDP dual" begin
    maxit = 50
    mpts = DualSDDP.MPTS([:FRA, :GER], 10, 5)
    sddpprimal = DualSDDP.initprimal(mpts)
    ubp, stdp = DualSDDP.runprimal!(sddpprimal, maxiterations=maxit)
    lb = sddpprimal.stats.lowerbound
    @testset "Vanilla dual" begin
        sddpdual = DualSDDP.initdual(mpts, sddpprimal)
        lbdual, timedual = DualSDDP.rundual!(sddpdual, sddpprimal,
                                            maxiterations=maxit)
        ub = lbdual[end]
        @test lb <= ub
        # Test that gap is below 1%
        @test (ub - lb) / max(1.0, abs(lb)) <= 0.01
    end

    @testset "Joint primal/dual" begin
        sddpdual = DualSDDP.initdual(mpts, sddpprimal)
        lbdual, timedual, ubp, stdp = DualSDDP.runjoint!(sddpprimal, sddpdual, maxiterations=maxit)
        ub = lbdual[end]
        @test lb <= ub
        # Test that gap is below 1%
        @test (ub - lb) / max(1.0, abs(lb)) <= 0.01
    end
end
