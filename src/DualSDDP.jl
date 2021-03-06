################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################

module DualSDDP

using Printf, DelimitedFiles, Random
using Statistics, LinearAlgebra
# optimization packages
using JuMP, StochDynamicProgramming
# use Gurobi as LP solver
using Xpress
#using Gurobi
# import Clustering for kmeans quantization of noises
using Clustering

const SDDP = StochDynamicProgramming

# utilities
include("utils.jl")
include("dualutils.jl")
# problem's data
include("data.jl")
# definition of noises process
include("noises.jl")
# production transport problem
include("mpts.jl")
# SDDP code
include("sddp.jl")
# inner approximation
include("innerapprox.jl")
# Philpott's lower bound
include("dp.jl")

end
