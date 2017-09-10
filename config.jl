################################################################################
################################################################################
# Configuration file
################################################################################
################################################################################
# Specify the valley to consider:
# "dam2" "dam3" "dam4" "dam5" "dam6" "dam8" "damX" "damXII" "damRC1" "damRC2"
const VALLEY = "dam1"

# Directory storing the data:
DATA_SOURCE = "data/$VALLEY"

# Specify if we want to consider linear costs:
const LINEAR_COST = true

################################################################################
#------------------------------
# Model parameters
#------------------------------
################################################################################
# Number of timesteps (as we manage the dams over a year, it is equal to the
# number of months):
const TF = 3

# Capacity of dams:
const STOCK_MAX = 80
const STOCK_MIN = 0

# Target for dam3, dam4, dam5, dam6:
const STOCK0 = 40

# Specify the maximum flow that could be turbined:
const CONTROL_MAX = 40
const CONTROL_MIN = 0

const EFFICIENCY = 0.066::Float64
const COST_CONV = 1000.::Float64
const COST_COEFF = 500.::Float64


if LINEAR_COST
    const EPSILON_U = .0
    const EPSILON_Z = .0
else
    const EPSILON_U = .1
    const EPSILON_Z = .1
end

################################################################################
#------------------------------
# SDDP parameters
#------------------------------
################################################################################
const MAX_ITER = 50
const FORWARD_PASS = 1
# Admissible gap between lower-bound and upper-bound:
const EPSILON = 0.001
# Specify size of Monte-Carlo simulation to estimate upper-bound:
const MONTE_CARLO_SIZE = 1000
const FINAL_MONTE_CARLO_SIZE = 10000
const CONFIDENCE_LEVEL = .975
# Compute upper bound every %% iterations:
const UPPER_BOUND = 0
# Prune cuts every %% iterations:
const PRUNE_CUTS = 100


