# ==============================================================================
# Generate XML file
# ==============================================================================
#
# PURPOSE
# Generate XML files for BEASTX/2 from alignment
#
# APPROACH
# --------
#
#
# OUTPUTS
# -------
#
# ============================================================================== 

library(beautier)

# NOTE: this doesn't work for MASCOT
create_beast2_input_file(
  "output/simulation/1/random_subsample_150.fasta",
  "output/simulation/1/random_subsample_150.xml",
  clock_model = create_clock_model_strict(),
  tree_prior = create_bd_tree_prior()
)
