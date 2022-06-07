# This code provides an example for the analytical computations of 
# "projection contribution" and "system-level long-term effects" 
# as described in the manuscript "Understanding the variability of 
# pairwise coexistence within multispecies systems" by 
# Deng, Taylor and Saavedra.

### load functions in toolbox
source("toolbox.R")

# interaction matrix of a 3-species system
A <- matrix(data = c(-1,-0.215188447823731,-0.563672454021047,
                     -0.215188447823731,-1,-0.372965112033583,
                     -0.563672454021047,-0.372965112033583,-1),
            nrow = 3, ncol = 3)

# index of pairs, i.e., 3 pairs
index_pairs <- combn(3, 2, simplify = FALSE)

# calculations of projection contributions
# and system-level long-term effects for each pair
for (pair in index_pairs){
  # system-level projection contribution
  projection <- Omega_proj(A, pair)
  print(paste0("The system-level projection contribution to pair {", paste0(pair,collapse = ",") ,"} is ", projection))
  
  # probability of pairwise coexistence in isolation
  feasibility_isolation <- Omega(A[pair, pair])
  print(paste0("The feasibility of pair {", paste0(pair,collapse = ",") ,"} in isolation is ", feasibility_isolation))
  
  # probability of pairwise coexistence in system
  feasibility_system <- Omega_comm_analytical(A, pair)
  print(paste0("The feasibility of pair {", paste0(pair,collapse = ",") ,"} in system is ", feasibility_system))
  
  # system-level long-term effect
  long_term_effect <- feasibility_system/feasibility_isolation
  print(paste0("The system-level long-term effect on pair {", paste0(pair,collapse = ",") ,"} is ", long_term_effect))
}
