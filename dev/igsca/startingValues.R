# Reusing single-group IGSCA weights as starting values for a multigroup fit
#
# For IGSCA (models containing both composites and common factors estimated
# by GSCA), completely specified starting values (i.e., the free weights of
# every construct are given) bypass the internal plain-GSCA initialization
# (see hasCompleteStartingValues() / calculateWeightsIGSCA()). A natural
# source of such values are the weights of a previously fitted single-group
# model.

library(cSEM)

model_igsca <- "
# Composite measurement models
OrgIden <~ ma1 + ma2 + ma3
AffJoy  <~ orgcmt1 + orgcmt2 + orgcmt3

# Common factor measurement model
OrgPres =~ cei1 + cei2 + cei3

# Structural model
OrgIden ~ OrgPres
AffJoy  ~ OrgIden
"

# Fit the single-group model
res_igsca_single <- csem(BergamiBagozzi2000, model_igsca,
                         .approach_weights = "GSCA",
                         .GSCA_modes = "NCMP"
                         )

# Convert the (construct x indicator) weight matrix into a named list of
# named vectors: for each construct, the weights of its indicators
W  <- res_igsca_single$Estimates$Weight_estimates
sv <- lapply(rownames(W), function(construct) {
  w <- W[construct, ]
  w[w != 0]
})
names(sv) <- rownames(W)

# Fit the multigroup model using the single-group weights as starting values
res_igsca_multi <- csem(BergamiBagozzi2000, model_igsca,
                        .approach_weights = "GSCA",
                        .GSCA_modes = "NCMP",
                        .id = "gender",
                        .starting_values = sv
                        )
verify(res_igsca_multi)
