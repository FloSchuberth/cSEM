library(cSEM)
model_Bergami_Bagozzi_Hwang = "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model 
OrgIden ~ OrgPres 
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

out_Hwang <- csem(
  .data = BergamiBagozzi2000,
  .model = model_Bergami_Bagozzi_Hwang,
  .approach_weights = "GSCA",
  .id = "gender",
  .tolerance = 1e-06
)

assess(
  out_Hwang,
  .quality_criterion = c(
    'dg',
    'dl',
    'dml',
    'df',
    'chi_square',
    'chi_square_df',
    'cfi',
    'gfi',
    'cn',
    'ifi',
    'nfi',
    'nnfi',
    'rmsea',
    'rms_theta',
    'srmr',
    'FIT',
    'FIT_m',
    'FIT_s',
    'gof'
  )
)