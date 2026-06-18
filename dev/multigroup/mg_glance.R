devtools::load_all()
model_GSCA = "
# Measurement models
OrgPres <~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
OrgIden <~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model 
OrgIden ~ OrgPres 
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

model_GSCAM = "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model 
OrgIden ~ OrgPres 
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

model_IGSCA = "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
OrgIden <~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model 
OrgIden ~ OrgPres 
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

out_Hwang_single_GSCAM <- csem(
  .data = BergamiBagozzi2000,
  .model = model_GSCAM,
  .approach_weights = "GSCA",
  .tolerance = 1e-06
)

calculateFIT(out_Hwang_single_GSCAM)
calculateAFIT(out_Hwang_single_GSCAM)
calculateFIT_m(out_Hwang_single_GSCAM)
calculateFIT_s(out_Hwang_single_GSCAM)

out_Hwang_mg_GSCA <- csem(
  .data = BergamiBagozzi2000,
  .model = model_GSCA,
  .approach_weights = "GSCA",
  .id = "gender",
  .tolerance = 1e-06
)

out_Hwang_mg_GSCAM <- csem(
  .data = BergamiBagozzi2000,
  .model = model_GSCAM,
  .approach_weights = "GSCA",
  .id = "gender",
  .tolerance = 1e-06
)

out_Hwang_mg_IGSCA <- csem(
  .data = BergamiBagozzi2000,
  .model = model_IGSCA,
  .approach_weights = "GSCA",
  .id = "gender",
  .tolerance = 1e-06
)

calculateFIT(out_Hwang_mg_GSCA)
calculateFIT(out_Hwang_mg_GSCAM)
calculateFIT(out_Hwang_mg_IGSCA)

calculateFIT_m(out_Hwang_mg_GSCA)
calculateFIT_m(out_Hwang_mg_GSCAM)
calculateFIT_m(out_Hwang_mg_IGSCA)

calculateFIT_s(out_Hwang_mg_GSCA)
calculateFIT_s(out_Hwang_mg_GSCAM)
calculateFIT_s(out_Hwang_mg_IGSCA)