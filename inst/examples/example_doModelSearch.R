# Perform model search for a linear model without second-order constructs 
model_Bergami_Bagozzi_Henseler="
 # Measurement models
 OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
 OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
 AffLove =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
 AffJoy  =~ orgcmt5 + orgcmt8
 Gender  <~ gender
 
 # Structural model 
 OrgIden ~ OrgPres
 AffLove ~ OrgPres + OrgIden + Gender 
 AffJoy  ~ OrgPres + OrgIden + Gender 
 "
 
 out <- csem(.data = BergamiBagozzi2000, 
             .model = model_Bergami_Bagozzi_Henseler,
             .PLS_weight_scheme_inner = 'factorial',
             .tolerance = 1e-06
 )

 
 outSearch <- doModelSearch(.object = out,
               .pop_size = 20,
               .n_generations = 20,
               .prob_mutation = 0.5,
               .prob_crossover = 0.8,
               .fbar = -100000,
               .ms_criterion = 'bic',
               .seed = 1234) 
 outSearch
 
#  Estimate the found model
 outNew <- csem(.data = outSearch$Inputcsem$data,
                .model = outSearch$Inputcsem$model[[1]])

 # summarize(outNew) 
 