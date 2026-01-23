\dontrun{
model_Bergami_int="
  # Common factor and composite models
  OrgPres <~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
  OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
  AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
  AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

  # Structural model 
  OrgIden ~ OrgPres 
  AffLove ~ OrgPres+OrgIden+OrgPres.OrgIden
  AffJoy  ~ OrgPres+OrgIden
  "
  
  outBergamiInt <- csem(.data = BergamiBagozzi2000,.model = model_Bergami_int,
                        .disattenuate = T,
                        .PLS_weight_scheme_inner = 'factorial',
                        .tolerance = 1e-6,
                        .resample_method = 'none')
  
  outPlot <- plot(outBergamiInt)
  outPlot
  savePlot(outPlot,.file='plot.pdf')
  savePlot(outPlot,.file='plot.png')
  savePlot(outPlot,.file='plot.svg')
  savePlot(outPlot,.file='plot.dot')
}
