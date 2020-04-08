\dontrun{
  # Example from Henseler (2020) about importance performance matrix analysis 
  # using the SQ dataset
  dat=cbind(SQ,SQ[,1:22],SQ[,1:22])
  # Give these indicators names
  colnames(dat)=c(colnames(SQ),paste0(colnames(SQ)[1:22],'rep'),paste0(colnames(SQ)[1:22],'reprep'))

  # Specify and estimate second-order model using the repeated indicators approach
    model<-'
  A01 <~ a01
  A02 <~ a02
  A03 <~ a03
  A04 <~ a04
  A05 <~ a05
  A06 <~ a06
  A07 <~ a07
  A08 <~ a08
  A09 <~ a09
  A10 <~ a10
  A11 <~ a11
  A12 <~ a12
  A13 <~ a13
  A14 <~ a14
  A15 <~ a15
  A16 <~ a16
  A17 <~ a17
  A18 <~ a18
  A19 <~ a19
  A20 <~ a20
  A21 <~ a21
  A22 <~ a22

  Tangibles <~ a01rep+a02rep+a03rep+a04rep
  Reliability<~ a05rep+a06rep+a07rep+a08rep+a09rep
  Responsiveness <~ a10rep+a11rep+a12rep+a13rep
  Assurance<~a14rep+a15rep+a16rep+a17rep
  Empathy<~a18rep+a19rep+a20rep+a21rep+a22rep
  Sat <~ sat
  SerQ<~a01reprep+a02reprep+a03reprep+a04reprep+a05reprep+a06reprep+a07reprep+a08reprep+a09reprep+a10reprep+a11reprep+a12reprep+a13reprep+a14reprep+a15reprep+a16reprep+a17reprep+a18reprep+a19reprep+a20reprep+a21reprep+a22reprep

  SerQ~Tangibles+Reliability+Responsiveness+Assurance+Empathy
  Sat~SerQ

  Tangibles~A01+A02+A03+A04
  Reliability~A05+A06+A07+A08+A09
  Responsiveness~A10+A11+A12+A13
  Assurance~A14+A15+A16+A17
  Empathy~A18+A19+A20+A21+A22
  '
    # Estimate the model
    out <- csem(.data = dat,.model = model,
                .PLS_weight_scheme_inner = 'factorial',
                .tolerance = 1e-06)
    
    # Apply doIPMA function to obtain the neccesary outcome to 
    # plot the importance performance matrix
    outIPA <- doIPMA(out)
    
    plot(x = outIPA,.dependent = 'Sat',.level = 'construct',
         .attributes = c("Tangibles","Reliability","Responsiveness", "Assurance","Empathy"))
}