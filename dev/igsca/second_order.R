model <- "
# Path model / Regressions 
c4   ~ eta1
eta2 ~ eta1 + c4

# Reflective measurement model
c1   <~ y11 + y12 
c2   <~ y21 + y22 + y23 + y24
c3   <~ y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38
eta1 =~ y41 + y42 + y43
eta2 =~ y51 + y52 + y53

# Composite model (second order)
c4   =~ c1 + c2 + c3
"

# debugonce(csem)
# debugonce(foreman)

res_2stage <- csem(dgp_2ndorder_cf_of_c, model, .approach_2ndorder = "2stage", .approach_weights = "GSCA")