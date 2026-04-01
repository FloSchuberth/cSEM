library(cSEM.DGP)

set.seed(3883)

uniF_uniF = 'xi1 =~ 1*x11

               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23

               xi2 ~ .5*xi1'

uniF_uniF_mod = 'xi1 =~ x11

               xi2=~x21 + x22 + x23

               xi2 ~ xi1'

csem(
  .data = cSEM.DGP::generateData(uniF_uniF, .empirical = TRUE),
  .model = uniF_uniF_mod,
  .approach_weights = 'GSCA',
  .disattenuate = TRUE,
  .conv_criterion = "sum_diff_absolute"
)
