
# Lavaan -----------------------------------------------------------------
library(lavaan)
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=TRUE,
              auto.cov.lv.x=TRUE)

summary(fit, fit.measures=TRUE)



# GSCAm ------------------------------------------------------------------
HS.model_csem <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 
              visual ~~ textual
              speed ~~ textual
              visual ~~ speed'


HolzingerSwineford1939 <- lavaan::HolzingerSwineford1939
debugonce(calculateWeightsGSCAm)
debugonce(checkConvergence)

sum_diff_absolute_mod <- cSEM::csem(
  .data = HolzingerSwineford1939,
  .model = HS.model_csem,
  .approach_weights = "GSCA",
  .disattenuate = TRUE,
  .conv_criterion = 'sum_diff_absolute'
)

mean_diff_absolute_mod <- cSEM::csem(
  .data = HolzingerSwineford1939,
  .model = HS.model_csem,
  .approach_weights = "GSCA",
  .disattenuate = TRUE,
  .conv_criterion = 'mean_diff_absolute'
)

waldo::compare(sum_diff_absolute_mod, mean_diff_absolute_mod)

cSEM::csem(
  .data = HolzingerSwineford1939,
  .model = HS.model_csem,
  .approach_weights = "GSCA",
  .disattenuate = TRUE
)