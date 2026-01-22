## ----eval=FALSE----------------------------------------------------------
#  csem(.data = data, .model = model, .approach_weights = "GSCA")

## ---- eval=FALSE---------------------------------------------------------
#  Error: The following error occured in the `calculateWeightsGSCAm()` function:
#  GSCAm only applicable to pure common factor models. Use `.disattenuate = FALSE`.

## ----message=FALSE-------------------------------------------------------
library(cSEM)
data(satisfaction)

model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT
VAL  ~ EXPE + QUAL

# Measurement model (pure common factor)

EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY  =~ loy1  + loy2  + loy3  + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT  =~ sat1  + sat2  + sat3  + sat4
VAL  =~ val1  + val2  + val3  + val4
"

## ----warning=FALSE-------------------------------------------------------
results1 <- csem(satisfaction, model, .approach_weights = "GSCA", .disattenuate = TRUE)

## ------------------------------------------------------------------------
results1

## ---- eval=FALSE---------------------------------------------------------
#  assess(results1)

## ------------------------------------------------------------------------
summarize(results1)

## ------------------------------------------------------------------------
verify(results1)

## ------------------------------------------------------------------------
results1$Estimates$Path_estimates

## ------------------------------------------------------------------------
results1$Estimates$Path_estimates["EXPE","IMAG"]

## ------------------------------------------------------------------------
results1$Estimates$Loading_estimates["EXPE", "expe1"]

## ------------------------------------------------------------------------
results2 <- csem(satisfaction, model, .approach_weights = "GSCA", .disattenuate = FALSE)

## ------------------------------------------------------------------------
results2$Estimates$Path_estimates["EXPE","IMAG"]
results2$Estimates$Loading_estimates["EXPE", "expe1"]

## ----message=FALSE-------------------------------------------------------
model2 <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT
VAL  ~ EXPE + QUAL

# Measurement model (common factors and composites)

EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
LOY  =~ loy1  + loy2  + loy3  + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT  <~ sat1  + sat2  + sat3  + sat4
VAL  <~ val1  + val2  + val3  + val4
"

## ----eval=FALSE----------------------------------------------------------
#  csem(satisfaction, model2, .approach_weights = "GSCA", .disattenuate = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  Error: The following error occured in the `calculateWeightsGSCAm()` function:
#  GSCAm only applicable to pure common factor models. Use `.disattenuate = FALSE`.

## ----warning=FALSE-------------------------------------------------------
results4 <- csem(satisfaction, model2, .approach_weights = "GSCA", .disattenuate = FALSE)

## ------------------------------------------------------------------------
summarize(results4)

