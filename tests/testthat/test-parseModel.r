context("csem")

### Linear models ==============================================================
## 1. One exogenous and one endogenous construct -------------------------------

# 1.1 Construct of the measurement model forgotten/misspelled
# 1.1 Everything is correctly specified
model1 <- "
# Structural model
EXPE ~ IMAG

# Composite model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"


# 1.2 One measurement equation is redundant
model2 <- "
# Structural model
EXPE ~ IMAG

# Composite model
EXPE <~ expe1 + expe2

# Measurement model
IMAG =~ imag1 + imag4
QUAL =~ qual1 + qual3
"

# 1.3. At least one indicator is connected to several constructs
model3 <- "
# Structural model
EXPE ~ IMAG

# Composite model
EXPE <~ expe1 + expe2

# Measurement model
IMAG =~ imag1 + imag4 + expe1
"
# 1.4 Construct of the measurement model forgotten/misspelled
model4 <- "
# Structural model
EXPE ~ IMAG

# Composite model
EXPE <~ expe1 + expe2
"

# 1.5 Measurement errors across blocks 
model5 <- "
# Structural model
EXPE ~ IMAG

# Composite model
EXPE <~ expe1 + expe2
IMAG =~ imag1 + imag2

imag1 ~~ expe1
"

## Tests
test_that("Linear model: incorrectly specified models provide an error", {
  
  expect_error(parseModel(model2))
  expect_error(parseModel(model3))
  expect_error(parseModel(model4))
  expect_error(parseModel(model5))
})

test_that("Linear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model1), "cSEMModel")
  expect_output(str(parseModel(model1)), "List of 13")
})

## 2. Several endogenous and exogenous constructs ------------------------------
# Including
# - Exogenous construct correlations
# - Measurement error correlations

model <- "
# Structural model
QUAL ~ IMAG + VAL + SAT
VAL  ~ IMAG + EXPE
SAT  ~ EXPE

# Composite Model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2

# Measurement model
SAT  =~ sat1 + sat4
QUAL =~ qual1 + qual2 + qual3
VAL  =~ val3 + val4

# Measurement correlation
qual1 ~~ qual2

# Construct correlation
IMAG ~~ EXPE
"

## Tests
test_that("Linear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model), "cSEMModel")
  expect_output(str(parseModel(model)), "List of 13")
})

## 3. Only measurement model (no structural model) -----------------------------
# 3.1 Only one measurement equation
model1 <- "
# Composite model
IMAG <~ imag1 + imag2 + imag3
"

# 3.2 Only measurement equations without construct correlations given
model2 <- "
# Measurement and composite model
EXPE =~ expe1 + expe2
IMAG <~ imag1 + imag2
"

# 3.3 Two measurement equations and correlation given
model3 <- "
# Construct correlations
EXPE ~~ IMAG

# Measurement and composite model
EXPE =~ expe1 + expe2
IMAG <~ imag1 + imag2
"

test_that("Linear model: incorrectly specified models provide an error", {
  expect_error(parseModel(model1))
  expect_error(parseModel(model2))
})

test_that("Linear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model3), "cSEMModel")
  expect_output(str(parseModel(model3)), "List of 13")
})

### Nonlinear models ===========================================================
## 1. One exogenous and one endogenous construct -------------------------------
## 1.1 Everything correctly specified
model1 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"

## 1.2 An interaction term does not appears individually in the structural model
model2 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG + QUAL.QUAL

# Measurement model
EXPE <~ expe1 + expe2
IMAG =~ imag1 + imag4
QUAL =~ qual1 + qual3
"

## 1.3 Interaction term as a dependent variable in a measurement equation
model3 <- "
# Structural model
QUAL ~ IMAG + IMAG.IMAG + EXPE + EXPE.EXPE + EXPE.IMAG

# Measurement model
EXPE      <~ expe1 + expe2
IMAG      =~ imag1 + imag2
IMAG.IMAG =~ imag3 + imag4
QUAL      =~ qual1 + qual2
"

## 1.4 Interaction term as a dependent variable in a structural equation
model4 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG + QUAL.QUAL
QUAL.QUAL ~ IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG =~ imag1 + imag4
QUAL =~ qual1 + qual3
"

# 1.5 Construct of the measurement model forgotten/misspelled
model5 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Measurement model
EXPE <~ expe1 + expe2
"

# 1.6 Correlation between exogenous construct and interaction term is specified
model6 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

EXPE ~~ IMAG.IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"

## Tests

test_that("Nonlinear model: incorrectly specified models provide an error", {
  
  expect_error(parseModel(model2))
  expect_error(parseModel(model3))
  expect_error(parseModel(model4))
  expect_error(parseModel(model5))
  expect_error(parseModel(model6))
})

## Tests
test_that("Nonlinear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model1), "cSEMModel")
  expect_output(str(parseModel(model1)), "List of 13")
})

## 2. Several endogenous and exogenous constructs ------------------------------
model <- "
# Structural model
QUAL ~ IMAG + VAL + SAT + SAT.SAT
VAL  ~ IMAG + EXPE + EXPE.IMAG
SAT  ~ EXPE

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
SAT  =~ sat1 + sat4
QUAL =~ qual1 + qual2
VAL  =~ val3 + val4
"

## Tests
test_that("Nonlinear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model), "cSEMModel")
  expect_output(str(parseModel(model)), "List of 13")
})

### Second-order model =========================================================
## 1. One second order construct -----------------------------------------------
# 1.1 One measurement equation for a construct used to define the second order 
#     construct got forgotten but appears in the structural model (IMAG)

model1 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT + IMAG

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
QUAL =~ IMAG + EXPE
"

# 1.2 One measurement equation for a construct used to define the second order 
#     construct got forgotten and does not (!) appear in the structural model

model2 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
QUAL =~ IMAG + EXPE
" # not detected by parseModel() but csem() throws an error

# 1.3. One second-order construct has also indicators attached
model3 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
QUAL =~ IMAG + EXPE + imag3
IMAG <~ imag1 + imag2
" # not detected by parseModel() but csem() throws an error

# 1.4 Second order construct appears as an interaction term with itself 
#     in the structural model

model4 <- "
# Structural model
SAT ~ QUAL + QUAL.QUAL
VAL ~ SAT + SAT.SAT

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
QUAL =~ IMAG + EXPE
IMAG <~ imag1 + imag2
"

# 1.5 Second order construct appears as an interaction term with another 
#     first order constructs in the structural model

model5 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT
LOY ~ QUAL + VAL + QUAL.VAL

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
LOY =~ loy1 + loy2
QUAL =~ IMAG + EXPE
IMAG <~ imag1 + imag2
"

# 1.6 "Normal" second order model including a nonlinear term

model6 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT

# Measurement model
EXPE <~ expe1 + expe2

SAT =~ sat1 + sat2
VAL =~ val1 + val2
QUAL =~ IMAG + EXPE
IMAG <~ imag1 + imag2
"

## Tests

test_that("Second-order model: incorrectly specified models provide an error", {
  
  expect_error(parseModel(model1))
  expect_error(csem(satisfaction, model2, .normality = FALSE))
  expect_error(csem(satisfaction, model3, .normality = FALSE))
})

## Tests
test_that("Second-order model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model4), "cSEMModel")
  expect_output(str(parseModel(model4)), "List of 13")

  expect_s3_class(parseModel(model5), "cSEMModel")
  expect_output(str(parseModel(model5)), "List of 13")

  expect_s3_class(parseModel(model6), "cSEMModel")
  expect_output(str(parseModel(model6)), "List of 13")
})


## 2. Two second-order construct-----------------------------------------------
## 2.1 Everything correctly specified
model1 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT + LOY

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
HELP <~ val3 + val4
QUAL =~ IMAG + EXPE
LOY =~ HELP
IMAG <~ imag1 + imag2
"

## 2.2 Two second-order construct with at least one 2nd order construct attatched
##     to a nonlinear term
model2 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT + LOY + LOY.SAT

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
HELP <~ val3 + val4
QUAL =~ IMAG + EXPE + IMAG.IMAG + EXPE.IMAG
LOY =~ HELP
IMAG <~ imag1 + imag2
"

# 2.3 Construct order > 2 

model3 <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + SAT.SAT + LOY

# Measurement model
EXPE <~ expe1 + expe2
SAT =~ sat1 + sat2
VAL =~ val1 + val2
HELP =~ val3 + val4
QUAL =~ IMAG + EXPE + LOY
LOY =~ HELP
IMAG <~ imag1 + imag2
"

## Tests

test_that("Second-order model: incorrectly specified models provide an error", {
  expect_error(parseModel(model2))
  expect_error(parseModel(model3))
})

## Tests
test_that("Second-order model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model1), "cSEMModel")
  expect_output(str(parseModel(model1)), "List of 13")
})

## .check errors argument ======================================================
models <- list(
  'C ~ A', # single path
  "B =~ x21", # single loading
  "D <~ x41", # single weight
  "A ~ C", # wrong single path, both constructs exist
  "D ~ E", # wrong path as construct E does not exist
  "B =~ x31",# wrongly assigned indicator
  "D <~ x51" # Indicator that does not exist
  )

for(i in seq_along(models)) {
  test_that(paste0("Model: ", models[[i]], " does not throw an error when .check_errors = FALSE."), {
    expect_s3_class(parseModel(models[[i]], .check_errors = FALSE), "cSEMModel")
    expect_output(str(parseModel(models[[i]], .check_errors = FALSE)), "List of 13")
  })
}
