context("csem")

### Linear models ==============================================================
## 1. One exogenous and one endogenous construct -------------------------------
# 1.1 Construct of the measurement model forgotten/misspelled
model1 <- "
# Structural model
EXPE ~ IMAG

# Measurement model
EXPE <~ expe1 + expe2
"

# 1.2 One measurement equation is redundant
model2 <- "
# Structural model
EXPE ~ IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG =~ imag1 + imag4
QUAL =~ qual1 + qual3
"

# 1.3. At least one indicator is connected to several constructs
model3 <- "
# Structural model
EXPE ~ IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG =~ imag1 + imag4 + expe1
"

# 1.4 Everything is correctly specified
model4 <- "
# Structural model
EXPE ~ IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"

## Tests

test_that("Linear model: incorrectly specified models provide an error", {
  
  expect_error(parseModel(model1))
  expect_error(parseModel(model2))
  expect_error(parseModel(model3))
})

test_that("Linear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model4), "cSEMModel")
  expect_output(str(parseModel(model4)), "List of 6")
  expect_equal(names(parseModel(model4)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
})

## 2. Several endogenous and exogenous constructs -----------------------------
model <- "
# Structural model
QUAL ~ IMAG + VAL + SAT
VAL  ~ IMAG + EXPE
SAT  ~ EXPE

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
SAT  =~ sat1 + sat4
QUAL =~ qual1 + qual2
VAL  =~ val3 + val4
"

## Tests
test_that("Linear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model), "cSEMModel")
  expect_output(str(parseModel(model)), "List of 6")
  expect_equal(names(parseModel(model)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
})

### Nonlinear models ===========================================================
## 1. One exogenous and one endogenous construct -------------------------------
# 1.1 Construct of the measurement model forgotten/misspelled
model1 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Measurement model
EXPE <~ expe1 + expe2
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

## 1.5 Everything correctly specified
model5 <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"

## Tests

test_that("Nonlinear model: incorrectly specified models provide an error", {
  
  expect_error(parseModel(model1))
  expect_error(parseModel(model2))
  expect_error(parseModel(model3))
  expect_error(parseModel(model4))
})

## Tests
test_that("Nonlinear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model5), "cSEMModel")
  expect_output(str(parseModel(model5)), "List of 6")
  expect_equal(names(parseModel(model5)), c("structural", "measurement", 
                                           "error_cor", "construct_type", 
                                           "construct_order","model_type"))
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
  expect_output(str(parseModel(model)), "List of 6")
  expect_equal(names(parseModel(model)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
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
  expect_output(str(parseModel(model4)), "List of 6")
  expect_equal(names(parseModel(model4)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
  expect_s3_class(parseModel(model5), "cSEMModel")
  expect_output(str(parseModel(model5)), "List of 6")
  expect_equal(names(parseModel(model5)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
  expect_s3_class(parseModel(model6), "cSEMModel")
  expect_output(str(parseModel(model6)), "List of 6")
  expect_equal(names(parseModel(model6)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
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
  
  expect_error(parseModel(model3))
})

## Tests
test_that("Second-order model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model1), "cSEMModel")
  expect_output(str(parseModel(model1)), "List of 6")
  expect_equal(names(parseModel(model1)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
  expect_s3_class(parseModel(model2), "cSEMModel")
  expect_output(str(parseModel(model2)), "List of 6")
  expect_equal(names(parseModel(model2)), c("structural", "measurement", 
                                            "error_cor", "construct_type", 
                                            "construct_order","model_type"))
})