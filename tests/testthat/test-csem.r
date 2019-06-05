context("csem")

model <- "
# Structural model
EXPE ~ IMAG

# Measurement model

IMAG <~ imag1 + imag2
EXPE <~ expe1 + expe2 + expe3
"

test_that("No data provided should give an error", {
  expect_error(
    csem(.model = model)
  )
})