x <- rnorm(1000, 2, 300)
scaled_x <-scale(x)
normed_x <- scaled_x / sqrt(length(x-1))

sum(scaled_x^2) == (length(x) - 1)

all.equal(sum(normed_x^2), 1)


# With help from Claude

## Demonstration: Standardized Parameter Estimates from Different Data Scalings
## This shows how standardized coefficients can be obtained regardless of 
## whether data is raw, standardized (mean=0, sd=1), or normalized (length=1)

set.seed(123)
n <- 100

library(broom)
# Generate raw data
X1 <- rnorm(n, mean = 50, sd = 10)
X2 <- rnorm(n, mean = 100, sd = 20)
y <- 2*X1 + 3*X2 + rnorm(n, mean = 0, sd = 5)

# Create a data frame
data_raw <- data.frame(y = y, X1 = X1, X2 = X2)

# ============================================
# Method 1: Raw data with manual standardization
# ============================================
model_raw <- lm(y ~ X1 + X2, data = data_raw)
coef(model_raw)

# Manually compute standardized coefficients
sd_y <- sd(y)
sd_X1 <- sd(X1)
sd_X2 <- sd(X2)

beta_std_X1 <- coef(model_raw)["X1"] * (sd_X1 / sd_y)
beta_std_X2 <- coef(model_raw)["X2"] * (sd_X2 / sd_y)

cat("\nStandardized coefficients (manual calculation):\n")
cat(sprintf("X1: %.6f\n", beta_std_X1))
cat(sprintf("X2: %.6f\n", beta_std_X2))

# ============================================
# Method 2: Standardized data (mean=0, sd=1)
# ============================================
cat("\n\nMETHOD 2: Standardized Data (mean=0, sd=1)\n")
cat("-------------------------------------------\n")

# Standardize: subtract mean, divide by sd
y_std <- scale(y)
X1_std <- scale(X1)
X2_std <- scale(X2)

data_standardized <- data.frame(y = y_std, X1 = X1_std, X2 = X2_std)

cat("Check standardization:\n")
cat(sprintf("Mean of X1_std: %.10f, SD: %.6f\n", mean(X1_std), sd(X1_std)))
cat(sprintf("Mean of X2_std: %.10f, SD: %.6f\n", mean(X2_std), sd(X2_std)))

model_std <- lm(y ~ X1 + X2, data = data_standardized)
cat("\nCoefficients from standardized data:\n")
# print(coef(model_std))
tidy(model_std)

# ============================================
# Method 3: Normalized data (length=1)
# ============================================
cat("\n\nMETHOD 3: Normalized Data (length=1)\n")
cat("-------------------------------------\n")

# Normalize: divide by L2 norm (length)
# y_norm <- y / sqrt(sum(y^2))
# X1_norm <- X1 / sqrt(sum(X1^2))
# X2_norm <- X2 / sqrt(sum(X2^2))

y_norm <- scale(y) / sqrt(length(y) - 1)
X1_norm <- scale(X1) / sqrt(length(X1) - 1)
X2_norm <- scale(X2) / sqrt(length(X2) - 1)

data_normalized <- data.frame(y = y_norm, X1 = X1_norm, X2 = X2_norm)

cat("Check normalization:\n")
cat(sprintf("Length of X1_norm: %.10f\n", sum(X1_norm^2)))
cat(sprintf("Length of X2_norm: %.10f\n", sum(X2_norm^2)))

model_norm <- lm(y ~ X1 + X2, data = data_normalized)
cat("\nCoefficients from normalized data:\n")
# print(coef(model_norm))
tidy(model_norm)



tidy(model_raw)
tidy(model_std)
tidy(model_norm)



coefs <- replicate(
  10000,
  lm(y ~ X1 + X2, data = data_standardized * runif(1, min = -1000, max = 1000)) |>
    coef()
)

plot(coefs[1,])
plot(coefs[2,])
plot(coefs[3,])

coefs_list <- list()
for (i in seq(10000)) {
  mult_factor <- runif(1, min = -1000, max = 1000)
  if (identical(data_standardized * mult_factor, data_standardized)) {
    stop("Multiplying the data by a factor does not work the way I think it does")
  }
  coef_vec <- lm(y ~ X1 + X2, data = data_standardized * mult_factor) |>
    coef()

  coefs_list[[i]] <- c(coef_vec, "mult_factor" = mult_factor)
}

# As the reader may note, regardless of what factor we multiply the data by, the regression coefficients stay the same.
Reduce(rbind, coefs_list)
