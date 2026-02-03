### Without higher order constructs --------------------------------------------
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"
  
# Estimate
out <- csem(threecommonfactors, model)
  
# Check admissibility
verify(out) # ok!

## Examine the structure of a cSEMVerify object
str(verify(out))

### With higher order constructs -----------------------------------------------
# If the model containes higher order constructs both the first and the second-
# stage estimates estimates are checked for admissibility

\dontrun{
require(cSEM.DGP) # download from https://m-e-rademaker.github.io/cSEM.DGP/
  
# Create DGP with 2nd order construct. Loading for indicator y51 is set to 1.1
# to produce a failing first stage model
  
dgp_2ndorder <- "
## Path model / Regressions
eta2 ~ 0.5*eta1
eta3 ~ 0.35*eta1 + 0.4*eta2

## Composite model
eta1 =~ 0.8*y41 + 0.6*y42 + 0.6*y43
eta2 =~ 1.1*y51 + 0.7*y52 + 0.7*y53
c1   =~ 0.8*y11 + 0.4*y12
c2   =~ 0.5*y21 + 0.3*y22

## Higher order composite
eta3 =~ 0.4*c1 + 0.4*c2
"
  
dat <- generateData(dgp_2ndorder) # requires the cSEM.DGP package
out <- csem(dat, .model = dgp_2ndorder)

verify(out) # not ok
}