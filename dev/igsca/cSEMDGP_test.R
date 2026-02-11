# Test IGSCA
library(cSEM.DGP)
library(lavaan)
# library(cSEM)
devtools::load_all()


# General Population Values ----------------------------------------------


# TODO: This assumes canonical weights. How so?

# These are the actual weights used to construct the population correlation matrix,
# from which the data is sampled
Sxi1 <- matrix(c(1, 0.4, -0.3,
                 0.4, 1, 0.4,
                 -0.3, 0.4, 1),
                  3, 3)
w1 <- c(0.4, 0.3, 0.2)

w1 %*% Sxi1 %*% w1

cmp_canonical_weights <- w1 / c(sqrt(w1 %*% Sxi1 %*% w1))

# GSCA -------------------------------------------------------------------
modelpop <- ' 
xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
x11~~0.4*x12 + -0.3*x13
x12~~0.4*x13


xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
x21~~0.4*x22 + -0.3*x23
x22~~0.4*x23



xi2~ 0.3*xi1

'
datapop=generateData(.model = modelpop,.empirical = T)

modelest <- ' 
xi1<~x11 + x12 + x13

xi2<~x21 + x22 + x23


xi2~ xi1
'

out_ccmp = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .disattenuate = FALSE,
    .GSCA_modes = "CCMP"
)
tidy(out_ccmp)

out_ncmp = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .disattenuate = FALSE,
    .GSCA_modes = "NCMP"
)
tidy(out_ncmp)


# GSCA M -----------------------------------------------------------------
modelpop <- ' 
xi1=~0.4*x11 + 0.3*x12 + 0.2*x13

xi2=~0.4*x21 + 0.3*x22 + 0.2*x23

xi2~ 0.3*xi1

'

datapop=generateData(.model = modelpop,.empirical = T)

modelest <- ' 
xi1=~x11 + x12 + x13

xi2=~x21 + x22 + x23


xi2~ xi1
'

out = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute"
)
tidy(out)


# IGSCA ------------------------------------------------------------------
modelpop <- ' 
xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
x11~~0.4*x12 + -0.3*x13
x12~~0.4*x13


xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
x21~~0.4*x22 + -0.3*x23
x22~~0.4*x23

xi3=~0.8*x31 + 0.6*x32 + 0.7*x33

xi3~ .5*xi1 + .5*xi2

xi1~~0.3*xi2
'

datapop = generateData(.model = modelpop, .empirical = T)

modelest <- ' 
xi1<~x11 + x12 + x13

xi2<~x21 + x22 + x23
xi3 =~ x31 + x32 + x33

xi3~ xi2+xi1
'

out = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute"
)
tidy(out)


# TODO: Add single indicator composite
# TODO: Add single indicator factor


# Mixture composite and common factors
modelpop <- ' 
xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
x11~~0.4*x12 + -0.3*x13
x12~~0.4*x13


xi2=~0.8*x21 + 0.6*x22 + 0.7*x23

xi3 =~ 1*x31

xi4 <~ 1*x41


xi3~ -0.3*xi1 + 0.4*xi2
xi1~~0.3*xi2
xi4 ~ 0.3*xi3
'

datapop = generateData(.model = modelpop, .empirical = T)


# library(cSEM)
modelest <- ' 
xi1<~x11 + x12 + x13

xi2=~x21 + x22 + x23

xi3 =~ x31

xi4 <~ x41


xi3~ xi1 + xi2
xi4 ~ xi3

'


# debugonce(igsca)
out = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .tolerance = 0.00001,
    .conv_criterion = "sum_diff_absolute"
)
summarize(out)

csa <- lavaan::sem(
    modelest,
    datapop,
    optim.gradient = "numerical"
)
broom::tidy(csa)

out1=csem(.data = datapop,.model = modelest,.approach_weights = 'PLS-PM', .disattenuate = TRUE)
summarize(out1)


# Other tests ------------------------------------------------------------

# Common factors only
modelpop <- ' 
xi1=~0.4*x11 + 0.3*x12 + 0.2*x13



xi2=~0.8*x21 + 0.6*x22 + 0.7*x23

xi3 =~ 1*x31

xi4 =~ 1*x41


xi3~ -0.3*xi1 + 0.4*xi2
xi1~~0.3*xi2
xi4 ~ 0.3*xi3
'

datapop = generateData(.model = modelpop, .empirical = T)

out = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .tolerance = 0.0001,
    .conv_criterion = "sum_diff_absolute"
)
summarize(out)

# Mixture composite and common factors
modelpop <- ' 
xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
x11~~0.4*x12 + -0.3*x13
x12~~0.4*x13


xi2<~0.8*x21 + 0.6*x22 + 0.7*x23

x21~~-.3*x22 + 0.2*x23
x22~~-0.3*x23

xi3 <~ 1*x31

xi4 <~ 1*x41


xi3~ -0.3*xi1 + 0.4*xi2
xi1~~0.3*xi2
xi4 ~ 0.3*xi3
'

datapop=generateData(.model = modelpop,.empirical = T)


out=csem(.data = datapop,.model = modelest,.approach_weights = 'GSCA')
summarize(out)

# --- without single indicator constructs
modelpop <- ' 
xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
x11~~0.4*x12 + -0.3*x13
x12~~0.4*x13


xi2=~0.8*x21 + 0.6*x22 + 0.7*x23

xi3 =~ 0.7*x31 + 0.8*x32



xi3~ -0.3*xi1 + 0.4*xi2
xi1~~0.3*xi2

'

datapop=generateData(.model = modelpop,.empirical = T)

modelest <- ' 
xi1<~x11 + x12 + x13

xi2=~x21 + x22 + x23

xi3 =~ x31 + x32


xi3~ xi1 + xi2


'

out=csem(.data = datapop,.model = modelest,.approach_weights = 'GSCA', .disattenuate = TRUE)
summarize(out)




