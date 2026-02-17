# Test IGSCA
library(cSEM.DGP)
library(lavaan)
library(dplyr)
library(purrr)
# library(cSEM)
devtools::load_all()


# General Population Values ----------------------------------------------

# Note: cSEM.DGP assumes canonical weights

# These are the actual weights used to construct the population correlation matrix,
# from which the data is sampled
Sxi1 <- matrix(c(1, 0.4, -0.3,
                 0.4, 1, 0.4,
                 -0.3, 0.4, 1),
                  3, 3)
w1 <- c(0.4, 0.3, 0.2)

w1 %*% Sxi1 %*% w1

(cmp_canonical_weights <- w1 / c(sqrt(w1 %*% Sxi1 %*% w1)))

(single_indicator_weight <- .4 / sqrt(.4 %*% 1 %*% .4))

# GSCA -------------------------------------------------------------------

gsca_pops <- list(
    uni_uni = 'xi1 <~ .4*x11

               xi2 <~ .3*x21

               xi2 ~ .5*xi1',
    uni_tri = 'xi1 <~ .4*x11

               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23

               xi2 ~ .5*xi1',
    tri_uni = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2 <~ .3*x21
               
               xi2 ~ .5*xi1',
    tri_tri = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23
               
               xi2 ~ 0.5*xi1'
)


gsca_datapop <- lapply(gsca_pops, cSEM.DGP::generateData, .empirical = TRUE)

gsca_model_spec <- list(
    uni_uni = 'xi1 <~ x11

               xi2 <~ x21

               xi2 ~ xi1',
    uni_tri = 'xi1 <~ x11

               xi2<~x21 + x22 + x23

               xi2 ~ xi1',
    tri_uni = 'xi1<~x11 + x12 + x13
               
               xi2 <~ x21
               
               xi2 ~ xi1',
    tri_tri = 'xi1<~x11 + x12 + x13
               
               xi2<~x21 + x22 + x23
               
               xi2 ~ xi1'
)

gsca_mods <- mapply(
    cSEM::csem,
    .data = gsca_datapop,
    .model = gsca_model_spec,
    SIMPLIFY = FALSE,
    .approach_weights = 'GSCA',
    .disattenuate = FALSE,
    .GSCA_modes = "CCMP"
)

tidied_gsca_mods <- lapply(gsca_mods, function(x) {
    tidy(x) |>
        dplyr::filter(op %in%  c('<~', '~')) |>
        dplyr::select(term, estimate)
}) |> 
    list_rbind(names_to = 'mod')


# GSCA M -----------------------------------------------------------------
set.seed(1234)
gscam_pop <- list(
    uni_uni = 'xi1 =~ 1*x11

               xi2 =~ 1*x21

               xi2 ~ .5*xi1',
    uni_tri = 'xi1 =~ 1*x11

               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23

               xi2 ~ .5*xi1',
    tri_uni = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2 =~ 1*x21
               
               xi2 ~ .5*xi1',
    tri_tri = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23
               
               xi2 ~ 0.5*xi1'
)


gscam_datapop <- lapply(gscam_pop, cSEM.DGP::generateData, .empirical = TRUE)

gscam_model_spec <- list(
    uni_uni = 'xi1 =~ x11

               xi2 =~ x21

               xi2 ~ xi1',
    uni_tri = 'xi1 =~ x11

               xi2=~x21 + x22 + x23

               xi2 ~ xi1',
    tri_uni = 'xi1=~x11 + x12 + x13
               
               xi2 =~ x21
               
               xi2 ~ xi1',
    tri_tri = 'xi1=~x11 + x12 + x13
               
               xi2=~x21 + x22 + x23
               
               xi2 ~ xi1'
)

gscam_mods <- mapply(
    cSEM::csem,
    .data = gscam_datapop,
    .model = gscam_model_spec,
    SIMPLIFY = FALSE,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute"
)

tidied_gscam_mods <- lapply(gscam_mods, function(x) {
    tidy(x) |>
        dplyr::filter(op %in%  c('=~', '~')) |>
        dplyr::select(term, estimate)
}) |> 
    list_rbind(names_to = 'mod')
    
# Note: The current GSCA_M may sometimes have trouble in computing D2 and U, especially when there's only one indicator for a common factor with a small loading


# IGSCA ------------------------------------------------------------------

igsca_pop <- list(
    uniC_uniF = 'xi1 <~ .4*x11

               xi2 =~ 1*x21

               xi2 ~ .5*xi1',
    uniC_triF = 'xi1 <~ .4*x11

               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23

               xi2 ~ .5*xi1',
    triC_uniF = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2 =~ 1*x21
               
               xi2 ~ .5*xi1',
    triC_triF = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23
               
               xi2 ~ 0.5*xi1',
    uniF_uniC = 'xi1 =~ 1*x11

               xi2 <~ .4*x21

               xi2 ~ .5*xi1',
    uniF_triC = 'xi1 =~ 1*x11

               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23

               xi2 ~ .5*xi1',
    triF_uniC = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2 =~ 1*x21
               
               xi2 ~ .5*xi1',
    triF_triC = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23
               
               xi2 ~ 0.5*xi1'
)


igsca_datapop <- lapply(igsca_pop, cSEM.DGP::generateData, .empirical = TRUE)

igsca_model_spec <- list(
    uniC_uniF = 'xi1 <~ x11

               xi2 =~ x21

               xi2 ~ xi1',
    uniC_triF = 'xi1 <~ x11

               xi2=~x21 + x22 + x23

               xi2 ~ xi1',
    triC_uniF = 'xi1<~x11 + x12 + x13
               
               xi2 =~ x21
               
               xi2 ~ xi1',
    triC_triF = 'xi1<~x11 + x12 + x13
               
               xi2=~x21 + x22 + x23
               
               xi2 ~ xi1',
    uniF_uniC = 'xi1 =~ x11

               xi2 <~ x21

               xi2 ~ xi1',
    uniF_triC = 'xi1 =~ x11

               xi2<~x21 + x22 + x23

               xi2 ~ xi1',
    triF_uniC = 'xi1=~x11 + x12 + x13
               
               xi2 <~ x21
               
               xi2 ~ xi1',
    triF_triC = 'xi1=~x11 + x12 + x13
               
               xi2<~x21 + x22 + x23
               
               xi2 ~ xi1'
)

igsca_mods <- mapply(
    cSEM::csem,
    .data = igsca_datapop,
    .model = igsca_model_spec,
    SIMPLIFY = FALSE,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.001,
    .iter_max = 1000
)

tidied_igsca_mods <- lapply(igsca_mods, function(x) {
    tidy(x) |>
        dplyr::filter(op %in%  c('=~', '~', '<~')) |>
        dplyr::select(term, estimate)
}) |> 
    list_rbind(names_to = 'mod')
    




# Old Test Code ----------------------------------------------------------
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



out_CCMP = csem(
    .data = datapop,
    .model = modelest,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP"
)
tidy(out_CCMP) |> View()

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
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP"
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




