# Postestimation: Assessing a model

## Introduction

As indicated by the name,
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md) is
used to assess a model estimated using the
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
function.

In **cSEM** model assessment is considered to be any task that in some
way or another seeks to assess the quality of the estimated model
*without conducting* *a statistical test* (tests are covered by the
`test_*` family of functions). Quality in this case is taken to be a
catch-all term for all common aspects of model assessment. This mainly
comprises fit indices, model selection criteria, reliability estimates,
common validity assessment criteria, effect sizes, and other related
quality measures/indices that do not rely on a formal test procedure.
Hereinafter, we will refer to a generic (fit) index, quality or
assessment measure as a **quality criterion**.

Currently the following quality criteria are implemented:

- Convergent and discriminant validity assessment:
  - The **average variance extracted** (AVE)
  - The **Fornell-Larcker** criterion
  - The **heterotrait-monotrait ratio of correlations** (HTMT)
- **Congeneric reliability** ($\rho_{C}$), also known as e.g.: composite
  reliability, construct reliability, (unidimensional) omega, Jöreskog’s
  $\rho$, $\rho_{A}$, or $\rho_{B}$.
- **Tau-equivalent reliability** ($\rho_{T}$), also known as e.g.:
  Cronbach alpha, alpha, $\alpha$, coefficient alpha, Guttman’s
  $\lambda_{3}$, KR-20.
- Distance measures
  - The **standardized root mean square residual** (SRMR)
  - The **geodesic distance** (DG)
  - The **squared Euclidian distance** (DL)
  - The **maximum-likelihood distance** (DML)
- Fit indices
  - The $\chi^{2}$-**statistic**
  - The $\chi^{2}/df$-**statistic**
  - The **comparative fit index** (CFI)
  - The **goodness-of-fit index** (GFI)
  - The **standardized root mean square residual** (SRMR)
  - The **root mean square error of approximation** (RMSEA)
  - The **normed fit index** (NFI)
  - The **non-normed fit index** (NNFI)
  - The **comparative fit index** (CFI)
  - The **incremental fit index** (IFI)
  - The **root mean square outer residual covariance**
    ($\text{RMS}_{\theta}$)
- The **Goodness-of-Fit** (GoF) proposed by Tenenhaus, Amanto, and Vinzi
  (2004).
- The **variance inflation factors** (VIF) for the structural equations
  as well as for Mode B regression equations (if
  `.approach_weights = "PLS-PM"`).
- The coefficient of determination and the adjusted coefficient of
  determination ($R^{2}$ and $R_{adj}^{2}$)
- A measure of the **effect size** (Cohen’s $f^{2}$).
- Direct, indirect and total effect assessment.
- Several model selection criteria as described in Sharma et al. (2019).

For implementation details see the [Methods & Formulae](#methods)
section.

## Syntax & Options

``` r
assess(
  .object              = NULL, 
  .only_common_factors = TRUE, 
  .quality_criterion   = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
                           "hqc", "mallows_cp", "ave",
                           "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "dg", "dl", "dml", "df",
                           "effects", "f2", "fl_criterion", "chi_square", "chi_square_df",
                           "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB"),
  ...
)
```

- `.object`:

  An object of class `cSEMResults` resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- `.quality_criterion`:

  A character string or a vector of character strings naming the quality
  criterion to compute. By default all quality criteria are computed
  (`"all"`). See
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  for a list of possible candidates.

- `.only_common_factors`:

  Logical. Should only concepts modeled as common factors be included
  when calculating one of the following quality criteria: AVE, the
  Fornell-Larcker criterion, HTMT, and all reliability estimates.
  Defaults to `TRUE`.

- `...`:

  Further arguments passed to functions called by
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).
  See
  [args_assess_dotdotdot](https://floschuberth.github.io/cSEM/reference/args_assess_dotdotdot.html)
  for a complete list of available arguments.

Like all postestimation functions
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
can be called on any object of class `cSEMResults`. The output is a
named list of the quality criteria given to `.quality_criterion`. By
default all possible quality criteria are calculated
(`.quality_criterion = "all"`).

### Details

In line with all of **cSEM**’s postestimation functions,
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md) is
a generic function with methods for objects of class
`cSEMResults_default`, `cSEMResults_multi`, `cSEMResults_2ndorder`. In
**cSEM** every `cSEMResults_*` object must also have class `cSEMResults`
for internal reasons. When using one of the major postestimation
functions, method dispatch is therefore technically done on one of the
`cSEMResults_*` class attributes, ignoring the `cSEMResults` class
attribute. As long as
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md) is
used directly, method dispatch is not of any practical concern to the
end-users.

#### Composite models vs. common factor models

Some assessment measures are inherently tied to the common factor model.
It is therefore unclear how to interpret their results in the context of
a composite model. Consequently, their computation is suppressed by
default for constructs modeled as composites. Currently, this applies to
the following quality criteria:

- AVE and validity assessment based thereon (i.e., the Fornell-Larcker
  criterion)
- HTMT and validity assessment based thereon
- All reliability measures

It is possible to force computation of all quality criteria for
constructs modeled as composites using `.only_common_factors = FALSE`,
however, we explicitly warn to interpret results, as they may not even
have a conceptual meaning.

All quality criteria assume that the estimated loadings, construct
correlations and path coefficients involved in the computation of a
specific quality measure are consistent estimates for their theoretical
population counterpart. If the user deliberately chooses an approach
that yields inconsistent estimates (by setting `.disattenuate = FALSE`
in [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
when the estimated model contains constructs modeled as common factors)
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
will still estimate all quantities, however, quantities such as the AVE
or the congeneric reliability $\rho_{C}$ inherit inconsistency.

## Methods & Formulae

This section provides technical details and relevant formulae. For the
relevant notation and terminology used in this section, see the
[Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
and the
[Termionology](https://floschuberth.github.io/cSEM/articles/Terminology.html)
help files.

### Average Variance Extracted (AVE)

#### Definition

The average variance extracted (AVE) was first proposed by Fornell and
Larcker (1981). Several definitions exist. For ease of comparison to
extant literature the most common definitions are given below:

- The AVE for a generic construct/latent variable $\eta$ is an estimate
  of how much of the variation of its indicators is due to the assumed
  latent variable. Consequently, the share of unexplained, i.e. error
  variation is 1 - AVE.
- The AVE for a generic construct/latent variable $\eta$ is the share of
  the total indicator variance (i.e., the sum of the indicator variances
  of all indicators connected to the construct), that is captured by the
  (indicator) true scores.
- The AVE for a generic construct/latent variable $\eta$ is the ratio of
  the sum of the (indicator) true score variances (explained variation)
  relative to the sum of the total indicator variances (total variation,
  i.e., the sum of the indicator variances of all indicators connected
  to the construct).
- Since for the regression of $x_{k}$ on $\eta_{k}$, the R squared
  ($R_{k}^{2})$ is equal to the share of variation of $x_{k}$ explained
  by $\eta_{k}$ relative to the total variation of $x_{k}$, the AVE for
  a generic construct/latent variable $\eta$ is equal to the average
  over all $R_{k}^{2}$.
- The AVE for a generic construct/latent variable $\eta$ is the sum of
  the squared correlation between indicator $x_{k}$ and the (indicator)
  true score $\eta_{k}$ relative to the sum of the indicator variances
  of all indicators connected to the construct in question.

It is important to stress that, although different in wording, all
definitions are synonymous!

The AVE is inherently tied to the common factor model. It is therefore
unclear how to interpret the AVE for constructs modeled as composites.
Consequently, the computation is suppressed by default for constructs
modeled as common factors. It is possible to force computation of the
AVE for constructs modeled as composites using
`.only_common_factors = FALSE`, however, we explicitly warn to interpret
results, as they may not even have a conceptual meaning.

#### Formulae

Using the results and notation derived and defined in the
[Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
help file, the AVE for a generic construct is:
$$AVE = \frac{\text{Sum indicator true score variances}}{\text{Sum indicator variances}} = \frac{\sum Var\left( \eta_{k} \right)}{\sum Var\left( x_{k} \right)} = \frac{\sum\lambda_{k}^{2}}{\sum\left( \lambda_{k}^{2} + Var\left( \varepsilon_{k} \right) \right)}$$
If $x_{k}$ is standardized (i.e., $Var\left( x_{k} \right) = 1$) the
denominator reduces to $K$ and the AVE for a generic construct is:
$$AVE = \frac{1}{K}\sum\lambda_{k}^{2} = \frac{1}{K}\sum\rho_{x_{k},\eta}^{2}$$
As an important consequence, the AVE is closely tied to the communality.
**Communality** ($COM_{k}$) is definied as the proportion of variation
in an indicator that is explained by its common factor. Empirically, it
is the square of the standardized loading of the $k$’th indicator
($\lambda_{k}^{2}$). Since indicators, scores/proxies and subsequently
loadings are always standardized in **cSEM**, the squared loading is
simply the squared correlation between the indicator and its related
construct/common factor. The AVE is also directly related to the
**indicator reliability**, defined as the squared correlation between an
indicator $k$ and its related proxy true score (see section
[Reliability](#reliability) below), which is again simply
$\lambda_{k}^{2}$. Therefore in **cSEM** we always have:

$$AVE = \frac{1}{K}\sum COM_{k} = \frac{1}{K}\sum\text{Indicator reliability}_{k} = \frac{1}{K}\sum\lambda_{k}^{2} = \frac{1}{K}\sum R_{k}^{2}$$

#### Implementation

The function is implemented as:
[`calculateAVE()`](https://floschuberth.github.io/cSEM/reference/calculateAVE.md).

#### See also

The AVE is the basis for the Fornell-Larcker criterion.

### Degrees of freedom

#### Definition

Degrees of freedom are calculated as the difference between the number
of non-redundant free elements of the empirical indicator correlation
matrix $\mathbf{S}$ and the model parameters.

Although, composite-based estimators retrieve parameters of the
postulated models by forming composites, which involves the estimation
of weights the computation of the degrees of freedom eventually depends
on the postulated model and the parameters implied by the model. Most
notably, a common factor model estimated by a composite-based approach
such as PLS has the same degrees of freedom compared to e.g., classical
maximum likelihood estimation of the same model.

#### Formulae

$$\begin{aligned}
\text{df} & {= {{\text{\# non-redundant off-diagonal elements of the empirical indicator correlation matrix}\mspace{6mu}}\mathbf{S}}} \\
 & {- \text{\# model parameters}}
\end{aligned}$$

If the model contains only linear terms the model parameters are:

- \# free correlations between exogenous constructs
- \# specified correlations between endogenous constructs
- \# structural parameters

In addition, for each construct $\eta_{j}$:

- \# of loadings if $\eta_{j}$ is modeled as a common factor
- \# of specified measurement error correlations between items of
  constructs $\eta_{j}$ if $\eta_{j}$ is modeled as a common factor
- \# of weights of $\eta_{j}$**minus 1** if $\eta_{j}$ is modeled as a
  composite. One weight per block is fixed and hence not counted as a
  model parameter since the variance of the composite is scaled to be
  unity.
- \# of non-redundant off-diagonal elements of $\mathbf{\Sigma}_{j}$ if
  $\eta_{j}$ is modeled as a composite.

If the model contains second-order terms the model parameters are
similar:

- \# free correlations between exogenous constructs
- \# specified correlations between endogenous constructs
- \# structural parameters. Note: relations between constructs
  measuring/forming the second-order construct are not path!

In addition, for each construct $\eta_{j}$ (including the second-order
constructs):

- \# of loadings if $\eta_{j}$ is modeled as a common factor
- \# of specified measurement error correlations between items of
  constructs $\eta_{j}$ if $\eta_{j}$ is modeled as a common factor
- \# of weights of $\eta_{j}$**minus 1** if $\eta_{j}$ is modeled as a
  composite. One weight per block is fixed and hence not counted as a
  model parameter since the variance of the composite is scaled to be
  unity.
- \# of non-redundant off-diagonal elements of $\mathbf{\Sigma}_{j}$ if
  $\eta_{j}$ is modeled as a composite.

##### Notes

1.  If all constructs are allowed to freely covary, i.e., there is no
    structural model and no structural parameters, all constructs are
    considered exogenous.
2.  If the structural model contains nonlinear terms (e.g.,
    $\eta_{1}^{2}$ or $\eta_{1}\eta_{2}$), degrees of freedom
    computation is currently unclear (at least to us). A warning is
    printed to inform the user that the calculation is may not be
    correct.

#### Implementation

The function is implemented as
[`calculateDf()`](https://floschuberth.github.io/cSEM/reference/calculateDf.md).

#### See also:

Degrees of freedom are required for several [fit
measures](#fit_indices).

### Fit Indices

#### Definition

Fit indices for confirmatory factor analysis (CFA) were first introduced
by Bentler and Bonett (1980). Since then a large number of indices has
been defined. Contrary to exact tests of model fit, the purpose of fit
indices is to measure the fit of a structural equation model on a
continuous scale. For normed fit indices this scale is between 0 and 1.
Fit indices can be divided into two classes:

- ‘badness of fit’ (resp. ‘lack of fit’) indices; a smaller value
  indicates a better fit.
- ‘goodness of fit’ indices; a higher value represents a better fit.

Several studies have analyzed the empirical and theoretical properties
of fit indices in the context of CFA where concepts are expressed by
latent variables. only little is known about the properties and the
performance of fit indices in composite models and for models estimated
using a composite-based approach. **cSEM** offers a number of fit
indices that are known from factor-based SEM. However, applied users
should be aware that only little is known about their applicability,
intuition, and interpretability in the context of models containing
constructs modeled as composites or for models estimated using a
composite-based approach.

Independent of the approach and model used, a particularly controversial
issue are cutoff values for fit indices (e.g., Marsh, Hau, and Wen
2004). In factor-based SEM cutoff values are rather popular. The basis
for these are numerous simulation studies, most notably Hu and Bentler
(1999). In contrast for composite models - for better or worse - no
cutoff values have been suggested.[¹](#fn1) Using
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md) to
calculate fit indices, the user should always keep in mind that the
value of a fit index is just *some* indication of good or bad fit. Other
aspects related to model fit must be considered as well. It is
unreasonable to make a binary decision about rejection or non-rejection
of a model by solely comparing the value of a fit index with a (more or
less) arbitrary cutoff value.

The definitions of fit indices calculated by
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
are given in the following:

- The $\chi^{2}$-**statistic** is the value of the fitting function
  times the sample size minus 1.
- The $\chi^{2}/df$-ratio is the $\chi^{2}$-statistic divided by its
  degrees of freedom.
- The **goodness-of-fit index** (GFI) measures the relative increase in
  fit of the specified model compared to no model at all.
- The **standardized root mean square residual** (SRMR) is the square
  root of the mean of squared residual correlations.
- The **root mean square error of approximation** (RMSEA) is the square
  root of the discrepancy due to approximation per degree of freedom.
- The **normed fit index** (NFI) measures the increase in fit when
  specifying the model under consideration relative to the fit of a
  certain baseline model called the “null model”.
- The **non-normed fit index** (NNFI) accounts for the degrees of
  freedom of the involved models. It is the ratio of the distance
  between the fit of the baseline model and the fit of the specified
  model (each per degree of freedom) and the distance between the fit of
  the baseline model and the expected fit of the specified model (each
  per degree of freedom).
- The **comparative fit index** (CFI) estimates the relative decrease in
  non-centrality when specifying the model under consideration instead
  of the baseline model.
- The **incremental fit index** (IFI) is the ratio of the distance
  between the fit of the baseline model and the fit of the specified
  model and the distance between the fit of the baseline model and the
  expected fit of the specified model. Its definition differs only
  marginally from the definition of the NNFI.
- The **root mean square outer residual covariance**
  ($\text{RMS}_{\theta}$) is defined as the square root of the mean
  squared covariances of the residuals of the outer model. The
  calculation of the indicator’s residual covariance matrix involves the
  calculation of the construct’s covariance matrix. See Lohmöller
  (1989).

It should be stressed again that (with the possible exception of the
$\text{RMS}_{\theta}$) none of the above mentioned fit indices were
originally designed for composite models. The indices RMSEA and CFI are
non-centrality based and require specific assumptions on model and data
typically made in CFA. The same applies for IFI and NNFI since their
calculation relies on the properties (primarily the expectation) of the
test statistic when data follows a normal distribution. In general,
those assumptions are not made in composite models and composite-based
estimators, respectively. For this reason, the intuition behind these
indices does not hold for composite-based SEM. Nevertheless, calculation
of these indices is also possible in this case. Whether the values of
these indices are still meaningful in a sense that they can be used for
assessment of model fit is an open question. Furthermore, values of fit
indices for composite-based estimators and factor-based estimators may
not be compared. Users should always keep this aspect and the general
limitations of fit indices in mind.

#### Formulae

The exact formulae of the fit indices as implemented in **cSEM** are
given in the following. The term
$F = F\left( \mathbf{S},\mathbf{\Sigma}\left( \widehat{\mathbf{θ}} \right) \right) = F\left( \mathbf{S},\widehat{\mathbf{\Sigma}} \right)$
stands for the value of the maximum likelihood fitting function
evaluated at $\mathbf{S}$ (the empirical covariance matrix of the
indicators) and $\widehat{\mathbf{\Sigma}}$ (the estimated model-implied
covariance matrix of the indicators). The value of the maximum
likelihood fitting function is computed by
[`calculateDML()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md).

##### The $\chi^{2}$-statistic

The $\chi^{2}$-**statistic** is defined as:
$$\chi^{2} = (N - 1) \cdot F$$ where $N$ is the sample size.

Main reference: K. G. Jöreskog (1969)

##### The $\chi^{2}/\text{df}$-ratio

The $\chi^{2}/\text{df}$-**statistic** is defined as:
$$\chi^{2} = (N - 1) \cdot F/\text{df}_{M}$$ where $N$ the sample size
and $\text{df}_{M}$ the degrees of freedom of the estimated model.

Main reference: K. G. Jöreskog (1969)

##### The goodness-of-fit index (GFI)

The GFI is generally defined in analogy to the coefficient of
determination ($R^{2}$) known from regression analysis as 1 minus the
share of the weighted unexplained variance (SSE; the difference between
$\mathbf{S}$ and $\widehat{\mathbf{\Sigma}}$) relative to the weighted
total variance (SST; the variance of $\mathbf{S}$):
$$GFI = 1 - \frac{\text{trace}\left\{ \left( \mathbf{W}^{- \frac{1}{2}}\lbrack\mathbf{S} - \widehat{\mathbf{\Sigma}}\rbrack\mathbf{W}^{- \frac{1}{2}} \right)^{2} \right\}}{\text{trace}\left\{ \left( \mathbf{W}^{- \frac{1}{2}}\mathbf{S}\mathbf{W}^{- \frac{1}{2}} \right)^{2} \right\}}$$
The matrix $\mathbf{W}$ is a weight matrix. Depending on the estimation
technique used to obtain $\widehat{\mathbf{θ}}$ different types of GFI
may be computed by choosing a particular weight.

1.  If $\mathbf{W} = \widehat{\mathbf{\Sigma}}$, the GFI is based on the
    SSE and the SST from a maximum likelihood estimation.
2.  If $\mathbf{W} = \widehat{\mathbf{S}}$, the GFI is based on SSE and
    the SST from a generalized least squares (GLS) estimation.
3.  If $\mathbf{W} = \widehat{\mathbf{I}}$, the GFI is based on SSE and
    the SST from a unweighted least squares (ULS) estimation.

Note that for any quadratic matrix , we have:
$\text{trace}\left( \mathbf{X}^{2} \right) = \sum_{i,j}x_{i}^{2}$.

Main references: Karl G. Jöreskog and Sörbom (1982), Mulaik et al.
(1989) and Tanaka and Huba (1985)

##### The standardized root mean square residual (SRMR)

The SRMR is defined as
$$\text{SRMR} = \sqrt{2\sum\limits_{j = 1}^{K}\sum\limits_{i = 1}^{j}\frac{\lbrack\left( s_{ij} - {\widehat{\sigma}}_{ij} \right)/\left( s_{ii}s_{jj} \right)^{1/2}\rbrack^{2}}{K(K + 1)}}$$
where $K$ stands for the number of indicators, $s_{ij}$ for the
empirical covariance between indicators $i$ and $j$, and
${\widehat{\sigma}}_{ij}$ for the estimated model-implied counterpart.
The SRMR describes with which distance the observed correlations are
reproduced on average by the model. Therefore, smaller values are
associated with a better fit. If data is standardized,
$s_{ii} = s_{jj} = 1$ holds, and the formula reduces to:
$$\text{SRMR} = \sqrt{2\sum\limits_{j = 1}^{K}\sum\limits_{i = 1}^{j}\frac{\left( s_{ij} - {\widehat{\sigma}}_{ij} \right)^{2}}{K(K + 1)}}$$

Main reference: Bentler (2006)

##### The root mean square error of approximation (RMSEA)

The RMSEA is defined as \$\$ \hat{\epsilon} =
\sqrt{\frac{\hat{F}\_0}{\text{df}\_{M}}} \quad \text{where} \quad
\hat{F}\_{0} = \max \Bigl( 0, F - \frac{\text{df}\_{M}}{N-1} \Bigr) \$\$
In this formula, $\text{df}_{M}$ stands for the degrees of freedom of
the specified model (see the [Degrees of Freedom](#df) section for
details on how the degrees of freedom are calculated). The term
${\widehat{F}}_{0}$ is an estimator for the discrepancy due to
approximation. Thus, the RMSEA measures the discrepancy due to
approximation per degree of freedom.

Main reference: Browne and Cudeck (1992)

##### The normed and non-normed fit index (NFI and NNFI)

The fit indices NFI and NNFI were among the first fit indices to be
introduced (Bentler and Bonett 1980). They are defined as:
$$\text{NFI} = \frac{F_{B} - F_{M}}{F_{B}}\quad\text{and}\quad\text{NNFI} = \frac{F_{B}/\text{df}_{B} - F_{M}/\text{df}_{M}}{F_{B}/\text{df}_{B} - 1/(N - 1)}$$
The term $F_{B}$ refers to the value of the fitting function in the null
model, $F_{M}$ to the value of the fitting function in the model under
consideration. Thus, the NFI measures the increase in fit relative to
the fit of the null model when specifying the model. The intuition of
NNFI is that (in factor-based methods) the expectation of
$F_{M}/\text{df}_{M}$ is equal to $1/N - 1$. This does not automatically
hold for composite-based estimators.

The NNFI measures the relative departure of the numerator’s term from
it’s expectation (in the denominator). That is why, the NNFI is not
normed and can take values larger than $1$.

Main reference: Bentler and Bonett (1980)

##### The comparative fit index (CFI)

The CFI is defined as:
$$\text{CFI} = 1 - \frac{\max\left( 0,(N - 1)F_{M} - \text{df}_{M} \right)}{\max\left( 0,(N - 1)F_{M} - \text{df}_{M},(N - 1)F_{B} - \text{df}_{B} \right)}$$
Like the RMSEA, the CFI is a non-centrality based index. It measures the
increase in fit (that is to say the reduction in non-centrality) when
specifying the model under consideration relative to the fit of the null
model. The CFI is a normed index with a value of $1$ indicating the best
fit. Since it makes use of the assumptions in factor-based methods, its
intuition does not apply to composite-based estimators.

Main reference: Bentler (1990).

##### The incremental fit index (IFI)

The IFI is defined as:
$$\text{IFI} = \frac{F_{B} - F_{M}}{F_{B} - df_{M}/(N - 1)}$$ The
rationale underlying the IFI is that the term $F_{B} - F_{M}$ (in the
numerator) is compared with its expectation
$F_{B} - \text{df}_{M}/(N - 1)$ (in the denominator).

Main reference: Bollen (1989)

##### The root mean square outer residual covariance

#### Implementation

The functions are implemented as:
[`calculateChiSquare()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateChiSquareDf()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateCFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateNFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateNNFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateIFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateGFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateRMSEA()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateRMSTheta()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateSRMR()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md).

#### See also

Several fit indices require a fitting function, i.e., a distance measure
like the geodesic distance, the squared Euclidean distance or the
maximum likelihood distance. These are implemented as:
[`calculateDG()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
[`calculateDL()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
and
[`calculateDML()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md).

### Reliability

#### Definition

Reliability is the **consistency of measurement**, i.e., the degree to
which a hypothetical repetition of the same measure would yield the same
results. As such, reliability is the closeness of a measure to an error
free measure. It is not to be confused with validity as a perfectly
reliable measure may be invalid.

Practically, reliability must be empirically assessed based on a
theoretical framework. The dominant theoretical framework against which
to compare empirical reliability results to is the well-known [true
score](https://floschuberth.github.io/cSEM/articles/Terminology.html)
framework which provides the foundation for the measurement model
described in the
[Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
help file. Based on the true score framework and using the terminology
and notation of the
[Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
and
[Termniology](https://floschuberth.github.io/cSEM/articles/Terminology.html)
help files, reliability of a generic measurement is defined as:

1.  The amount of proxy true score variance,
    $Var\left( \bar{\eta} \right)$, relative to the the proxy or test
    score variance, $Var\left( \widehat{\eta} \right)$.
2.  This is identical to the squared correlation between the common
    factor and its proxy/composite or test score:
    $\rho_{\eta,\widehat{\eta}}^{2} = Cor\left( \eta,\widehat{\eta} \right)^{2}$.

This “kind” of reliability is commonly referred to as **internal
consistency reliability**.

Based on the true score theory three major types of measurement models
are distinguished. Each type implies different assumptions which give
rise to the formulae written below. The well-established names for the
different types of measurement model provide natural naming candidates
for their corresponding (internal consistency) reliabilities measure:

1.  **Parallel** – Assumption:
    $\left. \eta_{kj} = \eta_{j}\rightarrow\lambda_{kj} = \lambda_{j} \right.$
    and
    $Var\left( \varepsilon_{kj} \right) = Var\left( \varepsilon_{j} \right)$.
2.  **Tau-equivalent** – Assumption:
    $\left. \eta_{kj} = \eta_{j}\rightarrow\lambda_{kj} = \lambda_{j} \right.$
    and
    $Var\left( \varepsilon_{kj} \right) \neq Var\left( \varepsilon_{lj} \right)$.
3.  **Congeneric** – Assumption: $\eta_{kj} = \lambda_{kj}\eta_{j}$ and
    $Var\left( \varepsilon_{kj} \right) \neq Var\left( \varepsilon_{lj} \right)$.

In principal the test score $\widehat{\eta}$ is a weighted linear
combinations of the indicators, i.e., a proxy or stand-in for the true
score/common factor. Historically, however, the test score is generally
assumed to be a simple sum score, i.e., a weighted sum of indicators
with all weights assumed to be equal to one. Hence, well-known
reliability measures such as Jöreskog’s $\rho$ or Cronbach’s $\alpha$
are defined with respect to a test score that indeed represents a simple
sum score. Yet, all reliability measures originally developed assuming a
sum score may equally well be computed with respect to a composite,
i.e., a weighted score with weights not necessarily equal to one.

Apart form the distinction between congeneric (i.e., Jöreskog’s $\rho$)
and tau-equivalent reliability (i.e., Cronbach’s $\alpha$) we therefore
distinguish between reliability estimates based on a test score
(composite) that uses the weights of the weight approach used to obtain
`.object` and a test score (proxy) based on unit weights. The former is
indicated by adding “**weighted**” to the original name.

#### Formulae

The most general formula for reliability is the **(weighted) congeneric
reliability**:

$$\rho_{C;\text{weighted}} = \frac{Var\left( \bar{\eta} \right)}{Var\left( {\widehat{\eta}}_{k} \right)} = \frac{(\mathbf{w}\prime{\mathbf{λ}})^{2}}{\mathbf{w}\prime\mathbf{\Sigma}\mathbf{w}}$$
Assuming $\mathbf{w} = {\mathbf{ι}}$, i.e., unit weights, the
“classical” formula for congeneric reliability (i.e., Jöreskog’s
$\rho$), follows:
$$\rho_{C} = \frac{Var\left( \bar{\eta} \right)}{Var\left( {\widehat{\eta}}_{k} \right)} = \frac{\left( \sum\lambda_{k} \right)^{2}}{\left( \sum\lambda_{k} \right)^{2} + Var\left( \bar{\varepsilon} \right)}$$
Using the assumptions imposed by the tau-equivalent measurement model we
obtain the **(weighted) tau-equivalent reliability, i.e., (weighted)
Cronbach’s alpha)**:

$$\rho_{T;\text{weighted}} = \frac{\lambda^{2}\left( \sum w_{k} \right)^{2}}{\lambda^{2}\left( \sum w_{k} \right)^{2} + \sum w_{k}^{2}Var\left( \varepsilon_{k} \right)} = \frac{{\bar{\sigma}}_{x}\left( \sum w_{k} \right)^{2}}{{\bar{\sigma}}_{x}\left\lbrack \left( \sum w_{k} \right)^{2} - \sum w_{k}^{2} \right\rbrack + \sum w_{k}^{2}Var\left( x_{k} \right)}$$
where we used the fact that if $\lambda_{k} = \lambda$
(tau-equivalence), $\lambda^{2}$ equals the average covariance between
indicators:
$${\bar{\sigma}}_{x} = \frac{1}{K(K - 1)}\sum\limits_{k = 1}^{K}\sum\limits_{l = 1}^{K}\sigma_{kl}$$
Again, assuming $w_{k} = 1$, i.e., unit weights, the “classical” formula
for tau-equivalent reliability (Cronbach’s $\alpha$) follows:
$$\rho_{T} = \frac{\lambda^{2}K^{2}}{\lambda^{2}K^{2} + \sum Var\left( {\bar{\varepsilon}}_{k} \right)} = \frac{{\bar{\sigma}}_{x}K^{2}}{{\bar{\sigma}}_{x}\left\lbrack K^{2} - K \right\rbrack + KVar\left( x_{k} \right)}$$
Using the assumptions imposed by the parallel measurement model we
obtain the **parallel reliability**:

$$\rho_{P} = \frac{\lambda^{2}\left( \sum w_{k} \right)^{2}}{\lambda^{2}\left( \sum w_{k} \right)^{2} + Var(\varepsilon)\sum w_{k}^{2}} = \frac{{\bar{\sigma}}_{x}\left( \sum w_{k} \right)^{2}}{{\bar{\sigma}}_{x}\left\lbrack \left( \sum w_{k} \right)^{2} - \sum w_{k}^{2} \right\rbrack + Var(x)\sum w_{k}^{2}}$$

In **cSEM** indicators are always standardized and weights are chosen
such that $Var\left( {\widehat{\eta}}_{k} \right) = 1$. This is done by
scaling the weight vector $\mathbf{w}$ by
$(\mathbf{w}\prime\mathbf{\Sigma}\mathbf{w})^{- \frac{1}{2}}$. This
simplifies the formulae: $$\begin{aligned}
\rho_{C;\text{weighted}} & {= \left( \sum w_{k}\lambda_{k} \right)^{2} = (\mathbf{w}\prime{\mathbf{λ}})^{2}} \\
{\rho_{T;\text{weighted}} = \rho_{P;\text{weighted}}} & {= {\bar{\rho}}_{x}\left( \sum w_{k} \right)^{2}} \\
 & 
\end{aligned}$$ where ${\bar{\rho}}_{x} = {\bar{\sigma}}_{x}$ is the
average correlation between indicators. Consequently, parallel and
tau-equivalent reliability are always identical in **cSEM**.

So far formulae have been motivated theoretically. Since
$\mathbf{\Sigma}$ is unknown it can be replaced by $\mathbf{S}$ (the
empirical indicator correlation matrix) or $\widehat{\mathbf{\Sigma}}$
(the model-implied indicator correlation matrix), however, $\mathbf{S}$
and $\widehat{\mathbf{\Sigma}}$ are generally not equal. The practical
implication is that if $\rho_{C}$ is computed as
$(\mathbf{w}\prime{\mathbf{λ}})^{2}$ using unit weights the weights can
in fact be scaled by both
$(\mathbf{w}\prime\mathbf{S}\mathbf{w})^{- \frac{1}{2}}$ or
$\left( \mathbf{w}\prime\widehat{\mathbf{\Sigma}}\mathbf{w} \right)^{- \frac{1}{2}}$!
Similarly, $\rho_{C;\text{weighted}}$ can be computed using weights
scaled using either $\mathbf{S}$ or $\widehat{\mathbf{\Sigma}}$.
Consequently there are in fact four types of congeneric reliability
depending the type of weight and the type of scaling for the weights.
Hence, the calculation is of “the” congeneric reliability is always:
$$(\mathbf{w}\prime{\mathbf{λ}})^{2}$$ where $\mathbf{w}$ can be:

1.  a vector of unit weights scaled by
    $\left( \mathbf{w}\prime\widehat{\mathbf{\Sigma}}\mathbf{w} \right)^{- \frac{1}{2}}$.
    This is typically what people refer to as *the* congeneric
    reliability (Jöreskog’s $\rho$). We label this type of reliability
    estimate $\rho_{C}$.
2.  a vector of unit weights scaled by
    $(\mathbf{w}\prime\mathbf{S}\mathbf{w})^{- \frac{1}{2}}$. This has
    no known name. Its usefulness is an open question. We label this
    type of reliability estimate $\rho_{C;mm}$.
3.  a vector of weights obtained using a composite-based estimator
    (e.g. PLS-PM) scaled by
    $(\mathbf{w}\prime\mathbf{S}\mathbf{w})^{- \frac{1}{2}}$. This is
    Dijkstra Henseler’s $\rho_{A}$. We label this type of reliability
    estimate $\rho_{C;\text{weighted}}$.
4.  a vector of weights obtained using a composite-based estimator
    (e.g. PLS-PM) scaled by
    $\left( \mathbf{w}\prime\widehat{\mathbf{\Sigma}}\mathbf{w} \right)^{- \frac{1}{2}}$.
    This has no known name. Its usefulness is an open question. We label
    this type of reliability estimate $\rho_{C;\text{weighted};mm}$

##### A note on the terminology

A vast bulk of literature dating back to seminal work by Spearman (e.g.,
Spearman (1904)) has been written on the subject of reliability.
Inevitably, definitions, formulae, notation and terminology conventions
are unsystematic and confusing. This is particularly true for newcomers
to structural equation modeling or applied users whose primary concern
is to apply the appropriate method to the appropriate case without
poring over books and research papers to understand each intricate
detail.

In **cSEM** we seek to make working with reliabilities as consistent as
possible by relying on a paper by Cho (2016) who proposed uniform
formula-generating methods and a systematic naming conventions for all
common reliability measures. Naturally, some of the conventional
terminology is deeply entrenched within the nomenclatura of a particular
filed (e.g., coefficient alpha alias Cronbach’s alpha in pychometrics)
such that a new, albeit consistent, naming scheme seems superfluous at
best. However, we belief the merit of a “standardized” naming pattern
will eventually be helpful to all users as it helps clarify potential
misconceptions thus preventing potential misuse, such as the (ab)use of
Cronbach alpha as a reliability measure for congeneric measurement
models.

Apart from these considerations, this package takes a pragmatic stance
in a sense that we use consistent naming because it naturally provides a
consistent naming scheme for the functions and the systematic formula
generating methods because they make code maintenance easier.
Eventually, what matters is the formula and more so its correct
application. To facilitate the translation between different naming
systems and conventions we provide a “translation table” below:

|          Systematic names           |        Mathematical        |                                                       Synonymous terms                                                        |
|:-----------------------------------:|:--------------------------:|:-----------------------------------------------------------------------------------------------------------------------------:|
|        Parallel reliability         |         $\rho_{P}$         |                  Spearman-Brown formula, Spearman-Brown prophecy, Standardized alpha, Split-half reliability                  |
|     Tau-equivalent reliability      |         $\rho_{T}$         |                         Cronbach’s alpha, $\alpha$, Coefficient alpha Guttman’s $\lambda_{3}$, KR-20                          |
| Tau-equivalent reliability weighted | $\rho_{T;\text{weighted}}$ |                                                               –                                                               |
|       Congeneric reliability        |         $\rho_{C}$         | Composite reliability, Jöreskog’s $\rho$, Construct reliability, $\omega$, reliability coefficient, Dillon-Goldstein’s $\rho$ |
|   Congeneric reliability weighted   | $\rho_{C;\text{weighted}}$ |                                                Dijkstra-Henseler’s $\rho_{A}$                                                 |

Systematic names and common synonymous names for the reliability
estimates found in the literature

##### Closed-form confidence interval

Trinchera, Marie, and Marcoulides (2018) proposed a closed-form
confidence interval (CI) for the tau-equivalent reliability (Cronbach’s
alpha). To compute the CI, set `.closed_form_ci = TRUE` when calling
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md) or
invoke `calculateRhoT(..., .closed_form_ci = TRUE)` directly. The level
of the CI can be changed by supplying a single value or a vector of
values to `.alpha`.

#### Implementation

The functions are implemented as
[`calculateRhoC()`](https://floschuberth.github.io/cSEM/reference/reliability.md)
and
[`calculateRhoT()`](https://floschuberth.github.io/cSEM/reference/reliability.md).

### The Goodness of Fit (GoF)

#### Definition

Calculate the Goodness of Fit (GoF) proposed by Tenenhaus, Amanto, and
Vinzi (2004). Note that, contrary to what the name suggests, the GoF is
**not** a measure of (overall) model fit in a $\chi^{2}$-fit test sense.
See e.g. Henseler and Sarstedt (2012) for a discussion.

#### Formulae

The GoF is defined as:

$$\text{GoF} = \sqrt{\varnothing\text{COM}_{k} \times \varnothing R_{structural}^{2}} = \sqrt{\frac{1}{k}\sum\limits_{k = 1}^{K}\lambda_{k}^{2} + \frac{1}{M}\sum\limits_{m = 1}^{M}R_{m;structural}^{2}}$$
where $COM_{k}$ is the communality of indicator $k$, i.e. the variance
in the indicator that is explained by its connected latent variable and
$R_{m;structural}^{2}$ the R squared of the $m$’th equation of the
structural model.

#### Implementation

The function is implemented as:
[`calculateGoF()`](https://floschuberth.github.io/cSEM/reference/calculateGoF.md).

### The Heterotrait-Monotrait-Ratio of Correlations (HTMT)

#### Definition

The heterotrait-monotrait ratio of correlations (HTMT) was first
proposed by  
Henseler, Ringle, and Sarstedt (2015) to assess convergent and
discriminant validity.

#### Formulae

See: Henseler, Ringle, and Sarstedt (2015) on page 121 (equation (6))

#### Implementation

The function is implemented as:
[`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md).

## Literature

Bentler, Peter M. 1990. “Comparative Fit Indexes in Structural Models.”
*Psychological Bulletin* 107 (2): 238–46.

———. 2006. *EQS 6 Structural Equations Program Manual* (version 6).
Encino, CA: Multivariate Software, Inc.

Bentler, Peter M., and Douglas G. Bonett. 1980. “Significance Tests and
Goodness of Fit in the Analysis of Covariance Structures.”
*Psychological Bulletin* 88 (3): 588–606.

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley-Interscience.

Browne, Michael W., and Robert Cudeck. 1992. “Alternative Ways of
Assessing Model Fit.” *Sociological Methods & Research* 21 (2): 230–58.

Cho, Eunseong. 2016. “Making Reliability Reliable.” *Organizational
Research Methods* 19 (4): 651–82.
<https://doi.org/10.1177/1094428116656239>.

Fornell, C., and D. F. Larcker. 1981. “Evaluating Structural Equation
Models with Unobservable Variables and Measurement Error.” *Journal of
Marketing Research* XVIII: 39–50.

Henseler, Jörg, Christian M. Ringle, and Marko Sarstedt. 2015. “A New
Criterion for Assessing Discriminant Validity in Variance-Based
Structural Equation Modeling.” *Journal of the Academy of Marketing
Science* 43 (1): 115–35. <https://doi.org/10.1007/s11747-014-0403-8>.

Henseler, Jörg, and Marko Sarstedt. 2012. “Goodness-of-Fit Indices for
Partial Least Squares Path Modeling.” *Computational Statistics* 28 (2):
565–80. <https://doi.org/10.1007/s00180-012-0317-1>.

Hu, Li-tze, and Peter M. Bentler. 1999. “Cutoff Criteria for Fit Indexes
in Covariance Structure Analysis: Conventional Criteria Versus New
Alternatives.” *Structural Equation Modeling* 6 (1): 1–55.

Jöreskog, K. G. 1969. “A General Approach to Confirmatory Maximum
Likelihood Factor Analysis.” *Psychometrika* 34 (2): 183–202.
<https://doi.org/10.1007/bf02289343>.

Jöreskog, Karl G., and Dag Sörbom. 1982. “Recent Developments in
Structural Equation Modeling.” *Journal of Marketing Research* 19 (4):
404–16.

Lohmöller, Jan-Bernd. 1989. *Latent Variable Path Modeling with Partial
Least Squares*. Physica, Heidelberg.

Marsh, Herbert W., Kit-Tai Hau, and Zhonglin Wen. 2004. “In Search of
Golden Rules: Comment on Hypothesis-Testing Approaches to Setting Cutoff
Values for Fit Indexes and Dangers in Overgeneralizing Hu and Bentler’s
(1999) Findings.” *Structural Equation Modeling: A Multidisciplinary
Journal* 11 (3): 320–41. <https://doi.org/10.1207/s15328007sem1103_2>.

Mulaik, Stanley A., Larry R. James, Judith Van Alstine, Nathan Bennett,
Sherri Lind, and C. Dean Stilwell. 1989. “Evaluation of Goodness-of-Fit
Indices for Structural Equation Models.” *Psychological Bulletin* 105
(3): 430–45. <https://doi.org/10.1037/0033-2909.105.3.430>.

Sharma, Pratyush, Marko Sarstedt, Galit Shmueli, Kevin H. Kim, and Kai
O. Thiele. 2019. “PLS-Based Model Selection: The Role of Alternative
Explanations in Information Systems Research.” *Journal of the
Association for Information Systems* 20 (4).

Tanaka, J. S., and G. J. Huba. 1985. “A Fit Index for Covariance
Structure Models Under Arbitrary GLS Estimation.” *British Journal of
Mathematical and Statistical Psychology* 38 (2): 197–201.
<https://doi.org/10.1111/j.2044-8317.1985.tb00834.x>.

Tenenhaus, Michel, Silvano Amanto, and Vincenzo Esposito Vinzi. 2004. “A
Global Goodness-of-Fit Index for PLS Structural Equation Modelling.” In
*Proceedings of the XLII SIS Scientific Meeting*, 739–42.

Trinchera, Laura, Nicolas Marie, and George A. Marcoulides. 2018. “A
Distribution Free Interval Estimate for Coefficient Alpha.” *Structural
Equation Modeling: A Multidisciplinary Journal* 25 (6): 876–87.
<https://doi.org/10.1080/10705511.2018.1431544>.

------------------------------------------------------------------------

1.  There are some cutoffs such as e.g., the SRMR should be less than
    0.08 or 0.1, however, these values are essentially arbitrary as they
    have never been formally motivated. Reference is usually done to Hu
    and Bentler (1999) which based the cut-off on a simulation using
    factor-based SEM.
