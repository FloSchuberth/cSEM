# Introduction to cSEM

## Preface

Structural equation modeling (SEM) has been used and developed for
decades across a variety of research fields, including psychology,
sociology, and business research.

As an almost inevitable consequence, a different terminology system and,
to some extent, mathematical notation has evolved within each field over
the years. This “terminological mess” is one of the major obstacles in
interdisciplinary research and (scientific) debate, as it hinders a
broader understanding of methodological issues and, even worse, promotes
systematic misuse (e.g., the use of Cronbach’s alpha as an estimator for
congeneric reliability). This is especially true for users who are new
to SEM, or for practitioners who are overwhelmed by the terminology of
their own field, only to find that the term they thought to have finally
understood is defined differently in another field, adding to the
confusion. A prime example is the term “formative” (measurement) which
has been used to describe both a causal-formative and [a composite
model](#cmodel), see e.g., Henseler (2017) for a clarification.

Ultimately, this is a matter of (mis)communication, which we believe can
only be satisfactorily resolved by providing a clear, unambiguous
definition for each term and symbol used in the package. We emphasize
that we are not trying to impose “our” conventions, nor are we claiming
they are the “correct” conventions, but merely trying to make
communication between us (the authors of the package) and you (users of
the package) as unambiguous and error-free as possible.

Therefore, we provide a
[Terminology](https://floschuberth.github.io/cSEM/articles/Terminology.html)
file and a
[Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
file, which contain key terms and mathematical notation/symbols that we
consider important, along with our definition. Users are encouraged to
read these files carefully to avoid potential misunderstandings.

1.  The
    [Terminology](https://floschuberth.github.io/cSEM/articles/Terminology.html)
    file contains any term that we feel should be defined and explained
    to ensure that package users understand the supplementary help files
    and vignettes in the way the package authors intended.
2.  The
    [Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
    file contains some fundamental mathematical notation/symbols used in
    the package documentation, along with a definition. With a few
    exceptions, we mostly follow standard notation as laid out in e.g.,
    Bollen (1989).

The package has been designed according to a set of principles and
terminology that are in part different from those of other commonly used
open source or commercial software packages with similar content (e.g.,
SmartPLS). Together with the
[Terminology](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and
[Notation](https://floschuberth.github.io/cSEM/articles/Notation.html)
files this introduction explains these principles. Finally, to use
**cSEM** effectively, it is helpful to understand its design. Therefore,
the architecture and design of the package and what we call the “cSEM
workflow” will be discussed.

## Composite-based structural equation modeling

### What is structural equation modeling (SEM)

[Structural equation modeling
(SEM)](https://floschuberth.github.io/cSEM/articles/Terminology.html) is
about analyzing, i.e., modeling, estimating, assessing, and testing, the
(causal) relationships between
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html) -
an entity defined by a conceptual definition - with other
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and/or observable quantities generally referred to as [indicators,
manifest variables or
items](https://floschuberth.github.io/cSEM/articles/Terminology.html).
Broadly speaking, two modeling approaches for the
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and their relationship exist. We refer to the first as [the latent
variable or common factor model](#lvmodel) and to the second as [the
composite model](#cmodel). Each approach entails a set of methods, test,
and evaluation criteria as well as a specific terminology that may or
may not be adequate within the realm of the other approach.

#### The classical latent variable or common factor model

Assuming a researcher identifies
$J$[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and $K$ indicators, the fundamental feature of the latent variable model
is the assumption of the existence of a set of $J$[latent variables
(common
factors)](https://floschuberth.github.io/cSEM/articles/Terminology.html)
that each serve as a representation of one of the
$J$[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
to be studied in a sense that each [latent
variable](https://floschuberth.github.io/cSEM/articles/Terminology.html)
is causally responsible for the manifestations of a set of
$K_{j}$[indicators](https://floschuberth.github.io/cSEM/articles/Terminology.html)
which are supposed to measure the
[concept](https://floschuberth.github.io/cSEM/articles/Terminology.html)
in question. The entirety of these measurement relations is captured by
the **measurement model** which relates
[indicators](https://floschuberth.github.io/cSEM/articles/Terminology.html)
to [latent
variables](https://floschuberth.github.io/cSEM/articles/Terminology.html)
according to the researchers theory of how observables are related to
the
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
in question. The entirety of the relationships between
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
(i.e., its representation in a statistical model, the
[construct](https://floschuberth.github.io/cSEM/articles/Terminology.html))
is captured by the **structural model** whose parameters are usually at
the center of the researchers interest. Caution is warranted though as
the [common
factor](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and its respective
[concept](https://floschuberth.github.io/cSEM/articles/Terminology.html)
are not the same thing. Within the “classical” [covariance-based or
factor-based](https://floschuberth.github.io/cSEM/articles/Terminology.html)
literature
[concept](https://floschuberth.github.io/cSEM/articles/Terminology.html),
[construct](https://floschuberth.github.io/cSEM/articles/Terminology.html),
[latent
variable](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and its representation the [common
factor](https://floschuberth.github.io/cSEM/articles/Terminology.html)
have often been used interchangeably (Rigdon 2012, 2016; Rigdon, Becker,
and Sarstedt 2019). This will not be the case in **cSEM** and readers
are explicitly made aware of the fact that
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
are the abstract entity which *may* be modeled by a [common
factor](https://floschuberth.github.io/cSEM/articles/Terminology.html),
however, no assertion as to the correctness of this approach in terms of
“closeness of the [common
factor](https://floschuberth.github.io/cSEM/articles/Terminology.html)
and its related
[concept](https://floschuberth.github.io/cSEM/articles/Terminology.html)”
are made.

Parameters in latent variable models are usually retrieved by maximum
likelihood (ML). The basic idea of ML is to find parameters such that
the difference between the model-implied and the empirical indicator
covariance matrix is minimized. Such estimation methods are therefore
often referred to as [covariance-based
methods](https://floschuberth.github.io/cSEM/articles/Terminology.html).

#### The composite model

The second approach is known as the composite model. As opposed to the
latent variable or common factor model, composites do not presuppose the
existence of a [latent
variable](https://floschuberth.github.io/cSEM/articles/Terminology.html).
Hence, designed entities (artifacts) such as the “OECD Better Life
Index” that arguably have no latent counterpart may be adequately
described by a
[composite](https://floschuberth.github.io/cSEM/articles/Terminology.html),
i.e., the linear combination of observables defining the composite.
Composites may also be formed to represent [latent variables/common
factors](https://floschuberth.github.io/cSEM/articles/Terminology.html)
(or more precisely
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
modeled as [common
factors](https://floschuberth.github.io/cSEM/articles/Terminology.html))
in which case the
[composite](https://floschuberth.github.io/cSEM/articles/Terminology.html)
serves as a [proxy or
stand-in](https://floschuberth.github.io/cSEM/articles/Terminology.html)
for the [latent
variable](https://floschuberth.github.io/cSEM/articles/Terminology.html).
However, in **cSEM**, the term “composite ***model***” is only used to
refer to a model in a the former sense, i.e., a model in which the
[composite](https://floschuberth.github.io/cSEM/articles/Terminology.html)
is a direct representation of
[concept/construct](https://floschuberth.github.io/cSEM/articles/Terminology.html)!

Parameters in composite models are retrieved by a composite-based
approach such as partial least squares path modeling (PLS-PM),
generalized structured component analysis (GSCA) or dimension reduction
techniques such as principal component analysis (PCA). The basic idea of
any composite-based approach is to build scores/composites for each
concept and subsequently retrieve structural model parameters by a
series of (linear) regressions. Such estimation methods are therefore
often referred to as [variance-based
methods](https://floschuberth.github.io/cSEM/articles/Terminology.html)
as regression maximizes the explained variance of the dependent
variable.

### What is composite-based SEM?

Composite-based SEM is the entirety of methods, approaches, procedures,
algorithms that in some way or another involve linear compounds
([composites/proxies/scores](https://floschuberth.github.io/cSEM/articles/Terminology.html)),
i.e., linear combinations of observables when retrieving (estimating)
quantities of interest such as the coefficients of the structural model.
It is crucial to clearly distinguish between the **composite model** and
**composite-based SEM**. They are not the same. While the former is
“only” ***a*** *statistical model* relating
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
to observables, the latter simply states that
[composites](https://floschuberth.github.io/cSEM/articles/Terminology.html) -
linear compounds, i.e., weighted linear combinations of observables -
are used to retrieve quantities of interest! Hence, composite-based SEM
as a way of obtaining/estimating parameters of interest may thus be used
for [the latent variable or common factor model](#lvmodel) as well as
[the composite model](#cmodel). However, interpretation of the parameter
estimates is fundamentally different since the underlying models differ!

As sketched above, common factor and composite models fundamentally
differ in how the relation between observables and
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)
is modeled. Naturally, results and, most notably, their
(correct/meaningful) interpretation critically hinge on what type of
model the user specifies. Across the package we therefore strictly
distinguish between

- [Concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)/[Constructs](https://floschuberth.github.io/cSEM/articles/Terminology.html)
  modeled as common factors (or alternatively latent variables)
- [Concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)/[Constructs](https://floschuberth.github.io/cSEM/articles/Terminology.html)
  modeled as composites

These phrases will repeatedly appear in help files and other
complementary files. It is therefore crucial to remember what they
supposed to convey.

## Using cSEM

The idea of cSEM is twofold:

1.  Provide a unified framework to all common composite-based SEM
    approaches, including typical postestimation procedures. In
    principal, it should be similar to what
    [lavaan](https://lavaan.ugent.be/) does for
    covariance-bases/factor-based SEM.
2.  Make the user experience as hassle-free as possible.

The first point is an always ongoing task since approaches are
constantly evolving with new developments appearing at a pace that we,
the package authors, will not be able to keep up with. The second point,
however, was particularly important to us as we have been frustrated
ourselves by how technical, unfriendly packages in R can be. Hence, from
the very start we envisioned a workflow that essentially only comprises
three steps:

1.  **Get the essential**: no estimator or approach works without data
    and a description of what parameters are to be estimated and how
    data is related to these parameters, i.e. a model. Hence, we always
    need a data set and model. Since, model specification in [lavaan
    model syntax](https://lavaan.ugent.be/tutorial/syntax1.html) is
    probably unbeatable in its ease and well known to R users that have
    an interest in SEM, to us, [lavaan model
    syntax](https://lavaan.ugent.be/tutorial/syntax1.html) is the
    obvious tool for users to specify their model. Experience tells,
    that for R beginners the biggest obstacle has been to get the data
    into R. However, largely thanks to [RStudio](https://posit.co/),
    data import and data transformation are nowadays relatively easy to
    handle. See the [Preparing the data](#preparedata) and the
    [Specifying a model](#specifyingamodel) sections below.

2.  **Estimate**: no matter the model and type of data, estimation is
    always done using one central function with the data as its first
    and the model as its second argument:

    ``` r
    csem(.data = my_data, .model = my_model)
    ```

    Naturally, the
    [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
    function has a number of additional arguments to fine-tune the
    estimation, however, since
    [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
    automatically recognizes, for instance, whether a concept was
    modeled as a common factor or composite and automatically applies
    appropriate correction for attenuation, default arguments are often
    sufficient. See the [Estimate using csem()](#estimate) section
    below.

3.  **Postestimate**: Inspired by the [grammar of data
    manipulation](https://dplyr.tidyverse.org/) underlying the
    [dplyr](https://dplyr.tidyverse.org/) package, cSEM provides 5
    postestimation verbs that concisely cover all common postestimation
    tasks as well as 4 additional test commands and 2 general do
    commands:

    1.  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
    2.  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
    3.  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
    4.  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
    5.  [`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md)
    6.  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
    7.  [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md)
    8.  [`testHausman()`](https://floschuberth.github.io/cSEM/reference/testHausman.md)
    9.  [`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md)
    10. [`doIPMA()`](https://floschuberth.github.io/cSEM/reference/doIPMA.md)
    11. `doNonlinearRedundancyAnalysis()`
    12. [`doRedundancyAnalysis()`](https://floschuberth.github.io/cSEM/reference/doRedundancyAnalysis.md)

    All verbs accept the result of a call to
    [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) as
    input which makes working with these function extremely simple. You
    only need to remember the word, not any specific syntax or
    arguments. Of course, all functions have a number of additional
    arguments to fine-tune the postestimation. See the [Apply
    postestimation functions](#postestimate) sections below. For details
    on the arguments consult the individual help files.

The price we pay for an increase in flexibility is primarily a, mostly
minor, loss in computational speed, in particular, when intense
resampling is involved (i.e., 5000 bootstrap run for a complex model
with, say, 1000 observations). Users looking for the most efficient
implementation of common resampling routines may find faster
implementations. That said, we believe, the time saved when using a
standardized estimate-postestimate workflow, no matter the model or data
used, well outweighs the potential loss in computational efficiency.

The following sections describe the workflow in more detail.

### The cSEM-Workflow

As described in the previous section, working with **cSEM** consists of
3-4 steps:

1.  Prepare/load the data to analyze
2.  Specify the model to estimate
3.  Estimate using the
    [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
    function
4.  Apply postestimation functions to the result of 3.

#### Prepare the data

Technically, preparing the data does not require the **cSEM** and is
therefore better considered a preparation task, i.e., a “pre-cSEM” task.
The reason why this step is nevertheless considered an explicit part of
the cSEM-workflow is motivated by the experience that applied/causal
users tend to shy away from software like R because “just getting the
data in” and understanding how to show, manipulate and work with data
can be frustrating if one is not aware of R’s rich and easy to learn
data import and data processing capabilities. While these topics may
have been overwhelming for newcomers several years ago, data import and
data transformation have become extremely simple and user-friendly if
the right tools packages are used. The best place to start is the
[Rstudio Cheat sheet webpage](https://posit.co/resources/cheatsheets/),
especially the *Data Import* and the *Data Transformation* cheat sheets.

**cSEM** is relatively flexible as to the type of data accepted.
Currently the following data types/structures are accepted:

1.  A `data.frame` or `tibble` with column names matching the indicator
    names used in the lavaan model description of the measurement or
    composite model. Possible column types or classes of the data
    provided are: `"logical"` (`TRUE`/`FALSE`), `"numeric"` (`"double"`
    or `"integer"`), `"factor"` (`"ordered"` and/or `"unordered"`) or a
    mix of several types. Additionally, the data may also include
    **one** character column whose column name must be given to `.id`.
    Values of this column are interpreted as group identifiers and
    [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
    will split the data by levels of that column and run the estimation
    for each level separately.

    **Example:**

    Assuming the following simple model is to be estimated:

    ``` r
    model <- "
    # Structural model
    EXPE ~ IMAG

    # Reflective measurement model
    EXPE =~ expe1 + expe2
    IMAG =~ imag1 + imag2
    "
    ```

    To estimate the model a data frame with $N$ rows (the observations)
    and $K = 4$ columns with column names `expe1`, `imag1`, `expe2`,
    `imag2` is required. The order of the columns in the dataset is
    irrelevant. In **cSEM** the order is defined by the order in which
    the names appear in the measurement or composite model equations in
    the model description. In this case any resulting matrix or vector
    whose (row/column) names contain the indicator names would have the
    order `expe1`, `expe2`, `imag1`, `imag2`. More one model
    specification below.

2.  A `matrix` with column names matching the indicator names used in
    the lavaan model description of the measurement model or composite
    model description.

3.  A list of data frames or matrices. In this case estimation is
    repeated for each data frame or matrix separately.

The current version 0.5.0 available on CRAN does not provide any tools
to handle missing values. Future versions are likely to include at least
the basic approaches for handling missing values. Regularly check
<https://github.com/FloSchuberth/cSEM/> to get the latest updates.

#### Specify a model

Models are defined using [lavaan model
syntax](https://lavaan.ugent.be/tutorial/syntax1.html). Currently, only
the “standard” lavaan model syntax is supported. This comprises:

1.  The definition of a latent variable/common factor (or more
    precisely: the definition of a concept modeled as a common factor)
    by the “`=~`” operator.
2.  The definition of a composite (or more precisely: the definition of
    a concept modeled as a composite) by the “`<~`” operator.
3.  The specification of regression equations by the “`~`” operator.
4.  The definition of error (co)variances, indicator correlations, or
    correlations between exogenous constructs using the “`~~`” operator.

**cSEM** handles linear, nonlinear and hierarchical models. Syntax for
each model is illustrated below using variables of the build-in
`satisfaction` dataset. For more information see the [lavaan syntax
tutorial](https://lavaan.ugent.be/tutorial/syntax1.html).

##### Linear models

A typical linear model would look like this:

``` r
model <- "
# Structural model
EXPE ~ IMAG
QUAL ~ EXPE
VAL  ~ EXPE + QUAL
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT

# Composite model
IMAG <~ imag1 + imag2 + imag3                  # composite
EXPE <~ expe1 + expe2 + expe3                  # composite
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5  # composite
VAL  <~ val1  + val2  + val3                   # composite

# Reflective measurement model
SAT  =~ sat1  + sat2  + sat3  + sat4           # common factor
LOY  =~ loy1  + loy2  + loy3  + loy4           # common factor

# Measurement error correlation
sat1 ~~ sat2
"
```

Note that the operator `<~` tells **cSEM** that the concept to its left
is modeled as a composite; the operator `=~` tells **cSEM** that the
concept to its left is modeled as a common factor. `~~` tells **cSEM**
that the measurement errors of `sat1` and `sat2` are assumed to
correlate.

##### Nonlinear models

Nonlinear terms are specified as interactions using the dot operator
`"."`. Nonlinear terms include interactions and exponential terms. The
latter is described in model syntax as an “interaction with itself”,
e.g., `x_1^3 = x1.x1.x1`. Currently the following terms are allowed

- Single, e.g., `eta1`
- Quadratic, e.g., `eta1.eta1`
- Cubic, e.g., `eta1.eta1.eta1`
- Two-way interaction, e.g., `eta1.eta2`
- Three-way interaction, e.g., `eta1.eta2.eta3`
- Quadratic and two-way interaction, e.g., `eta1.eta1.eta3`

A simple example would look like this:

``` r
model <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Composite model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"
```

##### Hierarchical (second order) models

Currently only second-order models are supported. Specification of the
second-order construct takes place in the measurement/composite model.

``` r
model <- "
# Structural model
SAT ~ QUAL
VAL ~ SAT + QUAL

# Reflective measurement model
SAT  =~ sat1 + sat2
VAL  =~ val1 + val2

# Composite model
IMAG <~ imag1 + imag2
EXPE <~ expe1 + expe2

# Second-order term
QUAL =~ IMAG + EXPE
"
```

In this case `QUAL` is modeled as a second-order common factor measured
by `IMAG` and `EXPE`, where both `IMAG` and `EXPE` are modeled as
composites.

#### Estimate using `csem()`

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) is the
central function of the package. Although it is possible to estimate a
model using individual functions called by
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) (such
as
[`parseModel()`](https://floschuberth.github.io/cSEM/reference/parseModel.md),
[`processData()`](https://floschuberth.github.io/cSEM/reference/processData.md),
[`calculateWeightsPLS()`](https://floschuberth.github.io/cSEM/reference/calculateWeightsPLS.md),
[`estimatePath()`](https://floschuberth.github.io/cSEM/reference/estimatePath.md)
etc.) using R’s `:::`mechanism for non-exported functions, it is
virtually always easier, safer and quicker to use
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
instead (this is why these functions are not exported).

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
accepts all models and data types described above. The result of a call
to [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)is
always an object of class `cSEMResults`. Technically, the resulting
object has an additional class attribute, namely `cSEMResults_default`,
`cSEMResults_multi` or `cSEMResults_2ndorder` that depends on the type
of model and/or data provided, however, users usually do not need to
worry since postestimation functions automatically work on all classes.

The simplest possible call to
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
involves a data set and a model:

``` r
require(cSEM)

model <- "
# Path model / Regressions
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Reflective measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

a <- csem(.data = threecommonfactors, .model = model)
a
```

    ## ________________________________________________________________________________
    ## ----------------------------------- Overview -----------------------------------
    ## 
    ## Estimation was successful.
    ## 
    ## The result is a list of class cSEMResults with list elements:
    ## 
    ##  - Estimates
    ##  - Information
    ## 
    ## To get an overview or help type:
    ## 
    ##  - ?cSEMResults
    ##  - str(<object-name>)
    ##  - listviewer::jsondedit(<object-name>, mode = 'view')
    ## 
    ## If you wish to access the list elements directly type e.g. 
    ## 
    ##  - <object-name>$Estimates
    ## 
    ## Available postestimation commands:
    ## 
    ##  - assess(<object-name>)
    ##  - infer(<object-name)
    ##  - predict(<object-name>)
    ##  - summarize(<object-name>)
    ##  - verify(<object-name>)
    ## ________________________________________________________________________________

This is equivalent to:

``` r
csem(
   .data                        = threecommonfactors,
   .model                       = model,
   .approach_cor_robust         = "none",
   .approach_nl                 = "sequential",
   .approach_paths              = "OLS",
   .approach_weights            = "PLS-PM",
   .conv_criterion              = "diff_absolute",
   .disattenuate                = TRUE,
   .dominant_indicators         = NULL,
   .estimate_structural         = TRUE,
   .id                          = NULL,
   .iter_max                    = 100,
   .normality                   = FALSE,
   .PLS_approach_cf             = "dist_squared_euclid",
   .PLS_ignore_structural_model = FALSE,
   .PLS_modes                   = NULL,
   .PLS_weight_scheme_inner     = "path",
   .reliabilities               = NULL,
   .starting_values             = NULL,
   .tolerance                   = 1e-05,
   .resample_method             = "none",
   .resample_method2            = "none",
   .R                           = 499,
   .R2                          = 199,
   .handle_inadmissibles        = "drop",
   .user_funs                   = NULL,
   .eval_plan                   = "sequential",
   .seed                        = NULL,
   .sign_change_option          = "no"
    )
```

See the
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
documentation for details on the arguments.

##### Inference

By default, no inferential quantities are calculated since
composite-based approaches, generally, do not have closed-form solutions
for standard errors. **cSEM** relies on the `bootstrap` or `jackknife`
to estimate standard errors, test statistics, critical quantiles, and
confidence intervals.

**cSEM** offers two ways to compute resamples:

1.  Inference can be done by first setting argument `.resample_method`
    to `"jackkinfe"` or `"bootstrap"` to perform resampling and
    subsequently use
    [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
    (or more conveniently
    [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
    which internally calls
    [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md))
    to compute the actual inferential quantities of interest.
2.  The same result is achieved by passing a `cSEMResults` object to
    [`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)
    and subsequently using
    [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
    or
    [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).

``` r
b1 <- csem(.data = threecommonfactors, .model = model, .resample_method = "bootstrap")
b2 <- resamplecSEMResults(a)
```

Several confidence intervals are implemented, see `?infer()`:

``` r
summarize(b1)
```

    ## ________________________________________________________________________________
    ## ----------------------------------- Overview -----------------------------------
    ## 
    ##  General information:
    ##  ------------------------
    ##  Estimation status                  = Ok
    ##  Number of observations             = 500
    ##  Weight estimator                   = PLS-PM
    ##  Inner weighting scheme             = "path"
    ##  Type of indicator correlation      = Pearson
    ##  Path model estimator               = OLS
    ##  Second-order approach              = NA
    ##  Type of path model                 = Linear
    ##  Disattenuated                      = Yes (PLSc)
    ## 
    ##  Resample information:
    ##  ---------------------
    ##  Resample method                    = "bootstrap"
    ##  Number of resamples                = 499
    ##  Number of admissible results       = 499
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 39337545
    ## 
    ##  Construct details:
    ##  ------------------
    ##  Name  Modeled as     Order         Mode      
    ## 
    ##  eta1  Common factor  First order   "modeA"   
    ##  eta2  Common factor  First order   "modeA"   
    ##  eta3  Common factor  First order   "modeA"   
    ## 
    ## ----------------------------------- Estimates ----------------------------------
    ## 
    ## Estimated path coefficients:
    ## ============================
    ##                                                              CI_percentile   
    ##   Path           Estimate  Std. error   t-stat.   p-value         95%        
    ##   eta2 ~ eta1      0.6713      0.0457   14.6890    0.0000 [ 0.5826; 0.7559 ] 
    ##   eta3 ~ eta1      0.4585      0.0806    5.6880    0.0000 [ 0.3037; 0.6182 ] 
    ##   eta3 ~ eta2      0.3052      0.0847    3.6007    0.0003 [ 0.1382; 0.4668 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                              CI_percentile   
    ##   Loading        Estimate  Std. error   t-stat.   p-value         95%        
    ##   eta1 =~ y11      0.6631      0.0415   15.9909    0.0000 [ 0.5746; 0.7423 ] 
    ##   eta1 =~ y12      0.6493      0.0377   17.2163    0.0000 [ 0.5759; 0.7164 ] 
    ##   eta1 =~ y13      0.7613      0.0344   22.1211    0.0000 [ 0.6924; 0.8222 ] 
    ##   eta2 =~ y21      0.5165      0.0547    9.4455    0.0000 [ 0.4037; 0.6144 ] 
    ##   eta2 =~ y22      0.7554      0.0379   19.9465    0.0000 [ 0.6784; 0.8220 ] 
    ##   eta2 =~ y23      0.7997      0.0371   21.5773    0.0000 [ 0.7244; 0.8742 ] 
    ##   eta3 =~ y31      0.8223      0.0344   23.9028    0.0000 [ 0.7538; 0.8817 ] 
    ##   eta3 =~ y32      0.6581      0.0385   17.1070    0.0000 [ 0.5774; 0.7277 ] 
    ##   eta3 =~ y33      0.7474      0.0398   18.7894    0.0000 [ 0.6718; 0.8295 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                              CI_percentile   
    ##   Weight         Estimate  Std. error   t-stat.   p-value         95%        
    ##   eta1 <~ y11      0.3956      0.0216   18.3009    0.0000 [ 0.3511; 0.4393 ] 
    ##   eta1 <~ y12      0.3873      0.0192   20.2160    0.0000 [ 0.3498; 0.4257 ] 
    ##   eta1 <~ y13      0.4542      0.0204   22.2910    0.0000 [ 0.4168; 0.4964 ] 
    ##   eta2 <~ y21      0.3058      0.0277   11.0436    0.0000 [ 0.2472; 0.3555 ] 
    ##   eta2 <~ y22      0.4473      0.0223   20.0708    0.0000 [ 0.4067; 0.4944 ] 
    ##   eta2 <~ y23      0.4735      0.0224   21.1640    0.0000 [ 0.4343; 0.5211 ] 
    ##   eta3 <~ y31      0.4400      0.0185   23.7537    0.0000 [ 0.4035; 0.4751 ] 
    ##   eta3 <~ y32      0.3521      0.0176   20.0198    0.0000 [ 0.3165; 0.3879 ] 
    ##   eta3 <~ y33      0.3999      0.0200   20.0445    0.0000 [ 0.3650; 0.4434 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   eta2 ~ eta1       0.6713      0.0457   14.6890    0.0000 [ 0.5826; 0.7559 ] 
    ##   eta3 ~ eta1       0.6634      0.0396   16.7646    0.0000 [ 0.5852; 0.7394 ] 
    ##   eta3 ~ eta2       0.3052      0.0847    3.6007    0.0003 [ 0.1382; 0.4668 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   eta3 ~ eta1          0.2049      0.0574    3.5709    0.0004 [ 0.0990; 0.3160 ] 
    ## ________________________________________________________________________________

Or directly via
[`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)

``` r
ii <- infer(b1, .quantity = c("CI_standard_z", "CI_percentile"), .alpha = c(0.01, 0.05))
ii$Path_estimates
```

    ## $CI_standard_z
    ##      eta2 ~ eta1 eta3 ~ eta1 eta3 ~ eta2
    ## 99%L   0.5563689   0.2510462  0.08354111
    ## 99%U   0.7918162   0.6663207  0.52013015
    ## 95%L   0.5845159   0.3006910  0.13573402
    ## 95%U   0.7636692   0.6166759  0.46793724
    ## 
    ## $CI_percentile
    ##      eta2 ~ eta1 eta3 ~ eta1 eta3 ~ eta2
    ## 99%L   0.5549906   0.2490596   0.1013221
    ## 99%U   0.7765964   0.6529602   0.5263406
    ## 95%L   0.5825700   0.3037477   0.1381978
    ## 95%U   0.7559078   0.6182413   0.4667838

Both bootstrap and jackknife resampling support platform-independent
multiprocessing as well as random seeds via the [future
framework](https://github.com/futureverse/future/). For multiprocessing
simply set `.eval_plan = "multiprocess"` in which case the maximum
number of available cores is used if not on Windows. On Windows as many
separate R instances are opened in the background as there are cores
available instead. Note that this naturally has some overhead.
Consequently, for a small number of resamples multiprocessing will
generally not be faster compared to sequential (single core) processing
(the default). Seeds are set via the `.seed` argument. A typical call
would look like this:

``` r
b <- csem(
  .data            = satisfaction,
  .model           = model,
  .resample_method = "bootstrap",
  .R               = 999,
  .seed            = 98234,
  .eval_plan       = "multiprocess")

# Output omitted
```

#### Apply postestimation functions

There are 5 major postestimation function and 4 test-family functions:

- [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md):

  Assess the quality of the estimated model *without conducting a
  statistical test*. Quality in this case is taken to be a catch-all
  term for all common aspects of model assessment. This mainly comprises
  fit indices, reliability estimates, common validity assessment
  criteria and other related quality measures/indices that do not rely
  on a formal test procedure. In **cSEM** a generic (fit) index or
  quality/assessment measure is referred to as a **quality criterion**.

- [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md):

  Calculate common inferential quantities. For users interested in the
  estimated standard errors and/or confidences intervals
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
  will usually be more helpful as it has a much more user-friendly print
  method.

- [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md):

  Predict indicator scores of endogenous constructs based on the
  procedure introduced by Shmueli et al. (2016).

- [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md):

  Summarize a model. The function is mainly called for its side effect,
  the printing of a structured summary of the estimates. It also
  provides most estimates in user-friendly data frames. The data frame
  format is usually much more convenient if users intend to present the
  results in e.g., a paper or a presentation.

- [`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md):

  Verify admissibility of the estimated quantities for a given model.
  Results based on an estimated model exhibiting one of the following
  defects are deemed inadmissible: non-convergence, loadings and/or
  (congeneric) reliabilities larger than 1, a construct VCV and/or a
  model-implied VCV matrix that is not positive (semi-)definite.

##### The `test_*` family of postestimation functions

- [`testHausman()`](https://floschuberth.github.io/cSEM/reference/testHausman.md):

  The regression-based Hausman test for SEM.

- [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md):

  Test for overall model fit based on Beran and Srivastava (1985). See
  also Dijkstra and Henseler (2015).

- [`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md):

  Test for group differences using several different approaches such as
  e.g., the one described in Klesel et al. (2019).

- [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md):

  Test of measurement invariance of composites proposed by Henseler,
  Ringle, and Sarstedt (2016)

##### The `do_*` family of postestimation functions

- [`doIPMA()`](https://floschuberth.github.io/cSEM/reference/doIPMA.md):

  Performs an importance-performance matrix analysis (IPMA).

- [`doNonlinearEffectsAnalysis()`](https://floschuberth.github.io/cSEM/reference/doNonlinearEffectsAnalysis.md):

  Performs nonlinear effects analysis such as floodlight and surface
  analysis as described in e.g., Spiller et al. (2013).

- [`doRedundancyAnalysis()`](https://floschuberth.github.io/cSEM/reference/doRedundancyAnalysis.md):

  Performs a redundancy analysis (RA) as proposed by Hair et al. (2016)
  with reference to Chin (1998).

Technically, postestimation functions are generic function with methods
for objects of class `cSEMResults_default`, `cSEMResults_multi`,
`cSEMResults_2ndorder`. In **cSEM** every `cSEMResults_*` object must
also have class `cSEMResults` for internal reasons. When using one of
the major postestimation functions, method dispatch is therefore
technically done on one of the `cSEMResults_*` class attributes,
ignoring the `cSEMResults` class attribute. As long as a postestimation
function is used directly method dispatch is not of any practical
concern to the end-user. The difference, however, becomes important if a
user seeks to directly invoke an internal function which is called by
one of the postestimation functions (e.g.,
[`calculateAVE()`](https://floschuberth.github.io/cSEM/reference/calculateAVE.md)
or
[`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)
as called by
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)).
In this case, only objects of class `cSEMResults_default` are accepted
as this ensures a specific structure. Therefore, it is important to
remember that *internal functions are generally **not** generic.*

### Principles underlying cSEM

**cSEM** is based on a number of principles, that have shaped its
design, terminology and scope. These principles are discussed below

#### Model vs. Estimation

The way different concepts and their relationship are *modeled* is
strictly distinct from how they are *estimated*. Hence we strictly
distinguish between concepts modeled as common factors (or composites)
and the actual estimation **for a given model**. In our opinion, these
differences are fundamental to understanding the scope and limits of a
certain approach. The most notable consequence is that approaches such
as partial least squares and everything related to it (e.g., the modes)
or generalized structured component analysis are “only” considered as
estimators/estimation approaches **for a given model**.

#### Composites in a composite model vs. composites in a common factor model and disattenuation

By virtue of the package, **cSEM** uses composite-based
estimators/approaches only. Depending on the postulated model, linear
compounds may therefore either serve as a composite as part of the
composite model or as a proxy/stand-in for a common factor. If a concept
is modeled as a common factor, proxy correlations, proxy-indicator
correlations and path coefficients are inconsistent estimates for their
supposed construct level counterparts (construct correlations, loadings
and path coefficients) unless the proxy is a perfect representation of
its construct level counterpart. This is commonly referred to as
**attenuation bias**. Several approaches have been suggested to correct
for these biases. In **cSEM** estimates are correctly disattenuated by
default if any of the concepts involved is modeled as a common factor!
Disattenuation is controlled by the `.disattenuate` argument of
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

**Example**

``` r
model <- "
## Structural model
eta2 ~ eta1

## Measurement model
eta1 <~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
"

# Identical
csem(threecommonfactors, model)
csem(threecommonfactors, model, .disattenuate = TRUE)

# To supress automatic disattenuation
csem(threecommonfactors, model, .disattenuate = FALSE)
```

Note that since `.approach_weights = "PLS-PM"` and
`.disattentuate = TRUE` by default (see for [The role of the weighting
scheme and partial least squares (PLS)](#roleofpls) below) and one of
the concepts in the model above is modeled as a common factor, composite
(proxy) correlations, loadings and path coefficients are adequately
disattenuated using the correction approach commonly known as
**consistent partial least squares (PLSc)**. If `.disattenuate = FALSE`
or all concepts are modeled as composites “proper” PLS values are
returned.

#### The role of the weighting scheme and partial least squares (PLS)

In principal, any weighted combination of appropriately chosen
observables can be used to estimate structural relationships between
these compounds. Hence, any conceptual or methodological issue discussed
based on a composite build by a given (weighting) approach may equally
well be discussed for any other potential weighting scheme. The
appropriateness or potential superiority of a specific weighting
approach such as “partial least squares path modeling” (PLS-PM) over
another such as “unit weights” (sum scores) or generalized structured
component analysis (GSCA) is therefore to some extent a question of
*relative* appropriateness and *relative* superiority.

As a notable consequence, we believe that well known approaches such
partial least squares path modeling (PLS-PM) and generalized structured
component analysis (GSCA) are - contrary to common belief - best
exclusively understood as prescriptions for forming linear compounds
based on observables, i.e., as weighting approaches. Not more, not
less.[¹](#fn1) In **cSEM** this is reflected by the fact that `"PLS"`
and `"GSCA"` are choices of the `.approach_weights` argument.

``` r
model <- "
## Structural model
eta2 ~ eta1

## Composite model
eta1 <~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
"

### Currently the following weight approaches are implemented
# Partial least squares path modeling (PLS)
csem(threecommonfactors, model, .approach_weights = "PLS-PM") # default

# Generalized canonical correlation analysis (Kettenring approaches)
csem(threecommonfactors, model, .approach_weights = "SUMCORR")
csem(threecommonfactors, model, .approach_weights = "MAXVAR")
csem(threecommonfactors, model, .approach_weights = "SSQCORR")
csem(threecommonfactors, model, .approach_weights = "MINVAR")
csem(threecommonfactors, model, .approach_weights = "GENVAR")

# Generalized structured component analysis (GSCA)
csem(threecommonfactors, model, .approach_weights = "GSCA")

# Principal component analysis (PCA)
csem(threecommonfactors, model, .approach_weights = "PCA")

# Factor score regression (FSR) using "unit", "bartlett" or "regression" weights
csem(threecommonfactors, model, .approach_weights = "unit")
csem(threecommonfactors, model, .approach_weights = "bartlett")
csem(threecommonfactors, model, .approach_weights = "regression")
```

## Literature

Beran, Rudolf, and Muni S. Srivastava. 1985. “Bootstrap Tests and
Confidence Regions for Functions of a Covariance Matrix.” *The Annals of
Statistics* 13 (1): 95–115. <https://doi.org/10.1214/aos/1176346579>.

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley-Interscience.

Chin, W. W. 1998. “Modern Methods for Business Research.” In, edited by
G. A. Marcoulides, 295–358. Mahwah, NJ: Lawrence Erlbaum.

Dijkstra, Theo K., and Jörg Henseler. 2015. “Consistent and
Asymptotically Normal PLS Estimators for Linear Structural Equations.”
*Computational Statistics & Data Analysis* 81: 10–23.

Hair, Joseph F, G Tomas M Hult, Christian Ringle, and Marko Sarstedt.
2016. *A Primer on Partial Least Squares Structural Equation Modeling
(PLS-SEM)*. Sage publications.

Henseler, Jörg. 2017. “Bridging Design and Behavioral Research with
Variance-Based Structural Equation Modeling.” *Journal of Advertising*
46 (1): 178–92. <https://doi.org/10.1080/00913367.2017.1281780>.

Henseler, Jörg, Christian M. Ringle, and Marko Sarstedt. 2016. “Testing
Measurement Invariance of Composites Using Partial Least Squares.”
*International Marketing Review* 33 (3): 405–31.
<https://doi.org/10.1108/imr-09-2014-0304>.

Klesel, Michael, Florian Schuberth, Jörg Henseler, and Bjoern Niehaves.
2019. “A Test for Multigroup Comparison Using Partial Least Squares Path
Modeling.” *Internet Research* 29 (3): 464–77.
<https://doi.org/10.1108/intr-11-2017-0418>.

Rigdon, Edward E. 2012. “Rethinking Partial Least Squares Path Modeling:
In Praise of Simple Methods.” *Long Range Planning* 45 (5-6): 341–58.
<https://doi.org/10.1016/j.lrp.2012.09.010>.

———. 2016. “Choosing PLS Path Modeling as Analytical Method in European
Management Research: A Realist Perspective.” *European Management
Journal* 34 (6). <https://doi.org/10.1016/j.emj.2016.05.006>.

Rigdon, Edward E., Jan-Michael Becker, and Marko Sarstedt. 2019. “Factor
Indeterminacy as Metrological Uncertainty: Implications for Advancing
Psychological Measurement.” *Multivariate Behavioral Research*, 1–15.
<https://doi.org/10.1080/00273171.2018.1535420>.

Shmueli, Galit, Soumya Ray, Juan Manuel Velasquez Estrada, and Suneel
Babu Chatla. 2016. “The Elephant in the Room: Predictive Performance of
PLS Models.” *Journal of Business Research* 69 (10): 4552–64.
<https://doi.org/10.1016/j.jbusres.2016.03.049>.

Spiller, Stephen A., Gavan J. Fitzsimons, John G. Lynch, and Gary H.
Mcclelland. 2013. “Spotlights, Floodlights, and the Magic Number Zero:
Simple Effects Tests in Moderated Regression.” *Journal of Marketing
Research* 50 (2): 277–88. <https://doi.org/10.1509/jmr.12.0420>.

------------------------------------------------------------------------

1.  In fact, labels such as PLS-PM and even more so PLS-SEM are
    misleading as they create the impression that PLS(-PM) is somehow
    capable of more than other composite-based approaches. While among
    composite-based approaches, methodological research surrounding
    *composites formed using weights obtained by the PLS(-PM) algorithm*
    is most advanced, the PLS algorithm remains a weighting scheme in
    its core.
