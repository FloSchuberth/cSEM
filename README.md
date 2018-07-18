
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM

WARNING: THIS IS WORK IN PROGRESS. BREAKING ARE VERY LIKELY. Do not use
the package before the first stable relase (which will be 0.0.1, towards
the end of 2018).

## Purpose

Estimate, analyse, test, and study linear and nonlinear structural
equation models using composite based approaches, procedures, and tests
including e.g. PLS, GSCA, 2SLS estimation, several tests and typical
postestimation procedures.

## Installation

``` r
# Currently only a development version from GitHub is available:
# install.packages("devtools")
devtools::install_github("M-E-Steiner/cSEM")
```

## Philopsophy/Goals/Ideas

  - Easy to use by non-R experts:
      - Functions `csem()` and `cca()` provide default choices for most
        of its arguments (similarity to the `sem` and `cfa` functions of
        the [lavaan](http://lavaan.ugent.be/) package is intended).
      - Well documented (Vignettes, HTML output, a website, intro
        course(s)). Of course this may take some time\!
      - There will be an extensive (non-expert) visually and
        didactically appealing documentation designed to make the
        learning curve of both the methods involved and the package as
        flat as possible.
      - Structured output/results that aims to be “easy”" in a sense
        that it is - … descriptive/verbose - … easy to export to other
        environments such as MS Word, Latex files etc. (exportability) -
        … easy to migrate from/to/between other PLS/VB/CB-based systems
        (lavaan, semPLS, ADANCO, SmartPLS) (this will also take a lot of
        time\!)
      - (In the future) Intro courses, accompaning website, cheatsheets.
  - The package is designed to be flexible/modular enough so that
    researchers developing new methods can take specific function
    provided by the package and alter them according to their need
    without working their way through a chain of other functions
    (naturally this will not always be possible). Modularity is largly
    inspired by the `matrixpls` package.
  - Modern in a sense that the package integrates modern developments
    within the R community. This mainly includes
    ideas/recommendations/design choices that fead into the packages of
    the [tidyverse](https://github.com/tidyverse/tidyverse).

## Basic usage

    require(cSEM)
    data(satisfaction)
    
    model <- "
    # Structural model
    EXPE ~ IMAG
    QUAL ~ EXPE
    VAL  ~ EXPE + QUAL + QUAL.EXPE
    SAT  ~ IMAG + EXPE + QUAL + VAL + IMAG.IMAG
    LOY  ~ IMAG + SAT
    
    # Measurement model
    
    IMAG <~ imag1 + imag2 + imag3
    EXPE <~ expe1 + expe2 + expe3 
    QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
    VAL  <~ val1  + val2  + val3
    SAT  =~ sat1  + sat2  + sat3  + sat4
    LOY  =~ loy1  + loy2  + loy3  + loy4
    "
    
    a <- csem(.data = satisfaction, .model = model)
    a
    
    summary(a) # summary is work in progress
