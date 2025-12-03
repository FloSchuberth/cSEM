# Internal: Classify structural model terms by type

Classify terms of the structural model according to their type.

## Usage

``` r
classifyConstructs(.terms = args_default()$.terms)
```

## Arguments

- .terms:

  A vector of construct names to be classified.

## Value

A named list of length equal to the number of terms provided containing
a data frame with columns "*Term_class*", "*Component*",
"*Component_type*", and "*Component_freq*".

## Details

Classification is required to estimate nonlinear structural
relationships. Currently the following terms are supported

- Single, e.g., `eta1`

- Quadratic, e.g., `eta1.eta1`

- Cubic, e.g., `eta1.eta1.eta1`

- Two-way interaction, e.g., `eta1.eta2`

- Three-way interaction, e.g., `eta1.eta2.eta3`

- Quadratic and two-way interaction, e.g., `eta1.eta1.eta3`

Note that exponential terms are modeled as "interactions with itself" as
in i.e., `eta1^3 = eta1.eta1.eta1`.
