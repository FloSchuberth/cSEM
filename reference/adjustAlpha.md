# Internal: Multiple testing correction

Adjust a given significance level `.alpha` to accommodate multiple
testing. The following corrections are implemented:

- `none`:

  (Default) No correction is done.

- `bonferroni`:

  A Bonferroni correction is done, i.e., alpha is divided by the number
  of comparisons `.nr_comparisons`.

## Usage

``` r
adjustAlpha(
 .alpha                 = args_default()$.alpha,
 .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
 .nr_comparisons        = args_default()$.nr_comparisons
)
```

## Arguments

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .approach_alpha_adjust:

  Character string. Approach used to adjust the significance level to
  accommodate multiple testing. One of "*none*" or "*bonferroni*".
  Defaults to "*none*".

- .nr_comparisons:

  Integer. The number of comparisons. Defaults to `NULL`.

## Value

A vector of (possibly adjusted) significance levels.
