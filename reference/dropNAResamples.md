# Internal: Drop permutation runs containing missing values

Removes those elements of a list of permutation (or resample) test
statistics that contain at least one `NA` or `NaN` anywhere in their
structure, including inside nested lists.

## Usage

``` r
dropNAResamples(.ref_dist = NULL)
```

## Arguments

- .ref_dist:

  A list of test statistics, one element per permutation run.

## Value

A list containing only those elements of `.ref_dist` that are free of
`NA` and `NaN`. Names and relative order are preserved.

## Details

Elements that are the literal `NA` sentinel (written when a permutation
run is inadmissible and `.handle_inadmissibles = "drop"`) are removed as
well.
