# Internal: Handle arguments

Internal helper function to handle arguments passed to any function
within `cSEM`.

## Usage

``` r
handleArgs(.args_used)
```

## Arguments

- .args_used:

  A list of argument names and user picked values.

## Value

The
[args_default](https://floschuberth.github.io/cSEM/reference/args_default.md)
list, with default values changed to the values given by the user.
