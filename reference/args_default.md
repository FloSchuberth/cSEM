# Show argument defaults or candidates

Show all arguments used by package functions including default or
candidate values. For argument descriptions see:
[csem_arguments](https://floschuberth.github.io/cSEM/reference/csem_arguments.md).

## Usage

``` r
args_default(.choices = FALSE)
```

## Arguments

- .choices:

  Logical. Should candidate values for the arguments be returned?
  Defaults to `FALSE`.

## Value

A named list of argument names and defaults or accepted candidates.

## Details

By default `args_default()`returns a list of default values by argument
name. If the list of accepted candidate values is required instead, use
`.choices = TRUE`.

## See also

[`handleArgs()`](https://floschuberth.github.io/cSEM/reference/handleArgs.md),
[csem_arguments](https://floschuberth.github.io/cSEM/reference/csem_arguments.md),
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[`foreman()`](https://floschuberth.github.io/cSEM/reference/foreman.md)
