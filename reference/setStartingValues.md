# Internal: Set starting values

Set the starting values.

## Usage

``` r
setStartingValues(
  .W               = args_default()$.W,
  .starting_values = args_default()$.starting_values
  )
```

## Arguments

- .W:

  A (J x K) matrix of weights.

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

## Value

The (J x K) matrix of starting values.
