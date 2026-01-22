# Internal: Function that mutates a vector

Function that randomly mutates one of the elements of the vectorized
structural model by flipping it from 0 to 1 or vice versa) taking into
account that not every element is allowed to be mutated, The elements in
the rows of exogenous constructs and the elements on the diagonal are
not mutated. Similarly, elements are not mutated when this would result
in feedback loops. This function adheres to the necessary requirements
for a function to be supplied to
[`ga`](https://github.com/luca-scr/GA/reference/ga.html).

## Usage

``` r
mutateVector(.object,
                   .parent,
                   .model_org)
```

## Value

A numeric vector; Vectorized version of the structural model in which
one element has been flipped from 0 to 1 or vice versa
