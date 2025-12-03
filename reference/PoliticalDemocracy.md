# Data: political democracy

The Industrialization and Political Democracy dataset. This dataset is
used throughout Bollen's 1989 book (see pages 12, 17, 36 in chapter 2,
pages 228 and following in chapter 7, pages 321 and following in chapter
8; Bollen (1989) ). The dataset contains various measures of political
democracy and industrialization in developing countries.

## Usage

``` r
PoliticalDemocracy
```

## Format

A data frame of 75 observations of 11 variables.

- `y1`:

  Expert ratings of the freedom of the press in 1960

- `y2`:

  The freedom of political opposition in 1960

- `y3`:

  The fairness of elections in 1960

- `y4`:

  The effectiveness of the elected legislature in 1960

- `y5`:

  Expert ratings of the freedom of the press in 1965

- `y6`:

  The freedom of political opposition in 1965

- `y7`:

  The fairness of elections in 1965

- `y8`:

  The effectiveness of the elected legislature in 1965

- `x1`:

  The gross national product (GNP) per capita in 1960

- `x2`:

  The inanimate energy consumption per capita in 1960

- `x3`:

  The percentage of the labor force in industry in 1960

## Source

The [lavaan](https://lavaan.ugent.be/) package (version 0.6-3).

## References

Bollen KA (1989). *Structural Equations with Latent Variables*.
Wiley-Interscience. ISBN 978-0471011712.

## Examples

``` r
#============================================================================
# Example is taken from the lavaan website
#============================================================================
# Note: example is modified. Across-block correlations are removed
model <- "
# Measurement model
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  
# Regressions / Path model
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  
# residual correlations
  y2 ~~ y4
  y6 ~~ y8
"

aa <- csem(PoliticalDemocracy, model)
```
