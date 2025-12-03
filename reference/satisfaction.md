# Data: satisfaction

A data frame with 250 observations and 27 variables. Variables from 1 to
27 refer to six latent concepts: `IMAG`=Image, `EXPE`=Expectations,
`QUAL`=Quality, `VAL`=Value, `SAT`=Satisfaction, and `LOY`=Loyalty.

- imag1-imag5:

  Indicators attached to concept `IMAG` which is supposed to capture
  aspects such as the institutions reputation, trustworthiness,
  seriousness, solidness, and caring about customer.

- expe1-expe5:

  Indicators attached to concept `EXPE` which is supposed to capture
  aspects concerning products and services provided, customer service,
  providing solutions, and expectations for the overall quality.

- qual1-qual5:

  Indicators attached to concept `QUAL` which is supposed to capture
  aspects concerning reliability of products and services, the range of
  products and services, personal advice, and overall perceived quality.

- val1-val4:

  Indicators attached to concept `VAL` which is supposed to capture
  aspects related to beneficial services and products, valuable
  investments, quality relative to price, and price relative to quality.

- sat1-sat4:

  Indicators attached to concept `SAT` which is supposed to capture
  aspects concerning overall rating of satisfaction, fulfillment of
  expectations, satisfaction relative to other banks, and performance
  relative to customer's ideal bank.

- loy1-loy4:

  Indicators attached to concept `LOY` which is supposed to capture
  aspects concerning propensity to choose the same bank again,
  propensity to switch to other bank, intention to recommend the bank to
  friends, and the sense of loyalty.

## Usage

``` r
satisfaction
```

## Format

An object of class `data.frame` with 250 rows and 27 columns.

## Source

The [plspm](https://github.com/gastonstat/plspm/) package (version
0.4.9). Original source according to plspm: "Laboratory of Information
Analysis and Modeling (LIAM). Facultat d'Informatica de Barcelona,
Universitat Politecnica de Catalunya".

## Details

This dataset contains the variables from a customer satisfaction study
of a Spanish credit institution on 250 customers. The data is identical
to the dataset provided by the
[plspm](https://github.com/gastonstat/plspm/) package but with the last
column (`gender`) removed. If you are looking for the original dataset
use the
[satisfaction_gender](https://floschuberth.github.io/cSEM/reference/satisfaction_gender.md)
dataset.
