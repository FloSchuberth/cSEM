# Data: Second order common factor of composites

A dataset containing 500 standardized observations on 19 indicator
generated from a population model with 6 concepts, three of which
(`c1-c3`) are composites forming a second order common factor (`c4`).
The remaining two (`eta1`, `eta2`) are concepts modeled as common
factors .

## Usage

``` r
dgp_2ndorder_cf_of_c
```

## Format

A matrix with 500 rows and 19 variables:

- y11-y12:

  Indicators attached to `c1`. Population weights are: 0.8; 0.4.
  Population loadings are: 0.925; 0.65

- y21-y24:

  Indicators attached to `c2`. Population weights are: 0.5; 0.3; 0.2;
  0.4. Population loadings are: 0.804; 0.68; 0.554; 0.708

- y31-y38:

  Indicators attached to `c3`. Population weights are: 0.3; 0.3; 0.1;
  0.1; 0.2; 0.3; 0.4; 0.2. Population loadings are: 0.496; 0.61; 0.535;
  0.391; 0.391; 0.6; 0.5285; 0.53

- y41-y43:

  Indicators attached to `eta1`. Population loadings are: 0.8; 0.7; 0.7

- y51-y53:

  Indicators attached to `eta1`. Population loadings are: 0.8; 0.8; 0.7

The model is: \$\$\`c4\` = gamma1 \* \`eta1\` + zeta1\$\$ \$\$\`eta2\` =
gamma2 \* \`eta1\` + beta \* \`c4\` + zeta2\$\$

with population values `gamma1` = 0.6, `gamma2` = 0.4 and `beta` = 0.35.
The second order common factor is \$\$\`c4\` = lambdac1 \* \`c1\` +
lambdac2 \* \`c2\` + lambdac3 \* \`c3\` + epsilon\$\$
