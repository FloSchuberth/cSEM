#set page(paper: "us-letter", margin: 1in)
#set text(font: "Times New Roman", size: 12pt)
#set par(leading: 1em, first-line-indent: 0.5in, spacing: 1.5em)
#show heading: it => { it + par("") }
#set heading(numbering: "1.")
#set quote(block: true, quotes: true)
#show quote: set pad(x: 5em)
#show quote: it => emph(it)
#set math.equation(numbering: "(1)")
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#show: codly-init.with()
#codly(languages: codly-languages, stroke: 1pt + red)
#import "@preview/iridis:0.1.0" // For colored parentheses
// iridis-show also installs a `show math.equation` rule that is broken on Typst
// 0.14 (it reads a removed `accent.size` field, dying on accents like overline()).
// Apply only its code-parenthesis colorizer, leaving math untouched.
#show raw: iridis.internals.colorize-code(
  state("parenthesis", 0),
  ("(", "[", "{"), (")", "]", "}"),
  iridis.iridis-palette,
)

#align(center)[
  *Different Ways of Using Variance Explained in Model Comparison*

  #v(0.5em)
  Michael S. Truong

  York University

]


This document discusses different ways of using variance explained $(R^2)$ as a means of model comparison. The primary aim is to reproduce and evaluate the model comparison procedure described in Hwang and Takane (2014) using the toy example of multiple linear regression (MLR). The secondary aim is to compare the replicated procedure against other potential approaches in terms of Type I error rates and statistical. 

The direct quotes from Hwang & Takane (2014) are as follows: #quote(attribution: [p.27])[Similar to bootstrapping R-squared in linear regression (Ohtani 2000), we can compute the bootstrapped standard errors or confidence intervals of the difference in FIT between two models so as to examine whether there is a statistically significant difference in the FIT values. This procedure can be regarded as a nonparametric version of the paired $t$ test for two groups of FIT values.] 

and 

#quote(attribution: [p. 41])[Furthermore, we conducted a paired $t$-test for the FIT values between the original and simpler models. As described earlier, this test was equivalent to testing whether the difference in the FIT values calculated from 100 bootstrap samples for each model was equal to zero. We found the $t$-statistic nonsignificant statistically [$t$(99) = −0.59, $p$ = 0.55, 95% CI = −0.00−0.00], indicating that there was no statistically significant difference in FIT between the original and simpler models.]

The above procedure is problematic because it is dependent on the number of bootstrap samples $(B)$ and does not appropriately sample from the null distribution of no diffference. As a reminder, the formula for a one-sample _t_-test is#footnote[Note that the paired _t_-test is equivalent to the formula for the one-sample _t_-test $()$. See Hogg et al. (2019, page 280), Wackerly et al. (2008, page 521, 646 to 648) and #raw("getAnywhere('t.test.default')", lang: "r") for more information.]  

$
 t = (overline(X) - mu_0)/(S"/"sqrt(n))
$

and is the same formula used for the paired _t_-test.#footnote[I note the equivalence between the formulas because something that confused me was whether the paired _t_-test comparing the pairs of bootstrap differences in FIT statistic would be equivalent to the one-sample _t_-test of the differences against 0. They are equivalent.] Specifically, the paired _t_-test compares the mean of the paired differences $(overline(X)=(x_(1n) - x_(2n))/n)$ against $mu_0$ divided by the standard error of the paired differences $S/sqrt(n)$. Additionally, in this context, $n = B$ where $B$ is the number of bootstraps, therefore the (paired) _t_-statistic becomes arbitrarily large as you increase the number of bootstrap samples. This is a highly undesirable property of the procedure described by Hwang and Takane (2014). Also, the _df_ (degrees of freedom) of the paired _t_-test is dependent on the number of bootstrap samples.

= General Method

This is a two part simulation study. The first part uses the simple example of MLR in OLS to evaluate the original procedure by Hwang and Takane (2014) against other alternatives. The second part extends the simulation to the case of GSCA after a general literature review.

Claude Opus 4.8 helped suggest what conditions to compare for Study 1. Both simulations were coded with the assistance of Claude Opus 4.8. "All AI-assisted code was thoroughly reviewed, modified and independently verified by the author prior to execution."#footnote[Direct quote of #link("http://arxiv.org/abs/2603.04003")]

= Study 1: Multiple Linear Regression

== Method

=== Conditions

The fixed conditions are

```R
sim_fixed <- list(
  beta0 = 0, beta1 = 0.3, beta2 = 0.3, sigma = 1, N = 200, alpha = 0.05,
  fit_compare = fit_compare                                  
)
```

which corresponds to @tbl-fixed-1.

#figure(
  table(
    columns: 2,
    align: (left, center),
    table.header([*Parameter*], [*Value*]),
    $beta_0$,  $0$,
    $beta_1$,  $0.3$,
    $beta_2$,  $0.3$,
    $sigma$,   $1$,
    $N$,       $200$,
    $alpha$,   $0.05$,
  ),
  caption: [Fixed parameters of the data-generating process.],
) <tbl-fixed-1>


The varying conditions are

```R
SimDesign::createDesign(
  delta = c(0, 0.1, 0.2, 0.4),
  B     = c(100, 1000),   # 100 reproduces H&T's df = 99
  prop1 = c(0.5, 0.7, 0.9) # fraction of N in part 1: 100/100, 140/60, 180/20
)
```

which corresponds to @tbl-varying-1.

#figure(
  table(
    columns: 2,
    align: (left, left),
    table.header([*Factor*], [*Levels*]),
    $delta$,   [$0$, $0.1$, $0.2$, $0.4$],
    $B$,       [$100$, $1000$],
    `prop1`,   [$0.5$, $0.7$, $0.9$],
  ),
  caption: [Simulation design grid. $B = 100$ reproduces Hwang and Takane (2014)'s $italic("df") = 99$; `prop1` is the fraction of $N$ in part 1 (splits of 100/100, 140/60, 180/20).],
) <tbl-varying-1>


=== Data Generating Process

A simple two-predictor data for a normally distributed response variable is generated. The coefficients $beta_1$ and $beta_2$ are incremented by $delta$ based on whether case $i$ is equal to group 2 or not. 

```R
Generate <- function(condition, fixed_objects) {
  N     <- fixed_objects$N
  # `prop1` = fraction of N in part 1 (0.5 = balanced). Clamp so both parts stay
  # estimable (need n >= the 3 regression coefficients).
  n1    <- min(N - 4L, max(4L, round(N * condition$prop1)))
  group <- factor(rep(c(1L, 2L), c(n1, N - n1)))
  x1    <- rnorm(N)
  x2    <- rnorm(N)
  d2    <- as.integer(group == 2L)                           # 1 in part 2, else 0
  b1    <- fixed_objects$beta1 + condition$delta * d2        # slopes differ by delta
  b2    <- fixed_objects$beta2 + condition$delta * d2
  mu    <- fixed_objects$beta0 + b1 * x1 + b2 * x2
  y     <- rnorm(N, mean = mu, sd = fixed_objects$sigma)
  tibble::tibble(x1 = x1, x2 = x2, y = y, group = group)
}
```

=== Model Comparison Technique


// TODO: How does stratified bootstrap work? 

```R
# TODO: How does the stratification work? How do I know it's correct?
boot_mat <- boot::boot(dat, function(d, i) fit_compare(d[i, ]), 
                         R = B, strata = dat$group)$t
```

==== Hwang and Takane (2014) Procedure

// TODO: Two sided or one sided?

```R
# TODO: Do we want the two-sided or one-sided p-value? How do we get that? What does it mean regarding the sampling distribution?
p_ttest_R2    <- t.test(boot_mat[, "R2_split"],    boot_mat[, "R2_pool"],    paired = TRUE)$p.value
```

==== Wald Statistic Procedure 

// TODO: Is this studentized test the two-sided or one-sided p-value?
```R
p_se_R2    <- 2 * pnorm(-abs(T_R2  / sd(dR2)))
p_se_adjR2 <- 2 * pnorm(-abs(T_adj / sd(dAdj)))
```

==== Percentile / ASL bootstrap CI of difference

```R
# (3a) Percentile / ASL bootstrap CI of the difference.
p_pctl_R2    <- min(1, 2 * min(mean(dR2  <= 0), mean(dR2  >= 0)))
p_pctl_adjR2 <- min(1, 2 * min(mean(dAdj <= 0), mean(dAdj >= 0)))
```

==== Label Permutation Test

// TODO: How does permutation bootstrap work? 
// TODO: Get the exact citations for the formulae for the permutation test and double check its correctness

```R
# TODO: How does the permutation work? How do I know it's correct?
perm_stat <- function(d, i) {
    d$group <- d$group[i] # randomly reassign group labels
    fc <- fit_compare(d)
    c(R2  = unname(fc["R2_split"]    - fc["R2_pool"]),
      adj = unname(fc["adjR2_split"] - fc["adjR2_pool"]))
  }
perm <- boot::boot(dat, perm_stat, R = B, sim = "permutation")
p_null_R2    <- (1 + sum(perm$t[, 1] >= perm$t0[1])) / (B + 1)
p_null_adjR2 <- (1 + sum(perm$t[, 2] >= perm$t0[2])) / (B + 1)
```

==== How $R^2$ was computed
For a single model, the R code is
```R
SSE_pool <- sum(.lm.fit(X, y)$residuals^2)
SST <- sum((y - mean(y))^2)
R2_pool = 1 - SSE_pool  / SST
```
This corresponds to the standard mathematical fomulation (James et al. 2021, Page 70, Eq. 3.17):
$
  R^2_("pooled") &= 1 - "SS"_("E")/"SS"_("T")  \
  &= 1 - (sum_(i=1)^(N)(hat(y)_i - y_i)^2)/(sum_(i=1)^N (y_i - overline(y)_i)^2)
$
where $R^2_("pooled") in [0, 1]$.

For two separate models, the R code for $R^2_("split")$ is like so

```R
SSE_split <- sum(.lm.fit(X[g1, , drop = FALSE], y[g1])$residuals^2) + 
               sum(.lm.fit(X[g2, , drop = FALSE], y[g2])$residuals^2)
SST <- sum((y - mean(y))^2)
R2_split    = 1 - SSE_split / SST
```

This `R2_split` approach is appropriate because it's analogous to what is done for GSCA in a multigroup model.#footnote[The aim of this study is to guide a simulation study for the context of GSCA, not to investigate what the optimal procedure for a OLS-only context is.] Mathematically, this is equivalent to

$
  R^2_("split") &= 1 - ("SS"_("E1") + "SS"_("E2"))/"SS"_("T") \
  &= 1 - (sum_(i=1)^(n_1)(hat(y)_i - y_i)^2 + sum_(i=n_1+1)^(N)(hat(y)_i - y_i)^2)/(sum_(i=1)^N (y_i - overline(y)_i)^2)
$

where the data for model one ranges from case _i_ to $n_1$ and the data for model two ranges from case $n_1 + 1$ to $N$. 

We also compute adjusted-$R^2$ to better control for the number of additional parameters estimated when using split versus pooled models (James et al. 2021, Page 234, Eq. 6.4):#footnote[The formula here computes $N-p$ as opposed to the $N-k-1$ shown in James et al. (2021). The implicit context of James et al. (2021) is that the extra $-1$ term is needed to remove the intercept and that $k$---as it is used in James et al. (2021)---does not include the intercept. A proof of this is found in @adjR2Proof.]

$
  "Adjusted" R^2_("pooled") &= 1 - ("SS"_E "/" (N - p))/("SS"_T "/" (N-1)) \
  &= (sum_(i=1)^(N)(hat(y)_i - y_i)^2 "/" (N-p))/(sum_(i=1)^N (y_i - overline(y)_i)^2 "/" (N-1))
$

where _p_ is the number of $beta$-coefficients estimated in a model. In R, this can be done like so:

```R
adjR2_pool = 1 - (SSE_pool / (n - p)) / (SST / (n - 1))
```
using the same `SSE_pool` and `SST` as computed above.

Mathematically, this adjusted-$R^2$ can be extended to the split model, like so:

$
  "Adjusted" R^2_("split") &= 1 - (("SS"_("E1") + "SS"_("E2")) "/" (N - (p_1 + p_2)))/("SS"_T "/" (N-1)) \
  &= 1 - ((sum_(i=1)^(n_1)(hat(y)_i - y_i)^2 + sum_(i=n_1+1)^(N)(hat(y)_i - y_i)^2)"/" (N - 2p))/(sum_(i=1)^N (y_i - overline(y)_i)^2 "/" (N-1))
$ <splitAdjR2>

There are two assertions for proving @splitAdjR2. The first is that for our specific case

$ p_1 +p_2 = 2p $ <2pAssertion>

where $p$ is the number of parameters in the pooled model ($p = 3$). @2pAssertion is proved by the fact that 

```R
y ~ group * (x1 + x2)
```

is equal to

#figure(
  caption: "Expanded Interaction with 2p parameters"
)[```R
y ~ 1  + group + x1 + x2 + group:x1 + group:x2
```] <InteractionSplitCode>
which has _2p_ parameters. The interaction model is relevant because @InteractionSplitCode is equivalent to fitting two separate `lm` models in `R` like so:#footnote[Due to dummy coding, the intercept in @InteractionSplitCode is the baseline of the dummy variable and group represents its contrast.]

#figure(caption: [Two split `lm` models with 3 parameters each])[```R
lm(y ~ 1 + x1 + x2, subset = group == 1)
lm(y ~ 1 + x1 + x2, subset = group == 2)
```]

Therefore, @splitAdjR2 is the correct formula for the adjusted-$R^2_("split")$  for our purposes.

#strong([Importantly, the above proof shows how @splitAdjR2 can  be adapted to compute the adjusted FIT index for a multigroup GSCA type model. It shows how to sum the naive FIT values and it also shows why you can count the number of parameters, as if every parameter was interacting with 'group'. Note that things might be more complex with constrained parameter estimates, but that is outside of the scope of my dissertation.])


=== Outcome

Looking at Type I error rate and statistical power.

// TODO: How does Phil's EDR function work? How does it need to be modified to account for Type I error versus power? One-sided or two-sided?


== Results

Ignore:

+ pctl method
  - Over reject
+ Chow method
  - Seems to make assumptions about the parameters that aren't made in GSCA. So won't pursue further. 

== Discussion

= Study 2: Generalized Structured Component Analysis

== Method

== Results

== Discussion

= General Conclusion


#pagebreak()

#counter(heading).update(0)
#set heading(numbering: "A.1")

= Appendix

== Proofs

=== Adjusted $R^2$ Formula <adjR2Proof>

Here is a LLM-assisted numerical proof of why the $-1$ term is needed for James et al. (2021)'s version of the adjusted $R^2$ formula, but not for ours in @splitAdjR2. 

#figure(caption:[R proof of altered formula.])[```R
fit  <- lm(mpg ~ hp + wt, data = mtcars)
y    <- mtcars$mpg
X    <- model.matrix(fit)          # [1, hp, wt]
n    <- length(y)
p    <- ncol(X)                    
SSE  <- sum(residuals(fit)^2)
SST  <- sum((y - mean(y))^2)

lm_adj   <- summary(fit)$adj.r.squared
manual_np   <- 1 - (SSE/(n - p))     / (SST/(n - 1))   # our formula
manual_np1  <- 1 - (SSE/(n - p - 1)) / (SST/(n - 1))   # the "extra -1" variant

paste0("n = ", n, ", p = ncol(X) = ", p)
paste0("lm() adj.r.squared      : ", round(lm_adj, 3))
paste0("manual  1-(SSE/(n-p))   : ", round(manual_np, 3),  "   identical? ", isTRUE(all.equal(lm_adj, manual_np)))
paste0("manual  1-(SSE/(n-p-1)) : ", round(manual_np1, 3), "   identical? ", isTRUE(all.equal(lm_adj, manual_np1)))
```]
which returns
```
"n = 32, p = ncol(X) = 3"
"lm() adj.r.squared      : 0.815"
"manual  1-(SSE/(n-p))   : 0.815   identical? TRUE"
"manual  1-(SSE/(n-p-1)) : 0.808   identical? FALSE"
```