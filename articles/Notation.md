# Notation

## The structural model

The structural model specifies the relationships between
[constructs](https://floschuberth.github.io/cSEM/articles/Terminology.html)
(i.e., the statistical representation of a
[concept](https://floschuberth.github.io/cSEM/articles/Terminology.html))
via paths (arrows) and associated path coefficients. The path
coefficients - sometimes also called structural coefficients - express
the magnitude of the influence exerted by the construct at the start of
the arrow on the variable at the arrow’s end. In [composite-based
SEM](https://floschuberth.github.io/cSEM/articles/Terminology.html)
constructs are always operationalized (not modeled!!) as
[composites](https://floschuberth.github.io/cSEM/articles/Terminology.html),
i.e., weighted linear combinations of its respective
[indicators](https://floschuberth.github.io/cSEM/articles/Terminology.html).
Consequently, depending on how a given construct is modeled, such a
composite may either serve as a
[proxy](https://floschuberth.github.io/cSEM/articles/Terminology.html)
for an underlying [latent
variable](https://floschuberth.github.io/cSEM/articles/Terminology.html)
([common
factor](https://floschuberth.github.io/cSEM/articles/Terminology.html))
or as a composite in its own right. Despite this crucial difference, we
stick with the common - although somewhat ambivalent - notation and
represent both the construct and the latent variable (which is only
**a** possible construct) by $`\eta`$. Let
$`x_{kj}`$$`(k = 1,\dots, K_j)`$ be an indicator (observable) belonging
to construct $`\eta_j`$$`(j = 1\dots, J)`$ and $`w_{kj}`$ be a weight. A
composite is definied as:
``` math
\hat{\eta}_j = \sum^{K_j}_{k = 1} w_{kj} x_{kj} 
```
Again, $`\hat{\eta}_j`$ may represent a latent variable $`\eta_j`$ but
may also serve as composite in its own right in which case we would
essentially say that  
$`\hat{\eta}_j = \eta_j`$ and refer to $`\eta_j`$ as a construct instead
of a latent variable. Since $`\hat{\eta}_j`$ generally does not have a
natural scale, weights are usually chosen such that $`\hat{\eta}_j`$ is
standardized. Therefore, unless otherwise stated:

``` math
E(\hat\eta_j) = 0\quad\quad \text{and}\quad\quad Var(\hat\eta_j) = E(\hat\eta^2_j) = 1
```

Since the relations between
[concepts](https://floschuberth.github.io/cSEM/articles/Terminology.html)(or
its statistical sibling the constructs) are a product of the
researcher’s theory and assumptions to be analyzed, some constructs are
typically not directly connected by a path. Technically this implies a
restriction of the path between construct a path we call the structural
model
[saturated](https://floschuberth.github.io/cSEM/articles/Terminology.html).
If at least one path is restricted to zero, the structural model is
called
[non-saturated](https://floschuberth.github.io/cSEM/articles/Terminology.html).

## The reflective measurement model

Define the general reflective (congeneric) measurement model as:
``` math
 x_{kj} = \eta_{kj} + \varepsilon_{kj} = \lambda_{kj}\eta_j +  \varepsilon_{kj}\quad\text{for}\quad k = 1, \dots, K_j\quad\text{and}\quad j = 1, \dots, J
```

Call $`\eta_{kj} = \lambda_{kj}\eta_j`$ the (indicator) true/population
score and $`\eta_j`$ the underlying latent variable supposed to be the
common factor or cause of the $`K_j`$ indicators connected to latent
variable $`\eta_j`$. Call $`\lambda_{kj}`$ the loading or direct effect
of the latent variable on its indicator. Let $`x_{kj}`$ be an indicator
(observable), $`\varepsilon_{kj}`$ be a measurement error and  
``` math
\hat{\eta}_j = \sum^{K_j}_{k = 1} w_{kj} x_{kj} = \sum^{K_j}_{k = 1} w_{kj} \eta_{kj} + \sum^{K_j}_{k = 1} w_{kj} \varepsilon_{kj}
= \bar\eta_{j} + \bar\varepsilon_{j} = \eta_j\sum_{k=1}^{K_J}w_{kj}\lambda_{kj} + \bar\varepsilon_{kj}, 
```
be a proxy/test score/composite/stand-in for/of $`\eta_j`$ based on a
weighted sum of observables, where $`w_{kj}`$ is a weight to be
determined and $`\bar\eta_j`$ the proxy true score, i.e., a weighted sum
of (indicator) true scores. Note the distinction between what we refer
to as the **indicator true score** $`\eta_{kj}`$ and the **proxy true
score** which is the true score for $`\hat\eta_j`$ (i.e, the true score
of a score that is in fact a linear combination of (indicator) scores!).

We will usually refer to $`\hat\eta_j`$ as a proxy for $`\eta_j`$ as it
stresses the fact that $`\hat\eta_j`$ is generally not the same as
$`\eta_j`$ unless $`\bar\varepsilon_{j} = 0`$ and
$`\sum_{k=1}^{K_J}w_{kj}\lambda_{kj} = 1`$.

Assume that
$`E(\varepsilon_{kj}) = E(\eta_j) = Cov(\eta_j, \varepsilon_{kj}) = 0`$.
Further assume that $`Var(\eta_j) = E(\eta^2_j) = 1`$ to determine the
scale.

It often suffices to look at a generic test score/latent variable. For
the sake of clarity the index $`j`$ is therefore dropped unless it is
necessary to avoid confusion.

Note that most of the classical literature on quality criteria such as
reliability is centered around the idea that the proxy $`\hat\eta`$ is a
in fact a simple sum score which implies that all weighs are set to one.
Treatment is more general here since $`\hat{\eta}`$ is allowed to be
*any* weighted sum of related indicators. Readers familiar with the
“classical treatment” may simply set weights to one (unit weights) to
“translate” results to known formulae.

Based on the assumptions and definitions above the following quantities
necessarily follow:

\$\$
``` math
\begin{align}
Cov(x_k, \eta) &= \lambda_k \\
Var(\eta_k) &= \lambda^2_k \\
Var(x_k)    &= \lambda^2_k + Var(\varepsilon_k) \\
Cor(x_k, \eta) &= \rho_{x_k, \eta} = \frac{\lambda_k}{\sqrt{Var(x_k)}} \\

Cov(\eta_k, \eta_l) &= Cor(\eta_k, \eta_l) = E(\eta_k\eta_l) = \lambda_k\lambda_lE(\eta^2) = \lambda_k\lambda_l \\

Cov(x_k, x_l) &= \lambda_k\lambda_lE(\eta^2) + \lambda_kE(\eta\varepsilon_k) + \lambda_lE(\eta\varepsilon_l) + E(\varepsilon_k\varepsilon_l) = \lambda_k\lambda_l + \delta_{kl} \\

Cor(x_k, x_l) &= \frac{\lambda_k\lambda_l + \delta_{kl}}{\sqrt{Var(x_k)Var(x_l)}} \\

Var(\bar\eta) &= E(\bar\eta^2) = \sum w_k^2\lambda^2_k + 2\sum_{k < l} w_k w_l \lambda_k\lambda_l = \left(\sum w_k\lambda_k \right)^2 = (\boldsymbol{\mathbf{w}}'\boldsymbol{\mathbf{\lambda}})^2 \\

Var(\bar\varepsilon) &= E(\bar\varepsilon^2) = \sum w_k^2E(\varepsilon_k^2) + 2\sum_{k < l} w_k w_lE(\varepsilon_k\varepsilon_l)\\

Var(\hat\eta) &= E(\hat\eta^2) = \sum w_k^2(\lambda^2_k + Var(\varepsilon_k)) + 2\sum_{k < l} w_k w_l (\lambda_k\lambda_l + \delta_{kl}) \\
&= \sum w_k^2\lambda^2_k + 2\sum_{k < l} w_k w_l \lambda_k\lambda_l + \sum w_k^2Var(\varepsilon_k)  + 2\sum_{k < l} w_k w_l \delta_{kl} \\
&=Var(\bar\eta) + Var(\bar\varepsilon) = (\boldsymbol{\mathbf{w}}'\boldsymbol{\mathbf{\lambda}})^2 + Var(\bar\varepsilon) = \boldsymbol{\mathbf{w}}'\boldsymbol{\mathbf{\Sigma}}\boldsymbol{\mathbf{w}} \\

Cov(\eta, \hat\eta) &= E\left(\sum w_k \lambda_k \eta^2\right) = \sum w_k\lambda_k = \boldsymbol{\mathbf{w}}'\boldsymbol{\mathbf{\lambda}}= \sqrt{Var(\bar\eta)}
\end{align}
```
\$\$

where $`\delta_{kl} = Cov(\varepsilon_{k}, \varepsilon_{l})`$ for
$`k \neq l`$ is the measurement error covariance and
$`\boldsymbol{\mathbf{\Sigma}}`$ is the indicator variance-covariance
matrix implied by the measurement model:

``` math
\boldsymbol{\mathbf{\Sigma }}= \begin{pmatrix}
\lambda^2_1 + Var(\varepsilon_1) & \lambda_1\lambda_2 + \delta_{12}  & \dots & \lambda_1\lambda_K + \delta_{1K} \\
\lambda_2\lambda_ 1 + \delta_{21} & \lambda^2_2 + Var(\varepsilon_2) & \dots & \lambda_2\lambda_K +\delta_{1K} \\
 \vdots & \vdots & \ddots & \vdots \\
\lambda_{K}\lambda_1 + \delta_{K1} & \lambda_K\lambda_2 + \delta_{K2} &\dots &\lambda^2_K + Var(\varepsilon_K)
\end{pmatrix}
```

In **cSEM** indicators are always standardized and weights are always
appropriately scaled such that the variance of $`\hat\eta`$ is equal to
one. Furthermore, unless explicitly specified measurement error
covariance is restricted to zero. As a consequence, it necessarily
follows that:

``` math
\begin{align}
Var(x_k) &= 1 \\
Cov(x_k, \eta) &= Cor(x_k, \eta) \\
Cov(x_k, x_l) &= Cor(x_k, x_l) \\
Var(\hat\eta) &= \boldsymbol{\mathbf{w}}'\boldsymbol{\mathbf{\Sigma}}\boldsymbol{\mathbf{w}} = 1 \\
Var(\varepsilon_k) &= 1 - Var(\eta_k) = 1 - \lambda^2_k \\
Cov(\varepsilon_k, \varepsilon_l) &= 0 \\
Var(\bar\varepsilon) &= \sum w_k^2 (1 - \lambda_k^2)
\end{align}
```
For most formulae this implies a significant simplification, however,
for ease of comparison to extant literature formulae we stick with the
“general form” here but mention the “simplified form” or “cSEM form” in
the Methods and Formula sections.

## Notation table

| Symbol | Dimension | Description |
|:---|:---|:---|
| $`x_{kj}`$ | $`(1 \times 1)`$ | The $`k`$’th indicator of construct $`j`$ |
| $`\eta_{kj}`$ | $`(1 \times 1)`$ | The $`k`$’th (indicator) true score related to construct $`j`$ |
| $`\eta_j`$ | $`(1 \times 1)`$ | The $`j`$’th common factor/latent variable |
| $`\lambda_{kj}`$ | $`(1 \times 1)`$ | The $`k`$’th (standardized) loading or direct effect of $`\eta_j`$ on $`x_{kj}`$ |
| $`\varepsilon_{kj}`$ | $`(1 \times 1)`$ | The $`k`$’th measurement error or error score |
| $`\hat\eta_j`$ | $`(1 \times 1)`$ | The $`j`$’th test score/composite/proxy for $`\eta_j`$ |
| $`w_{kj}`$ | $`(1 \times 1)`$ | The $`k`$’th weight |
| $`\bar\eta_j`$ | $`(1 \times 1)`$ | The $`j`$’th (proxy) true score, i.e. the weighted sum of (indicator) true scores |
| $`\delta_{kl}`$ | $`(1 \times 1)`$ | The covariance between the $`k`$’th and the $`l`$’th measurement error |
| $`\boldsymbol{\mathbf{w}}`$ | $`(K \times 1)`$ | A vector of weights |
| $`\boldsymbol{\mathbf{\lambda}}`$ | $`(K \times 1)`$ | A vector of loadings |
