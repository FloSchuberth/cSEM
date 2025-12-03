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
**a** possible construct) by $\eta$. Let
$x_{kj}$$\left( k = 1,\ldots,K_{j} \right)$ be an indicator (observable)
belonging to construct $\eta_{j}$$(j = 1\ldots,J)$ and $w_{kj}$ be a
weight. A composite is definied as:
$${\widehat{\eta}}_{j} = \sum\limits_{k = 1}^{K_{j}}w_{kj}x_{kj}$$
Again, ${\widehat{\eta}}_{j}$ may represent a latent variable $\eta_{j}$
but may also serve as composite in its own right in which case we would
essentially say that  
${\widehat{\eta}}_{j} = \eta_{j}$ and refer to $\eta_{j}$ as a construct
instead of a latent variable. Since ${\widehat{\eta}}_{j}$ generally
does not have a natural scale, weights are usually chosen such that
${\widehat{\eta}}_{j}$ is standardized. Therefore, unless otherwise
stated:

$$E\left( {\widehat{\eta}}_{j} \right) = 0\quad\quad\text{and}\quad\quad Var\left( {\widehat{\eta}}_{j} \right) = E\left( {\widehat{\eta}}_{j}^{2} \right) = 1$$

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
$$x_{kj} = \eta_{kj} + \varepsilon_{kj} = \lambda_{kj}\eta_{j} + \varepsilon_{kj}\quad\text{for}\quad k = 1,\ldots,K_{j}\quad\text{and}\quad j = 1,\ldots,J$$

Call $\eta_{kj} = \lambda_{kj}\eta_{j}$ the (indicator) true/population
score and $\eta_{j}$ the underlying latent variable supposed to be the
common factor or cause of the $K_{j}$ indicators connected to latent
variable $\eta_{j}$. Call $\lambda_{kj}$ the loading or direct effect of
the latent variable on its indicator. Let $x_{kj}$ be an indicator
(observable), $\varepsilon_{kj}$ be a measurement error and  
$${\widehat{\eta}}_{j} = \sum\limits_{k = 1}^{K_{j}}w_{kj}x_{kj} = \sum\limits_{k = 1}^{K_{j}}w_{kj}\eta_{kj} + \sum\limits_{k = 1}^{K_{j}}w_{kj}\varepsilon_{kj} = {\bar{\eta}}_{j} + {\bar{\varepsilon}}_{j} = \eta_{j}\sum\limits_{k = 1}^{K_{J}}w_{kj}\lambda_{kj} + {\bar{\varepsilon}}_{kj},$$
be a proxy/test score/composite/stand-in for/of $\eta_{j}$ based on a
weighted sum of observables, where $w_{kj}$ is a weight to be determined
and ${\bar{\eta}}_{j}$ the proxy true score, i.e., a weighted sum of
(indicator) true scores. Note the distinction between what we refer to
as the **indicator true score** $\eta_{kj}$ and the **proxy true score**
which is the true score for ${\widehat{\eta}}_{j}$ (i.e, the true score
of a score that is in fact a linear combination of (indicator) scores!).

We will usually refer to ${\widehat{\eta}}_{j}$ as a proxy for
$\eta_{j}$ as it stresses the fact that ${\widehat{\eta}}_{j}$ is
generally not the same as $\eta_{j}$ unless
${\bar{\varepsilon}}_{j} = 0$ and
$\sum_{k = 1}^{K_{J}}w_{kj}\lambda_{kj} = 1$.

Assume that
$E\left( \varepsilon_{kj} \right) = E\left( \eta_{j} \right) = Cov\left( \eta_{j},\varepsilon_{kj} \right) = 0$.
Further assume that
$Var\left( \eta_{j} \right) = E\left( \eta_{j}^{2} \right) = 1$ to
determine the scale.

It often suffices to look at a generic test score/latent variable. For
the sake of clarity the index $j$ is therefore dropped unless it is
necessary to avoid confusion.

Note that most of the classical literature on quality criteria such as
reliability is centered around the idea that the proxy $\widehat{\eta}$
is a in fact a simple sum score which implies that all weighs are set to
one. Treatment is more general here since $\widehat{\eta}$ is allowed to
be *any* weighted sum of related indicators. Readers familiar with the
“classical treatment” may simply set weights to one (unit weights) to
“translate” results to known formulae.

Based on the assumptions and definitions above the following quantities
necessarily follow:

\$\$ $$\begin{aligned}
{Cov\left( x_{k},\eta \right)} & {= \lambda_{k}} \\
{Var\left( \eta_{k} \right)} & {= \lambda_{k}^{2}} \\
{Var\left( x_{k} \right)} & {= \lambda_{k}^{2} + Var\left( \varepsilon_{k} \right)} \\
{Cor\left( x_{k},\eta \right)} & {= \rho_{x_{k},\eta} = \frac{\lambda_{k}}{\sqrt{Var\left( x_{k} \right)}}} \\
{Cov\left( \eta_{k},\eta_{l} \right)} & {= Cor\left( \eta_{k},\eta_{l} \right) = E\left( \eta_{k}\eta_{l} \right) = \lambda_{k}\lambda_{l}E\left( \eta^{2} \right) = \lambda_{k}\lambda_{l}} \\
{Cov\left( x_{k},x_{l} \right)} & {= \lambda_{k}\lambda_{l}E\left( \eta^{2} \right) + \lambda_{k}E\left( \eta\varepsilon_{k} \right) + \lambda_{l}E\left( \eta\varepsilon_{l} \right) + E\left( \varepsilon_{k}\varepsilon_{l} \right) = \lambda_{k}\lambda_{l} + \delta_{kl}} \\
{Cor\left( x_{k},x_{l} \right)} & {= \frac{\lambda_{k}\lambda_{l} + \delta_{kl}}{\sqrt{Var\left( x_{k} \right)Var\left( x_{l} \right)}}} \\
{Var\left( \bar{\eta} \right)} & {= E\left( {\bar{\eta}}^{2} \right) = \sum w_{k}^{2}\lambda_{k}^{2} + 2\sum\limits_{k < l}w_{k}w_{l}\lambda_{k}\lambda_{l} = \left( \sum w_{k}\lambda_{k} \right)^{2} = (\mathbf{w}\prime{\mathbf{λ}})^{2}} \\
{Var\left( \bar{\varepsilon} \right)} & {= E\left( {\bar{\varepsilon}}^{2} \right) = \sum w_{k}^{2}E\left( \varepsilon_{k}^{2} \right) + 2\sum\limits_{k < l}w_{k}w_{l}E\left( \varepsilon_{k}\varepsilon_{l} \right)} \\
{Var\left( \widehat{\eta} \right)} & {= E\left( {\widehat{\eta}}^{2} \right) = \sum w_{k}^{2}\left( \lambda_{k}^{2} + Var\left( \varepsilon_{k} \right) \right) + 2\sum\limits_{k < l}w_{k}w_{l}\left( \lambda_{k}\lambda_{l} + \delta_{kl} \right)} \\
 & {= \sum w_{k}^{2}\lambda_{k}^{2} + 2\sum\limits_{k < l}w_{k}w_{l}\lambda_{k}\lambda_{l} + \sum w_{k}^{2}Var\left( \varepsilon_{k} \right) + 2\sum\limits_{k < l}w_{k}w_{l}\delta_{kl}} \\
 & {= Var\left( \bar{\eta} \right) + Var\left( \bar{\varepsilon} \right) = (\mathbf{w}\prime{\mathbf{λ}})^{2} + Var\left( \bar{\varepsilon} \right) = \mathbf{w}\prime\mathbf{\Sigma}\mathbf{w}} \\
{Cov\left( \eta,\widehat{\eta} \right)} & {= E\left( \sum w_{k}\lambda_{k}\eta^{2} \right) = \sum w_{k}\lambda_{k} = \mathbf{w}\prime{\mathbf{λ}} = \sqrt{Var\left( \bar{\eta} \right)}}
\end{aligned}$$ \$\$

where $\delta_{kl} = Cov\left( \varepsilon_{k},\varepsilon_{l} \right)$
for $k \neq l$ is the measurement error covariance and $\mathbf{\Sigma}$
is the indicator variance-covariance matrix implied by the measurement
model:

$$\mathbf{\Sigma} = \begin{pmatrix}
{\lambda_{1}^{2} + Var\left( \varepsilon_{1} \right)} & {\lambda_{1}\lambda_{2} + \delta_{12}} & \ldots & {\lambda_{1}\lambda_{K} + \delta_{1K}} \\
{\lambda_{2}\lambda_{1} + \delta_{21}} & {\lambda_{2}^{2} + Var\left( \varepsilon_{2} \right)} & \ldots & {\lambda_{2}\lambda_{K} + \delta_{1K}} \\
\vdots & \vdots & \ddots & \vdots \\
{\lambda_{K}\lambda_{1} + \delta_{K1}} & {\lambda_{K}\lambda_{2} + \delta_{K2}} & \ldots & {\lambda_{K}^{2} + Var\left( \varepsilon_{K} \right)}
\end{pmatrix}$$

In **cSEM** indicators are always standardized and weights are always
appropriately scaled such that the variance of $\widehat{\eta}$ is equal
to one. Furthermore, unless explicitly specified measurement error
covariance is restricted to zero. As a consequence, it necessarily
follows that:

$$\begin{aligned}
{Var\left( x_{k} \right)} & {= 1} \\
{Cov\left( x_{k},\eta \right)} & {= Cor\left( x_{k},\eta \right)} \\
{Cov\left( x_{k},x_{l} \right)} & {= Cor\left( x_{k},x_{l} \right)} \\
{Var\left( \widehat{\eta} \right)} & {= \mathbf{w}\prime\mathbf{\Sigma}\mathbf{w} = 1} \\
{Var\left( \varepsilon_{k} \right)} & {= 1 - Var\left( \eta_{k} \right) = 1 - \lambda_{k}^{2}} \\
{Cov\left( \varepsilon_{k},\varepsilon_{l} \right)} & {= 0} \\
{Var\left( \bar{\varepsilon} \right)} & {= \sum w_{k}^{2}\left( 1 - \lambda_{k}^{2} \right)}
\end{aligned}$$ For most formulae this implies a significant
simplification, however, for ease of comparison to extant literature
formulae we stick with the “general form” here but mention the
“simplified form” or “cSEM form” in the Methods and Formula sections.

## Notation table

| Symbol                 | Dimension      | Description                                                                     |
|:-----------------------|:---------------|:--------------------------------------------------------------------------------|
| $x_{kj}$               | $(1 \times 1)$ | The $k$’th indicator of construct $j$                                           |
| $\eta_{kj}$            | $(1 \times 1)$ | The $k$’th (indicator) true score related to construct $j$                      |
| $\eta_{j}$             | $(1 \times 1)$ | The $j$’th common factor/latent variable                                        |
| $\lambda_{kj}$         | $(1 \times 1)$ | The $k$’th (standardized) loading or direct effect of $\eta_{j}$ on $x_{kj}$    |
| $\varepsilon_{kj}$     | $(1 \times 1)$ | The $k$’th measurement error or error score                                     |
| ${\widehat{\eta}}_{j}$ | $(1 \times 1)$ | The $j$’th test score/composite/proxy for $\eta_{j}$                            |
| $w_{kj}$               | $(1 \times 1)$ | The $k$’th weight                                                               |
| ${\bar{\eta}}_{j}$     | $(1 \times 1)$ | The $j$’th (proxy) true score, i.e. the weighted sum of (indicator) true scores |
| $\delta_{kl}$          | $(1 \times 1)$ | The covariance between the $k$’th and the $l$’th measurement error              |
| $\mathbf{w}$           | $(K \times 1)$ | A vector of weights                                                             |
| $\mathbf{λ}$           | $(K \times 1)$ | A vector of loadings                                                            |
