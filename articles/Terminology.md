# Terminology

### Common factor

A common factor or latent variable is a type of [construct](#construct).
The name common factor is motivated by the theorized relationship to its
[indicators](#indicator). This relation is commonly referred to as the
[measurement model](#mm). Two kinds of measurement models exist: the
reflective and the (causal-)formative measurement model. The defining
feature of a common factor in a reflective measurement model is the idea
that the common factor is the common underlying (latent) cause of the
realizations of a set of indicators that are said to measure the common
factor. The idea of the reflective measurement model is closely related
to the [true score](#truescore) theory. Accordingly, indicators related
to the common factor are modeled as measurement error-prone
manifestations of a common variable, cause or factor. Although there are
subtle conceptional differences, the terms common factor, true score,
and latent variable are mostly used synonymously in **cSEM**.

The common factor is also the central entity in the causal-formative
measurement model. Here the indicators are modeled as causing the/one
common factor. This type of measurement looks similar to the [composite
(measurement) model](#composite), however, both models are, in fact,
quite different. While the causal-formative measurement model assumes
that the common factor is imperfectly measured by its indicators (i.e.,
there is an error term with non-zero variance) the composite in a
composite model is build/defined by its indicators, i.e., error-free by
definition.

### Composite

A composite is a weighted sum of [indicators](#indicator). Composites
may either serve as [constructs](#construct) in their own right or as
[proxies](#proxy) for a [latent variable](#latentvariable) which, in
turn, serves as a “statistical proxy” for a [concept](#concept) under
study. The nature of the composite is therefore defined by the type of
(measurement) model. If composites are error-free representations of a
concept we refer to the measurement model as the composite (measurement)
model. If composites are used as [stand-ins](#standin) for a [latent
variable](#latentvariable), the measurement model is called
causal-formative (Henseler 2017).

Note that, although we sometimes use it this way as well, the term
“measurement” is actually rather inadequate for the composite model
since in a composite model the construct is build/formed by its related
indicators. Hence no measurement in the actual sense of the word takes
place.

### Composite-based methods

Composite-based methods or composite-based SEM refers to the entirety of
methods centered around the use as of [composites](#composite) (linear
compounds of observables) as [stand-ins](#standin) or error-free
representations for the [concepts](#concept) under investigation.
Composite-based methods are often also called variance-based methods as
focal parameters are usually retrieved such that explained variance of
the dependent constructs in a structural model is maximized.

### Composite-based SEM

See: [Composite-based methods](#cbased)

### Concept

An entity defined by a conceptual/theoretical definition. In line with
Rigdon (2016) variables representing or subsuming a concept are called
conceptual variables. The precise nature/meaning of a conceptual
variable depends on “different assumptions, philosophies or worldviews
\[of the researcher\]” (Rigdon 2016, 2). Unless otherwise stated, in
**cSEM**, it is sufficient to think of concepts as entities that exist
simply because they have been defined. Hence, abstract terms such as
“loyalty” or “depression” as well as designed entities (artifacts) such
as the “OECD Better Life Index” are covered by the definition.

### Construct

Construct refers to a representation of a [concept](#concept) within a
given statistical model. While a concept is defined by a conceptual
(theoretical) definition, a construct for a concept is defined/created
by the researcher’s operationalization of a concept within a statistical
model. Concepts are either modeled as common factors/latent variables or
as composites. Both operationalizations - the common factor and the
composite - are called constructs in **cSEM**. As opposed to concepts,
constructs therefore exist because they arise as the result of the act
of modeling their relation to the observable variables (indicators)
based on a specific set of assumptions. Constructs may therefore best be
understood as [stand-ins](#standin), i.e. statistical proxies for
concepts. Consequently, constructs do not necessarily represent the
concept they seek to represent, i.e., there may be a validity gap.

### Covariance-based SEM

See: [Factor-based methods](#fbmethods)

### Factor-based methods

Factor-based methods or factor-based SEM refers to the entirety of
methods centered around the use of [common factors](#commonfactor) as
statistical [proxies](#proxy) for the [concepts](#concept) under
investigation. Factor-based methods are also called covariance-based
methods as focal parameters are retrieved such that the difference
between the model-implied $\mathbf{\Sigma}(\theta)$ and the empirical
indicator covariance matrix $\mathbf{S}$ is minimized.

### Indicator

An observable variable. In **cSEM** observable variables are generally
referred to as indicators, however, terms such as item, manifest
variable, or observable (variable) are sometimes used synonymously.

### Latent variable

See: [Common factor](#commonfactor)

### Measurement model

The measurement model is a statistical model relating
[indicators](#indicator) to [constructs](#construct) (the statistical
representation of a [concept](#concept)). If the concept under study is
modeled as a [common factor](#commonfactor) two measurement models
exist:

- The reflective measurement model
- The causal-formative measurement model

If the concept under study is modeled as a [composite](#composite) we
call the measurement model:

- composite (measurement) model

Note that, although we sometimes use it this way as well, the term
“measurement” is actually rather inadequate for the composite model
since in a composite model the construct is build/formed by its related
indicators. Hence no measurement in the actual sense of the word takes
place.

### Model

There are many ways to define the term “model”. In **cSEM** we use the
term model to refer to both the theoretical and the statistical model.

Simply speaking theoretical models are a formalized set of hypotheses
stating if and how entities (observable or unobservable) are related.
Since theoretical models are by definition theoretical any statistical
analysis inevitably mandates a statistical model.

A statistical model is typically defined by a set of (testable)
restrictions. Statistical models are best understood as the
operationalized version of the theoretical model. Note that the act of
operationalizing a given theoretical model always entails the
possibility for error. Using a [construct](#construct) modeled as a
[composite](#composite) or a [common factor](#commonfactor), for
instance, is an attempt to map a theoretical entity (the
[concept](#concept)) from the theoretical space into the statistical
space. If this mapping is not one-to-one the construct is only an
imperfect representation of the concept. Consequently, there is a
validity gap.

NOTE: Note that we try to refrain from using the term “model” when
describing an estimation approach or algorithm such as partial least
squares (PLS) as this helps to clearly distinguish between the model and
the approach used to estimate *a given model*.

### Test score

A [proxy](#proxy) for a [true score](#truescore). Usually, the test
score is a simple (unweighted) sum score of
[observables/indicators](#indicators), i.e. unit weights are assumed
when building the test score. More generally, the test score can be any
weighted sum of observables (i.e. a [composite](#composite)), however,
the term “test score” is historically closely tied to the idea that it
is indeed a simple sum score. Hence, whenever it is important to
distinguish between a true score as representing a sum score and a true
score as representing a weighted sum of indicators (where indicator
weights a not necessarily one) we will explicitly state what kind of
test score we mean.

### True score

The term true score derives from the true score theory which
theorizes/models an [indicator](#indicator) or outcome as the sum of a
true score and an error score. The term is closely linked to the [latent
variable/common factor](#commonfactor) model in that the true score of a
set of indicators are linear functions of some underlying common factor.
Mathematically speaking, the correspondence is
$\eta_{jk} = \lambda_{jk}\eta_{j}$  
where $\eta_{jk}$ is the true score, $\eta_{j}$ the underlying latent
variable and $\lambda_{jk}$ the loading. Despite some differences, the
term true score can generally be used synonymously to the terms common
factor and latent variable in **cSEM** without risking a
misunderstanding.

### Proxy

Any quantity that functions as a representation of some other quantity.
Prominent proxies are the [test score](#testscore) - which serves as a
stand-in/proxy for the [true score](#truescore) - and the
[composite](#composite) - which serves as a stand-in/proxy for a common
factor if not used as a composite in its own right. Proxies are
usually - but not necessarily - error-prone representations of the
quantity they seek to represent.

### Saturated and non-saturated models

A structural model is called “saturated” if all [constructs](#construct)
of the model are allowed to freely covary. This is equivalent to saying
that none of the path of the structural model are restricted to zero. A
saturated model has zero degrees of freedom and hence carries no
testable restrictions. If at least one path is restricted to zero, the
structural model is called “non-saturated”.

### Stand-in

See: [Proxy](#proxy)

### Structural Equation Modeling (SEM)

The entirety of a set of related theories, mathematical models, methods,
algorithms and terminologies related to analyzing the relationships
between [concepts](#concept) and/or observables.

### Variance-based methods

See: [Composite-based methods](#cbased)

## Literature

Henseler, Jörg. 2017. “Bridging Design and Behavioral Research with
Variance-Based Structural Equation Modeling.” *Journal of Advertising*
46 (1): 178–92. <https://doi.org/10.1080/00913367.2017.1281780>.

Rigdon, Edward E. 2016. “Choosing PLS Path Modeling as Analytical Method
in European Management Research: A Realist Perspective.” *European
Management Journal* 34 (6). <https://doi.org/10.1016/j.emj.2016.05.006>.
