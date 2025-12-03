# Internal: Calculate consistent moments of a nonlinear model

Collection of various moment estimators. See
[classifyConstructs](https://floschuberth.github.io/cSEM/reference/classifyConstructs.md)
for a list of possible moments.

## Usage

``` r
SingleSingle(.i, .j, .Q, .H)

SingleQuadratic(.i, .j, .Q, .H)

SingleCubic(.i, .j, .Q, .H)

SingleTwInter(.i, .j, .Q, .H)

SingleThrwInter(.i, .j, .Q, .H)

SingleQuadTwInter(.i, .j, .Q, .H)

QuadraticQuadratic(.i, .j, .Q, .H)

QuadraticCubic(.i, .j, .Q, .H)

QuadraticTwInter(.i, .j, .Q, .H)

QuadraticThrwInter(.i, .j, .Q, .H)

QuadraticQuadTwInter(.i, .j, .Q, .H)

CubicCubic(.i, .j, .Q, .H)

CubicTwInter(.i, .j, .Q, .H)

CubicThrwInter(.i, .j, .Q, .H)

CubicQuadTwInter(.i, .j, .Q, .H)

TwInterTwInter(.i, .j, .Q, .H)

TwInterThrwInter(.i, .j, .Q, .H)

TwInterQuadTwInter(.i, .j, .Q, .H)

ThrwInterThrwInter(.i, .j, .Q, .H)

ThrwInterQuadTwInter(.i, .j, .Q, .H)

QuadTwInercQuadTwInter(.i, .j, .Q, .H)
```

## Arguments

- .i:

  Row index

- .j:

  Column index

- .Q:

  A vector of composite-construct correlations with element names equal
  to the names of the J construct names used in the measurement model.
  Note Q^2 is also called the reliability coefficient.

- .H:

  The (N x J) matrix of construct scores.

## Details

M is the matrix of the sample counterparts (estimates) of the left-hand
side terms in Equation (21) - (24) (Dijkstra and Schermelleh-Engel 2014)
. The label "M" did not appear in the paper and is only used in the
package. Similar is suggested by Wall and Amemiya (2000) using classical
factor scores.

## References

Dijkstra TK, Schermelleh-Engel K (2014). “Consistent Partial Least
Squares For Nonlinear Structural Equation Models.” *Psychometrika*,
**79**(4), 585–604.  
  
Wall MM, Amemiya Y (2000). “Estimation for polynomial structural
equation models.” *Journal of the American Statistical Association*,
**95**(451), 929–940.
