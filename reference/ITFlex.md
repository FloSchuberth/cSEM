# Data: ITFlex

A data frame containing 16 variables with 100 observations.

## Usage

``` r
ITFlex
```

## Format

A data frame containing the following variables:

- `ITCOMP1`:

  Software applications can be easily transported and used across
  multiple platforms.

- `ITCOMP2`:

  Our firm provides multiple interfaces or entry points (e.g., web
  access) for external end users.

- `ITCOMP3`:

  Our firm establishes corporate rules and standards for hardware and
  operating systems to ensure platform compatibility.

- `ITCOMP4`:

  Data captured in one part of our organization are immediately
  available to everyone in the firm.

- `ITCONN1`:

  Our organization has electronic links and connections throughout the
  entire firm.

- `ITCONN2`:

  Our firm is linked to business partners through electronic channels
  (e.g., websites, e-mail, wireless devices, electronic data
  interchange).

- `ITCONN3`:

  All remote, branch, and mobile offices are connected to the central
  office.

- `ITCONN4`:

  There are very few identifiable communications bottlenecks within our
  firm.

- `MOD1`:

  Our firm possesses a great speed in developing new business
  applications or modifying existing applications.

- `MOD2`:

  Our corporate database is able to communicate in several different
  protocols.

- `MOD3`:

  Reusable software modules are widely used in new systems development.

- `MOD4`:

  IT personnel use object-oriented and prepackaged modular tools to
  create software applications.

- `ITPSF1`:

  Our IT personnel have the ability to work effectively in
  cross-functional teams.

- `ITPSF2`:

  Our IT personnel are able to interpret business problems and develop
  appropriate technical solutions.

- `ITPSF3`:

  Our IT personnel are self-directed and proactive.

- `ITPSF4`:

  Our IT personnel are knowledgeable about the key success factors in
  our firm.

## Source

The data was collected through a survey by Benitez et al. (2018) .

## Details

The dataset was studied by Benitez et al. (2018) and is used in Henseler
(2021) for demonstration purposes, see the corresponding tutorial. All
questionnaire items are measured on a 5-point scale.

## References

Benitez J, Ray G, Henseler J (2018). “Impact of Information Technology
Infrastructure Flexibility on Mergers and Acquisitions.” *MIS
Quarterly*, **42**(1), 25–43.  
  
Henseler J (2021). *Composite-Based Structural Equation Modeling:
Analyzing Latent and Emergent Variables*. Guilford Press, New York.

## Examples

``` r
#============================================================================
# Example is taken from Henseler (2020)
#============================================================================
model_IT_Fex="
# Composite models
ITComp  <~ ITCOMP1 + ITCOMP2 + ITCOMP3 + ITCOMP4
Modul   <~ MOD1 + MOD2 + MOD3 + MOD4
ITConn  <~ ITCONN1 + ITCONN2 + ITCONN3 + ITCONN4
ITPers  <~ ITPSF1 + ITPSF2 + ITPSF3 + ITPSF4

# Saturated structural model
ITPers ~ ITComp + Modul + ITConn
Modul  ~ ITComp + ITConn 
ITConn ~ ITComp 
"

out <- csem(.data = ITFlex, .model = model_IT_Fex,
           .PLS_weight_scheme_inner = 'factorial',
           .tolerance = 1e-06,
           .PLS_ignore_structural_model = TRUE)
```
