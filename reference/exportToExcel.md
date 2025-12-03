# Export to Excel (.xlsx)

**\[experimental\]**

## Usage

``` r
exportToExcel(
  .postestimation_object = NULL, 
  .filename              = "results.xlsx",
  .path                  = NULL
  )
```

## Arguments

- .postestimation_object:

  An object resulting from a call to one of cSEM's postestimation
  functions (e.g.
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)).

- .filename:

  Character string. The file name.

- .path:

  Character string. Path of the directory to save the file to. Defaults
  to `NULL`.

## Details

Export results from postestimation functions
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
and
[`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
to an .xlsx (Excel) file. The function uses the
[openxlsx](https://rdrr.io/pkg/openxlsx/man/openxlsx.html) package which
does not depend on Java!

The function is deliberately kept simple: it takes all the relevant
elements in `.postestimation_object` and writes them (worksheet by
worksheet) into an .xlsx file named `.filename` in the directory given
by `.path` (the current working directory by default).

If `.postestimation_object` has class attribute `_2ndorder` two .xlsx
files named `".filename_first_stage.xlsx"` and
`".filename_second_stage.xlsx"` are created. If `.postestimation_object`
is a list of appropriate objects, one file for each list elements is
created.

Note: rerunning `exportToExcel()` without changing `.filename` and
`.path` overwrites the file!

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md),
[`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md),
[`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
