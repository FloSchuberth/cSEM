## Submission 
## Notes

- This is version 0.7.0

## Tests

Tested using GitHub Actions on 

*macos-latest (release)
*windows-latest (release)
*windows-latest (4.1)
*ubuntu-latest (devel)
*ubuntu-latest  (release)

There were no errors or warnings.


I also tested the package using devtools::check() with default arguments.
There were no errors or warnings but the following note:

❯ checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'Rdpack'
    All declared Imports should be used.
     
The Rdpack package is required for referencing.
