## Submission 
## Notes

- This is version 0.6.1
- Added NeedsCompilation: no to the description file as asked by the CRAN team.

## Tests

Tested using GitHub Actions on 

*macos-latest (r: 'release')
*windows-latest (r: 'release')
*windows-latest r: '4.1'}
*ubuntu-latest (r: 'devel')
*ubuntu-latest  (r: 'release)

There were no errors or warnings.


I also tested the package using devtools::check() with default arguments.
There were no errors or warnings but the following 2 notes:

❯ checking for future file timestamps ... NOTE
  unable to verify current time
   
This is a known issue with devtools::check() and does not affect the package.


❯ checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'Rdpack'
    All declared Imports should be used
     
The Rdpack package is required for referencing.
