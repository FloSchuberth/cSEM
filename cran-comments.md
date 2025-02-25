## Resubmission 
Package CITATION file contains call(s) to old-style personList() or
as.personList().  Please use c() on person objects instead.
Package CITATION file contains call(s) to old-style citEntry().  Please
use bibentry() instead.

Fix issues in the citation file. 


Found the following (possibly) invalid URLs:
  URL: https://github.com/HenrikBengtsson/future.apply/ (moved to https://github.com/futureverse/future.apply)
    From: man/resampleData.Rd
          man/resamplecSEMResults.Rd
    Status: 301
    Message: Moved Permanently
  URL: https://github.com/HenrikBengtsson/future/ (moved to https://github.com/futureverse/future)
    From: man/resamplecSEMResults.Rd
          inst/doc/cSEM.html
          README.md
          
Replaced the two links.


checkRd: (-1) calculateEffects.Rd:24: Lost braces; missing escapes or markup?
    24 | equals (I-B)^{(-1)}Gamma. The indirect effect equals the difference between
       |  
       
Removed brackets around -1.

## Notes

- This is version 0.6.0
- Add a couple of bug fixes, new features, as described in the NEWS file.

## Tests

Tested using GitHub Actions on 

*macos-latest (r: 'release')
*windows-latest (r: 'release')
*windows-latest r: '4.1'}
*ubuntu-latest (r: 'devel')
*ubuntu-latest  (r: 'release)

There were no errors or warnings.


I also tested the package using devtools::check() with default arguments.
There were no errors or warnings but the following 4 notes:

checking package dependencies ... NOTE
  Packages suggested but not available for checking:
    'nnls', 'plotly', 'listviewer'
  
  Imports includes 21 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.


I checked the dependencies and all listed packages are required.


 checking for future file timestamps ... NOTE
  unable to verify current time
   
This is a known issue with devtools::check() and does not affect the package.


 checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'Rdpack'
    All declared Imports should be used.
     
The Rdpack package is required for referencing.


 checking Rd files ... NOTE
  checkRd: (-1) calculateEffects.Rd:24: Lost braces; missing escapes or markup?
      24 | equals (I-B)^{(-1)}Gamma. The indirect effect equals the difference between
         |              ^
         
I checked the brackets and they are okay. 