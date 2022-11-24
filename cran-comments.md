## Notes

- This is version 0.5.0
- Add a couple of bug fixes, new features, as described in the NEWS file.

## Tests

Tested using GitHub Actions on 

* Windows-latest (x86_64-w64-mingw32 (64-bit)); R-current (R 4.2.2)
* macOS-latest (x86_64-apple-darwin17.0 (64-bit)); R-current (R 4.2.2)
* ubuntu-20.04 (x86_64-pc-linux-gnu (64-bit)); R-current (R 4.2.2)
* ubuntu 20.04 R-devel

There were no errors or warnings.


I also tested the package using devtools::check() with default arguments.
There were no errors or warnings but the following 2 notes:

N  checking dependencies in R code (3.5s) Namespace in Imports field not imported from: 'Rdpack'
     All declared Imports should be used.
     
The Rdpack package is required for referencing.

N  checking package dependencies (3.7s)
   Imports includes 21 non-default packages.
   Importing from so many packages makes the package vulnerable to any of
   them becoming unavailable.  Move as many as possible to Suggests and
   use conditionally.

I checked the dependencies and all listed packages are required.

## Other issues

- None
