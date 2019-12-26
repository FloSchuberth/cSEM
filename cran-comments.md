## Release summary

This is the initial release

## Test environments

* local: Windows 10 x64 (build 17763)
* ubuntu xenial 16.04 (on travis-ci), R-devel, R-current (R 3.6.1)
* mac OS 10.13 (on travis-ci) R-current (R 3.6.1)
* windows x86_64-w64-mingw32/x64 (on appveyor) R-current (R 3.6.1)

### Note:

- All builds passed except for the mac os 10.13 build which causes a segfault
  ("address 0x0, cause 'memory not mapped'") when trying to install the package 
  on Travis. Not sure if this is a problem on Travis side. If yes, this issue can
  be ignored.

## R CMD check results (for windows 10 x64 and ubuntu xenial 16.04)

0 ERRORs | 0 WARNINGs | 0 NOTE

## Other issues

- None
