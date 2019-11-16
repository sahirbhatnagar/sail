## Test environments
* local R installation, R 3.6.1
* ubuntu 16.04 (on travis-ci), R 3.6.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

### Win-builder results 
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'doMC'
* checking Rd cross-references ... NOTE
Package unavailable to check Rd xrefs: 'doMC'

I realize doMC is not available for Windows, and therefore have removed all examples, and tests 
which used the doMC package in order to avoid all errors on Win-builder. However, I'm not sure
how to address these notes, which are caused by putting doMC in Suggests. 
