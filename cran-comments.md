## Test environments

* ubuntu 24.04.3 (GitHub Actions), R 4.5.1
* macOS 15.5 (local + GitHub Actions), R 4.5.1
* Windows 2022 10.0.20348 (GitHub Actions), R 4.5.1

## R CMD check results

* There were 0 ERRORs.
* There were 0 WARNINGs.
* There were 0 NOTEs.

## Downstream dependencies

* The `txshift` and `survML` packages list import this package.

## Additional notes

* This package was recently identified as containing a documentation error (see
  below), which is now resolved.

  > Specifically, please see the NOTEs about Rd file(s) with Rd \link{}
  > targets missing package anchors in the "Rd cross-references" check.

  > CRAN is currently changing its package web pages to providing (static)
  > HTML refmans in addition to PDF refmans, which needs Rd cross-references
  > to Rd \link{} targets not in the package itself nor in the base packages
  > to use package anchors, i.e., use \link[PKG]{FOO} (see section
  > "Cross-references" in "Writing R Extensions"): otherwise these links
  > will not work.
