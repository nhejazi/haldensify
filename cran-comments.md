## Test environments
* ubuntu 20.04 (local + GitHub Actions), R 4.1.1
* macOS 10.15 (local + GitHub Actions), R 4.1.1
* windows 2019 (on GitHub Actions), R 4.1.1

## R CMD check results
* There were 0 ERRORs.
* There were 0 WARNINGs.
* There were 0 NOTEs.

## Downstream dependencies
* There is one downstream dependency on CRAN: `txshift`.

## Additional notes
* This is an update to an existing CRAN package, submitted to correct build
  failures introduced by a very recent update to a dependency `hal9001`, per
  correspondence with CRAN maintainers.
* The time-intensive nature of the unit testing has been reduced so that the
  full suite of tests runs to completion in under the 10min allowance. Note
  that CRAN Windows machines seem to take longer than GitHub Actions instances
  for running both examples and tests longer than expected; moreover, the time
  taken on CRAN Windows instances exceeds that on CRAN Debian by about 3 times.
  U. Ligges noted in e-mail that "We can disable regular tests on Windows if
  the difference is that drastic, just let us know when you submitted with
  shorter examples." (05 October)
* This re-submission reduces the time-intensivity of three examples.
