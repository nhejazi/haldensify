# Print: Highly Adaptive Lasso Conditional Density Estimates

Print: Highly Adaptive Lasso Conditional Density Estimates

## Usage

``` r
# S3 method for class 'haldensify'
print(x, ...)
```

## Arguments

- x:

  An object of class `haldensify`.

- ...:

  Other options (not currently used).

## Value

None. Called for the side effect of printing an informative summary of
slots of objects of class `haldensify`.

## Details

The `print` method for objects of class `haldensify`

## Examples

``` r
# simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
set.seed(11249)
n_train <- 50
w <- runif(n_train, -4, 4)
a <- rnorm(n_train, w, 0.5)

# learn relationship A|W using HAL-based density estimation procedure
haldensify_fit <- haldensify(
  A = a, W = w, n_bins = c(3, 5),
  lambda_seq = exp(seq(-1, -15, length = 50L)),
  max_degree = 2, smoothness_orders = 0, reduce_basis = 0.1
)
print(haldensify_fit)
#> HAL Conditional Density Estimation
#> Number of bins over support of A: 5
#> CV-selected lambda: 0.0038
#> Summary of fitted HAL:
#>          coef               term
#>         <num>             <char>
#>  1:  13.66335        (Intercept)
#>  2: -18.73701 [ I(W >= -2.086) ]
#>  3: -15.56878 [ I(W >= -1.105) ]
#>  4:  15.39282 [ I(W >= -1.612) ]
#>  5:  15.17157 [ I(bin_id >= 2) ]
#>  6:  13.70079 [ I(bin_id >= 3) ]
#>  7: -13.35268  [ I(W >= 1.006) ]
#>  8:  12.79074 [ I(bin_id >= 4) ]
#>  9: -11.35484  [ I(W >= 1.884) ]
#> 10:  11.30173  [ I(W >= 1.487) ]
```
