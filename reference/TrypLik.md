# Calculate Tryptase Variant Likelihoods

Calculate Tryptase Variant Likelihoods

## Usage

``` r
TrypLik(wt, fs, areads, breads, dreads, pop = 0)
```

## Arguments

- wt:

  wildtype

- fs:

  frameshift

- areads:

  alpha counts

- breads:

  beta counts

- dreads:

  delta counts

- pop:

  population to use (0=EUR, 1=AFR, 2=EAS, 3=SAS) as either a number or
  character string

## Value

a `data.frame` containing probabilities of the relevant likelihoods

## Examples

``` r
# as individual arguments
TrypLik(26, 0, 347, 316, 304, pop = 0)
#>   Alpha_count Beta_count Beta_FS_count Posterior_likelihood
#> 1           2          2             0              0.99890
#> 2           2          3             0              0.00008
#> 3           3          2             0              0.00082
#> 4           3          3             0              0.00020

# as a `counts` vector
TrypLik(c(26, 0, 347, 316, 304), pop = 0)
#>   Alpha_count Beta_count Beta_FS_count Posterior_likelihood
#> 1           2          2             0              0.99890
#> 2           2          3             0              0.00008
#> 3           3          2             0              0.00082
#> 4           3          3             0              0.00020

# as the output of `preprocess`
TrypLik(c(24, 0, 22, 0, 26, 0, 347, 316, 304), pop = 0)
#>   Alpha_count Beta_count Beta_FS_count Posterior_likelihood
#> 1           2          2             0              0.99890
#> 2           2          3             0              0.00008
#> 3           3          2             0              0.00082
#> 4           3          3             0              0.00020
```
