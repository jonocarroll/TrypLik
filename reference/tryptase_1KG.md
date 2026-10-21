# Retrieve the 1000 Genomes data, processed for tryptase allele counts

See the background vignette of `TrypLik` for additional details.

## Usage

``` r
tryptase_1KG()
```

## Value

A `data.frame` of 2A/D, 2B/D, and the inferred genotype group for the
503 relevant samples

## Examples

``` r
head(tryptase_1KG())
#>   2A/D   2B/D Geno Group
#> 1    0 2.9885      0α:4β
#> 2    0 3.2629      0α:4β
#> 3    0 3.4384      0α:4β
#> 4    0 3.4791      0α:4β
#> 5    0 3.4792      0α:4β
#> 6    0 3.5294      0α:4β
```
