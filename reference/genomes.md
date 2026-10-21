# 1000 Genomes Tryptase Data

A subset of 503 samples from the 1000 Genomes Project, processed for
tryptase allele counts. See the background vignette of `TrypLik` for
additional details.

## Usage

``` r
data(genomes)
```

## Format

### `genomes`

A data frame with 503 rows and 3 columns:

- 2A/D:

  double ratio of alpha to delta alleles

- 2B/D:

  double ratio of beta to delta alleles

- Geno Group:

  Genotype classification group

## Source

<https://www.internationalgenome.org/1000-genomes-summary/>

## Examples

``` r
data(genomes)
head(genomes)
#>   2A/D   2B/D Geno Group
#> 1    0 2.9885      0α:4β
#> 2    0 3.2629      0α:4β
#> 3    0 3.4384      0α:4β
#> 4    0 3.4791      0α:4β
#> 5    0 3.4792      0α:4β
#> 6    0 3.5294      0α:4β
```
