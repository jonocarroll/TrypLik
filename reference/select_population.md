# Select a genetic-ancestry population

Select a genetic-ancestry population

## Usage

``` r
select_population(pop = NULL)
```

## Arguments

- pop:

  population to use (0=EUR, 1=AFR, 2=EAS, 3=SAS) as either a number or
  character string

## Value

a named integer mapping the population name to the index

## Examples

``` r
# by index
select_population(pop = 0)
#> EUR 
#>   0 

# by name
select_population(pop = "EUR")
#> EUR 
#>   0 
```
