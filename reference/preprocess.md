# Run preprocessing script

Run preprocessing script

## Usage

``` r
preprocess(
  cram_file,
  output_dir = dirname(cram_file),
  consensus_dir = system.file(package = "TrypLik"),
  nthreads = 8
)
```

## Arguments

- cram_file:

  file path of .cram file to be processed

## Value

produces a temp.sam file

## Examples

``` r
tryCatch(preprocess("missing_file.cram"), error = function(e) e)
#> <simpleError in preprocess("missing_file.cram"): file.exists(cram_file) is not TRUE>

if (FALSE) { # \dontrun{
preprocess("HG00100.final.cram")
} # }
```
