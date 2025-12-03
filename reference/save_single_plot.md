# Internal: save_single_plot Helper function to save a single DiagrammeR plot based on the file extension

`diagrammer_obj` DiagrammeR plot object to be saved. `out_file` The name
of the file to save the plot to (supports 'pdf', 'png', 'svg', and 'dot'
formats).

## Usage

``` r
save_single_plot(
diagrammer_obj, 
out_file,
.path = NULL)
```

## Arguments

- .path:

  Character string. Path of the directory to save the file to. Defaults
  to `NULL`.

## Value

NULL.
