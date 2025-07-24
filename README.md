# spring-migration-paper
[![DOI](https://zenodo.org/badge/983015429.svg)](https://doi.org/10.5281/zenodo.16411086)

Code and data for the spring paper Tal et al.

## Additional data
### Filtered count data
The filtered count data used in the analysis scripts below can be found in the `data/Count data` directory. For both the spring and autumn seasons, `Annual-totals-per-distance-zone_FILTERED.csv` contains the annual totals, while `Daily-totals-per-distance-zone_FILTERED.csv` provides the daily totals. Before using these files, please consult the `Meta data.md` for explanations of each column.

### Count data analysis
The script used to process data, create tables, and compare the significance of autumn and spring count data can be found in `1. BRC - Processing count data.R`. The code for visualizing the output from this processing step is in `2. BRC - Visualizing count data.R`. An overview of the additional non-raptor species counted (as reported in the supplements of Tal et al.) is provided in `3. BRC - Summarizing non-raptor count data.R`.

### Digital Elevation Model
[30-m SRTM Digital Elevation Model](https://lpdaac.usgs.gov/products/srtmgl1v003/) data can easily be downloaded via [this tool](https://dwtkns.com/srtm30m/). The file used for the generation of least cost paths in Georgia can be found in `data/dem/Georgia_DEM_1200x800.tif`.

### Least cost paths
Procedure to generate least cost paths can be found in `leastcostpaths.R`. Paths are stored in `data/lcps/`. Only the shortest routes, files suffixed by `_shortest`, are presented in the paper.
