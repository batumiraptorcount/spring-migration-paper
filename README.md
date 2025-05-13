# spring-migration-paper
Code for the spring paper Tal et al.

## Additional data
### Digital Elevation Model
[30-m SRTM Digital Elevation Model](https://lpdaac.usgs.gov/products/srtmgl1v003/) data can easily be downloaded via [this tool](https://dwtkns.com/srtm30m/). The file used for the generation of least cost paths in Georgia can be found in `data/dem/Georgia_DEM_1200x800.tif`.

### Least cost paths
Procedure to generate least cost paths can be found in `leastcostpaths.R`. Paths paths are stored in `data/lcps/`. Only the shortest routes, files suffixed by `_shortest`, are presented in the paper.