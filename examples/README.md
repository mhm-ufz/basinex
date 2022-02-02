# Example for basinex

This is the river Elbe with one gauge at Decin.
Needed files are `facc.asc` with the flow accumulation and `fdir.asc` with the flow direction.
In addition, we cut out a NetCDF file containing precipitation.

Just execute the basin extractor in this directory:

```bash
basinex
```

## Data sources

### DEM related

The two DEM derivatives (`fdir` and `facc`) are based on upscaled terrain elevation data,
that was collected from USGS EROS Archive—Digital Elevation—Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010),
available at: https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation.

### Precipitation

Sample ERA5 meteo file, available at: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

> Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., Nicolas, J., Peubey, C., Radu, R., Schepers, D. and Simmons, A., 2020.
> The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), pp.1999-2049

### Gauge

Gauge IDs are obtained from the GRDC database: https://www.bafg.de/GRDC/EN/Home/homepage_node.html
