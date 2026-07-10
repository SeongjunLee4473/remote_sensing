# MODIS MOD13Q1 NDVI/EVI time series & phenology

Builds a multi-year vegetation index time series from NASA's MODIS MOD13Q1.061 (250 m, 16-day composite NDVI/EVI) over the Honam plain rice-paddy region, South Korea, and walks through extracting simple phenological metrics (start/peak/end of growing season) from the seasonal curve.

Unlike the [`hls_l30_s30`](../hls_l30_s30/) tutorial (single-scene, per-pixel spatial comparison, UTM grid), this one is about a **time axis**: one fixed AOI observed roughly every 16 days for 9 years, so the interesting questions are about the shape of the seasonal curve rather than spatial harmonization.

## Data source

NASA LP DAAC, MOD13Q1.061, requested via [AppEEARS](https://appeears.earthdatacloud.nasa.gov/) as an **Area** sample request:

- **Layers**: `_250m_16_days_NDVI`, `_250m_16_days_EVI`, `_250m_16_days_pixel_reliability`, `_250m_16_days_composite_day_of_the_year`
- **Projection**: Geographic (WGS84) — AppEEARS reprojects MODIS's native Sinusoidal grid for you, so this notebook works directly in longitude/latitude with no reprojection step
- **Output format**: NetCDF (one file, every date bundled with a shared `time`/`lat`/`lon` grid — see the root README's note on why NetCDF was chosen over per-date GeoTIFFs for this time-series use case)

Place the downloaded file at:

```
tutorials/modis_vi_timeseries/data/MOD13Q1.061_250m_aid0001.nc
```

or point the notebook at a different location via the `MODIS_NC_PATH` environment variable. Raw data is not committed to the repo (see root `.gitignore`).

## What the notebook covers

1. Load the NetCDF export via `MODIS.VITimeSeriesReader`
2. Decode Pixel Reliability via `MODIS.is_reliable` and mask unreliable (cloudy/snow/fill) composites
3. Compute an AOI-mean NDVI/EVI time series across the full multi-year record
4. Zoom into one representative year and compare the NDVI and EVI seasonal curves
5. Extract simple phenology metrics (start/peak/end of season) via a threshold-crossing method
6. Plot a single-date NDVI snapshot as a map, directly in longitude/latitude (no UTM tick-relabeling trick needed here, since the data is already in a geographic grid)

**Caveat**: the phenology thresholding here is a simple, illustrative amplitude-threshold method, not a validated phenology-detection algorithm — treat the extracted dates as approximate.
