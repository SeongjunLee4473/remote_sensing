# Tutorials

Worked examples showing raw-to-analysis-ready preprocessing for individual satellite products, as a companion to the reusable modules in [`GIS/`](../GIS) and [`HLS.py`](../HLS.py).

## Environment

Same environment as the rest of the repo (see the root [README](../README.md#installation)) — `rasterio`, `GDAL`, `cartopy`, `geopandas`.

## Contents

| Notebook | Focus | Data source |
|---|---|---|
| [`hls_l30_s30/preprocessing_hls_l30_s30.ipynb`](hls_l30_s30/preprocessing_hls_l30_s30.ipynb) | Cross-sensor comparison of NASA HLS L30 (Landsat) and S30 (Sentinel-2): valid-pixel ROI intersection, RGB composites, cross-sensor NDVI check, and each product's unique bands (L30 thermal, S30 red edge) | NASA LP DAAC HLS v2.0 |
| [`modis_vi_timeseries/preprocessing_modis_vi_timeseries.ipynb`](modis_vi_timeseries/preprocessing_modis_vi_timeseries.ipynb) | Multi-year MODIS NDVI/EVI time series over rice-paddy cropland: Pixel Reliability masking, NDVI-vs-EVI seasonal comparison, start/peak/end-of-season phenology extraction | NASA LP DAAC MOD13Q1.061 (via AppEEARS) |

Each tutorial folder has its own `README.md` with the exact data source and how to point the notebook at your own download (raw data is not committed to this repo).
