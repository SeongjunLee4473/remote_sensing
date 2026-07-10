# HLS L30 vs S30 preprocessing & comparison

Compares NASA's Harmonized Landsat Sentinel-2 (HLS) v2.0 products for the same MGRS tile: **L30** (derived from Landsat 8/9 OLI) and **S30** (derived from Sentinel-2 MSI). Both are distributed on the same 30 m UTM grid, but a scene pair still differs in acquisition date, cloud cover, and swath coverage — this notebook walks through reconciling that before any pixel-level comparison.

## Data source

NASA LP DAAC, HLS v2.0, tile `T52SCG`:
- `HLS.L30.T52SCG.<DOY>.v2.0.*.tif` — search/download via [NASA Earthdata Search](https://search.earthdata.nasa.gov/) (product: `HLSL30`)
- `HLS.S30.T52SCG.<DOY>.v2.0.*.tif` — same portal (product: `HLSS30`)

Download all per-band `.tif` files (B01…, Fmask, angle bands) for one L30 granule and one S30 granule into:

```
tutorials/hls_l30_s30/data/hls_l30/
tutorials/hls_l30_s30/data/hls_s30/
```

or point the notebook at a different location via the `HLS_L30_DIR` / `HLS_S30_DIR` environment variables. Raw bands are not committed to the repo (see root `.gitignore`).

## What the notebook covers

1. Load bands + tags (`scale_factor`, `unit`) for both products
2. Cloud/shadow masking via `HLS.is_cloud`
3. Intersect valid-pixel footprints (nodata ∪ cloud/shadow) and crop both products to the shared ROI
4. RGB true-color composites, with maps overlaid with a longitude/latitude grid
5. Cross-sensor NDVI comparison (shared red/NIR bands)
6. Sensor-unique bands: L30 thermal (B10 → surface brightness temperature) and S30 red edge (B05 + B08 → NDRE)

**Caveat:** the two granules here are ~11 days apart with different sun/view geometry, so this demonstrates grid/format harmonization — not a same-moment cross-calibration.
