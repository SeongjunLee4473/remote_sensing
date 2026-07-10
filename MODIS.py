"""Reader, QA decoder, and plotter for MODIS MOD13Q1/MYD13Q1 vegetation index time series.

Expects a single NetCDF file as exported by NASA AppEEARS (an "Area" sample request,
Geographic projection) for MOD13Q1.061 or MYD13Q1.061: one file bundling every requested
16-day composite's NDVI, EVI, Pixel Reliability, and Composite Day of the Year on a
shared time/lat/lon grid. Variable names are matched by suffix (``_NDVI``, ``_EVI``, ...)
rather than hardcoded, so the same reader works for the 250m or 500m product, Terra or
Aqua, without changes.
"""

from __future__ import annotations

import datetime

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# Pixel Reliability rank codes (MOD13Q1/MYD13Q1, both collections; see product ATBD).
RELIABILITY_FILL = -1
RELIABILITY_GOOD = 0
RELIABILITY_MARGINAL = 1
RELIABILITY_SNOW_ICE = 2
RELIABILITY_CLOUDY = 3

#: Ranks kept by default: usable data, whether pristine or merely marginal.
DEFAULT_RELIABLE_RANKS = frozenset({RELIABILITY_GOOD, RELIABILITY_MARGINAL})


def is_reliable(pixel_reliability, valid_ranks: frozenset[int] = DEFAULT_RELIABLE_RANKS) -> np.ndarray:
    """Return a boolean array, True where the Pixel Reliability rank is in ``valid_ranks``."""
    return np.isin(pixel_reliability, tuple(valid_ranks))


def _find_var(ds: xr.Dataset, suffix: str) -> str:
    matches = [name for name in ds.data_vars if str(name).endswith(suffix)]
    if not matches:
        raise KeyError(f"No variable ending in {suffix!r} found (have: {list(ds.data_vars)})")
    if len(matches) > 1:
        raise KeyError(f"Multiple variables end in {suffix!r}: {matches}")
    return matches[0]


class VITimeSeriesReader:
    """Reads an AppEEARS MOD13Q1/MYD13Q1 NetCDF export into masked NDVI/EVI DataArrays."""

    def __init__(self, nc_path: str):
        self.dataset = xr.open_dataset(nc_path)
        self._ndvi_var = _find_var(self.dataset, "_NDVI")
        self._evi_var = _find_var(self.dataset, "_EVI")
        self._reliability_var = _find_var(self.dataset, "_pixel_reliability")
        self._doy_var = _find_var(self.dataset, "_composite_day_of_the_year")

    @property
    def time(self) -> np.ndarray:
        return self.dataset.time.values

    @property
    def lat(self) -> np.ndarray:
        return self.dataset.lat.values

    @property
    def lon(self) -> np.ndarray:
        return self.dataset.lon.values

    def reliability(self) -> xr.DataArray:
        return self.dataset[self._reliability_var]

    def composite_day_of_year(self) -> xr.DataArray:
        return self.dataset[self._doy_var]

    def ndvi(
        self, mask_unreliable: bool = True, valid_ranks: frozenset[int] = DEFAULT_RELIABLE_RANKS
    ) -> xr.DataArray:
        return self._masked(self._ndvi_var, mask_unreliable, valid_ranks)

    def evi(
        self, mask_unreliable: bool = True, valid_ranks: frozenset[int] = DEFAULT_RELIABLE_RANKS
    ) -> xr.DataArray:
        return self._masked(self._evi_var, mask_unreliable, valid_ranks)

    def _masked(self, var_name: str, mask_unreliable: bool, valid_ranks: frozenset[int]) -> xr.DataArray:
        data = self.dataset[var_name]
        if not mask_unreliable:
            return data
        valid = is_reliable(self.dataset[self._reliability_var].values, valid_ranks)
        return data.where(valid)

    def roi_mean(
        self, data: xr.DataArray, bbox: tuple[float, float, float, float] | None = None
    ) -> xr.DataArray:
        """Spatial mean of ``data`` over ``bbox`` (west, south, east, north); whole grid if ``None``."""
        if bbox is not None:
            west, south, east, north = bbox
            data = data.sel(lon=slice(west, east), lat=slice(north, south))  # lat is stored descending
        return data.mean(dim=("lat", "lon"), skipna=True)


class VIPlotter:
    """Matplotlib helpers for ``VITimeSeriesReader`` output: map snapshots and time series."""

    def __init__(self, reader: VITimeSeriesReader):
        self.reader = reader

    def plot_map(
        self,
        data_2d,
        title: str | None = None,
        cmap: str = "RdYlGn",
        vmin: float = -1,
        vmax: float = 1,
        colorbar_label: str | None = None,
    ):
        """Plot one lat/lon snapshot (e.g. ``reader.ndvi().sel(time=...)``). Already in WGS84,
        so unlike UTM-gridded products (see HLS.py) no reprojection or tick relabeling is needed."""
        extent = [self.reader.lon.min(), self.reader.lon.max(), self.reader.lat.min(), self.reader.lat.max()]
        fig, ax = plt.subplots(figsize=(7, 7))
        im = ax.imshow(data_2d, extent=extent, origin="upper", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, linestyle="--", color="gray", alpha=0.4)
        if title:
            ax.set_title(title)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=colorbar_label)
        return fig, ax

    def plot_timeseries(self, series: xr.DataArray, label: str | None = None, ax=None, **plot_kwargs):
        """Plot a 1D (time,) series, e.g. ``reader.roi_mean(reader.ndvi(), bbox=...)``.

        AppEEARS NetCDF exports encode ``time`` with cftime objects (non-Gregorian
        calendar attribute), which matplotlib cannot plot directly without the optional
        ``nc-time-axis`` package -- converted here to plain ``datetime.datetime`` by
        year/month/day instead, since that's all a 16-day-composite date needs.
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 4))
        dates = [datetime.datetime(t.year, t.month, t.day) for t in series.time.values]
        ax.plot(dates, series.values, label=label, **plot_kwargs)
        ax.set_xlabel("Date")
        ax.set_ylabel("Index value")
        if label:
            ax.legend()
        return ax
