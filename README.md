# pmwDPR

A Python package for collocating passive microwave (PMW) brightness temperatures with Global Precipitation Measurement (GPM) Dual-frequency Precipitation Radar (DPR) observations.

## Installation

Requires Python 3.10 or newer. Clone or extract the repository, then run:

```bash
python -m pip install -e '.[maps,notebook]'
```

## Quick start

```python
from pmwDPR import PMWDPR

results = PMWDPR(
    sensor="ATMS",
    satellite="NPP",
    start_date="2015-01-01",
    end_date="2015-01-01",
    version="V08A",
    delta_t=7.5,
    max_distance_km=45,
    track_distance_km=15,
    output_dir="./NPP_ATMS_DPR",
    plot=True,
    show=True,
)
```

The package prompts for your NASA PPS registered email. Set `PPS_PASSWORD` if your password differs from your email. Use sensor and satellite identifiers matching PPS filenames. Dates are inclusive (UTC).

## Method and outputs

- Discovers DPR and PMW granules from the PPS archive and identifies temporal candidate pairs.
- Maps each PMW channel independently to the DPR grid using the nearest valid observation within the specified time (minutes) and distance (km) limits, without interpolating brightness temperatures.
- Retains full DPR geometry and precipitation; unmatched PMW values remain NaN. Scan times are used when pixel times are unavailable.
- Saves compressed NPZ files with contributing PMW orbit numbers in the filename, JSON metadata, orbit catalogs (CSV), and optional diagnostic figures.

Use `download=False` for catalog discovery only, `dpr_orbits=[4793]` to select a DPR orbit, or `show=False` for batch processing. Track proximity is diagnostic, not a prerequisite for footprint matching.

## Local files

```python
from pmwDPR.automatic import process

process("DPR.HDF5", ["PMW.HDF5"], "matched.npz",
        delta_t=7.5, max_distance_km=45, plot=True)
```

The automatic reader supports standard FS DPR and top-level PMW swaths with ScanTime fields. Other layouts may require adapters. Live PPS access and sensor-specific layouts should be verified before large-scale processing.

## Command line and tests

```bash
pmwDPR --help
pmwDPR auto
pmwDPR inspect PMW.HDF5
pmwDPR run examples/gmi_dpr.yaml
python -m pip install -e '.[test]'
python -m pytest -q
```

The YAML-based explicit-file workflow retains its original spatial-nearest-first matching behavior; the automatic workflow applies time and distance screening before nearest-valid selection.
