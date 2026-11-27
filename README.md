# pmwDPR

**Passive microwave–precipitation radar swath collocation in Python.**

[![PyPI](https://img.shields.io/pypi/v/pmwDPR)](https://pypi.org/project/pmwDPR/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

`pmwDPR` automates the discovery, downloading and spatiotemporal collocation of passive microwave (PMW) brightness temperatures with Global Precipitation Measurement (GPM) Dual-frequency Precipitation Radar (DPR) observations. It preserves the native DPR grid and provides channel-wise collocation information for precipitation retrieval, satellite intercomparison and machine-learning dataset development.

## Installation

Python 3.10 or newer is required.

```bash
pip install pmwDPR
```

For geographical visualization and notebook support:

```bash
pip install "pmwDPR[maps,notebook]"
```

## Quick start

```python
from pmwDPR import PMWDPR

results = PMWDPR(
    sensor="SSMIS",
    satellite="F16",
    version="V08A",
    start_date="2015-01-01",
    end_date="2015-01-01",
    delta_t=7.5,
    max_distance_km=45,
    track_distance_km=15,
    output_dir="./F16_SSMIS_DPR_V08A",
    download=True,
    plot=True,
    show=True,
    n_plot_channels=3,
)
```

The package prompts for your NASA PPS registered email unless `email` or `PPS_USERNAME` is supplied. Set `PPS_PASSWORD` if your archive password differs from your email. Register for archive access at [NASA PPS](https://pps.gsfc.nasa.gov/).

## Sensor selection

Example sensor/satellite filename tokens:

| Sensor | Satellite | Visualization |
|---|---|---|
| SSMIS | F16, F17, F18 | pcolormesh |
| ATMS | NPP | pcolormesh |
| AMSR2 | GCOMW1 | pcolormesh |
| MHS | PPS satellite token | pcolormesh |
| GMI | GPM | Automatically disabled |

Archive availability, product version and file layout vary by sensor and date. The automatic reader supports standard FS DPR and top-level PMW swaths with latitude, longitude, Tc or Tb, and ScanTime fields. Other layouts may require adapters.

## Matching methodology

1. Discover DPR and PMW granules in the NASA PPS archive.
2. Identify temporally overlapping candidate orbits, including the specified time tolerance.
3. Search spatial neighbors using a spherical-coordinate KD-tree.
4. Screen candidates by absolute time difference and maximum great-circle distance.
5. Select the nearest valid PMW footprint independently for each channel and DPR pixel.
6. Merge contributing PMW granules and save all available channels on the native DPR grid.

Brightness temperatures are not spatially interpolated or averaged. Unmatched values remain NaN. Scan times are broadcast across pixels when beam-level timing is unavailable. The full DPR geometry and precipitation field are retained.

## Parameters

| Parameter | Meaning |
|---|---|
| `sensor`, `satellite` | PPS filename tokens |
| `version` | Product version filter, e.g. V08A |
| `start_date`, `end_date` | Inclusive UTC dates |
| `delta_t` | Maximum absolute time difference (minutes) |
| `max_distance_km` | Maximum footprint separation (km) |
| `track_distance_km` | Track proximity threshold for diagnostics, not footprint matching |
| `output_dir` | Output directory |
| `download` | False for catalog discovery only |
| `plot` | Generate diagnostic figures |
| `show` | Display diagnostic figures |
| `n_plot_channels` | Number of channels requested for visualization only |
| `dpr_orbits` | Optional list of selected DPR orbit numbers |

## Output files

```text
F16_SSMIS_DPR_V08A/
├── DPR/                         # Original DPR granules
├── PMW/                         # Original PMW granules
├── matches/
│   ├── <DPR>_F16_SSMIS_<PMW_orbits>.npz
│   ├── <DPR>_F16_SSMIS_<PMW_orbits>.json
│   └── Diagnostic intersection figures (.png)
├── DPR_available_orbits.csv
├── PMW_available_orbits.csv
└── temporal_pairs.csv
```

### Filename convention

```text
<DPR_filename_stem>_<satellite>_<sensor>_<PMW_orbit_1>_<PMW_orbit_2>....npz
```

PMW orbit numbers are sorted and formatted with at least six digits. Only orbits contributing to the final matched channels are included. If no valid footprints are matched, the suffix is `no_match`. A JSON sidecar records source granules, contributing orbits, timing modes and matching thresholds. The returned `results` list contains saved file paths (or None when processing yields no merged swaths).

### NPZ variables

N = DPR scans, M = DPR cross-track pixels, C = total available PMW channels.

| Variable | Shape | Description |
|---|---|---|
| `latitude` | (N, M) | DPR latitude (degrees) |
| `longitude` | (N, M) | DPR longitude (degrees) |
| `dpr_precipitation` | (N, M) | DPR near-surface precipitation rate (mm/h) |
| `pmw_tb` | (N, M, C) | Mapped brightness temperatures (K) |
| `channel_names` | (C,) | Channel identifiers |
| `dpr_time` | (N, M) | DPR time (Unix seconds, UTC) |
| `pmw_time` | (N, M, C) | Matched PMW time (Unix seconds, UTC) |
| `distance_km` | (N, M, C) | Spatial separation (km) |
| `delta_t_minutes` | (N, M, C) | Absolute temporal separation (minutes) |
| `coincidence_mask` | (N, M, C) | Boolean valid-match mask |

**All available channels are saved.** The visualization channel count does not affect the output data. Channels are concatenated in sorted swath-group order. Generic channel names identify the swath and channel index; the recognized 11-channel SSMIS layout uses frequency/polarization labels.

### Reading the results

```python
import numpy as np

for filename in results:
    if filename is None:
        continue
    with np.load(filename) as data:
        print(filename)
        print(data["pmw_tb"].shape)
        print(data["channel_names"])
        print("Matched fraction:", data["coincidence_mask"].mean())
```

## Visualization and batch processing

Diagnostic intersection plots compare original PMW swaths and mapped TB fields using `pcolormesh`. GMI automatically skips visualization while retaining collocation and file output. For bulk processing of other sensors, use `plot=False, show=False`.

## Processing local files

```python
from pmwDPR.automatic import process

saved_file = process(
    "DPR.HDF5", ["PMW.HDF5"], "matched.npz",
    delta_t=7.5, max_distance_km=45, plot=False,
)
```

Use actual PPS granule filenames with orbit identifiers. The output name receives the contributing-orbit suffix.

## Command-line interface

```bash
pmwDPR --help
pmwDPR auto
pmwDPR inspect PMW.HDF5
```

The explicit-file YAML workflow (`pmwDPR run`) retains its original spatial-nearest-first behavior; the automatic workflow screens both time and distance before nearest-valid selection.

## Development

```bash
git clone https://github.com/Buddha-subedi/pmwDPR_python_package.git
cd pmwDPR_python_package
pip install -e ".[maps,notebook,test]"
python -m pytest -q
```

## Author and support

Developed by Buddha Subedi. Report issues or request features through [GitHub Issues](https://github.com/Buddha-subedi/pmwDPR_python_package/issues).

Data products are provided by NASA GPM and PPS. Please acknowledge the relevant satellite data providers when using collocated datasets in publications.
