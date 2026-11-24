# pmwmatch v0.1.0 — MSI research starter package

`PMWDPR` maps **original PMW brightness temperatures onto the DPR grid**. `PMWMRMS` maps **original PMW brightness temperatures onto the MRMS grid**. The source PMW pixel is selected by **spatial nearest neighbor**, then rejected if its per-beam observation time differs by more than **7.5 minutes** or its distance exceeds **214 km** (configurable). Missing matches remain NaN. **No temporal interpolation of TBs** is performed.

## What works now

- Configurable HDF5 swath and channel readers (with per-scan time components or per-pixel Unix timestamps).
- HDF5 TB/swath inventory: `pmwmatch inspect FILE.HDF5`.
- DPR–PMW centerline crossing candidate detection for different satellites, default ±7 minutes.
- Same-platform GMI–DPR option that bypasses crossing requirement but still applies pixel-level screening.
- PMW to DPR/MRMS nearest-neighbor mapping with spatial and temporal diagnostics.
- NetCDF or GRIB2 MRMS readers (GRIB2 needs optional cfgrib/eccodes).
- PNG plots of reference precipitation, mapped TBs, masks, distance, time difference, and orbit centerlines.
- HTTPS file downloader with `.part` resume, environment-based PPS credentials, and explicit URL lists.
- Unit tests and Slurm job example.

## Important limitations — read before research use

**This is a functional starter implementation, not a validated production archive downloader.** Automated PPS orbit discovery, listing historical DPR/PMW products, identifying available swath groups by sensor, and automated historical MRMS catalog lookup are **not implemented**. You must provide direct verified file URLs or local files. The download helper does not discover file URLs. PPS login method/endpoint can change; verify against your current PPS account. Authentication is HTTPS Basic; a provider requiring a different login method needs an adapter. Do not assume any file URL format or HDF5 group path is universal.

**Swath selection is manual** via the reader configuration; `inspect` prints candidate datasets but does not prompt to select. Some sensors have different TB dataset names, dimensions, scaling conventions, and time representations. The supplied GMI/DPR paths are *examples only*, not verified against your particular product version. Confirm calibration, fill values, units, and time conventions for every instrument. The code reads raw numerical values and applies fill flags; it does **not** automatically apply instrument-specific scale/offset calibration or quality screening.

The crossing detector uses local tangent-plane segment intersections after candidate pruning. Its fixed candidate radius is not robust for very long scan gaps or extreme polar cases; verify candidate crossings visually and against orbit times. A 214 km radius can be far larger than a physical PMW footprint; reduce for science-grade footprint matching. The algorithm intentionally chooses the nearest spatial PMW pixel *before* checking time; it does not search for a farther time-valid pixel. The MRMS workflow maps onto the entire input grid; crop to the relevant region for performance. It does not implement antenna-footprint averaging or MRMS quality flags. `plot` uses scatter geolocation, not a cartographic coastline projection.

## Install on MSI

```bash
unzip pmwmatch_msi.zip
cd pmwmatch_package
# Activate an existing Python >=3.10 environment, or create a new one:
conda create -n pmwmatch python=3.11 -y
conda activate pmwmatch
python -m pip install -e '.[mrms,test]'
pytest -q
pmwmatch --help
```

For MSI OnDemand Jupyter, use `python -m ipykernel install --user --name pmwmatch --display-name 'Python (pmwmatch)'` after installing `ipykernel`.

## Discover swaths and set the target

```bash
pmwmatch inspect /path/to/GMI_FILE.HDF5
```

Inspect the output and set `pmw_reader.lat_path`, `lon_path`, `tb_path` and time fields in YAML. Example GMI S1 and S2 are separate groups; map each separately by running two configurations. For SSMIS and other sensors, check the *actual* swath names in the chosen product. DPR always supplies target geometry in `PMWDPR`; MRMS supplies target geometry in `PMWMRMS`.

## Run

```bash
# Edit paths and HDF5 dataset paths first!
pmwmatch run examples/gmi_dpr.yaml
pmwmatch run examples/ssmis_dpr.yaml
pmwmatch run examples/gmi_mrms.yaml
pmwmatch plot /path/to/output.npz --channel 0 --step 4
```

Or in Python:

```python
from pmwmatch import PMWDPR, PMWMRMS
result = PMWDPR(pmw_file='gmi.h5', dpr_file='dpr.h5',
                pmw_config={'lat_path':'S1/Latitude','lon_path':'S1/Longitude',
                            'tb_path':'S1/Tc','time_path':'S1/UnixTime'},
                dpr_config={'lat_path':'FS/Latitude','lon_path':'FS/Longitude',
                            'precip_path':'FS/SLV/precipRateNearSurface','time_path':'FS/UnixTime'},
                output_file='gmi_dpr.npz', same_platform=True)
```

The example `UnixTime` fields above are illustrative; PPS may store scan time as separate year/month/day/hour/minute/second arrays. Use YAML time_components as shown in `examples/gmi_dpr.yaml`.

## Download known file URLs

Create `pps_urls.txt` containing one **verified direct HTTPS file URL per line**. Never put passwords in this file.

```bash
export PPS_USERNAME='your_registered_username'
read -rs -p 'PPS password: ' PPS_PASSWORD; echo
export PPS_PASSWORD
pmwmatch download pps_urls.txt --pps --out /scratch.global/1subed042/pmwmatch/raw
unset PPS_PASSWORD
```

For public MRMS direct HTTPS URLs, omit `--pps`. For NOAA historical MRMS archives, first determine the correct source and product (2-minute instantaneous precipitation rate versus accumulation); then supply direct URLs or local NetCDF/GRIB2 files. The downloader skips complete existing files and attempts HTTP Range resume for partial downloads; it does not verify checksums.

## Output

`npz` arrays: `latitude`, `longitude`, `time_unix_s`, `dpr_precipitation` or `mrms_precipitation`, `pmw_tb`, `valid_mask`, `distance_km`, `delta_t_minutes`. Sidecar JSON has match counts and accepted centerline crossings. `pmw_tb` shape is `[target_scan, target_pixel, channel]` or `[MRMS_y, MRMS_x, channel]`.

## Science checks before publication

1. Validate PPS HDF5 paths, sensor calibration and quality flags with a real file.
2. Compare actual PMW/DPR beam times, especially GMI S1/S2.
3. Verify longitude wrapping and crossing behavior at poles/dateline.
4. Confirm MRMS timestamp, precipitation variable, units, quality flags, and coverage.
5. Decide if 214 km is scientifically appropriate (default is intentionally user-requested, not a footprint size).
6. Compare matched pixel statistics and quicklooks with independently known coincidences.

## Version 0.2: automated PPS–DPR workflow (development)

Install `pip install -e '.[maps,test]'`, then run `pmwmatch auto` or:

```python
from pmwmatch import PMWDPR
PMWDPR()
```

The first prompt requests the PPS registered email. The archive backend uses HTTPS Basic authentication (email as password by default; override with PPS_PASSWORD). It discovers filenames in PPS daily `1C` and `radar` directories rather than submitting STORM web orders. Register for PPS access first. Live archive authentication and real sensor products have not been validated in this development build.

```python
PMWDPR(sensor='SSMIS', satellite='F16', start_date='2015-01-01',
       end_date='2015-01-01', delta_t=7.5, max_distance_km=15,
       output_dir='./data', plot=True)
```

Use exact sensor and platform tokens in PPS filenames. Dates are inclusive UTC granule-overlap dates; full DPR granules are retained. The default 15 km radius is configurable, not a universal footprint radius. Every discovered Tc/Tb swath is mapped independently, selecting the nearest valid observation per channel after time and distance filtering. No interpolation. Unmatched TB, distance and PMW time are NaN. Scan timestamps are broadcast when beam timestamps are unavailable; the JSON records this. Channel labels are generic swath/channel indices, not verified frequency labels. Output includes full DPR precipitation and coordinates, but plots contain only TBs. Three matched channels are plotted when available. Cartopy coastlines may require a Natural Earth download.

For already downloaded files, use the new multi-swath workflow directly:

```python
from pmwmatch.automatic import process
process('DPR.HDF5', ['PMW1.HDF5', 'PMW2.HDF5'], 'matched.npz',
        delta_t=7.5, max_distance_km=15, plot=True)
```

Legacy explicit-file and MRMS interfaces remain available and retain their original nearest-first matching semantics. Use the automatic workflow or `process` for the new time-constrained nearest-valid algorithm. The automatic reader currently supports standard FS DPR and top-level PMW swaths with ScanTime calendar fields; other layouts fail explicitly and need adapters. Do not assume all PPS products are supported until tested. Avoid mixing product versions in one run (set `version` when appropriate).

Run tests: `python -m pytest -q`.

## v0.2.1: PPS text catalog and orbit diagnostics

```python
from pmwmatch import PMWDPR
results = PMWDPR(sensor='SSMIS', satellite='F16', start_date='2015-01-01',
                 end_date='2015-01-01', delta_t=7.5, max_distance_km=15,
                 output_dir='./SSMIS_DPR', plot=True, show=True)
```
Use `download=False` for catalog-only mode, `dpr_orbits=[4793]` for a single orbit,
`version='V08A'` to constrain the archive version, and `show=False` for batch jobs.
Install with `pip install -e '.[maps,notebook]'` for notebook figures.
Catalog CSVs list full granules and temporal candidate pairs. Track figures are
created for pairs with matched TB footprints; black annotations indicate separate
track proximity events (default 15 km, configurable via `track_distance_km`).
If spacecraft navigation is absent, the central footprint is explicitly labeled.
Track proximity is not a prerequisite for footprint matching. Scan timestamps
are broadcast when beam timestamps are absent. No inferred beam times are used.
V08A HDF5-backed NetCDF paths were confirmed from user inspection logs; live
end-to-end matching has not been tested by the package author in this environment.
