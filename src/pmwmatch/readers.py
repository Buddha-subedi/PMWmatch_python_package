"""Configurable HDF5/NetCDF readers. Inspect file paths before running."""
from pathlib import Path
import numpy as np
import h5py
import xarray as xr

def inspect_swaths(path):
    """Return HDF5 groups containing candidate brightness temperature arrays."""
    out = []
    with h5py.File(path, 'r') as h:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and any(k in name.lower() for k in ('tc','tb','brightness')):
                out.append((name, obj.shape, str(obj.dtype)))
        h.visititems(visit)
    return out

def _get(h, name):
    if name not in h:
        raise KeyError(f"Dataset {name!r} absent. Run `pmwmatch inspect FILE` to inspect candidates.")
    a = h[name][...]
    if np.ma.isMaskedArray(a): a = a.filled(np.nan)
    a = np.asarray(a)
    if np.issubdtype(a.dtype, np.number):
        a = a.astype('float64')
        attrs = h[name].attrs
        for key in ('_FillValue', 'missing_value'):
            if key in attrs:
                a[a == float(np.asarray(attrs[key]).flat[0])] = np.nan
    return a

def _timestamps(h, cfg, shape):
    """Produce Unix seconds. Either explicit per-pixel seconds or scan time components."""
    if cfg.get('time_path'):
        a = _get(h, cfg['time_path'])
        unit = cfg.get('time_unit', 'unix_s')
        if unit == 'unix_ms': a = a / 1000.
        elif unit != 'unix_s': raise ValueError('time_unit must be unix_s or unix_ms')
    elif cfg.get('time_components'):
        from datetime import datetime, timezone
        c = cfg['time_components']
        keys = ['year','month','day','hour','minute','second']
        vals = {k: _get(h, c[k]).astype(int) if k in c else np.zeros(shape[0], dtype=int) for k in keys}
        a = np.array([datetime(int(vals['year'][i]), int(vals['month'][i]), int(vals['day'][i]), int(vals['hour'][i]), int(vals['minute'][i]), int(vals['second'][i]), tzinfo=timezone.utc).timestamp() for i in range(shape[0])])
        if 'millisecond' in c: a += _get(h, c['millisecond']) / 1000.
    else:
        raise ValueError('Configure time_path or time_components for exact beam-level times')
    if a.shape == shape: return a
    if a.shape == (shape[0],): return np.broadcast_to(a[:,None], shape).copy()
    raise ValueError(f'Time shape {a.shape} cannot broadcast to geolocation shape {shape}')

def _read_hdf(path, cfg, value_key):
    with h5py.File(path, 'r') as h:
        lat = _get(h, cfg['lat_path']); lon = _get(h, cfg['lon_path'])
        if lat.shape != lon.shape or lat.ndim != 2: raise ValueError('Latitude/longitude must be matching 2D arrays')
        data = _get(h, cfg[value_key]); t = _timestamps(h, cfg, lat.shape)
        if data.shape[:2] != lat.shape: raise ValueError(f'Data shape {data.shape} incompatible with geo shape {lat.shape}')
        return dict(lat=lat, lon=lon, time=t, data=data)

def read_pmw(path, cfg): return _read_hdf(path, cfg, 'tb_path')
def read_dpr(path, cfg): return _read_hdf(path, cfg, 'precip_path')

def read_mrms(path, cfg):
    """Read gridded NetCDF or GRIB2; timestamps may come from the file or configuration."""
    path = Path(path)
    engine = 'cfgrib' if path.suffix.lower() in ('.grib2','.grb2','.grib') else None
    with xr.open_dataset(path, engine=engine) as ds:
        da = ds[cfg['variable']].squeeze(drop=True)
        if da.ndim != 2: raise ValueError('MRMS variable must be a single 2D time slice')
        latname=cfg.get('lat_name','latitude'); lonname=cfg.get('lon_name','longitude')
        lat=np.asarray(ds[latname]); lon=np.asarray(ds[lonname]); data=np.asarray(da.values, dtype=float)
        if lat.ndim == 1 and lon.ndim == 1: lon,lat=np.meshgrid(lon,lat)
        if lat.shape != data.shape: raise ValueError('MRMS coordinate shape mismatch')
        time_name=cfg.get('time_name','time')
        if 'unix_time' in cfg: t=float(cfg['unix_time'])
        elif time_name in ds:
            raw=np.asarray(ds[time_name].values).reshape(-1)[0]
            t=float(raw.astype('datetime64[ns]').astype('int64') / 1e9) if np.issubdtype(type(raw),np.datetime64) or np.issubdtype(np.asarray(raw).dtype,np.datetime64) else float(raw)
        else: raise ValueError('MRMS time missing; set unix_time in config')
        data[(~np.isfinite(data)) | (data < cfg.get('valid_min',0))] = np.nan
        return dict(lat=lat,lon=lon,time=np.full(lat.shape,t),data=data)
