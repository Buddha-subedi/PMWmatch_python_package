"""Public workflows, local files first; downloader is a separate explicit stage."""
from pathlib import Path
import json
import numpy as np
from .readers import read_pmw,read_dpr,read_mrms
from .matching import map_nearest
from .orbits import find_intersections

def _save(outpath,target,mapped,reference_name, crossings=None):
    outpath=Path(outpath);outpath.parent.mkdir(parents=True,exist_ok=True)
    payload={'latitude':target['lat'].astype('float32'),'longitude':target['lon'].astype('float32'),
             'time_unix_s':target['time'].astype('float64'),reference_name:target['data'].astype('float32'),
             'pmw_tb':mapped['tb'],'valid_mask':mapped['valid_mask'],
             'distance_km':mapped['distance_km'],'delta_t_minutes':mapped['delta_t_minutes']}
    np.savez_compressed(outpath,**payload)
    summary={'output':str(outpath),'reference':reference_name,'pixels':int(mapped['valid_mask'].size),
             'matched_pixels':int(mapped['valid_mask'].sum()),'matched_fraction':float(mapped['valid_mask'].mean()),
             'crossings':crossings or []}
    outpath.with_suffix('.json').write_text(json.dumps(summary,indent=2))
    return summary

def PMWDPR(pmw_file,dpr_file,pmw_config,dpr_config,output_file,
           coincidence_minutes=7.,pixel_time_minutes=7.5,max_distance_km=214.,
           enforce_distance=True,require_centerline=True,same_platform=False):
    """Map original PMW TBs onto DPR geometry. GMI may set same_platform=True.

    Coincidence filter applies to crossing platforms; same-platform bypasses it.
    Pixel screening always applies, including GMI.
    """
    pmw=read_pmw(pmw_file,pmw_config); dpr=read_dpr(dpr_file,dpr_config)
    crossings=[] if same_platform else find_intersections(pmw,dpr,coincidence_minutes)
    if require_centerline and not same_platform and not crossings:
        return {'skipped':True,'reason':'No centerline crossing within tolerance','pmw_file':str(pmw_file),'dpr_file':str(dpr_file)}
    mapped=map_nearest(pmw,dpr,max_distance_km,pixel_time_minutes,enforce_distance)
    return _save(output_file,dpr,mapped,'dpr_precipitation',crossings)

def PMWMRMS(pmw_file,mrms_file,pmw_config,mrms_config,output_file,
            pixel_time_minutes=7.5,max_distance_km=214.,enforce_distance=True):
    """Map original PMW TBs to a single MRMS time slice."""
    pmw=read_pmw(pmw_file,pmw_config);mrms=read_mrms(mrms_file,mrms_config)
    mapped=map_nearest(pmw,mrms,max_distance_km,pixel_time_minutes,enforce_distance)
    return _save(output_file,mrms,mapped,'mrms_precipitation')
