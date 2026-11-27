"""Swath meshes without scatter fallback."""
import numpy as np
from scipy.ndimage import distance_transform_edt

def is_gmi(path):
    import re
    from pathlib import Path
    return bool(re.search(r'(^|[._-])GMI([._-]|$)', Path(path).name, re.I))

def swath_mesh(ax, lon, lat, values, mask=None, **kwargs):
    x=np.asarray(lon,dtype=float); y=np.asarray(lat,dtype=float); z=np.asarray(values,dtype=float)
    geo=np.isfinite(x)&np.isfinite(y)&(abs(y)<=90)&(abs(x)<=360)
    valid=geo&np.isfinite(z)
    if mask is not None: valid &= mask
    if not valid.any(): return None
    # Fill invalid coordinates only for mesh geometry; hide adjacent cells
    # whose inferred corners would depend on missing coordinates.
    from scipy.ndimage import binary_dilation
    if not geo.all():
        idx=distance_transform_edt(~geo,return_distances=False,return_indices=True)
        x=x[tuple(idx)]; y=y[tuple(idx)]
        valid &= ~binary_dilation(~geo,structure=np.ones((3,3)))
    # Hide cells adjacent to longitude discontinuities, not bridge the seam.
    for axis in (0,1):
        jump=abs(np.diff(x,axis=axis))>180
        a=[slice(None)]*2; b=a.copy(); a[axis]=slice(None,-1); b[axis]=slice(1,None)
        valid[tuple(a)] &= ~jump; valid[tuple(b)] &= ~jump
    return ax.pcolormesh(x,y,np.ma.masked_where(~valid,z),shading='nearest',rasterized=True,**kwargs)
