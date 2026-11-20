"""Ground-track crossing discovery with longitude unwrapping and interpolated times."""
import numpy as np
from scipy.spatial import cKDTree
from .matching import _xyz, EARTH_KM

def _centerline(swath):
    mid=swath['lat'].shape[1]//2
    return swath['lat'][:,mid], swath['lon'][:,mid], swath['time'][:,mid]

def find_intersections(pmw, dpr, coincidence_minutes=7., search_radius_km=25.):
    """Approximate polyline intersection via closest centerline segment pairs.

    Uses 3D proximity to identify candidates, then local tangent-plane line segment
    intersections; avoids dateline discontinuity. Returns all unique crossings.
    """
    la,lo,ta=_centerline(pmw); lb,lp,tb=_centerline(dpr)
    good_a=np.isfinite(la)&np.isfinite(lo)&np.isfinite(ta)
    good_b=np.isfinite(lb)&np.isfinite(lp)&np.isfinite(tb)
    ia=np.flatnonzero(good_a[:-1]&good_a[1:]); ib=np.flatnonzero(good_b[:-1]&good_b[1:])
    if not len(ia) or not len(ib): return []
    mid_a=_xyz((la[ia]+la[ia+1])/2,(lo[ia]+lo[ia+1])/2)
    mid_b=_xyz((lb[ib]+lb[ib+1])/2,(lp[ib]+lp[ib+1])/2)
    tree=cKDTree(mid_b)
    # Expand by segment length; large gaps must not be treated as actual crossings.
    candidates=tree.query_ball_point(mid_a,r=2*np.sin((search_radius_km+40)/(2*EARTH_KM)))
    results=[]
    for ai,neighbors in enumerate(candidates):
        i=ia[ai]
        for bj in neighbors:
            j=ib[bj]
            lat0=np.deg2rad(np.mean([la[i],la[i+1],lb[j],lb[j+1]]))
            lon0=lo[i]
            def xy(lat,lon):
                dlon=((lon-lon0+180)%360)-180
                return np.array([EARTH_KM*np.deg2rad(dlon)*np.cos(lat0),EARTH_KM*np.deg2rad(lat)])
            p=xy(la[i],lo[i]); r=xy(la[i+1],lo[i+1])-p
            q=xy(lb[j],lp[j]); s=xy(lb[j+1],lp[j+1])-q
            cross=lambda u,v:u[0]*v[1]-u[1]*v[0]
            den=cross(r,s)
            if abs(den)<1e-12: continue
            u=cross(q-p,s)/den; v=cross(q-p,r)/den
            if not (0<=u<=1 and 0<=v<=1): continue
            tpmw=ta[i]+u*(ta[i+1]-ta[i]); tdpr=tb[j]+v*(tb[j+1]-tb[j])
            dt=abs(tpmw-tdpr)/60
            if dt>coincidence_minutes: continue
            lat=float(la[i]+u*(la[i+1]-la[i]))
            lon=float(((lo[i]+u*(((lo[i+1]-lo[i]+180)%360)-180)+180)%360)-180)
            if not any(abs(x['pmw_time']-tpmw)<1 and abs(x['dpr_time']-tdpr)<1 for x in results):
                results.append(dict(latitude=lat,longitude=lon,pmw_time=float(tpmw),dpr_time=float(tdpr),delta_minutes=float(dt),pmw_scan=int(i),dpr_scan=int(j)))
    return sorted(results,key=lambda x:x['dpr_time'])
