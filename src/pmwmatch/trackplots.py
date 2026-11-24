"""Orbit tracks are diagnostic only; footprint matching is independent."""
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import h5py
from scipy.spatial import cKDTree
from .matching import _xyz, EARTH_KM
from .readers import _get

def track(path,group,swath):
    with h5py.File(path) as h:
        for a,b in [('navigation/scLat','navigation/scLon'),('SCstatus/SClatitude','SCstatus/SClongitude')]:
            if f'{group}/{a}' in h and f'{group}/{b}' in h:
                lat=_get(h,f'{group}/{a}').reshape(-1); lon=_get(h,f'{group}/{b}').reshape(-1)
                if len(lat)==swath['lat'].shape[0] and len(lon)==len(lat): break
        else:
            c=swath['lat'].shape[1]//2; lat=swath['lat'][:,c]; lon=swath['lon'][:,c]
            a='central footprint'
    t=swath['time'][:,swath['time'].shape[1]//2]
    ok=np.isfinite(lat)&np.isfinite(lon)&np.isfinite(t)&(abs(lat)<=90)&(lon>=-180)&(lon<=360)
    return lat[ok],(lon[ok]+180)%360-180,t[ok],a

def events(a,b,minutes,km,separation=500):
    x=_xyz(a[0],a[1]); y=_xyz(b[0],b[1]); candidates=[]
    if not len(x) or not len(y): return []
    for i,js in enumerate(cKDTree(y).query_ball_point(x,2*np.sin(km/(2*EARTH_KM)))):
        for j in js:
            dt=abs(a[2][i]-b[2][j]); d=np.linalg.norm(x[i]-y[j])
            if dt<=minutes*60: candidates.append((d,dt,i,j))
    chosen=[]
    for d,dt,i,j in sorted(candidates):
        if all(np.linalg.norm(x[i]-x[k])>2*np.sin(separation/(2*EARTH_KM)) for k,l in chosen): chosen.append((i,j))
    return chosen

def _finish(fig,path,show):
    import matplotlib.pyplot as plt
    fig.savefig(path,dpi=250,bbox_inches='tight')
    if show:
        from IPython.display import display, Image
        display(Image(filename=str(path)))
    plt.close(fig)

def plot_pair(df,pf,dpr,swaths,matches,out,minutes,km,show):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cf
    g=next(iter(swaths)); a=track(df,'FS',dpr); b=track(pf,g,swaths[g]); ev=events(a,b,minutes,km)
    fig,ax=plt.subplots(figsize=(13,7),subplot_kw={'projection':ccrs.PlateCarree()}); ax.set_global(); ax.stock_img(); ax.coastlines(linewidth=.6); ax.add_feature(cf.BORDERS,linewidth=.3)
    gl=ax.gridlines(draw_labels=True,linestyle='--',alpha=.4); gl.top_labels=False; gl.right_labels=False
    for tr,color,label in [(a,'red','DPR'),(b,'blue','PMW')]:
        lon=tr[1].copy(); lon[1:][abs(np.diff(lon))>180]=np.nan
        ax.plot(lon,tr[0],color=color,lw=1.7,transform=ccrs.PlateCarree(),label=label+(' (central footprint)' if tr[3]=='central footprint' else ' spacecraft track'))
    mask=np.logical_or.reduce([np.isfinite(r['pmw_tb']).any(axis=-1) for r in matches]); ax.scatter(dpr['lon'][mask],dpr['lat'][mask],s=1,c='limegreen',transform=ccrs.PlateCarree(),label='Matched DPR footprints')
    fmt=lambda t:datetime.fromtimestamp(t,timezone.utc).strftime('%H:%M:%S')
    for i,j in ev:
        lat,lon=a[0][i],a[1][i]; ax.plot(lon,lat,'ko',transform=ccrs.PlateCarree()); ax.annotate(f'Δt: {abs(a[2][i]-b[2][j])/60:.1f} min\nDPR: {fmt(a[2][i])}\nPMW: {fmt(b[2][j])}',(lon,lat),xytext=(-12 if lon>130 else 12,-15 if lat>55 else 12),textcoords='offset points',ha='right' if lon>130 else 'left',va='top' if lat>55 else 'bottom',fontsize=8,bbox=dict(fc='white',alpha=.85),arrowprops=dict(arrowstyle='-'))
    ax.legend(loc='upper right',fontsize=8); ax.set_title(f'DPR: {Path(df).name}\nPMW: {Path(pf).name}\n{mask.sum():,} matched DPR pixels; {len(ev)} track proximity events',fontsize=10)
    _finish(fig,out.with_name(out.stem+'_'+Path(pf).stem+'_tracks.png'),show)

def plot_grid(p,out,minutes,km,show):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    fig,axes=plt.subplots(2,2,figsize=(14,8),subplot_kw={'projection':ccrs.PlateCarree()},layout='constrained')
    mask=p['coincidence_mask'].any(-1); lat=p['latitude']; lon=p['longitude']; ok=np.isfinite(lat)&np.isfinite(lon)&(abs(lat)<=90)
    for ax,v,title,cmap,lim in zip(axes.flat,[p['dpr_precipitation'],mask.astype(float),np.max(np.where(np.isfinite(p['delta_t_minutes']),p['delta_t_minutes'],-np.inf),axis=-1),np.max(np.where(np.isfinite(p['distance_km']),p['distance_km'],-np.inf),axis=-1)],['DPR precipitation (mm/h)','Coincidence mask','Maximum time difference (min)','Maximum distance (km)'],['turbo','Greens','viridis','plasma'],[10,1,minutes,km]):
        valid=ok&np.isfinite(v); ax.set_global(); ax.coastlines(); im=ax.scatter(lon[valid],lat[valid],c=v[valid],s=2,cmap=cmap,vmin=0,vmax=lim,transform=ccrs.PlateCarree()); ax.set_title(title); fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.7)
    _finish(fig,out.with_name(out.stem+'_grid.png'),show)
