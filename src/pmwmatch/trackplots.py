"""Orbit tracks are diagnostic only; footprint matching is independent."""
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import h5py
from scipy.spatial import cKDTree
from .matching import _xyz, EARTH_KM
from .readers import _get
from .orbits import find_intersections

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

def plot_grid(p, out, minutes, km, show):
    """Disabled legacy diagnostic plot; retained for API compatibility."""
    return None

# GPM Level 1C nominal center frequencies (GHz), in stored channel order.
# Reference: GPM Constellation L1C Product Format, sections 1, 3–6.
# Sideband channels are represented by their nominal center frequency.
SENSOR_LAYOUTS = {
    'ATMS': {'S1':[23.8], 'S2':[31.4], 'S3':[88.2],
             'S4':[165.5,183.31,183.31,183.31,183.31,183.31]},
    'AMSR2': {'S1':[10.65,10.65], 'S2':[18.7,18.7],
              'S3':[23.8,23.8], 'S4':[36.5,36.5],
              'S5':[89,89], 'S6':[89,89]},
    'GMI': {'S1':[10.65,10.65,18.7,18.7,23.8,36.5,36.5,89,89],
            'S2':[166,166,183.31,183.31]},
    'MHS': {'S1':[89,157,183.31,183.31,190.31]},
    'SSMIS': {'S1':[19.35,19.35,22.235], 'S2':[37,37],
              'S3':[150,183.31,183.31,183.31], 'S4':[91.665,91.665]},
}
FREQUENCY_OVERRIDES = {}

def channel_frequencies(path, sensor, swaths):
    result = {}
    with h5py.File(path, 'r') as h:
        for g, src in swaths.items():
            n = src['data'].shape[-1]
            group = h[g]
            ds = group['Tc'] if 'Tc' in group else group['Tb']
            for obj in (ds, group):
                for key, val in obj.attrs.items():
                    if 'frequen' not in key.lower(): continue
                    try: f = np.asarray(val, dtype=float).ravel()
                    except (ValueError, TypeError): continue
                    if len(f) == n and np.all((f > 0) & (f < 1000)):
                        result.update({(g,i): float(x) for i,x in enumerate(f)})
            for key, obj in group.items():
                if isinstance(obj,h5py.Dataset) and 'frequen' in key.lower() and obj.size == n:
                    try: f = np.asarray(obj[...],dtype=float).ravel()
                    except (ValueError,TypeError): continue
                    if np.all((f > 0)&(f < 1000)):
                        result.update({(g,i):float(x) for i,x in enumerate(f)})
    for g, fs in SENSOR_LAYOUTS.get(sensor, {}).items():
        if g not in swaths:
            continue
        n = swaths[g]['data'].shape[-1]
        if n == len(fs):
            for i, f in enumerate(fs):
                result.setdefault((g, i), f)
        else:
            print(f'[CHANNEL WARNING] {sensor} {g}: expected {len(fs)}, found {n}; using metadata only')
    result.update(FREQUENCY_OVERRIDES.get(sensor, {}))
    print('[CHANNELS]', sensor, {g: src['data'].shape[-1] for g,src in swaths.items()}, result)
    return result

def plot_pair(df, pf, dpr, swaths, matches, out, minutes, km, show):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from matplotlib.patches import Rectangle
    from matplotlib.colors import Normalize
    from matplotlib.ticker import MaxNLocator
    sensor = Path(pf).name.split('.')[2].upper()
    a = track(df,'FS',dpr)
    g0 = next(iter(swaths))
    b = track(pf,g0,swaths[g0])
    def as_swath(t):
        return {k: np.asarray(v)[:,None] for k,v in zip(('lat','lon','time'),t[:3])}
    crossings = find_intersections(as_swath(b),as_swath(a),coincidence_minutes=np.inf,search_radius_km=km)
    crossings = sorted((c for c in crossings if np.isfinite(c['delta_minutes']) and abs(c['delta_minutes']) <= minutes),key=lambda c:c['dpr_time'])
    if not crossings:
        print(f'[SKIP] {sensor}: no time-valid geometric intersection'); return
    frequencies = channel_frequencies(pf,sensor,swaths)
    pair = dict(zip(swaths,matches))
    for number,cross in enumerate(crossings,1):
        lon0,lat0 = cross['longitude'],cross['latitude']
        # Fixed regional extent; longitude relative to crossing.
        half_lat = 15.0; half_lon = 25.0
        south,north = max(-90.,lat0-half_lat),min(90.,lat0+half_lat)
        wrap = lambda x: (np.asarray(x)-lon0+180)%360-180
        def region(s):
            return (abs(wrap(s['lon']))<=half_lon)&(s['lat']>=south)&(s['lat']<=north)&np.isfinite(s['lat'])&np.isfinite(s['lon'])
        dm = region(dpr)
        eligible = []
        for (g,i),f in frequencies.items():
            if g not in pair or not 0 <= i < swaths[g]['data'].shape[-1]:
                continue
            mapped = pair[g]['pmw_tb'][...,i]
            native = swaths[g]['data'][...,i]
            count = np.count_nonzero(dm & np.isfinite(mapped))
            if count and np.any(region(swaths[g]) & np.isfinite(native)):
                eligible.append((count,g,i,f))
        choices = []
        for low, high in [(0,30),(30,90),(90,np.inf)]:
            options = [e for e in eligible if low < e[3] <= high]
            choices.append(max(options) if options else None)
        # If the instrument lacks a band, prefer a different nominal
        # frequency; then allow another polarization/sideband channel.
        for k, choice in enumerate(choices):
            if choice is not None:
                continue
            used = [e for e in choices if e is not None]
            remaining = [e for e in eligible if e not in used]
            if remaining:
                choices[k] = max(remaining, key=lambda e: (e[3] not in {u[3] for u in used}, e[0]))
                print(f'[BAND FALLBACK] {sensor}: column {k+1} uses {choices[k][3]:g} GHz')
        if any(choice is None for choice in choices):
            print(f'[SKIP FIGURE] {sensor} intersection {number}: fewer than three usable channels; choices={choices}')
            continue
        choices.sort(key=lambda e:e[3])
        print(f'[SELECTED CHANNELS] {sensor}: {choices}')
        geo=ccrs.PlateCarree(); local=ccrs.PlateCarree(central_longitude=lon0)
        with plt.rc_context({'font.family':'serif','font.size':10}):
            fig=plt.figure(figsize=(12,13.0))
            gs=fig.add_gridspec(4,3,height_ratios=[5.04,1.98,1.98,.12],hspace=.17,wspace=.05,left=.08,right=.92,top=.875,bottom=.085)
            ax=fig.add_subplot(gs[0,:],projection=local); ax.set_global(); ax.set_anchor('N'); ax.set_aspect('auto'); ax.set_box_aspect(0.5); ax.stock_img(); ax.coastlines('110m',linewidth=.55,color='black')
            for t,color,label in [(a,'blue','DPR track'),(b,'darkviolet',f'{sensor} track')]:
                ax.plot(t[1],t[0],transform=ccrs.Geodetic(),color=color,lw=1.3,label=label)
            ax.add_patch(Rectangle((-half_lon,south),2*half_lon,north-south,transform=local,fill=False,ec='black',lw=1.2,label='Coincidence region'))
            ax.plot(lon0,lat0,'ko',ms=5,transform=geo)
            ax.annotate(f"ΔT: {cross['delta_minutes']:.2f} min",xy=(lon0,lat0),xycoords=geo._as_mpl_transform(ax),xytext=(10,-12),textcoords='offset points',fontsize=9,fontweight='bold',bbox=dict(fc='white',alpha=.8,ec='none'))
            ax.legend(loc='lower center',ncol=3,fontsize=8); ax.text(.01,.97,'(a)',transform=ax.transAxes,va='top',fontweight='bold')
            gl=ax.gridlines(draw_labels=True,linewidth=.3,linestyle='--',alpha=.5); gl.right_labels=False; gl.bottom_labels=False; gl.xlabel_style={'weight':'bold'}; gl.ylabel_style={'weight':'bold'}
            def mesh(ax,s,z,mask,norm):
                # Crop to relevant scan rows and columns; scatter only if invalid coordinates prevent pcolormesh.
                ii,jj=np.where(mask & np.isfinite(z))
                if not len(ii): return None
                sl=(slice(max(0,ii.min()-1),ii.max()+2),slice(max(0,jj.min()-1),jj.max()+2))
                x=wrap(s['lon'][sl]); y=s['lat'][sl]; v=np.where(mask[sl],z[sl],np.nan)
                if min(x.shape)>=2 and np.isfinite(x).all() and np.isfinite(y).all() and np.all(abs(y)<=90) and np.all(abs(np.diff(x,axis=0))<180) and np.all(abs(np.diff(x,axis=1))<180):
                    return ax.pcolormesh(x,y,v,transform=local,shading='auto',cmap='jet',norm=norm,rasterized=True)
                return ax.scatter(wrap(s['lon'][ii,jj]),s['lat'][ii,jj],c=z[ii,jj],s=2,transform=local,cmap='jet',norm=norm,rasterized=True)
            for col,choice in enumerate(choices):
                axes=[fig.add_subplot(gs[row,col],projection=local) for row in (1,2)]
                for row,axr in enumerate(axes):
                    axr.set_extent([-half_lon,half_lon,south,north],crs=local); axr.set_aspect('auto'); axr.set_box_aspect(0.60); axr.set_anchor('C'); axr.coastlines('50m',linewidth=.4,color='gray')
                    gl=axr.gridlines(draw_labels=True,linewidth=.3,linestyle='--',alpha=.5); gl.top_labels=False; gl.right_labels=False; gl.left_labels=col==0; gl.bottom_labels=row==1; gl.xlocator=MaxNLocator(3); gl.ylocator=MaxNLocator(3); gl.xlabel_style={'weight':'bold'}; gl.ylabel_style={'weight':'bold'}
                    axr.text(.02,.97,f'({chr(98+row*3+col)})',transform=axr.transAxes,va='top',fontweight='bold',zorder=10)
                if choice is None:
                    for axr in axes: axr.text(.5,.5,'No verified channel\nwith regional matches',transform=axr.transAxes,ha='center',va='center',fontsize=8)
                    continue
                _,g,i,f=choice; src=swaths[g]; mapped=pair[g]['pmw_tb'][...,i]; native=src['data'][...,i]; sm=region(src)
                vals=np.concatenate([mapped[dm & np.isfinite(mapped)],native[sm & np.isfinite(native)]])
                lo,hi=np.percentile(vals,[2,98]); norm=Normalize(lo,hi if hi>lo else lo+1)
                im=mesh(axes[0],dpr,mapped,dm,norm); mesh(axes[1],src,native,sm,norm)
                for axr in axes:
                    for j in (0,-1): axr.plot(wrap(dpr['lon'][:,j]),dpr['lat'][:,j],transform=local,color='black',lw=.45)
                axes[0].set_title(f'{f:g} GHz ({g}, ch {i+1})',fontsize=10,fontweight='bold')
                cb=fig.colorbar(im,cax=fig.add_subplot(gs[3,col]),orientation='horizontal'); cb.set_label(f'TB {f:g} GHz (K)'); cb.locator=MaxNLocator(4); cb.update_ticks()
            for row_index, label in [(1, 'DPR swath'), (2, f'{sensor} swath')]:
                bounds = gs[row_index, 2].get_position(fig)
                fig.text(.955, (bounds.y0+bounds.y1)/2, label, rotation=90, ha='center', va='center', fontweight='bold')
            fig.suptitle(f'{sensor} – DPR Coincidence {number}', y=.985, fontsize=18, fontweight='bold')
            fig.text(.5, .951, f'DPR: {Path(df).name}', ha='center', va='top', fontsize=10, fontweight='bold')
            fig.text(.5, .930, f'{sensor}: {Path(pf).name}', ha='center', va='top', fontsize=10, fontweight='bold')
            path=Path(out).with_name(f'{Path(out).stem}_{Path(pf).stem}_intersection_{number:02d}.png')
            _finish(fig,path,show); print(f'[FIGURE] {path}')
