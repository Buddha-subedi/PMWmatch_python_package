"""Automatic PPS archive discovery and multi-swath DPR collocation."""
import os,re,json
from datetime import datetime,timedelta,timezone
from pathlib import Path
from urllib.parse import urljoin
from html.parser import HTMLParser
import numpy as np
import h5py,requests
from scipy.spatial import cKDTree
from .matching import _xyz,EARTH_KM
from .readers import _get
from .download import download_urls

BASE='https://arthurhouhttps.pps.eosdis.nasa.gov/gpmdata/'
PAT=re.compile(r'\.(\d{8})-S(\d{6})-E(\d{6})\.')
def interval(name):
    m=PAT.search(str(name))
    if not m: raise ValueError(f'Unrecognized granule name: {name}')
    a,b=[datetime.strptime(m[1]+v,'%Y%m%d%H%M%S').replace(tzinfo=timezone.utc) for v in (m[2],m[3])]
    return a,b+timedelta(days=b<a)

def orbit(name):
    return int(re.search(r'\.(\d+)\.V\d+[A-Z]?\.',str(name),re.I)[1])

def discover(session,start,end,folder,predicate,base=BASE):
    found=set(); day=start.date()
    while day<=end.date():
        url=urljoin(base,day.strftime('%Y/%m/%d/')+folder+'/')
        from urllib.parse import urlsplit
        u=urlsplit(url); listing=f'{u.scheme}://{u.netloc}/text{u.path}'
        r=session.get(listing,timeout=90)
        print(f'[PPS] {listing}: HTTP {r.status_code}',flush=True)
        if r.status_code==404: day+=timedelta(days=1); continue
        r.raise_for_status()
        if '<html' in r.text.lower(): raise RuntimeError('Expected plain-text PPS listing, received HTML')
        names=re.findall(r'[\w.+-]+\.(?:HDF5|H5|NC4|NC)(?![\w.])',r.text,re.I)
        for name in names:
            if predicate(name) and PAT.search(name):
                a,b=interval(name)
                if a<=end and b>=start: found.add(urljoin(url,name))
        print(f'[PPS] {len(names)} listed files',flush=True)
        day+=timedelta(days=1)
    unique={}
    for x in sorted(found):
        key=orbit(x); version=re.search(r'\.V(\d+)([A-Z]?)\.',x,re.I)
        rank=(int(version[1]),version[2],x)
        if key not in unique or rank>unique[key][0]: unique[key]=(rank,x)
    return sorted(v[1] for v in unique.values())

def read_group(h,g,value):
    lat,lon=(_get(h,f'{g}/{k}') for k in ('Latitude','Longitude'))
    data=_get(h,f'{g}/{value}'); shape=lat.shape
    if lat.ndim!=2 or lon.shape!=shape or data.shape[:2]!=shape: raise ValueError(f'Invalid geometry in {g}')
    s=h[f'{g}/ScanTime']; keys=['Year','Month','DayOfMonth','Hour','Minute','Second']
    arrays=np.broadcast_arrays(*[s[k][...] for k in keys]); ts=np.full(arrays[0].shape,np.nan)
    for i in np.ndindex(ts.shape):
        try: ts[i]=datetime(*[int(a[i]) for a in arrays],tzinfo=timezone.utc).timestamp()
        except (ValueError,OverflowError): pass
    if 'MilliSecond' in s: ts+=_get(s,'MilliSecond')/1000.
    mode='beam' if ts.shape==shape else 'scan_broadcast'
    if ts.shape==(shape[0],): ts=np.broadcast_to(ts[:,None],shape).copy()
    elif ts.shape==(shape[0],1): ts=np.broadcast_to(ts,shape).copy()
    if ts.shape!=shape: raise ValueError(f'Unsupported time shape {ts.shape} in {g}')
    # Optional explicitly documented beam offsets in seconds; never infer offsets.
    if 'beamTimeOffset' in s:
        ds=s['beamTimeOffset']; units=ds.attrs.get('units',''); units=units.decode() if isinstance(units,bytes) else str(units)
        if units.lower() in ('s','seconds','second') and mode=='scan_broadcast':
            ts=ts+np.broadcast_to(ds[...],shape); mode='beam_offset'
    data[(data<0)|(~np.isfinite(data))]=np.nan
    if value in ('Tc','Tb'): data[(data>400)|(data<20)]=np.nan
    return dict(lat=lat,lon=lon,time=ts,data=data,time_mode=mode)

def read_all(path):
    out={}
    with h5py.File(path,'r') as h:
        for g in sorted(h):
            if not isinstance(h[g],h5py.Group): continue
            for v in ('Tc','Tb'):
                if v in h[g] and 'Latitude' in h[g] and 'Longitude' in h[g]:
                    a=read_group(h,g,v)
                    if a['data'].ndim==2: a['data']=a['data'][...,None]
                    out[g]=a; break
    if not out: raise ValueError(f'No supported TB swaths in {path}')
    return out

def valid_geo(a):
    return np.isfinite(a['lat'])&np.isfinite(a['lon'])&np.isfinite(a['time'])&(abs(a['lat'])<=90)&(a['lon']>=-180)&(a['lon']<=360)

def map_group(src,tgt,delta_t,max_distance_km):
    """Nearest valid footprint per channel; search all spatial candidates."""
    shape=tgt['lat'].shape; nc=src['data'].shape[-1]; size=tgt['lat'].size
    tb=np.full((size,nc),np.nan,np.float32); dist=tb.copy(); dt=tb.copy(); times=np.full((size,nc),np.nan)
    ids=np.flatnonzero(valid_geo(src)); ti=np.flatnonzero(valid_geo(tgt))
    if len(ids) and len(ti):
        xyz=_xyz(src['lat'].ravel()[ids],src['lon'].ravel()[ids]); tree=cKDTree(xyz)
        st=src['time'].ravel()[ids]; vals=src['data'].reshape(-1,nc)[ids]; tt=tgt['time'].ravel()
        radius=2*np.sin(max_distance_km/(2*EARTH_KM))
        for chunk in np.array_split(ti,max(1,(len(ti)+4095)//4096)):
            q=_xyz(tgt['lat'].ravel()[chunk],tgt['lon'].ravel()[chunk]); neighbors=tree.query_ball_point(q,radius,workers=-1)
            for k,js in enumerate(neighbors):
                if not js: continue
                i=chunk[k]; js=np.asarray(js); js=js[abs(st[js]-tt[i])<=delta_t*60]
                if not len(js): continue
                chord=np.linalg.norm(xyz[js]-q[k],axis=1); order=np.argsort(chord,kind='stable'); js=js[order]; chord=chord[order]
                for c in range(nc):
                    good=np.flatnonzero(np.isfinite(vals[js,c]))
                    if not len(good): continue
                    j=good[0]; s=js[j]; tb[i,c]=vals[s,c]; times[i,c]=st[s]; dt[i,c]=abs(st[s]-tt[i])/60; dist[i,c]=2*EARTH_KM*np.arcsin(min(1,chord[j]/2))
    return {k:v.reshape(*shape,nc) for k,v in dict(pmw_tb=tb,distance_km=dist,delta_t_minutes=dt,pmw_time=times).items()}

def process(dpr_file,pmw_files,output_file,delta_t=7.5,max_distance_km=15,plot=True,n_plot_channels=3,show=False,track_distance_km=15):
    if not np.isfinite(delta_t) or delta_t<0 or not 0<max_distance_km<=20000: raise ValueError('Invalid tolerance')
    with h5py.File(dpr_file,'r') as h: dpr=read_group(h,'FS','SLV/precipRateNearSurface')
    from .trackplots import plot_pair, plot_grid
    merged={}; sources={}; modes={}
    pmw_files=list(pmw_files)
    out=Path(output_file); out.parent.mkdir(parents=True,exist_ok=True)
    for file in pmw_files:
        swaths=read_all(file); pair_matches=[]
        for g,src in swaths.items():
            r=map_group(src,dpr,delta_t,max_distance_km); modes.setdefault(g,set()).add(src['time_mode'])
            pair_matches.append(r)
            if g not in merged: merged[g]=r; sources[g]=[]
            else:
                old=merged[g]; take=np.isfinite(r['pmw_tb'])&(~np.isfinite(old['pmw_tb'])|(r['distance_km']<old['distance_km']))
                for k in r: old[k][take]=r[k][take]
            if plot: sources[g].append(src)
        if plot and any(np.isfinite(r['pmw_tb']).any() for r in pair_matches):
            plot_pair(dpr_file,file,dpr,swaths,pair_matches,out,delta_t,track_distance_km,show)
    if not merged: return None
    groups=sorted(merged); names=[f'{g}_channel_{i+1}' for g in groups for i in range(merged[g]['pmw_tb'].shape[-1])]
    if pmw_files and '.SSMIS.' in str(pmw_files[0]).upper() and groups==['S1','S2','S3','S4'] and [merged[g]['pmw_tb'].shape[-1] for g in groups]==[3,2,4,2]:
        names=['19.35V','19.35H','22.235V','37V','37H','150H','183.31±1H','183.31±3H','183.31±6.6H','91.665V','91.665H']
    payload={k:np.concatenate([merged[g][k] for g in groups],axis=-1) for k in next(iter(merged.values()))}
    payload.update(latitude=dpr['lat'],longitude=dpr['lon'],dpr_time=dpr['time'],dpr_precipitation=dpr['data'],channel_names=np.array(names),coincidence_mask=np.isfinite(payload['pmw_tb']))
    out=Path(output_file); out.parent.mkdir(parents=True,exist_ok=True); np.savez_compressed(out,**payload)
    out.with_suffix('.json').write_text(json.dumps(dict(dpr_file=str(dpr_file),pmw_files=list(map(str,pmw_files)),dpr_time_mode=dpr['time_mode'],pmw_time_modes={g:sorted(v) for g,v in modes.items()},delta_t_minutes=delta_t,max_distance_km=max_distance_km),indent=2))
    if plot and payload['coincidence_mask'].any(): plot_tb(dpr,merged,sources,out.with_suffix('.png'),n_plot_channels)
    if plot and payload['coincidence_mask'].any():
        plot_grid(payload,out,delta_t,max_distance_km,show)
        if show:
            from IPython.display import display, Image
            display(Image(filename=str(out.with_suffix('.png'))))
    return str(out)

def plot_tb(dpr,merged,sources,path,n=3):
    import matplotlib.pyplot as plt
    try: import cartopy.crs as ccrs
    except ImportError: ccrs=None
    choices=[(g,c) for g in sorted(merged) for c in range(merged[g]['pmw_tb'].shape[-1]) if np.isfinite(merged[g]['pmw_tb'][...,c]).any()]
    choices=[choices[i] for i in np.linspace(0,len(choices)-1,min(n,len(choices)),dtype=int)]
    fig,axes=plt.subplots(len(choices),2,figsize=(13,4*len(choices)),squeeze=False,layout='constrained',subplot_kw={'projection':ccrs.PlateCarree()} if ccrs else {})
    for row,(g,c) in enumerate(choices):
        mapped=merged[g]['pmw_tb'][...,c]; mask=np.isfinite(mapped); lon=dpr['lon'][mask]; lat=dpr['lat'][mask]; center=np.rad2deg(np.angle(np.mean(np.exp(1j*np.deg2rad(lon)))))
        wrap=lambda x:(x-center+180)%360-180+center
        x=wrap(lon); bounds=[max(-180,x.min()-3),min(180,x.max()+3),max(-90,lat.min()-3),min(90,lat.max()+3)]
        if bounds[0]>=bounds[1] or x.min()<-180 or x.max()>180: bounds=[-180,180,bounds[2],bounds[3]]
        for col,items in enumerate(([(s,s['data'][...,c]) for s in sources[g]],[(dpr,mapped)])):
            ax=axes[row,col]
            for s,v in items:
                ok=valid_geo(s)&np.isfinite(v); xx=(s['lon'][ok]+180)%360-180
                m=ax.scatter(xx,s['lat'][ok],c=v[ok],s=2,cmap='turbo',vmin=100,vmax=320,rasterized=True,**({'transform':ccrs.PlateCarree()} if ccrs else {}))
            if ccrs: ax.coastlines(linewidth=.5); ax.set_extent(bounds,crs=ccrs.PlateCarree()); ax.gridlines(draw_labels=True,linewidth=.3,alpha=.4)
            else: ax.set(xlim=bounds[:2],ylim=bounds[2:],xlabel='Longitude',ylabel='Latitude'); ax.grid(alpha=.2)
            ax.set_title(f'({chr(97+row*2+col)}) {g} channel {c+1}: '+('Original PMW' if col==0 else 'Mapped to DPR'))
        fig.colorbar(m,ax=axes[row,:].tolist(),orientation='horizontal',shrink=.7,label='Brightness temperature (K)')
    fig.savefig(path,dpi=200); plt.close(fig)

def automatic(sensor=None,satellite=None,start_date=None,end_date=None,delta_t=None,max_distance_km=15,output_dir='./PMWmatch_data',email=None,plot=True,n_plot_channels=3,base=BASE,version=None,show=True,download=True,dpr_orbits=None,track_distance_km=15):
    email=email or os.getenv('PPS_USERNAME') or input('PPS registered email: ').strip()
    sensor=(sensor or input('Radiometer (e.g. SSMIS, GMI, AMSR2, ATMS, MHS): ')).upper().strip()
    satellite=(satellite or input('Satellite filename token (e.g. F16, GPM, GCOMW1): ')).upper().strip()
    start=datetime.fromisoformat(start_date or input('Start date YYYY-MM-DD: ')).replace(tzinfo=timezone.utc)
    end=datetime.fromisoformat(end_date or input('End date YYYY-MM-DD (inclusive): ')).replace(tzinfo=timezone.utc)+timedelta(days=1)-timedelta(microseconds=1)
    delta_t=float(delta_t if delta_t is not None else input('Coincidence tolerance ± minutes [7.5]: ') or 7.5)
    if end<start or not np.isfinite(delta_t) or delta_t<0: raise ValueError('Invalid dates or tolerance')
    root=Path(output_dir); pad=timedelta(minutes=delta_t); password=os.getenv('PPS_PASSWORD') or email
    def pmw_pred(n):
        tokens=n.upper().split('.'); return tokens[0] in ('1C','1C-R') and sensor in tokens and satellite in tokens and (not version or version.upper() in tokens)
    with requests.Session() as s:
        s.auth=(email,password)
        dpr=discover(s,start-timedelta(days=1),end,'radar',lambda n:n.startswith('2A.GPM.DPR.') and 'SLH' not in n.upper() and (not version or version.upper() in n.upper()),base)
        dpr=[u for u in dpr if interval(u)[1]>=start and interval(u)[0]<=end]
        pmw=discover(s,start-pad-timedelta(days=1),end+pad+timedelta(days=1),'1C',pmw_pred,base)
    if not dpr or not pmw: raise RuntimeError(f'No matching archive products: DPR={len(dpr)}, PMW={len(pmw)}. Check dates, filename tokens and PPS access.')
    import csv
    root.mkdir(parents=True,exist_ok=True)
    for label,urls in [('DPR',dpr),('PMW',pmw)]:
        print(f'Available {label} orbits ({len(urls)}): {[orbit(u) for u in urls]}')
        with (root/f'{label}_available_orbits.csv').open('w',newline='') as f:
            w=csv.writer(f); w.writerow(['orbit','start_utc','end_utc','filename','url'])
            for u in urls: w.writerow([orbit(u),*interval(u),Path(u).name,u])
    pairs=[(u,v) for u in dpr for v in pmw if interval(v)[0]<=interval(u)[1]+pad and interval(v)[1]>=interval(u)[0]-pad]
    with (root/'temporal_pairs.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['DPR_orbit','PMW_orbit','DPR_url','PMW_url'])
        for u,v in pairs: w.writerow([orbit(u),orbit(v),u,v]); print(f'[Candidate] DPR {orbit(u)} / PMW {orbit(v)}')
    print(f'Temporal candidate pairs: {len(pairs)}')
    if not download: return []
    results=[]
    for u in dpr:
        if dpr_orbits is not None and orbit(u) not in dpr_orbits: continue
        a,b=interval(u); candidates=[v for v in pmw if interval(v)[0]<=b+pad and interval(v)[1]>=a-pad]
        if not candidates: continue
        df=download_urls([u],root/'DPR',email,password)[0]; pf=download_urls(candidates,root/'PMW',email,password)
        print(f'{Path(df).name}: {len(pf)} candidate granules',flush=True)
        results.append(process(df,pf,root/'matches'/f'{Path(df).stem}_{satellite}_{sensor}.npz',delta_t,max_distance_km,plot,n_plot_channels,show,track_distance_km))
    return results
