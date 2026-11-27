"""Publication-friendly quicklook figures; maps use lon/lat axes without cartopy."""
from pathlib import Path
import numpy as np
from .mesh import swath_mesh, is_gmi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_collocation(npz_file,output_dir=None,channel=0,step=1):
    if is_gmi(npz_file): return None
    z=np.load(npz_file)
    out=Path(output_dir or Path(npz_file).with_suffix(''))
    out.mkdir(parents=True,exist_ok=True)
    lon=z['longitude'][::step,::step];lat=z['latitude'][::step,::step]
    refkey='dpr_precipitation'
    tb=z['pmw_tb'][::step,::step,channel]
    ref=z[refkey][::step,::step]
    mask=(z['valid_mask'][::step,::step] if 'valid_mask' in z else z['coincidence_mask'][::step,::step,channel])
    dist=z['distance_km'][::step,::step];dt=z['delta_t_minutes'][::step,::step]
    if dist.ndim==3: dist=dist[...,channel]; dt=dt[...,channel]
    fields=[(ref,refkey.replace('_',' ').title(),'viridis'),(tb,f'PMW TB channel {channel}','turbo'),
            (mask.astype(float),'Valid collocation mask','gray_r'),(dist,'Nearest PMW distance (km)','magma'),
            (dt,'Absolute time difference (min)','plasma')]
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for ax,(a,title,cmap) in zip(axes.flat,fields):
        valid=np.isfinite(lon)&np.isfinite(lat)&np.isfinite(a)
        if valid.any():
            m=swath_mesh(ax,(lon+180)%360-180,lat,a,valid,cmap=cmap)
            fig.colorbar(m,ax=ax,shrink=.8)
        ax.set_title(title);ax.set_xlabel('Longitude');ax.set_ylabel('Latitude')
        ax.grid(alpha=.2)
    axes.flat[-1].axis('off')
    axes.flat[-1].text(.05,.85,f"Matched pixels: {int(mask.sum()):,} / {mask.size:,}\nFraction: {mask.mean():.2%}\nTB channel: {channel}",transform=axes.flat[-1].transAxes,va='top',fontsize=12)
    path=out/'collocation_overview.png';fig.savefig(path,dpi=180);plt.close(fig)
    return str(path)

def plot_orbit_tracks(pmw,dpr,intersections,output_path):
    fig,ax=plt.subplots(figsize=(13,6))
    for swath,name,color in [(pmw,'PMW','tab:blue'),(dpr,'DPR','tab:red')]:
        mid=swath['lat'].shape[1]//2
        lon=swath['lon'][:,mid].astype(float).copy();lat=swath['lat'][:,mid]
        lon[1:][np.abs(np.diff(lon))>180]=np.nan
        ax.plot(lon,lat,color=color,label=name,linewidth=1)
    for c in intersections:ax.plot(c['longitude'],c['latitude'],marker='x',color='black',markersize=7,linestyle='none')
    ax.set(xlabel='Longitude',ylabel='Latitude',title='Orbit centerlines and accepted crossings',xlim=(-180,180),ylim=(-90,90))
    ax.grid(alpha=.3);ax.legend();fig.savefig(output_path,dpi=180,bbox_inches='tight');plt.close(fig)
