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

def plot_pair(df, pf, dpr, swaths, matches, out, minutes, km, show):

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cf

    from pathlib import Path
    from datetime import datetime, timezone

    # ========================================================
    # IDENTIFY MICROWAVE SENSOR
    # ========================================================

    sensor = Path(pf).name.split(".")[2].upper()

    # ========================================================
    # EXTRACT SATELLITE TRACKS AND COINCIDENCE EVENTS
    # ========================================================

    g = next(iter(swaths))

    a = track(df, "FS", dpr)
    b = track(pf, g, swaths[g])

    def as_swath(tr):
        return {k: np.asarray(v)[:, None] for k, v in zip(('lat', 'lon', 'time'), tr[:3])}

    annotations = find_intersections(as_swath(b), as_swath(a),
                                     coincidence_minutes=np.inf, search_radius_km=km)
    for c in annotations:
        c['valid'] = c['delta_minutes'] <= minutes
        print(f"[CROSSING] lat={c['latitude']:.4f}, lon={c['longitude']:.4f}, "
              f"delta_t={c['delta_minutes']:.2f} min, valid={c['valid']}")
    print(f"[TRACK] DPR / {sensor}: {len(annotations)} geometric intersections; "
          f"{sum(c['valid'] for c in annotations)} within time tolerance")

    # ========================================================
    # INITIALIZE FIGURE
    # ========================================================

    fig = plt.figure(figsize=(10, 7))

    ax = fig.add_subplot(
        1, 1, 1,
        projection=ccrs.PlateCarree()
    )

    ax.set_extent(
        [-180, 180, -90, 90],
        crs=ccrs.PlateCarree()
    )

    ax.stock_img()

    ax.add_feature(
        cf.COASTLINE,
        linewidth=0.7
    )

    ax.add_feature(
        cf.BORDERS,
        linestyle=":",
        linewidth=0.5
    )

    # ========================================================
    # SATELLITE TRACKS
    # ========================================================

    dpr_line = ax.plot(
        a[1],
        a[0],
        "r-",
        transform=ccrs.Geodetic(),
        linewidth=2,
        zorder=2,
        label="DPR track"
    )

    pmw_line = ax.plot(
        b[1],
        b[0],
        "b-",
        transform=ccrs.Geodetic(),
        linewidth=2,
        zorder=2,
        label=f"{sensor} track"
    )

    # ========================================================
    # COINCIDENCE ANNOTATION LOCATIONS
    # ========================================================

    annotation_locations = [
        (c["longitude"], c["latitude"])
        for c in annotations
    ]

    # ========================================================
    # AUTOMATIC LEGEND PLACEMENT
    # ========================================================

    legend_corners = {
        "upper_left": {
            "pos": (0.02, 0.96),
            "anchor": "upper left",
            "penalty": 0
        },
        "upper_right": {
            "pos": (0.98, 0.96),
            "anchor": "upper right",
            "penalty": 0
        },
        "lower_left": {
            "pos": (0.02, 0.03),
            "anchor": "lower left",
            "penalty": 10
        },
        "lower_right": {
            "pos": (0.98, 0.03),
            "anchor": "lower right",
            "penalty": 10
        }
    }

    for corner in legend_corners.values():

        x_norm, y_norm = corner["pos"]

        count = 0

        for lon, lat in annotation_locations:

            ax_pos = ax.projection.transform_point(
                lon,
                lat,
                ccrs.PlateCarree()
            )

            ax_x = (
                (ax_pos[0] - ax.get_xlim()[0])
                / (ax.get_xlim()[1] - ax.get_xlim()[0])
            )

            ax_y = (
                (ax_pos[1] - ax.get_ylim()[0])
                / (ax.get_ylim()[1] - ax.get_ylim()[0])
            )

            if (
                abs(ax_x - x_norm) < 0.3
                and abs(ax_y - y_norm) < 0.3
            ):
                count += 1

        corner["score"] = (
            count * 10 + corner["penalty"]
        )

    best_corner = min(
        legend_corners.values(),
        key=lambda x: x["score"]
    )

    legend = ax.legend(
        [dpr_line[0], pmw_line[0]],
        ["DPR track", f"{sensor} track"],
        loc=best_corner["anchor"],
        bbox_to_anchor=best_corner["pos"],
        frameon=True,
        framealpha=0.7,
        facecolor="white",
        edgecolor="black",
        fontsize=9
    )

    for text in legend.get_texts():
        text.set_fontweight("bold")

    # ========================================================
    # GRIDLINES AND COORDINATE LABELS
    # ========================================================

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color="gray",
        alpha=0.6,
        linestyle="--"
    )

    gl.top_labels = False
    gl.right_labels = False

    gl.xlabel_style = {
        "size": 10,
        "weight": "bold"
    }

    gl.ylabel_style = {
        "size": 10,
        "weight": "bold"
    }

    # Explicit tick intervals matching the original figure.

    gl.xlocator = plt.FixedLocator(
        [-180, -120, -60, 0, 60, 120, 180]
    )

    gl.ylocator = plt.FixedLocator(
        [-60, -30, 0, 30, 60]
    )

    # ========================================================
    # FIGURE TITLE
    # ========================================================

    title = (
        f"DPR: {Path(df).name}\n"
        f"{sensor}: {Path(pf).name}"
    )

    ax.set_title(
        title,
        fontsize=12,
        fontweight="bold",
        pad=10
    )

    # Finalize map geometry before annotation positioning.

    fig.canvas.draw()

    # ========================================================
    # COINCIDENCE ANNOTATIONS
    # ========================================================

    prev_annotations = []

    collision_threshold = 1.5

    def fmt(t):
        return datetime.fromtimestamp(
            float(t),
            timezone.utc
        ).strftime("%H:%M:%S")

    for c in annotations:
        lat = c['latitude']
        lon = c['longitude']
        dpr_time = fmt(c['dpr_time'])
        pmw_time = fmt(c['pmw_time'])
        delta_t = c['delta_minutes']
        prefix = "" if c['valid'] else "Outside time tolerance\n"
        annotation_text = (
            f"{prefix}ΔT: {delta_t:.1f} min\n"
            f"DPR: {dpr_time}\n"
            f"{sensor}: {pmw_time}"
        )

        # Default annotation position.

        x_offset = 0
        y_offset = 12

        ha = "center"
        va = "bottom"

        bbox_props = dict(
            boxstyle="round,pad=0.3",
            fc="white",
            ec="black",
            lw=0.5,
            alpha=0.7
        )

        # Longitude boundary handling.

        if lon > 150:

            x_offset = -15
            ha = "right"

        elif lon < -150:

            x_offset = 15
            ha = "left"

        # Latitude boundary handling.

        if lat > 60:

            y_offset = -15
            va = "top"

        elif lat < -60:

            y_offset = 15
            va = "bottom"

        # Avoid overlapping annotations.

        for prev in prev_annotations:

            if (
                abs(lon - prev["lon"])
                < collision_threshold
                and abs(lat - prev["lat"])
                < collision_threshold
            ):

                y_offset = -prev["y_offset"]

                va = (
                    "top"
                    if va == "bottom"
                    else "bottom"
                )

                break

        # Convert geographic position to display coordinates.

        transform = (
            ccrs.PlateCarree()
            ._as_mpl_transform(ax)
        )

        x_ax, y_ax = transform.transform(
            (lon, lat)
        )

        # Convert annotation offset from points to pixels.

        points_to_pixels = fig.dpi / 72.0

        final_y = (
            y_ax + y_offset * points_to_pixels
        )

        final_y_norm = (
            (final_y - ax.bbox.y0)
            / ax.bbox.height
        )

        if final_y_norm > 0.96:

            y_offset -= (
                (final_y_norm - 0.96)
                * ax.bbox.height
                / points_to_pixels
            )

        elif final_y_norm < 0.03:

            y_offset += (
                (0.03 - final_y_norm)
                * ax.bbox.height
                / points_to_pixels
            )

        prev_annotations.append({
            "lon": lon,
            "lat": lat,
            "y_offset": y_offset
        })

        # Plot annotation.

        ax.annotate(
            annotation_text,
            xy=(lon, lat),
            xycoords=transform,
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=9,
            fontweight="bold",
            bbox=bbox_props,
            arrowprops=dict(
                arrowstyle="->",
                color="black",
                linewidth=0.7
            ),
            zorder=9,
            annotation_clip=True
        )

        # Black coincidence marker.

        ax.plot(
            lon,
            lat,
            "ko" if c["valid"] else "kx",
            markersize=8,
            transform=ccrs.Geodetic(),
            zorder=9
        )

    # ========================================================
    # SAVE AND DISPLAY FIGURE
    # ========================================================

    filename = out.with_name(
        out.stem
        + "_"
        + Path(pf).stem
        + "_tracks.png"
    )

    fig.savefig(
        filename,
        bbox_inches="tight",
        dpi=400
    )

    if show:

        from IPython.display import display, Image

        display(
            Image(filename=str(filename))
        )

    plt.close(fig)

    print(f"[FIGURE SAVED] {filename}")

def plot_grid(p,out,minutes,km,show):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    fig,axes=plt.subplots(2,2,figsize=(14,8),subplot_kw={'projection':ccrs.PlateCarree()},layout='constrained')
    mask=p['coincidence_mask'].any(-1); lat=p['latitude']; lon=p['longitude']; ok=np.isfinite(lat)&np.isfinite(lon)&(abs(lat)<=90)
    for ax,v,title,cmap,lim in zip(axes.flat,[p['dpr_precipitation'],mask.astype(float),np.max(np.where(np.isfinite(p['delta_t_minutes']),p['delta_t_minutes'],-np.inf),axis=-1),np.max(np.where(np.isfinite(p['distance_km']),p['distance_km'],-np.inf),axis=-1)],['DPR precipitation (mm/h)','Coincidence mask','Maximum time difference (min)','Maximum distance (km)'],['turbo','Greens','viridis','plasma'],[10,1,minutes,km]):
        valid=ok&np.isfinite(v); ax.set_global(); ax.coastlines(); im=ax.scatter(lon[valid],lat[valid],c=v[valid],s=2,cmap=cmap,vmin=0,vmax=lim,transform=ccrs.PlateCarree()); ax.set_title(title); fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.7)
    _finish(fig,out.with_name(out.stem+'_grid.png'),show)
