# -*- coding: utf-8 -*-
"""
Step 2 — Watershed Map, Any Basin Worldwide (open data)
=======================================================
Builds a three-panel watershed figure in the classic study-area layout:

  LEFT          Main map — basin boundary, drainage network in three classes
                (main river / tributaries / river branches), provincial
                boundaries with names, title, compass rose, scale bar and a
                decimal-degree graticule on the top and left edges.
  RIGHT, TOP    Location Map — host country among its neighbours, internal
                HydroBASINS (or province) lines, sea names, basin in grey
                with a red outline.
                Elevation Map — hillshaded DEM of the basin.
  RIGHT, BOTTOM Legend — map symbols, classed elevation swatches and the
                coordinate system.

Basin sources (BASIN_SOURCE)
    'hydrobasins'  HydroBASINS polygon containing the outlet point plus all
                   polygons upstream of it (NEXT_DOWN topology), or the
                   whole river system when HYBAS_WHOLE_BASIN = True
                   (MAIN_BAS). Best for medium–large basins.
    'dem'          D8 delineation from the Copernicus DEM (pysheds) with the
                   outlet snapped to a channel; EXPECTED_AREA_KM2 helps pick
                   the right river. Best for small basins / gauges.
    'file'         Your own boundary (.shp / .gpkg / .geojson).

Drainage network
    Always derived from the Copernicus DEM (D8, pysheds). The main river is
    traced upstream from the outlet along the largest contributing area;
    reaches of Strahler order ≥ TRIB_MIN_ORDER are tributaries, the rest
    are river branches. Grid resolution is chosen automatically: GLO-30
    for small windows, GLO-90 when the 30 m grid would be too large.

Open data sources
    DEM          Copernicus GLO-30 / GLO-90 (ESA, © DLR/Airbus), AWS Open Data
    Basins       HydroBASINS v1c (Lehner & Grill 2013), hydrosheds.org
    Countries    Natural Earth 1:50m countries + marine areas (public domain)
    Boundaries   geoBoundaries gbOpen (CC BY 4.0)

Inputs:
  none (everything is downloaded); optional WATERSHED_FILE

Outputs:
  outputs/<name>/<name>_watershed_map.png      ← 600 dpi raster
  outputs/<name>/<name>_watershed_map.pdf      ← vector for the journal
  outputs/<name>/<name>_watershed.gpkg         ← basin, streams, main river, outlet
  outputs/<name>/<name>_dem.tif                ← DEM clipped to the basin
  outputs/<name>/<name>_morphometry.csv

Notes
  Plotting is in WGS 84 lon/lat with a true-scale aspect at the map centre,
  so basins above ~70° latitude or crossing the antimeridian render less
  well. Basins > ~100 000 km² need a lot of RAM even at 90 m.

Dependencies: geopandas, rasterio, pysheds, matplotlib, requests
  conda create -n hydrology -c conda-forge python=3.11 geopandas rasterio \
        pysheds requests
"""

# ── Standard library ────────────────────────────────────────────────────────
import math                     # trigonometry for latitude corrections, log10 for scale bar
import os                       # file / folder handling
import sys                      # console encoding
import textwrap                 # wrap long legend and sea labels
import time                     # wait between download retries
import warnings                 # silence noisy library warnings

# ── Plotting ────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')           # non-interactive backend: render straight to file
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe                       # white halo behind labels
from matplotlib.colors import LinearSegmentedColormap, LightSource, Normalize
from matplotlib.patches import Ellipse, Polygon as MplPolygon, Rectangle

# ── Data / GIS ──────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import geopandas as gpd         # vector data (basins, countries, rivers)
import requests                 # HTTP downloads
import rasterio                 # raster I/O (DEM)
from rasterio.features import geometry_mask               # polygon → boolean raster mask
from rasterio.mask import mask as rio_mask                # clip raster to a polygon
from rasterio.merge import merge as rio_merge             # mosaic DEM tiles
from shapely.geometry import Point, LineString, MultiPolygon, Polygon, shape, box
from shapely.ops import polylabel, unary_union            # label placement, dissolve

from pysheds.grid import Grid   # D8 hydrology: fill sinks, flow direction, accumulation

warnings.filterwarnings('ignore')
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')   # ✓ ⚠ ✗ on Windows consoles

# ---------------------------------------------------------------------------
# Configuration — the only block to edit for a new basin
# ---------------------------------------------------------------------------
WATERSHED_NAME    = 'Ba River Basin'
TITLE             = None          # None → "Watershed Map of <name>, <country>"
BASIN_SOURCE      = 'hydrobasins' # 'hydrobasins' | 'dem' | 'file'
OUTLET_LAT        = 13.05         # any point on the river inside the basin (WGS 84)
OUTLET_LON        = 108.98

# 'hydrobasins' options
HYBAS_LEVEL       = 6             # 1–12; higher = finer sub-basins (used when not whole basin)
HYBAS_WHOLE_BASIN = True          # True = entire river system to the sea (MAIN_BAS)

# 'dem' options
EXPECTED_AREA_KM2 = None          # known area → robust outlet snapping; None = nearest strong channel
SNAP_RADIUS_M     = 1000          # outlet search radius
SNAP_MIN_KM2      = 0.5           # smallest upstream area accepted as a channel when snapping
SEARCH_BUFFER_DEG = None          # DEM window half-width; None = auto from EXPECTED_AREA_KM2

# 'file' option
WATERSHED_FILE    = None          # e.g. 'shp/my_basin.shp'

# DEM and drainage network
DEM_RES           = 'auto'        # 'auto' | 30 | 90 (metres)
DEM_MAX_CELLS_30M = 30e6          # above this many 30 m cells, 'auto' switches to 90 m
STREAM_DENSITY    = 600           # channel threshold = basin area / this (higher = denser network)
TRIB_MIN_ORDER    = None          # Strahler order counted as tributary; None = max order − 2
VERT_EXAG         = 2.0           # hillshade vertical exaggeration (Elevation Map)

# Context layers
COUNTRY_ISO3         = None       # e.g. 'VNM'; None = auto-detect
ADMIN_LEVEL          = 'ADM1'     # provincial boundaries on the main map ('ADM1', 'ADM2', None)
ADMIN_LABEL          = 'Provincial Boundary'
LOCATION_HYBAS_LEVEL = 6          # HydroBASINS lines inside the country on the Location Map; None = ADM1
LOCATION_PAD         = 0.08       # padding around the host country on the Location Map

# Output paths (derived from the basin name)
SLUG       = WATERSHED_NAME.lower().replace(' ', '_')   # 'Ba River Basin' → 'ba_river_basin'
OUTPUT_DIR = f'Study_area/outputs/{SLUG}'                          # results for this basin
CACHE_DIR  = 'outputs/_cache'                           # downloads shared by all basins
FIG_STEM   = f'{OUTPUT_DIR}/{SLUG}_watershed_map'       # figure file name without extension

# Typography — serif, as in the reference layout.
FONT_FAMILY   = ['Times New Roman', 'Times', 'Liberation Serif', 'DejaVu Serif']
FONT_SCALE    = 1.0               # multiply every font size (e.g. 1.2 = 20 % larger)
FIG_WIDTH_IN  = 7.2               # ≈183 mm double-column width
MAIN_FRAC     = 0.56              # share of the usable width given to the main map
TICK_DECIMALS = 2                 # 108.00°E style graticule labels (2 decimal places)
DPI           = 600               # PNG resolution

# Font sizes (pt, before FONT_SCALE) — one hierarchy for the whole figure.
FS_TITLE      = 10.5              # all three map titles + legend heading
FS_BODY       = 7.5               # legend entries, main-map place names
FS_TICK       = 7.0               # graticule labels, compass letters
FS_SCALE      = 6.5               # scale-bar numbers (all maps)
FS_SMALL      = 6.0               # country / sea names on the Location Map

# Palette
C_INK         = '#000000'         # frames, basin outline, text
C_RIVER_MAIN  = '#0B1F8A'         # main river — navy
C_RIVER_TRIB  = '#3D7DD8'         # tributaries — mid blue
C_RIVER_BR    = '#A9CDEB'         # river branches — light blue
C_SEA         = '#BEE3F5'         # sea on the Location Map
C_OTHER_LAND  = '#DCDCDC'         # neighbouring countries
C_OTHER_EDGE  = '#8C8C8C'         # neighbouring country borders
C_HOST        = '#FFF4C2'         # host country fill
C_HOST_LINES  = '#9C9277'         # HydroBASINS / province lines inside the host country
C_BASIN_FILL  = '#9A9A9A'         # basin fill on the Location Map
C_BASIN_EDGE  = '#E31A1C'         # basin outline on the Location Map (red)
C_SEA_LABEL   = '#1E4FD0'         # sea / gulf names
ELEV_COLORS   = ['#2E6B12', '#8DBF57', '#FFF9C9', '#F2C14E', '#C98B57', '#8A5A2B']  # low → high

# Download sources
NE_BASE   = 'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson'
HYBAS_URL = 'https://data.hydrosheds.org/file/HydroBASINS/standard/hybas_{reg}_lev{lev:02d}_v1c.zip'
HEADERS   = {'User-Agent': 'watershed-map/1.0 (research script)'}

# Create output folders if they don't exist yet.
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)


def fs(size):
    """Scale a font size by FONT_SCALE."""
    return size * FONT_SCALE


# Global matplotlib style for the whole figure.
plt.rcParams.update({
    'font.family':     'serif',
    'font.serif':      FONT_FAMILY,
    'font.size':       fs(7.5),
    'mathtext.fontset': 'stix',
    'axes.linewidth':  0.8,
    'axes.edgecolor':  C_INK,
    'pdf.fonttype':    42,        # embed TrueType → editable text in Illustrator
    'svg.fonttype':    'none',
})
HALO = [pe.withStroke(linewidth=2.0, foreground='white')]   # white outline that keeps labels readable


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------
def download(url, path, max_retries=3):
    """Download url to path once; reuse the cached copy afterwards."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path                                   # already cached → skip download
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(url, timeout=300, headers=HEADERS)
            r.raise_for_status()                      # fail on HTTP 4xx / 5xx
            with open(path, 'wb') as f:
                f.write(r.content)
            return path
        except Exception:
            if attempt == max_retries:
                raise                                 # give up after the last attempt
            time.sleep(2 * attempt)                   # back off: 2 s, 4 s, …


def upper_cols(gdf):
    """Upper-case all attribute column names (Natural Earth varies in case)."""
    gdf.columns = [c if c == 'geometry' else c.upper() for c in gdf.columns]
    return gdf


def get_countries():
    """Natural Earth 1:50m countries with 'Name' and 'ISO3' (EPSG:4326)."""
    path = download(f'{NE_BASE}/ne_50m_admin_0_countries.geojson',
                    f'{CACHE_DIR}/ne_50m_admin_0_countries.geojson')
    gdf = upper_cols(gpd.read_file(path).to_crs(4326))

    def pick_iso(row):
        # Some countries have '-99' in ISO_A3 → fall back to the other code columns.
        for col in ('ISO_A3', 'ISO_A3_EH', 'ADM0_A3'):
            val = str(row.get(col, ''))
            if len(val) == 3 and val != '-99':
                return val
        return None

    gdf['ISO3'] = gdf.apply(pick_iso, axis=1)
    gdf['Name'] = gdf['NAME' if 'NAME' in gdf.columns else 'ADMIN'].astype(str)
    return gdf[['Name', 'ISO3', 'geometry']]


def get_marine():
    """Natural Earth marine area names (seas, gulfs, bays); None if unavailable."""
    try:
        path = download(f'{NE_BASE}/ne_50m_geography_marine_polys.geojson',
                        f'{CACHE_DIR}/ne_50m_geography_marine_polys.geojson')
        gdf = upper_cols(gpd.read_file(path).to_crs(4326))
        gdf['Name'] = gdf['NAME'].astype(str)
        return gdf[gdf['Name'].str.lower() != 'none'][['Name', 'geometry']]   # drop unnamed areas
    except Exception as exc:
        print(f"  ⚠ Marine names unavailable ({exc})")
        return None


def get_geoboundaries(iso3, level, simplified=True):
    """
    Fetch administrative boundaries from geoBoundaries (gbOpen).

    Parameters
    ----------
    iso3 : str
        ISO 3166-1 alpha-3 code, e.g. 'VNM'.
    level : str
        'ADM0', 'ADM1', 'ADM2', …
    simplified : bool
        Use the simplified geometry (much smaller download).

    Returns
    -------
    gpd.GeoDataFrame (EPSG:4326, 'Name' column) or None if not published.
    """
    tag = 'simplified' if simplified else 'full'
    path = f'{CACHE_DIR}/{iso3}_{level}_{tag}.geojson'
    try:
        if not os.path.exists(path):
            # The API returns metadata JSON that holds the actual GeoJSON download link.
            meta = requests.get(f'https://www.geoboundaries.org/api/current/gbOpen/{iso3}/{level}/',
                                timeout=60, headers=HEADERS).json()
            key = 'simplifiedGeometryGeoJSON' if simplified else 'gjDownloadURL'
            download(meta[key], path)
        gdf = gpd.read_file(path).to_crs(4326)
        gdf['Name'] = gdf['shapeName'].astype(str)
        return gdf[['Name', 'geometry']]
    except Exception as exc:
        print(f"  ⚠ geoBoundaries {iso3} {level} unavailable ({exc})")
        return None


def guess_hybas_regions(lon, lat):
    """HydroBASINS region codes ordered by likelihood for a point."""
    # HydroBASINS is split into continental files; try the most likely one first.
    if lon < -25:
        first = ['sa'] if lat < 15 else ['na', 'ar', 'gr']     # Americas / Arctic / Greenland
    elif lon < 60 and lat < 37 and not (lon > 34 and lat > 12):
        first = ['af', 'eu']                                   # Africa
    elif lon < 60:
        first = ['eu', 'af', 'si']                             # Europe / Middle East
    elif lat > 55:
        first = ['si', 'as']                                   # Siberia
    elif (lon > 110 and lat < -10) or (lon > 140 and lat < 0):
        first = ['au', 'as']                                   # Australia / Oceania
    else:
        first = ['as', 'au', 'eu']                             # Asia
    rest = [r for r in ['as', 'eu', 'af', 'na', 'sa', 'au', 'si', 'ar', 'gr'] if r not in first]
    return first + rest


def get_hydrobasins(level, lon, lat):
    """
    Load the HydroBASINS region (standard v1c, given level) that contains
    the point and return (all polygons of that region, the hit polygon).
    A point just offshore (river mouth) takes the nearest polygon < 0.05°.
    """
    pt = Point(lon, lat)
    for reg in guess_hybas_regions(lon, lat):
        path = f'{CACHE_DIR}/hybas_{reg}_lev{level:02d}_v1c.zip'
        try:
            download(HYBAS_URL.format(reg=reg, lev=level), path)
            hb = gpd.read_file(path)                  # geopandas reads the zipped shapefile directly
        except Exception as exc:
            print(f"  ⚠ HydroBASINS '{reg}' level {level} unavailable ({exc})")
            continue
        hb = hb.set_crs(4326) if hb.crs is None else hb.to_crs(4326)
        hit = hb[hb.contains(pt)]                     # polygon containing the outlet point
        if hit.empty:
            d = hb.distance(pt)                       # point slightly offshore → nearest polygon
            if d.min() < 0.05:
                hit = hb.loc[[d.idxmin()]]
        if len(hit):
            print(f"✓ HydroBASINS region '{reg}', level {level} ({len(hb):,} polygons)")
            return hb, hit.iloc[0]
    raise ValueError(f"No HydroBASINS polygon found at {lat}, {lon} — check the coordinate")


def choose_dem_res(bounds):
    """30 m unless the window would exceed DEM_MAX_CELLS_30M cells."""
    if DEM_RES in (30, 90):
        return DEM_RES                                # user forced a resolution
    w, s, e, n = bounds
    cells_30m = (e - w) * 3600 * (n - s) * 3600       # 1 arc-second ≈ 30 m → 3600 cells per degree
    return 30 if cells_30m <= DEM_MAX_CELLS_30M else 90


def get_copernicus_dem(bounds, out_path, res=30):
    """
    Mosaic Copernicus GLO-30/GLO-90 tiles covering bounds (W, S, E, N).
    Tiles are Cloud-Optimised GeoTIFFs, so only the window is transferred.
    Missing tiles (open ocean) are filled with nodata.
    """
    if os.path.exists(out_path):
        return out_path                               # cached mosaic
    bucket, code = ('copernicus-dem-30m', '10') if res == 30 else ('copernicus-dem-90m', '30')
    west, south, east, north = bounds
    sources = []
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_HTTP_MAX_RETRY='3',
                      GDAL_HTTP_RETRY_DELAY='2'):
        # One tile per 1° × 1° cell, named by its south-west corner.
        for lat in range(math.floor(south), math.ceil(north)):
            for lon in range(math.floor(west), math.ceil(east)):
                ns = f'N{lat:02d}' if lat >= 0 else f'S{-lat:02d}'
                ew = f'E{lon:03d}' if lon >= 0 else f'W{-lon:03d}'
                tile = f'Copernicus_DSM_COG_{code}_{ns}_00_{ew}_00_DEM'
                url = f'/vsicurl/https://{bucket}.s3.amazonaws.com/{tile}/{tile}.tif'
                try:
                    sources.append(rasterio.open(url))   # opened remotely, read on demand
                    print(f"  ✓ Tile {tile}")
                except Exception:
                    print(f"  ⚠ Tile {tile} not available (ocean / no data)")
        if not sources:
            raise ValueError("No Copernicus DEM tiles found for the requested window")
        # Stitch the tiles and crop to the requested window.
        mosaic, transform = rio_merge(sources, bounds=bounds, nodata=-9999)
    profile = sources[0].profile.copy()
    profile.update(driver='GTiff', height=mosaic.shape[1], width=mosaic.shape[2],
                   transform=transform, compress='deflate', tiled=True,
                   blockxsize=256, blockysize=256, nodata=-9999, dtype=mosaic.dtype)
    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(mosaic)
    for s in sources:
        s.close()
    return out_path


# ---------------------------------------------------------------------------
# Geometry / terrain helpers
# ---------------------------------------------------------------------------
# pysheds D8 codes (N, NE, E, SE, S, SW, W, NW) → (row, col) step downstream.
DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)
D8_STEP = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1),
           4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1)}


def cell_size_m(transform, lat):
    """Return (dx, dy) in metres for a geographic raster at latitude lat."""
    dx = abs(transform.a) * 111_320 * math.cos(math.radians(lat))   # longitude shrinks with latitude
    dy = abs(transform.e) * 110_574                                  # metres per degree of latitude
    return dx, dy


def slope_deg(dem, dx, dy):
    """Slope in degrees from the elevation gradient."""
    gy, gx = np.gradient(dem, dy, dx)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def largest_polygon(geom):
    """Largest part of a (multi)polygon, without interior holes."""
    geom = geom.buffer(0)                             # repair invalid geometry
    if geom.geom_type == 'MultiPolygon':
        geom = max(geom.geoms, key=lambda g: g.area)
    return Polygon(geom.exterior)                     # exterior ring only → holes removed


def trace_main_stem(fdir_arr, acc_arr, row, col, min_acc):
    """
    Walk upstream from (row, col), always stepping into the inflowing
    neighbour with the largest accumulation, until it drops below min_acc.
    Returns a list of (row, col) from the outlet to the source.
    """
    n_r, n_c = acc_arr.shape
    path = [(row, col)]
    while True:
        best, best_acc = None, min_acc
        for code, (dr, dc) in D8_STEP.items():
            r, c = row - dr, col - dc          # this neighbour flows into (row, col)
            if 0 <= r < n_r and 0 <= c < n_c and fdir_arr[r, c] == code and acc_arr[r, c] >= best_acc:
                best, best_acc = (r, c), acc_arr[r, c]
        if best is None:
            return path                        # reached the source
        row, col = best
        path.append(best)


def padded_extent(bounds, pad=0.06):
    """Extent (W, E, S, N) around shapely bounds (W, S, E, N)."""
    w, s, e, n = bounds
    dx, dy = (e - w) * pad, (n - s) * pad
    return (w - dx, e + dx, s - dy, n + dy)


def true_ratio(ext):
    """Height / width of an extent drawn at true scale."""
    w, e, s, n = ext
    return (n - s) / ((e - w) * math.cos(math.radians((s + n) / 2)))


def fit_exact(ext, ratio):
    """Grow an extent until its true-scale height/width equals ratio."""
    w, e, s, n = ext
    for _ in range(3):                                 # a few passes (cos(lat) changes as it grows)
        k = math.cos(math.radians((s + n) / 2))
        r = (n - s) / ((e - w) * k)
        if r < ratio:                                  # too wide → add height
            g = ((e - w) * k * ratio - (n - s)) / 2
            s, n = max(s - g, -89.0), min(n + g, 89.0)
        else:                                          # too tall → add width
            g = ((n - s) / (ratio * k) - (e - w)) / 2
            w, e = w - g, e + g
    return (w, e, s, n)


def set_map(ax, ext):
    """Apply limits and a true-scale aspect without resizing the frame."""
    w, e, s, n = ext
    ax.set_xlim(w, e)
    ax.set_ylim(s, n)
    ax.set_aspect(1 / math.cos(math.radians((s + n) / 2)), adjustable='datalim')


def no_ticks(ax):
    """Hide coordinate ticks (small inset maps)."""
    ax.set_xticks([])
    ax.set_yticks([])


# ---------------------------------------------------------------------------
# Map furniture
# ---------------------------------------------------------------------------
def deg_label(v, axis):
    """Format a coordinate as e.g. 108.000°E / 13.000°N."""
    hemi = ('E' if v >= 0 else 'W') if axis == 'x' else ('N' if v >= 0 else 'S')
    return f"{abs(v):.{TICK_DECIMALS}f}°{hemi}"


def nice_deg_step(span, target=3):
    """Round graticule interval giving about `target` lines across `span` degrees."""
    for st in [0.005, 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20]:
        if span / st <= target + 0.5:
            return st
    return 30


def set_graticule(ax, ext, target=3):
    """Decimal-degree labels on the top and left edges (reference style)."""
    w, e, s, n = ext
    sx, sy = nice_deg_step(e - w, target), nice_deg_step(n - s, target)
    xt = np.round(np.arange(math.ceil(w / sx) * sx, e, sx), 6)   # longitudes inside the frame
    yt = np.round(np.arange(math.ceil(s / sy) * sy, n, sy), 6)   # latitudes inside the frame
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_xticklabels([deg_label(v, 'x') for v in xt])
    ax.set_yticklabels([deg_label(v, 'y') for v in yt], rotation=90, va='center')
    # Labels on top and left only; ticks point outward; increased padding to clear titles
    ax.tick_params(direction='out', length=3, width=0.6, labelsize=fs(FS_TICK), pad=6,
                   top=True, labeltop=True, bottom=False, labelbottom=False,
                   left=True, labelleft=True, right=False, labelright=False)


def scalebar_layout(ext, ax_in, loc=(0.04, 0.035), frac=0.3, segments=4, bar_h=0.018):
    """
    Geometry of a lower-left scale bar (true length at the centre latitude).
    ax_in = (width, height) of the axes in inches, used to size the text so
    the returned 'box' is the full footprint (bar + numbers + unit) in data
    coordinates — used to keep labels and the basin clear of it.
    """
    w, e, s, n = ext
    km_per_deg = 111.320 * math.cos(math.radians((s + n) / 2))    # km per degree of longitude here
    raw = (e - w) * km_per_deg * frac / segments                  # wanted segment length (km)
    mag = 10 ** math.floor(math.log10(raw))
    seg_km = min([1, 2, 2.5, 5, 10], key=lambda m: abs(m * mag - raw)) * mag   # round to 1/2/2.5/5/10
    seg = seg_km / km_per_deg                                     # segment length in degrees
    pt_x = fs(FS_SCALE) / 72 / ax_in[0] * (e - w)     # 1 em in data-x units
    pt_y = fs(FS_SCALE) / 72 / ax_in[1] * (n - s)     # 1 em in data-y units
    x0 = w + loc[0] * (e - w) + 0.3 * pt_x            # room for half of the "0"
    y0 = s + loc[1] * (n - s)
    h = (n - s) * bar_h                               # bar thickness
    x1 = x0 + segments * seg                          # right end of the bar
    last_half = 0.3 * pt_x * len(f"{seg_km * segments:g}")
    right = max(x1 + last_half, x1 + 0.4 * pt_x + 1.3 * pt_x)   # last number vs "km"
    # Whole area used by the bar, its numbers and the unit (plus a small margin).
    footprint = box(x0 - 0.3 * pt_x - 0.4 * pt_x, y0 - 0.5 * pt_y,
                    right + 0.4 * pt_x, y0 + h + 0.35 * pt_y + 1.2 * pt_y)
    return dict(x0=x0, y0=y0, h=h, seg=seg, seg_km=seg_km, segments=segments,
                pt_x=pt_x, pt_y=pt_y, box=footprint)


def add_scalebar(ax, sb, backing=True):
    """Alternating black/white bar; numbers above, unit to the right of the bar."""
    x0, y0, h, seg = sb['x0'], sb['y0'], sb['h'], sb['seg']
    if backing:
        # Semi-transparent white box so the bar stays readable over map features.
        bx = sb['box'].bounds
        ax.add_patch(Rectangle((bx[0], bx[1]), bx[2] - bx[0], bx[3] - bx[1],
                               facecolor='white', edgecolor='none', alpha=0.85, zorder=19))
    for i in range(sb['segments']):                   # alternating black / white segments
        ax.add_patch(Rectangle((x0 + i * seg, y0), seg, h,
                               facecolor=C_INK if i % 2 == 0 else 'white',
                               edgecolor=C_INK, linewidth=0.5, zorder=20))
    for i in range(sb['segments'] + 1):               # distance numbers above each tick
        ax.text(x0 + i * seg, y0 + h + 0.35 * sb['pt_y'], f"{sb['seg_km'] * i:g}",
                ha='center', va='bottom', fontsize=fs(FS_SCALE), zorder=21)
    ax.text(x0 + sb['segments'] * seg + 0.4 * sb['pt_x'], y0 + h / 2, 'km',   # unit beside the bar
            ha='left', va='center', fontsize=fs(FS_SCALE), zorder=21)


def clear_scalebar(ext, ratio, geom, ax_in, **kw):
    """Extend the map downward until the lower-left scale bar clears geom."""
    for _ in range(25):
        sb = scalebar_layout(ext, ax_in, **kw)
        if not sb['box'].intersects(geom):
            return ext, sb                            # no overlap with the basin → done
        w, e, s, n = ext
        ext = fit_exact((w, e, s - 0.03 * (n - s), n), ratio)   # add 3 % more space at the bottom
    return ext, scalebar_layout(ext, ax_in, **kw)


def add_compass(ax, xy=(0.88, 0.86), r=0.075, fontsize=FS_TICK):
    """Eight-point compass rose in axes coordinates (shape-corrected)."""
    pos = ax.get_position()
    fw, fh = ax.figure.get_size_inches()
    kx = (pos.height * fh) / (pos.width * fw)        # axes-x units per axes-y unit
    cx0, cy0 = xy

    def pt(dx, dy):
        # Offset from the centre, corrected so the rose stays round in a non-square frame.
        return (cx0 + dx * kx, cy0 + dy)

    # White disc behind the rose, then the inner ring.
    ax.add_patch(Ellipse((cx0, cy0), 2 * 1.45 * r * kx, 2 * 1.45 * r, transform=ax.transAxes,
                         facecolor='white', edgecolor='none', alpha=0.85, zorder=23))
    ax.add_patch(Ellipse((cx0, cy0), 2 * 0.62 * r * kx, 2 * 0.62 * r, transform=ax.transAxes,
                         fill=False, edgecolor=C_INK, linewidth=0.8, zorder=24))
    # (angle, length, half-width): four short diagonal arms + four long cardinal arms.
    arms = [(45, 0.62, 0.12), (135, 0.62, 0.12), (225, 0.62, 0.12), (315, 0.62, 0.12),
            (0, 1.0, 0.17), (90, 1.0, 0.17), (180, 1.0, 0.17), (270, 1.0, 0.17)]
    for ang, length, width in arms:
        a = math.radians(ang)
        d = (math.sin(a), math.cos(a))               # direction of the arm
        p = (math.cos(a), -math.sin(a))              # perpendicular to the arm
        tip = pt(d[0] * length * r, d[1] * length * r)
        b1 = pt(p[0] * width * r, p[1] * width * r)
        b2 = pt(-p[0] * width * r, -p[1] * width * r)
        # Each arm = one black and one white triangle (3-D look).
        for base, face in ((b1, C_INK), (b2, 'white')):
            ax.add_patch(MplPolygon([pt(0, 0), tip, base], closed=True, transform=ax.transAxes,
                                    facecolor=face, edgecolor=C_INK, linewidth=0.4, zorder=25))
    for lab, ang in (('N', 0), ('E', 90), ('S', 180), ('W', 270)):   # N / E / S / W letters
        a = math.radians(ang)
        x, y = pt(math.sin(a) * 1.25 * r, math.cos(a) * 1.25 * r)
        ax.text(x, y, lab, transform=ax.transAxes, ha='center', va='center',
                fontsize=fs(fontsize), zorder=26)


def polygons_of(geom):
    """Flatten any geometry into its Polygon parts."""
    if geom.geom_type == 'Polygon':
        return [geom]
    if hasattr(geom, 'geoms'):
        return [p for g in geom.geoms for p in polygons_of(g)]
    return []


def label_parts(ax, gdf, view, avoid, min_frac=0.03, **text_kw):
    """
    Label polygons at the most interior point (pole of inaccessibility) of
    their visible part, minus `avoid` — so names stay clear of the frame
    edge, coastlines, the basin, scale bar and compass.
    """
    for name, geom in zip(gdf['Name'], gdf.geometry):
        part = geom.intersection(view)               # only the part visible in the map
        if avoid is not None:
            part = part.difference(avoid)            # remove areas reserved for other features
        polys = polygons_of(part)
        if not polys or part.area < min_frac * view.area:
            continue                                 # too small to label
        big = max(polys, key=lambda p: p.area)
        try:
            rp = polylabel(big, tolerance=max(big.bounds[2] - big.bounds[0],
                                              big.bounds[3] - big.bounds[1]) / 200)
        except Exception:
            rp = big.representative_point()
        ax.text(rp.x, rp.y, name, **{'ha': 'center', 'va': 'center', 'zorder': 12,
                                     'clip_on': True, **text_kw})

# ============================================================================
# 1 — BASIN BOUNDARY
# ============================================================================
print("\n" + "=" * 80)
print(f"WATERSHED MAP — {WATERSHED_NAME.upper()}")
print("=" * 80)

hb_region = None                                     # full HydroBASINS region (reused later)
if BASIN_SOURCE == 'file':
    # User-supplied boundary: dissolve all features and keep the largest part.
    ws_geom = largest_polygon(unary_union(gpd.read_file(WATERSHED_FILE).to_crs(4326).geometry))
    print(f"✓ Boundary from file: {WATERSHED_FILE}")
elif BASIN_SOURCE == 'hydrobasins':
    print(f"Loading HydroBASINS level {HYBAS_LEVEL} …")
    hb_region, hit = get_hydrobasins(HYBAS_LEVEL, OUTLET_LON, OUTLET_LAT)
    if HYBAS_WHOLE_BASIN:
        # Every sub-basin of the same river system (same MAIN_BAS id).
        sel = hb_region[hb_region['MAIN_BAS'] == hit['MAIN_BAS']]
    else:
        # Only the hit polygon and everything upstream: walk the NEXT_DOWN links backwards.
        children = hb_region.groupby('NEXT_DOWN')['HYBAS_ID'].apply(list).to_dict()
        ids, stack = set(), [hit['HYBAS_ID']]
        while stack:
            h = stack.pop()
            if h not in ids:
                ids.add(h)
                stack.extend(children.get(h, []))
        sel = hb_region[hb_region['HYBAS_ID'].isin(ids)]
    ws_geom = largest_polygon(unary_union(sel.geometry))   # merge sub-basins into one boundary
    print(f"  HydroBASINS units merged:    {len(sel)}")
elif BASIN_SOURCE != 'dem':
    raise ValueError("BASIN_SOURCE must be 'hydrobasins', 'dem' or 'file'")

# ── DEM window ──────────────────────────────────────────────────────────────
if BASIN_SOURCE == 'dem':
    # Basin not known yet → download a square window around the outlet.
    if SEARCH_BUFFER_DEG:
        buf = SEARCH_BUFFER_DEG
    elif EXPECTED_AREA_KM2:
        # Elongated basins reach ~1.3 × √A from the outlet.
        buf = 1.3 * math.sqrt(EXPECTED_AREA_KM2) / (111.32 * math.cos(math.radians(OUTLET_LAT))) + 0.05
    else:
        buf = 0.30
    dem_bounds = (OUTLET_LON - buf, OUTLET_LAT - buf, OUTLET_LON + buf, OUTLET_LAT + buf)
else:
    # Basin known → its bounding box plus a 5 % margin.
    w0, s0, e0, n0 = ws_geom.bounds
    m = max(0.02, 0.05 * max(e0 - w0, n0 - s0))
    dem_bounds = (w0 - m, s0 - m, e0 + m, n0 + m)

dem_res = choose_dem_res(dem_bounds)                 # 30 m or 90 m
dem_raw = (f'{CACHE_DIR}/cop{dem_res}_' + '_'.join(f'{v:.3f}' for v in dem_bounds) + '.tif')
print(f"Fetching Copernicus GLO-{dem_res} DEM …")
get_copernicus_dem(dem_bounds, dem_raw, res=dem_res)
print(f"✓ DEM window: {dem_raw}")


# ============================================================================
# 2 — FLOW ROUTING AND DRAINAGE NETWORK
# ============================================================================
grid = Grid.from_raster(dem_raw)                     # pysheds grid in the DEM's geometry
dem = grid.read_raster(dem_raw)
transform = grid.affine                              # pixel ↔ lon/lat conversion
dx_m, dy_m = cell_size_m(transform, (dem_bounds[1] + dem_bounds[3]) / 2)
cell_km2 = dx_m * dy_m / 1e6                         # area of one cell in km²
print(f"  Grid:                        {dem.shape[0]:,} × {dem.shape[1]:,} cells, "
      f"{dx_m:.0f} × {dy_m:.0f} m")

print("Conditioning DEM and routing flow …")
# Hydrological conditioning: fill single-cell pits → fill depressions → slope the flats,
# so every cell has a downhill path to the edge.
inflated = grid.resolve_flats(grid.fill_depressions(grid.fill_pits(dem)))
fdir = grid.flowdir(inflated, dirmap=DIRMAP)         # D8 flow direction (one of 8 neighbours)
acc = grid.accumulation(fdir, dirmap=DIRMAP)         # number of upstream cells
fdir_arr = np.asarray(fdir)
acc_km2 = np.asarray(acc) * cell_km2                 # upstream area in km²
n_rows, n_cols = acc_km2.shape
print("✓ Flow direction and accumulation done")

# ── Outlet and (for 'dem') the basin itself ─────────────────────────────────
if BASIN_SOURCE == 'dem':
    # Snap the outlet to the best channel cell within SNAP_RADIUS_M.
    row0, col0 = rasterio.transform.rowcol(transform, OUTLET_LON, OUTLET_LAT)
    r_px = int(math.ceil(SNAP_RADIUS_M / min(dx_m, dy_m)))        # search radius in cells
    rr, cc = np.mgrid[row0 - r_px:row0 + r_px + 1, col0 - r_px:col0 + r_px + 1]
    valid = (rr >= 0) & (rr < n_rows) & (cc >= 0) & (cc < n_cols)
    rr_c, cc_c = np.clip(rr, 0, n_rows - 1), np.clip(cc, 0, n_cols - 1)
    dist = np.hypot((rr - row0) * dy_m, (cc - col0) * dx_m)      # distance from the given point (m)
    win_acc = acc_km2[rr_c, cc_c]
    cand = valid & (dist <= SNAP_RADIUS_M) & (win_acc >= SNAP_MIN_KM2)
    if not cand.any():
        raise ValueError(f"No channel cell (≥ {SNAP_MIN_KM2} km²) within {SNAP_RADIUS_M} m "
                         f"of the outlet — increase SNAP_RADIUS_M or check OUTLET_LAT/LON")
    if EXPECTED_AREA_KM2:
        # Known area: pick the cell whose upstream area is closest to it.
        score = np.where(cand, np.abs(np.log(win_acc / EXPECTED_AREA_KM2)), np.inf)
    else:
        # Nearest strong channel: distance penalised, larger rivers preferred nearby.
        score = np.where(cand, dist - 50 * np.log10(win_acc), np.inf)
    k = np.unravel_index(np.argmin(score), score.shape)
    out_row, out_col = int(rr_c[k]), int(cc_c[k])
    # Delineate everything draining to the snapped outlet.
    catch = grid.catchment(x=out_col, y=out_row, fdir=fdir, dirmap=DIRMAP,
                           xytype='index')   # exact snapped cell, no re-rounding
    polys = [shape(g) for g, v in grid.polygonize(catch.astype(np.uint8)) if v == 1]
    # Smooth the stair-stepped raster edge by half a cell; drop pit-fill holes.
    ws_geom = largest_polygon(max(polys, key=lambda p: p.area).simplify(abs(transform.a) * 0.5))
    wb = ws_geom.bounds
    tol = 2 * abs(transform.a)
    # Warn if the basin is cut off by the edge of the DEM window.
    if (wb[0] - dem_bounds[0] < tol or dem_bounds[2] - wb[2] < tol or
            wb[1] - dem_bounds[1] < tol or dem_bounds[3] - wb[3] < tol):
        print("  ⚠ Basin touches the DEM window edge — raise SEARCH_BUFFER_DEG")
    ws_mask = ~geometry_mask([ws_geom], out_shape=acc_km2.shape, transform=transform)
else:
    # Outlet = highest-accumulation cell inside the boundary.
    ws_mask = ~geometry_mask([ws_geom], out_shape=acc_km2.shape, transform=transform)
    out_row, out_col = np.unravel_index(np.argmax(np.where(ws_mask, acc_km2, -1)), acc_km2.shape)
    out_row, out_col = int(out_row), int(out_col)

# Basin area and perimeter, measured in metres in the local UTM zone.
out_lon, out_lat = rasterio.transform.xy(transform, out_row, out_col)
ws_gdf = gpd.GeoDataFrame({'Name': [WATERSHED_NAME]}, geometry=[ws_geom], crs=4326)
utm_crs = ws_gdf.estimate_utm_crs()
ws_utm = ws_gdf.to_crs(utm_crs)
area_km2 = ws_utm.area.iloc[0] / 1e6
perim_km = ws_utm.length.iloc[0] / 1e3
print(f"✓ Basin area:                  {area_km2:,.1f} km²")
print(f"  Outlet:                      {out_lat:.5f}, {out_lon:.5f}")

# ── Streams with Strahler order ─────────────────────────────────────────────
# The network is built directly from the D8 grid as one connected tree:
# only channel cells whose whole downstream path stays inside the basin and
# reaches the outlet are kept, so no reach is cut by the boundary. Every reach
# runs from a source/junction down to the next junction and shares that
# junction vertex with its neighbours → no gaps between line segments.
stream_thresh = max(20 * cell_km2, area_km2 / STREAM_DENSITY)   # min upstream area of a channel (km²)
cand = ws_mask & (acc_km2 >= stream_thresh)          # candidate channel cells inside the basin
cand[out_row, out_col] = True                        # the outlet is always part of the network
rows, cols = np.nonzero(cand)                        # row/col of every channel cell
n_cells = len(rows)
cell_id = np.full(acc_km2.shape, -1, dtype=np.int64)
cell_id[rows, cols] = np.arange(n_cells)             # grid position → channel-cell index
codes = fdir_arr[rows, cols]                         # flow direction of each channel cell
# Row/col of the cell each channel cell drains into.
d_r, d_c = rows.copy(), cols.copy()
for code, (dr, dc) in D8_STEP.items():
    sel = codes == code
    d_r[sel] += dr
    d_c[sel] += dc
in_grid = (d_r >= 0) & (d_r < n_rows) & (d_c >= 0) & (d_c < n_cols) & ((d_r != rows) | (d_c != cols))
down = np.full(n_cells, -1, dtype=np.int64)          # index of the downstream channel cell (−1 = none)
down[in_grid] = cell_id[d_r[in_grid], d_c[in_grid]]
i_out = int(cell_id[out_row, out_col])               # index of the outlet cell

# Downstream cells always carry more accumulation → descending order is topological.
topo = np.argsort(-acc_km2[rows, cols], kind='stable')
# Keep a cell only if its downstream neighbour is kept → connected tree ending at the outlet.
valid = np.zeros(n_cells, bool)
valid[i_out] = True
for i in topo:
    if i != i_out and down[i] >= 0 and valid[down[i]]:
        valid[i] = True

# Strahler order, upstream → downstream.
# Rule: a source has order 1; where two streams of the same (highest) order
# meet, the order goes up by 1; otherwise the highest incoming order continues.
strahler = np.zeros(n_cells, np.int32)
up_max = np.zeros(n_cells, np.int32)                 # highest order flowing into each cell
up_cnt = np.zeros(n_cells, np.int32)                 # how many inflows have that highest order
n_up = np.zeros(n_cells, np.int32)                   # number of inflowing channel cells
for i in topo[::-1]:                                 # sources first, outlet last
    if not valid[i]:
        continue
    s_ord = 1 if up_max[i] == 0 else up_max[i] + (1 if up_cnt[i] >= 2 else 0)
    strahler[i] = s_ord
    d = down[i]
    if i != i_out and d >= 0:
        # Pass this cell's order on to its downstream neighbour.
        n_up[d] += 1
        if s_ord > up_max[d]:
            up_max[d], up_cnt[d] = s_ord, 1
        elif s_ord == up_max[d]:
            up_cnt[d] += 1

# Cell-centre coordinates (lon/lat) of every channel cell.
cx = transform.c + (cols + 0.5) * transform.a
cy = transform.f + (rows + 0.5) * transform.e
tol = abs(transform.a) * 0.5                         # line smoothing tolerance: half a cell
# Cut the network into reaches: each starts at a source or junction and runs to the next junction.
# Build reach endpoints first so we can snap them precisely.
reach_data = []
for st in np.flatnonzero(valid & (n_up != 1)):        # sources and junctions
    path, j = [st], st
    while j != i_out:
        j = down[j]
        path.append(j)
        if n_up[j] != 1:                              # stop on the next junction
            break
    if len(path) < 2:
        continue
    reach_data.append((path, st, strahler[st]))

# Build reaches with precise junction snapping to ensure connectivity.
lines, orders = [], []
for path, start_idx, ord_val in reach_data:
    # Get coordinates for the full reach path
    coords = np.column_stack([cx[path], cy[path]])
    # Simplify the reach to reduce vertex count while preserving shape
    reach_line = LineString(coords).simplify(tol)
    # Extract and snap the endpoint coordinates to exact junction locations
    reach_coords = list(reach_line.coords)
    # Start point: snap to source/junction cell
    reach_coords[0] = (cx[path[0]], cy[path[0]])
    # End point: snap to downstream junction cell
    reach_coords[-1] = (cx[path[-1]], cy[path[-1]])
    # Reconstruct the reach with snapped endpoints
    snapped_line = LineString(reach_coords)
    lines.append(snapped_line)
    orders.append(int(ord_val))
streams = gpd.GeoDataFrame({'Strahler_Order': orders}, geometry=lines, crs=4326)
max_order = int(streams['Strahler_Order'].max()) if len(streams) else 1
trib_min = TRIB_MIN_ORDER or max(2, max_order - 2)   # lowest order drawn as a tributary
streams['River_Class'] = np.where(streams['Strahler_Order'] >= trib_min,
                                  'Tributary', 'River branch')

# ── Main river: upstream trace along the largest contributing area ──────────
net_mask = np.zeros(acc_km2.shape, bool)
net_mask[rows[valid], cols[valid]] = True            # restrict the trace to the connected network
stem = trace_main_stem(fdir_arr, np.where(net_mask, acc_km2, 0), out_row, out_col, stream_thresh)
xs, ys = rasterio.transform.xy(transform, [p[0] for p in stem], [p[1] for p in stem])
main_coords = list(zip(xs, ys))
main_line_simplified = LineString(main_coords).simplify(tol)
# Snap endpoints to exact outlet and source coordinates for connectivity
main_coords_final = list(main_line_simplified.coords)
main_coords_final[0] = main_coords[0]   # Outlet: snap to first point
main_coords_final[-1] = main_coords[-1] # Source: snap to last point
main_line = LineString(main_coords_final)
main_river = gpd.GeoDataFrame({'Name': ['Main river']}, geometry=[main_line], crs=4326)
main_len_km = main_river.to_crs(utm_crs).length.iloc[0] / 1e3     # length in km (UTM)
stream_len_km = streams.to_crs(utm_crs).length.sum() / 1e3        # total stream length in km
print(f"✓ Streams: {len(streams):,} reaches, Strahler 1–{max_order}, {stream_len_km:,.0f} km "
      f"(threshold {stream_thresh:.2f} km²)")
print(f"  Tributaries:                 order ≥ {trib_min}")
print(f"  Main river length:           {main_len_km:,.1f} km")


# ============================================================================
# 3 — BASIN DEM AND MORPHOMETRY
# ============================================================================
# Clip the DEM to the basin (with a 2-cell buffer so the hillshade edge is clean).
with rasterio.open(dem_raw) as src:
    clip, clip_tf = rio_mask(src, [ws_geom.buffer(abs(transform.a) * 2)], crop=True,
                             nodata=-9999, filled=True)
    clip_profile = src.profile.copy()
clip = clip[0].astype('float32')
clip[clip <= -9999] = np.nan                         # nodata → NaN
inside = ~geometry_mask([ws_geom], out_shape=clip.shape, transform=clip_tf)   # True inside the basin
dem_ws = np.where(inside, clip, np.nan)              # elevation inside the basin only
wb = ws_geom.bounds
cdx, cdy = cell_size_m(clip_tf, (wb[1] + wb[3]) / 2)
clip_filled = np.where(np.isnan(clip), np.nanmean(clip), clip)   # no NaN for gradient / hillshade
slope = slope_deg(clip_filled, cdx, cdy)
clip_extent = (clip_tf.c, clip_tf.c + clip_tf.a * clip.shape[1],   # (W, E, S, N) for imshow
               clip_tf.f + clip_tf.e * clip.shape[0], clip_tf.f)

# Save the basin DEM as GeoTIFF.
dem_out = f'{OUTPUT_DIR}/{SLUG}_dem.tif'
clip_profile.update(driver='GTiff', height=clip.shape[0], width=clip.shape[1],
                    transform=clip_tf, dtype='float32', nodata=np.nan, compress='deflate')
with rasterio.open(dem_out, 'w', **clip_profile) as dst:
    dst.write(dem_ws.astype('float32'), 1)
print(f"✓ Saved: {dem_out}")

# Basic morphometric parameters of the basin.
z = dem_ws[~np.isnan(dem_ws)]                        # all valid elevations
morph = {
    'Basin':                       WATERSHED_NAME,
    'Basin_Source':                BASIN_SOURCE,
    'DEM_Resolution_m':            dem_res,
    'Outlet_Lat':                  round(out_lat, 5),
    'Outlet_Lon':                  round(out_lon, 5),
    'Area_km2':                    round(area_km2, 1),
    'Perimeter_km':                round(perim_km, 1),
    'Elev_Min_m':                  round(float(z.min()), 0),
    'Elev_Max_m':                  round(float(z.max()), 0),
    'Elev_Mean_m':                 round(float(z.mean()), 0),
    'Slope_Mean_deg':              round(float(np.nanmean(slope[inside])), 1),
    'Main_River_Length_km':        round(main_len_km, 1),
    'Stream_Length_km':            round(stream_len_km, 1),
    'Drainage_Density_km_per_km2': round(stream_len_km / area_km2, 3),       # stream km per km²
    'Max_Strahler_Order':          max_order,
    'Circularity_Ratio':           round(4 * math.pi * area_km2 / perim_km ** 2, 3),   # 1 = circle
}
morph_csv = f'{OUTPUT_DIR}/{SLUG}_morphometry.csv'
pd.DataFrame([morph]).to_csv(morph_csv, index=False)
print(f"✓ Saved: {morph_csv}")
for key, val in morph.items():
    print(f"  {key:<30} {val}")

# Save all vector layers into one GeoPackage (open in QGIS / ArcGIS).
gpkg = f'{OUTPUT_DIR}/{SLUG}_watershed.gpkg'
ws_gdf.to_file(gpkg, layer='basin', driver='GPKG')
streams.to_file(gpkg, layer='streams', driver='GPKG')
main_river.to_file(gpkg, layer='main_river', driver='GPKG')
gpd.GeoDataFrame({'Name': ['Outlet']}, geometry=[Point(out_lon, out_lat)],
                 crs=4326).to_file(gpkg, layer='outlet', driver='GPKG')
print(f"✓ Saved: {gpkg}")


# ============================================================================
# 4 — COUNTRY, NEIGHBOURS AND ADMINISTRATIVE CONTEXT
# ============================================================================
print("\nResolving country and boundaries …")
world = get_countries()
marine = get_marine()
# Host country = the one covering most of the basin (or the nearest one).
hits = world[world.intersects(ws_geom)].copy()
if COUNTRY_ISO3:
    iso3 = COUNTRY_ISO3.upper()
elif len(hits):
    hits['Overlap'] = hits.intersection(ws_geom).area
    iso3 = hits.sort_values('Overlap', ascending=False)['ISO3'].iloc[0]
else:
    iso3 = world.loc[world.distance(ws_geom.centroid).idxmin(), 'ISO3']
country_ne = world[world['ISO3'] == iso3]
country_name = country_ne['Name'].iloc[0] if len(country_ne) else iso3
print(f"  Host country:                {country_name} ({iso3})")

# Detailed country outline from geoBoundaries (fallback: Natural Earth).
country = get_geoboundaries(iso3, 'ADM0')
if country is None or country.empty:
    country = country_ne[['Name', 'geometry']]
country_union = unary_union(country.geometry)
parts = list(country_union.geoms) if isinstance(country_union, MultiPolygon) else [country_union]
# The part holding the basin frames the Location Map (skips far-off islands/territories).
host_part = max(parts, key=lambda p: (p.intersection(ws_geom).area, p.area))

# Lines inside the host country on the Location Map.
loc_lines, loc_lines_label = None, None
if LOCATION_HYBAS_LEVEL:
    try:
        if hb_region is not None and LOCATION_HYBAS_LEVEL == HYBAS_LEVEL:
            hb_loc = hb_region                       # already loaded → reuse
        else:
            hb_loc, _ = get_hydrobasins(LOCATION_HYBAS_LEVEL, ws_geom.centroid.x, ws_geom.centroid.y)
        loc_lines = hb_loc[hb_loc.intersects(country_union)].clip(country_union)
        loc_lines_label = f"{country_name} – HydroBASINS (Level {LOCATION_HYBAS_LEVEL})"
    except Exception as exc:
        print(f"  ⚠ HydroBASINS for Location Map unavailable ({exc}); using ADM1")
if loc_lines is None:
    loc_lines = get_geoboundaries(iso3, 'ADM1')      # fallback: provinces
    loc_lines_label = f"{country_name} – Provinces"

# ── Page layout (all sizes in inches) ───────────────────────────────────────
W = FIG_WIDTH_IN
ML, MR, MB, GAP = 0.32, 0.06, 0.06, 0.3          # left / right / bottom margin, gap between columns
GAP_TITLE_TO_MAP = 0.5                  # gap between title baseline and map top (all three maps)
TITLE_BAND      = 0.24                   # height reserved for the title text
MW = (W - ML - MR - GAP) * MAIN_FRAC                 # main map width
base_ext = padded_extent(ws_geom.bounds, 0.06)       # basin + 6 % margin
main_ratio = min(max(true_ratio(base_ext), 0.95), 1.45)   # frame height/width, kept within limits
MH = MW * main_ratio                                 # main map height
# Fit the extent to the frame, then push the bottom down until the scale bar clears the basin.
main_ext, sb_main = clear_scalebar(fit_exact(base_ext, main_ratio), main_ratio, ws_geom,
                                   (MW, MH), frac=0.30, segments=4)
view_main = box(main_ext[0], main_ext[2], main_ext[1], main_ext[3])   # visible area as a polygon

# Provincial boundaries for every country visible on the main map.
admin_main = None
if ADMIN_LEVEL:
    frames = []
    for iso in world[world.intersects(view_main)]['ISO3'].dropna().unique():
        g = get_geoboundaries(iso, ADMIN_LEVEL)
        if g is not None:
            frames.append(g[g.intersects(view_main)])
    if frames:
        admin_main = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=4326)
        print(f"✓ {len(admin_main)} {ADMIN_LEVEL} units on the main map")


# ============================================================================
# 5 — FIGURE
# ============================================================================
print("\nRendering map …")
RX = ML + MW + GAP                       # right block, left edge (inches)
RW = W - MR - RX                         # right block width
SGAP = 0.3                              # gap between the two small maps
SW = (RW - SGAP) / 2                     # small map width
SH = min(MH * 0.6, SW * 1.6)           # small map height
# All three maps have their top edges aligned at the same height.
MAP_TOP = MB + MH                        # top edge of all three maps (common baseline)
s_bot = MAP_TOP - SH                     # bottom edge of the small maps
Y_TITLE = MAP_TOP + GAP_TITLE_TO_MAP     # title baseline (above all maps)
H = Y_TITLE + TITLE_BAND                 # total figure height
fig = plt.figure(figsize=(W, H), facecolor='white')

def add_axes_in(x, y, w, h):
    """Add an axes positioned in inches (converted to figure fractions)."""
    return fig.add_axes([x / W, y / H, w / W, h / H])


ax_m = add_axes_in(ML, MB, MW, MH)                   # main map (left)
ax_loc = add_axes_in(RX, s_bot, SW, SH)              # Location Map (top right)
ax_elev = add_axes_in(RX + SW + SGAP, s_bot, SW, SH) # Elevation Map (top right)
ax_leg = add_axes_in(RX, MB, RW, s_bot - 0.15 - MB)  # legend (bottom right)

# ── Main map ────────────────────────────────────────────────────────────────
ax_m.set_facecolor('white')
# Layers from bottom (low zorder) to top: provinces → countries → rivers → basin outline.
if admin_main is not None:
    admin_main.boundary.plot(ax=ax_m, color=C_INK, linewidth=0.55,
                             linestyle=(0, (1, 1.3)), zorder=2)        # dotted provincial lines
# world[world.intersects(view_main)].boundary.plot(ax=ax_m, color='#444444',
#                                                  linewidth=0.6, zorder=3)   # national borders / coast
branch = streams[streams['River_Class'] == 'River branch']
trib = streams[streams['River_Class'] == 'Tributary']
# Plot river network with snapped endpoints for seamless connections at junctions
branch.plot(ax=ax_m, color=C_RIVER_BR, linewidth=0.35, zorder=4)       # small streams, thin
trib.plot(ax=ax_m, color=C_RIVER_TRIB, linewidth=0.75, zorder=5)       # tributaries, medium
main_river.plot(ax=ax_m, color=C_RIVER_MAIN, linewidth=1.4, zorder=6)  # main river, thick
ws_gdf.boundary.plot(ax=ax_m, color=C_INK, linewidth=1.0, zorder=7)    # basin boundary on top
set_map(ax_m, main_ext)
set_graticule(ax_m, main_ext, target=3)
# Compass position (axes fractions) and its footprint in map coordinates.
COMPASS_XY, COMPASS_R = (0.89, 0.92), 0.05
mw_, me_, ms_, mn_ = main_ext
kx = MH / MW
compass_box = box(mw_ + (COMPASS_XY[0] - 1.6 * COMPASS_R * kx) * (me_ - mw_),
                  ms_ + (COMPASS_XY[1] - 1.6 * COMPASS_R) * (mn_ - ms_), me_, mn_)
if admin_main is not None:
    # Names outside the basin, clear of the compass and scale bar.
    pad = 0.015 * (me_ - mw_)
    label_parts(ax_m, admin_main, view_main,
                avoid=unary_union([ws_geom.buffer(0.01 * (me_ - mw_)), compass_box.buffer(pad),
                                   sb_main['box'].buffer(pad)]),
                min_frac=0.025, fontsize=fs(FS_BODY + 0.5), path_effects=HALO)
add_compass(ax_m, xy=COMPASS_XY, r=COMPASS_R)        # top right
add_scalebar(ax_m, sb_main)                          # bottom left

# ── Location Map ────────────────────────────────────────────────────────────
loc_ext = fit_exact(padded_extent(host_part.bounds, LOCATION_PAD), SH / SW)   # host country + padding
view_loc = box(loc_ext[0], loc_ext[2], loc_ext[1], loc_ext[3])
sb_loc = scalebar_layout(loc_ext, (SW, SH), frac=0.40, segments=2, bar_h=0.022)
sb_loc_zone = sb_loc['box'].buffer(0.03 * (loc_ext[1] - loc_ext[0]))         # keep labels out of here
ax_loc.set_facecolor(C_SEA)                          # background = sea
others = world[world.intersects(view_loc) & (world['ISO3'] != iso3)]
others.plot(ax=ax_loc, facecolor=C_OTHER_LAND, edgecolor=C_OTHER_EDGE, linewidth=0.3, zorder=1)   # neighbours
country.plot(ax=ax_loc, facecolor=C_HOST, edgecolor='none', zorder=2)                             # host country
if loc_lines is not None:
    loc_lines.boundary.plot(ax=ax_loc, color=C_HOST_LINES, linewidth=0.2, zorder=3)   # internal lines
country.boundary.plot(ax=ax_loc, color='#555555', linewidth=0.4, zorder=4)
ws_gdf.plot(ax=ax_loc, facecolor=C_BASIN_FILL, edgecolor=C_BASIN_EDGE, linewidth=0.8, zorder=5)   # the basin
if (wb[2] - wb[0]) < 0.03 * (loc_ext[1] - loc_ext[0]):
    # Basin too small to see at this scale → locator ring.
    ax_loc.plot(ws_geom.centroid.x, ws_geom.centroid.y, marker='o', markersize=7,
                markerfacecolor='none', markeredgecolor=C_BASIN_EDGE, markeredgewidth=1.0, zorder=6)
# Neighbouring country names.
label_parts(ax_loc, others, view_loc, avoid=sb_loc_zone, min_frac=0.03,
            fontsize=fs(FS_SMALL), zorder=8)
# Sea / gulf names in blue italics, wrapped to ~12 characters per line.
if marine is not None:
    sea_v = marine[marine.intersects(view_loc)].copy()
    sea_v['Name'] = [textwrap.fill(nm, 12) for nm in sea_v['Name']]
    label_parts(ax_loc, sea_v, view_loc, avoid=sb_loc_zone, min_frac=0.02,
                fontsize=fs(FS_SMALL - 0.3), style='italic', color=C_SEA_LABEL)
# Host country name in bold, kept away from the basin.
host_lbl = gpd.GeoDataFrame({'Name': [country_name]}, geometry=[host_part], crs=4326)
label_parts(ax_loc, host_lbl, view_loc,
            avoid=unary_union([ws_geom.buffer(0.05 * (loc_ext[1] - loc_ext[0])), sb_loc_zone]),
            min_frac=0.0, fontsize=fs(FS_SMALL + 1), fontweight='bold', path_effects=HALO)
set_map(ax_loc, loc_ext)
set_graticule(ax_loc, loc_ext, target=2)             # add graticule with fewer labels (target=2)
add_scalebar(ax_loc, sb_loc, backing=True)           # white box: corner is covered by land/sea

# ── Elevation Map ───────────────────────────────────────────────────────────
# Legend classes: a round step giving at most ~5 intervals between min and max elevation.
lo, hi = float(z.min()), float(z.max())
step = next(st for st in [10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000]
            if math.ceil(hi / st) - math.floor(lo / st) <= 5)
breaks = np.arange(math.floor(lo / step) * step, math.ceil(hi / step) * step + step / 2, step)
cmap = LinearSegmentedColormap.from_list('elev', ELEV_COLORS, N=256)   # green → yellow → brown
norm = Normalize(vmin=breaks[0], vmax=breaks[-1])
light = LightSource(azdeg=315, altdeg=45)            # sun from the north-west, 45° high
rgba = light.shade(clip_filled, cmap=cmap, norm=norm, blend_mode='soft',
                   vert_exag=VERT_EXAG, dx=cdx, dy=cdy)                # colour + hillshade
rgba[..., 3] = inside.astype(float)                  # transparent outside the basin
elev_ext, sb_elev = clear_scalebar(fit_exact(padded_extent(ws_geom.bounds, 0.04), SH / SW),
                                   SH / SW, ws_geom, (SW, SH), frac=0.40, segments=2,
                                   bar_h=0.022)
ax_elev.set_facecolor('white')
ax_elev.imshow(rgba, extent=clip_extent, interpolation='bilinear', zorder=1)
ws_gdf.boundary.plot(ax=ax_elev, color='#333333', linewidth=0.4, zorder=2)
set_map(ax_elev, elev_ext)
set_graticule(ax_elev, elev_ext, target=2)           # add graticule with fewer labels (target=2)
add_scalebar(ax_elev, sb_elev, backing=False)        # corner is empty → no white box needed

# ── Titles: all three outside their frames, on one common baseline ──────────
title = TITLE or f"Watershed Map of {WATERSHED_NAME}, {country_name}"
for ax, ttl in ((ax_m, title), (ax_loc, 'Location Map'), (ax_elev, 'Elevation Map')):
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)                        # same frame thickness on all maps
    pos = ax.get_position()
    fig.text(pos.x0 + pos.width / 2, Y_TITLE / H, ttl, ha='center', va='bottom',
             fontsize=fs(FS_TITLE))                  # centred above each frame

# ── Legend ──────────────────────────────────────────────────────────────────
ax_leg.axis('off')
ax_leg.set_xlim(0, 1)                                # legend drawn in 0–1 coordinates
ax_leg.set_ylim(0, 1)
ax_leg.text(0.0, 0.9, 'Legend', ha='left', va='top', fontsize=fs(FS_TITLE), fontweight='bold')

# (symbol type, style, label) for each map layer.
entries = [
    ('box',  dict(fc=C_HOST, ec=C_INK, lw=0.6, ls='-'),       loc_lines_label),
    ('box',  dict(fc=C_BASIN_FILL, ec=C_BASIN_EDGE, lw=0.9, ls='-'), WATERSHED_NAME),
    ('box',  dict(fc='white', ec=C_INK, lw=1.0, ls='-'),      'Basin Boundary'),
    ('box',  dict(fc='white', ec=C_INK, lw=0.6, ls=(0, (1, 1.3))), ADMIN_LABEL),
    ('line', dict(color=C_RIVER_MAIN, lw=1.4),                 'Main River'),
    ('line', dict(color=C_RIVER_TRIB, lw=0.9),                 'Tributaries'),
    ('line', dict(color=C_RIVER_BR, lw=0.6),                   'River Branches'),
]
if admin_main is None:
    entries = [e for e in entries if e[2] != ADMIN_LABEL]   # no provinces drawn → no legend entry
y_top, y_bot = 0.76, 0.20                            # vertical range used by the entries
dy = (y_top - y_bot) / max(len(entries) - 1, 1)      # uniform spacing to fill available space
sw_w, sw_h = 0.085, 0.055                            # swatch width / height
for i, (kind, sty, label) in enumerate(entries):     # left column: map symbols
    y = y_top - i * dy
    if kind == 'box':
        ax_leg.add_patch(Rectangle((0.0, y - sw_h / 2), sw_w, sw_h, facecolor=sty['fc'],
                                   edgecolor=sty['ec'], linewidth=sty['lw'], linestyle=sty['ls']))
    else:
        ax_leg.plot([0.0, sw_w], [y, y], color=sty['color'], linewidth=sty['lw'],
                    solid_capstyle='butt')
    ax_leg.text(sw_w + 0.04, y, textwrap.fill(label, 26), ha='left', va='center',
                fontsize=fs(FS_BODY), linespacing=1.1)

# Right column: elevation colour swatches.
ex = 0.66
ax_leg.text(ex, y_top + 0.02, 'Elevation (m)', ha='left', va='bottom', fontsize=fs(FS_BODY))
edy = (y_top - 0.08 - y_bot) / max(len(breaks) - 1, 1)  # uniform spacing for elevation swatches
for i, b in enumerate(breaks):
    y = y_top - 0.06 - i * edy
    ax_leg.add_patch(Rectangle((ex, y - sw_h / 2), sw_w * 0.9, sw_h,
                               facecolor=cmap(norm(b)), edgecolor='#555555', linewidth=0.3))
    ax_leg.text(ex + sw_w * 0.9 + 0.06, y, f"{b:,.0f}", ha='left', va='center',
                fontsize=fs(FS_BODY))

ax_leg.text(0.5, 0.02, "Coordinate System: WGS 84\nDatum: WGS 84", ha='center',
            va='bottom', fontsize=fs(FS_BODY), linespacing=1.3)

# Data credits in the bottom-right corner.
# fig.text(0.995, 0.002, "Data: Copernicus DEM (ESA); HydroBASINS; Natural Earth; geoBoundaries",
#          ha='right', va='bottom', fontsize=fs(4.5), color='#888888')

# Save as PNG (raster) and PDF (vector).
for ext in ('png', 'pdf'):
    out = f'{FIG_STEM}.{ext}'
    fig.savefig(out, dpi=DPI, facecolor='white', bbox_inches='tight', pad_inches=0.04)
    print(f"✓ Saved: {out}")

plt.show()
plt.close(fig)

print("\n✓ Watershed map complete")
for p in (f'{FIG_STEM}.png', f'{FIG_STEM}.pdf', gpkg, dem_out, morph_csv):
    print(f"  {p}")
