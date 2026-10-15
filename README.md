# Watershed Map Generator

A Python-based geospatial analysis tool that generates professional-quality watershed maps for any basin worldwide using open data sources. This tool is designed to produce publication-ready figures displaying drainage networks, elevation data, and comprehensive basin morphometry.

## Overview

This project automates the creation of multi-panel watershed maps in a classic study-area layout:

- **Main Map**: Basin boundary, drainage network (classified by Strahler order), provincial boundaries, compass rose, scale bar, and decimal-degree graticule
- **Location Map**: Host country context with internal administrative divisions and the basin highlighted in red
- **Elevation Map**: Hillshaded DEM of the basin showing topography
- **Legend**: Map symbols, elevation color scale, and coordinate system information

## Features

- **Global Coverage**: Works for any watershed worldwide using the HydroBASINS dataset
- **Automatic Data Retrieval**: Downloads all required data automatically from open sources
- **High-Resolution Output**: Generates both PNG (600 dpi) and PDF (vector) formats
- **GIS-Ready Outputs**: Exports basin geometry, streams, and DEM as GeoPackage for further analysis
- **Morphometric Analysis**: Calculates comprehensive basin statistics (area, perimeter, elevation, drainage density, etc.)
- **Flexible Basin Definition**:
  - HydroBASINS polygon matching (recommended for large/medium basins)
  - D8 delineation from Copernicus DEM using pysheds (best for small basins/gauges)
  - User-supplied watershed boundaries

## Data Sources

All data is retrieved automatically from open sources:

| Data Type | Source | Coverage |
|-----------|--------|----------|
| DEM | Copernicus GLO-30/GLO-90 (AWS Open Data) | Global, 30m or 90m resolution |
| Basin Boundaries | HydroBASINS v1c | Global, multiple hierarchy levels |
| Administrative Boundaries | geoBoundaries gbOpen | Countries and provinces |
| Country/Coast Data | Natural Earth 1:50m | Global reference data |

## Requirements

- Python 3.11+
- GDAL/geopandas stack

### Installation

```bash
conda create -n hydrology -c conda-forge python=3.11 geopandas rasterio pysheds requests
conda activate hydrology
```

### Dependencies
- **geopandas**: Vector geospatial data handling
- **rasterio**: Raster DEM processing
- **pysheds**: D8 flow routing and watershed delineation
- **matplotlib**: Publication-quality figure generation
- **requests**: Downloading data from web sources
- **shapely**: Geometry operations
- **pandas/numpy**: Data manipulation

## Usage

### Basic Configuration

Edit the configuration section in `Study_area.py` (lines 96–150):

```python
WATERSHED_NAME    = 'Ba River Basin'      # Basin name for output files
OUTLET_LAT        = 13.05                 # Outlet point (WGS 84)
OUTLET_LON        = 108.98

BASIN_SOURCE      = 'hydrobasins'        # 'hydrobasins' | 'dem' | 'file'
HYBAS_LEVEL       = 6                    # HydroBASINS detail level (1–12)
HYBAS_WHOLE_BASIN = True                 # True = entire river system to sea

DEM_RES           = 'auto'               # DEM resolution ('auto' | 30 | 90 metres)
STREAM_DENSITY    = 600                  # Channel threshold = basin area / this
COUNTRY_ISO3      = None                 # e.g., 'VNM' for Vietnam; None = auto-detect
ADMIN_LEVEL       = 'ADM1'               # Provincial boundaries ('ADM1', 'ADM2', None)
```

### Running the Script

```bash
python Study_area.py
```

The script will:
1. Download required data automatically (cached for reuse)
2. Delineate the basin and compute flow routing
3. Calculate stream networks and morphometry
4. Generate and save all outputs

## Output Files

For each basin, outputs are saved in `outputs/<basin_slug>/`:

### Maps
- `<name>_watershed_map.png` — 600 dpi raster image, ready for publication
- `<name>_watershed_map.pdf` — Vector format, editable in Illustrator/CorelDRAW

### Geospatial Data
- `<name>_watershed.gpkg` — GeoPackage containing:
  - `basin` — Watershed boundary polygon
  - `streams` — Stream network with Strahler order classification
  - `main_river` — Main stem traced upstream from outlet
  - `outlet` — Outlet point location

### DEM & Analysis
- `<name>_dem.tif` — DEM clipped to basin (GeoTIFF format)
- `<name>_morphometry.csv` — Basin statistics table

### Example Output Structure
```
outputs/
├── ba_river_basin/
│   ├── ba_river_basin_watershed_map.png
│   ├── ba_river_basin_watershed_map.pdf
│   ├── ba_river_basin_watershed.gpkg
│   ├── ba_river_basin_dem.tif
│   └── ba_river_basin_morphometry.csv
└── _cache/
    ├── cop30_27.671_85.344_0.30.tif
    ├── hybas_as_lev06_v1c.zip
    └── ... (downloaded data, reused across basins)
```

## Basin Source Options

### 1. HydroBASINS (Recommended for medium–large basins)

Uses the HydroBASINS v1c dataset:

```python
BASIN_SOURCE      = 'hydrobasins'
OUTLET_LAT        = 13.05
OUTLET_LON        = 108.98
HYBAS_LEVEL       = 6                # Higher = finer detail
HYBAS_WHOLE_BASIN = True             # True = entire river system to sea
```

**Pros**: Fast, accurate, consistent with global datasets  
**Cons**: Limited to pre-defined basin hierarchy

### 2. DEM Delineation (Best for small basins/gauges)

Uses D8 routing from Copernicus DEM with pysheds:

```python
BASIN_SOURCE      = 'dem'
OUTLET_LAT        = 13.05
OUTLET_LON        = 108.98
EXPECTED_AREA_KM2 = None              # Known area helps snapping
SNAP_RADIUS_M     = 1000              # Search radius for snapping
```

**Pros**: Works for any point, custom delineation  
**Cons**: Slower, requires more RAM for large basins (>100,000 km²)

### 3. User-Supplied Boundary (Maximum flexibility)

Provide your own watershed polygon:

```python
BASIN_SOURCE      = 'file'
WATERSHED_FILE    = 'shp/my_basin.shp'  # .shp, .gpkg, or .geojson
```

## Customization

### Typography & Style

```python
FONT_FAMILY       = ['Times New Roman', 'Times', 'Liberation Serif', 'DejaVu Serif']
FONT_SCALE        = 1.0               # Multiply all font sizes
FIG_WIDTH_IN      = 7.2               # Figure width (inches)
DPI               = 600               # PNG resolution
```

### Colors

Modify the color palette (lines 152–165):

```python
C_RIVER_MAIN      = '#0B1F8A'         # Main river (navy)
C_RIVER_TRIB      = '#3D7DD8'         # Tributaries (mid-blue)
C_RIVER_BR        = '#A9CDEB'         # River branches (light blue)
ELEV_COLORS       = ['#2E6B12', '#8DBF57', '#FFF9C9', '#F2C14E', '#C98B57', '#8A5A2B']
```

### Drainage Network Classification

```python
STREAM_DENSITY    = 600               # Higher = denser network
TRIB_MIN_ORDER    = None              # Strahler order threshold for tributaries
                                      # None = max_order − 2 (automatic)
```

## Morphometric Output

The `morphometry.csv` includes:

| Parameter | Unit | Description |
|-----------|------|-------------|
| Area_km2 | km² | Total watershed area |
| Perimeter_km | km | Watershed perimeter |
| Elev_Min_m | m | Minimum elevation |
| Elev_Max_m | m | Maximum elevation |
| Elev_Mean_m | m | Mean elevation |
| Slope_Mean_deg | ° | Mean slope gradient |
| Main_River_Length_km | km | Main river stem length |
| Stream_Length_km | km | Total stream length |
| Drainage_Density_km_per_km2 | km/km² | Stream density |
| Circularity_Ratio | – | 1 = perfect circle |
| Max_Strahler_Order | – | Maximum stream order |

## Technical Notes

### Limitations

- **Latitude Limits**: Plotting in WGS 84 lon/lat with true-scale aspect works best for basins <70° latitude and not crossing the antimeridian
- **Large Basins**: Basins >100,000 km² require significant RAM at 30m resolution; the script automatically switches to 90m for large windows
- **Ocean Coverage**: Copernicus DEM has gaps in some ocean regions; the script will warn if the basin touches these areas

### Performance

- **Small basins** (<1,000 km²): ~5–10 minutes (mostly downloads)
- **Medium basins** (1,000–50,000 km²): ~15–30 minutes
- **Large basins** (>50,000 km²): >1 hour (may require 8+ GB RAM at 30m resolution)

Execution time is dominated by DEM tile downloads. The `_cache/` directory preserves downloaded tiles, so subsequent analyses of nearby basins are much faster.

### Memory Requirements

- **30m DEM**: ~1 GB per degree of latitude × longitude
- **90m DEM**: ~100 MB per degree
- Automatic resolution selection keeps most analyses under 4 GB

## References

- Lehner, B., & Grill, G. (2013). Global river hydrography and network routing: baseline data and new approaches to study the world's large river systems. *Hydrological Processes*, 27(15), 2171–2186. [HydroBASINS](https://www.hydrosheds.org/)
- Kelso, N. V., et al. Natural Earth — free vector and raster geospatial data. [Natural Earth Data](https://www.naturalearthdata.com/)
- Copernicus Digital Elevation Model. European Space Agency. [AWS Open Data](https://registry.opendata.aws/copernicus-dem/)
- geoBoundaries. World Bank. [geoBoundaries](https://www.geoboundaries.org/)

## License

This project is provided as-is for research and educational use.

## Contact & Citation

For questions or improvements, refer to the repository issues.

If you use this tool in your research, please cite the underlying data sources:
- HydroBASINS (Lehner & Grill 2013)
- Copernicus DEM
- Natural Earth
- geoBoundaries

---

**Last Updated**: October 2026  
**Author**: Hanik Lakhe
