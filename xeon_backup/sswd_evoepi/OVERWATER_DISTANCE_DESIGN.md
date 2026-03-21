# Overwater Distance Calculation Design

## Executive Summary

This document specifies the design for computing shortest overwater distances between all pairs of sea star survey sites (~150-300 sites) across the Northeast Pacific (Aleutians to Baja California) for the SSWD-EvoEpi model. The approach will replace the current Haversine × tortuosity approximation with accurate pathfinding that respects coastline geometry, islands, fjords, and narrow passages.

**Recommended Approach**: Hybrid raster-Dijkstra with adaptive resolution and visibility graph fallback.

**Key Design Decisions**:
- **Primary**: 500m resolution raster with Dijkstra pathfinding on sea cells
- **Fallback**: Visibility graph for complex cases or performance optimization
- **CRS**: Alaska Albers (EPSG:3338) for northern sites, custom Lambert for southern sites
- **Antimeridian**: Split Aleutian region into separate processing zone
- **Memory target**: <8GB RAM, <30 minutes compute time for ~45,000 pairs

## Current Limitations

The existing SSWD-EvoEpi model uses:
```python
def compute_distance_matrix(lats, lons, tortuosity=1.5):
    # Haversine great-circle × uniform tortuosity
    return haversine_km(lat1, lon1, lat2, lon2) * tortuosity
```

**Problems**:
- Uniform tortuosity (1.5) doesn't capture coastline complexity variation
- Can't resolve narrow passages (Hood Canal, Johnstone Strait)
- No handling of island barriers
- Antimeridian issues for Aleutian Islands
- Doesn't account for fjord sill depths or complex inland waterways

## Literature Review

**Oceanographic Standards**:
- 1km resolution is standard for marine connectivity studies (Anadón et al. 2013)
- Rasterized approaches with land/sea binary grids are well-established
- Scientists commonly use Dijkstra/A* on sea cells for accurate distances

**Existing Libraries**:
- **searoute**: Uses MARNET maritime network, ~0.5s per route, good for shipping
- **scgraph**: Faster maritime network (~0.05s), less coastal detail
- **pyvisgraph**: Visibility graphs for polygon obstacles, O(n² log n)
- **extremitypathfinder**: Fast geometric pathfinding in 2D multi-polygon environments

**Assessment**: Existing maritime libraries focus on shipping routes and lack fine coastal detail needed for nearshore species. Custom raster approach provides optimal accuracy-performance balance.

## Recommended Approach: Hybrid Raster-Dijkstra

### 1. Primary Method: Adaptive Resolution Raster

**Algorithm**:
1. **Coastline Rasterization**: Convert coastline polygons to binary land/sea grid
2. **Sea Cell Graph**: Build 8-connected adjacency graph of sea cells
3. **Dijkstra Pathfinding**: Compute shortest paths between all site pairs
4. **Distance Conversion**: Convert cell paths to metric distances

**Resolution Strategy**:
- **Base resolution**: 500m (balance between accuracy and memory)
- **High-detail zones**: 250m for complex areas (Salish Sea, fjords)
- **Low-detail zones**: 1km for open coast segments
- **Memory estimate**: ~2GB for 500m, ~8GB for 250m full region

### 2. Fallback Method: Visibility Graph

For performance optimization or complex polygon handling:
- **Coastline simplification**: Douglas-Peucker algorithm to reduce vertices
- **Visibility graph construction**: Connect visible coastline vertices
- **Dijkstra on graph**: Find shortest paths avoiding land polygons
- **Use case**: High-level connectivity analysis, rapid prototyping

### 3. Coordinate Reference Systems

**Alaska Region (>55°N)**:
- **CRS**: Alaska Albers Equal Area (EPSG:3338)
- **Parameters**: Center longitude -154°, standard parallels 55°N and 65°N
- **Advantages**: Minimal distortion for Alaska, equal-area for distance calculations

**Pacific Coast Region (<55°N)**:
- **CRS**: Custom Lambert Conformal Conic
- **Parameters**: Center longitude -125°, standard parallels 35°N and 50°N
- **Coverage**: British Columbia to Baja California

**Antimeridian Handling**:
- **Problem**: Aleutian Islands span ±180° longitude
- **Solution**: 
  1. Split dataset at 170°W
  2. Process western Aleutians separately with longitude shift (+360°)
  3. Merge results with proper antimeridian distance calculation
  4. Use `--config CENTER_LONG 180` for GDAL operations

### 4. Edge Cases and Special Handling

**Narrow Passages**:
- **Hood Canal**: 1.6km width at sill, ensure 250m resolution captures passage
- **Johnstone Strait**: 2-3km width, critical for north-south connectivity
- **Solution**: Adaptive resolution increase near known narrow passages

**Fjord Sills**:
- **Current model**: Uses `sill_depth` for pathogen dispersal attenuation
- **Distance calculation**: Sills don't block larval transport, so overwater distance ignores depth
- **Implementation**: Treat sills as open water for distance calculation

**Island Barriers**:
- **Vancouver Island**: Major barrier requiring circumnavigation
- **Queen Charlotte Islands**: Separate northern/southern populations
- **Aleutian Chain**: Complex island-hopping paths

**Coastal Cells**:
- **Cost**: Base cost = 1.0 (1 cell = 500m distance)
- **Coastal penalty**: 1.2× cost for cells adjacent to land (reflects wave exposure, shallow water)
- **Deep water**: 0.9× cost for cells >20km from shore (faster pelagic transport)

### 5. Data Flow Pipeline

```
Coastline Shapefile
        ↓
    [Project to CRS]
        ↓
    [Rasterize to Binary Grid]
        ↓
    [Build Adjacency Graph]
        ↓
    [Load Site Coordinates]
        ↓
    [Map Sites to Grid Cells]
        ↓
    [Dijkstra All-Pairs Shortest Path]
        ↓
    [Convert Cell-Distance to Metric]
        ↓
    [Export Distance Matrix]
```

## Implementation Specification

### Required Dependencies

```python
# Core GIS and raster processing
pip install geopandas>=1.0 shapely>=2.0 fiona>=1.8
pip install rasterio>=1.3 pyproj>=3.5

# Pathfinding and numerical computation
pip install scipy>=1.10 scikit-image>=0.20 networkx>=3.0

# Optional: Maritime routing libraries for comparison
pip install searoute pyvisgraph extremitypathfinder
```

### Core Classes and Functions

```python
class OverwaterDistanceCalculator:
    """Compute shortest overwater distances between coastal sites."""
    
    def __init__(self, 
                 coastline_path: str,
                 resolution_m: float = 500.0,
                 crs: str = "EPSG:3338",
                 adaptive_resolution: bool = True):
        """Initialize with coastline data and resolution settings."""
        
    def load_coastline(self, path: str) -> None:
        """Load coastline shapefile/GeoJSON and reproject to working CRS."""
        
    def rasterize_coastline(self, bounds: tuple, resolution_m: float) -> np.ndarray:
        """Convert coastline polygons to binary land(0)/sea(1) raster."""
        
    def build_sea_graph(self, raster: np.ndarray) -> scipy.sparse.csr_matrix:
        """Build 8-connected adjacency graph of sea cells."""
        
    def map_sites_to_cells(self, site_coords: np.ndarray) -> np.ndarray:
        """Map lat/lon coordinates to raster cell indices."""
        
    def compute_all_pairs_distances(self, 
                                   site_cells: np.ndarray,
                                   sea_graph: scipy.sparse.csr_matrix,
                                   resolution_m: float) -> np.ndarray:
        """Use Dijkstra to compute all-pairs shortest distances."""
        
    def handle_antimeridian_split(self, sites: gpd.GeoDataFrame) -> np.ndarray:
        """Special handling for sites crossing ±180° longitude."""
        
    @classmethod
    def from_sites_file(cls, sites_csv: str, coastline_path: str) -> 'OverwaterDistanceCalculator':
        """Factory method to create calculator from site definitions."""

# Utility functions
def validate_distance_matrix(distances: np.ndarray, 
                           site_coords: np.ndarray) -> Dict[str, float]:
    """Validate distance matrix against Haversine baselines."""

def export_connectivity_matrices(distances: np.ndarray,
                               D_L: float, D_P: float,
                               output_dir: str) -> None:
    """Export C (larval) and D (pathogen) matrices for SSWD-EvoEpi."""
```

### Performance Specifications

**Memory Requirements**:
- **500m resolution**: ~2GB RAM for full NE Pacific region
- **250m high-detail**: ~8GB RAM (requires chunked processing)
- **Sparse graph storage**: CSR format reduces memory by ~80%

**Compute Time Estimates**:
- **150 sites (11,175 pairs)**: ~5-10 minutes on modern CPU
- **300 sites (44,850 pairs)**: ~15-30 minutes
- **Parallel processing**: Dijkstra can be parallelized across source nodes

**Accuracy Validation**:
- **Simple cases**: <5% difference from great-circle distance in open water
- **Complex cases**: Manual verification for key fjord/strait passages
- **Comparison benchmark**: Match literature values for well-studied routes

### Output Specifications

**Distance Matrix Format**:
```python
# distances.npz (compressed NumPy format)
{
    'distances': np.ndarray,      # (N, N) float32, kilometers
    'site_ids': np.ndarray,       # (N,) int32, original site IDs  
    'coordinates': np.ndarray,    # (N, 2) float64, lat/lon
    'metadata': dict              # Resolution, CRS, timestamp, version
}
```

**Integration with SSWD-EvoEpi**:
```python
# Replace current distance calculation in spatial.py
def compute_distance_matrix(lats, lons, 
                          coastline_path=None, 
                          use_overwater=True,
                          tortuosity=None):
    if use_overwater and coastline_path:
        calc = OverwaterDistanceCalculator.from_sites_file(
            sites_csv="data/nodes/node_definitions.csv",
            coastline_path=coastline_path
        )
        return calc.compute_all_pairs_distances()
    else:
        # Fallback to current Haversine approach
        return haversine_distance_matrix(lats, lons, tortuosity)
```

## Data Requirements

### Coastline Data Sources

**Primary**: Natural Earth Land polygons (1:10m resolution)
- **URL**: https://www.naturalearthdata.com/downloads/10m-physical-vectors/
- **Advantages**: Global coverage, well-maintained, appropriate detail level
- **Format**: Shapefile or GeoPackage

**Alternative**: GSHHG (Global Self-consistent Hierarchical High-resolution Geography)
- **URL**: https://www.soest.hawaii.edu/pwessel/gshhg/
- **Advantages**: Higher resolution available, optimized for ocean applications
- **Format**: Shapefile

**Regional Enhancement**: 
- **NOAA Coastal Change Analysis Program (C-CAP)**: High-resolution for US coasts
- **Canadian Hydrographic Service**: Detailed charts for British Columbia fjords

### Site Definitions

**Current Format**: `data/nodes/node_definitions.csv`
```csv
node_id,name,lat,lon,subregion,habitat_area,carrying_capacity,is_fjord
0,Sitka AK,57.06,-135.34,AK-SE,50000.0,1000,False
1,Howe Sound BC,49.52,-123.25,SS,20000.0,400,True
...
```

**Enhancement Needed**: Add complexity flags for adaptive resolution
```csv
...,complexity,narrow_passage,antimeridian_zone
...,high,False,False
...,very_high,True,False
```

## Validation and Testing

### Test Cases

**Unit Tests**:
1. **Rasterization accuracy**: Known coastline features preserved
2. **Graph connectivity**: No isolated sea regions
3. **Distance calculation**: Metric conversion accuracy
4. **CRS handling**: Correct projection transformations

**Integration Tests**:
1. **5-node network**: Match current model behavior in simple cases  
2. **Known routes**: Vancouver to Seattle, Sitka to Juneau
3. **Fjord passages**: Hood Canal, Howe Sound depth/width validation
4. **Antimeridian**: Attu Island to mainland Alaska

**Performance Tests**:
1. **Memory profiling**: Confirm <8GB for 300 sites
2. **Speed benchmarks**: <30min for full 300×300 matrix
3. **Scaling**: Linear growth with number of sites

### Error Handling

**Coastline Gaps**:
- **Detection**: Check for isolated site clusters
- **Resolution**: Manual gap filling or connection validation
- **Fallback**: Use Haversine × region-specific tortuosity

**Projection Issues**:
- **Antimeridian wrapping**: Validate longitude normalization
- **CRS distortion**: Check distance accuracy at projection edges
- **Site mapping**: Ensure all sites fall within valid raster bounds

**Memory Limitations**:
- **Chunked processing**: Process large matrices in blocks
- **Sparse methods**: Use scipy.sparse throughout
- **Early termination**: Stop on memory exhaustion, return partial results

## Future Enhancements

### Phase 2 (Optional)

**3D Bathymetry Integration**:
- **Data**: ETOPO1 or GEBCO bathymetry grids
- **Purpose**: Depth-weighted pathfinding for deep-water species
- **Cost function**: Penalize deep water based on species depth preferences

**Current-Driven Routing**:
- **Data**: Ocean circulation models (ROMS, HYCOM)
- **Purpose**: Current-assisted vs current-opposed pathfinding
- **Application**: Larval dispersal with realistic oceanography

**Seasonal Variation**:
- **Data**: Seasonal ice extent, storm patterns
- **Purpose**: Time-varying connectivity matrices
- **Application**: Match breeding/dispersal seasons to optimal routes

### Phase 3 (Research)

**Machine Learning Optimization**:
- **Training data**: GPS tracks from tagged sea stars (if available)
- **Algorithm**: Graph neural networks for pathway prediction
- **Validation**: Compare ML routes with observed dispersal patterns

**Uncertainty Quantification**:
- **Coastline uncertainty**: Monte Carlo with perturbed shorelines
- **Route confidence**: Multiple path analysis with confidence intervals
- **Parameter sensitivity**: Gaussian process surrogate models

## Implementation Timeline

**Week 1**: Core raster-Dijkstra implementation
- Load/project coastline data
- Rasterization and graph building
- Basic all-pairs pathfinding

**Week 2**: SSWD-EvoEpi integration  
- Replace distance matrix in spatial.py
- Validate against 5-node network
- Performance optimization

**Week 3**: Complex feature support
- Antimeridian handling for Aleutians
- Adaptive resolution implementation
- Edge case validation

**Week 4**: Full-scale testing
- 150-300 site distance matrix
- Connectivity matrix export (C, D)
- Documentation and deployment

**Dependencies**: 
- Coastline data procurement (Natural Earth): ~1 day
- Site definition enhancement: ~2 days
- HPC resource allocation (if needed): ~1 week lead time

---

## References

- Anadón et al. (2013). "Habitat-specific larval dispersal and marine connectivity." *Ecosphere* 4(7):1-32.
- Cowen & Sponaugle (2009). "Larval dispersal and marine population connectivity." *Ann Rev Mar Sci* 1:443-466.
- Lee (1978). "Proximity and reachability in the plane." *Technical Report*, Stanford University.
- Natural Earth. "1:10m Physical Vectors." https://www.naturalearthdata.com/
- White et al. (2010). "Connectivity, dispersal and recruitment." *Oceanography* 23(3):34-45.

**Code References**:
- extremitypathfinder: https://github.com/MrMinimal64/extremitypathfinder
- pyvisgraph: https://github.com/TaipanRex/pyvisgraph  
- searoute: https://github.com/eurostat/searoute
- scikit-image shortest_path: https://scikit-image.org/docs/dev/api/skimage.graph.html

---

*This design document will be read by the next cron job (overwater-2-implement) as implementation specification.*