//! Constants used throughout the plate solving library.

// =============================================================================
// Quad Matching Constants
// =============================================================================

/// Tolerance for matching star quads (normalized distance difference)
pub const MATCHED_TOLERANCE: f64 = 0.001;

/// Minimum number of matched quads required for a valid solution
pub const MIN_MATCHED_QUADS: usize = 3;

// =============================================================================
// Spiral Search Constants
// =============================================================================

/// Maximum number of spiral search iterations before giving up
pub const MAX_SPIRAL_ITERATIONS: usize = 500;

/// Number of positions to process in parallel during spiral search
pub const SPIRAL_BATCH_SIZE: usize = 16;

// =============================================================================
// Star Matching Refinement Constants
// =============================================================================

/// Maximum iterations for star matching refinement
pub const MAX_STAR_MATCH_ITERATIONS: usize = 10;

/// Initial tolerance for star matching in arcseconds
pub const INITIAL_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 60.0;

/// Minimum tolerance for star matching in arcseconds
pub const MIN_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 5.0;

/// Factor by which tolerance is reduced each iteration (tolerance /= factor^iteration)
pub const STAR_MATCH_TOLERANCE_REDUCTION_FACTOR: f64 = 1.5;

/// Convergence threshold for scale/rotation changes (radians)
pub const CONVERGENCE_THRESHOLD: f64 = 0.0001;

// =============================================================================
// Re-centering Constants
// =============================================================================

/// Threshold in degrees to trigger re-centering refinement
pub const RE_CENTER_THRESHOLD_DEG: f64 = 0.5;

/// Maximum iterations for re-centering refinement
pub const MAX_RE_CENTER_ITERATIONS: usize = 3;

/// Convergence threshold for position during re-centering (degrees)
pub const RE_CENTER_CONVERGENCE_DEG: f64 = 0.01;

// =============================================================================
// Star Detection Constants
// =============================================================================

/// Star detection threshold (fraction of max pixel value)
pub const STAR_THRESHOLD: f32 = 0.15;

/// Minimum number of connected pixels to consider a detection as a star
pub const MIN_STAR_PIXELS: u32 = 2;

/// Number of nearest neighbors for star graph construction
pub const STAR_GRAPH_K_NEIGHBORS: usize = 3;

// =============================================================================
// Camera/Telescope Constants
// =============================================================================

/// Camera pixel size in microns
pub const PIXEL_SIZE_MICRON: f64 = 6.0;

/// Telescope focal length in millimeters
pub const TELESCOPE_FOCAL_LENGTH: f64 = 1200.0;

/// Camera bit depth
pub const BITDEPTH: u16 = 14;

/// Conversion factor from focal length to pixel scale (arcsec/pixel = 206.265 * pixel_size_um / focal_length_mm)
pub const ARCSEC_PER_RADIAN: f64 = 206.265;

/// Multiplier for maximum stars to fetch from catalog (relative to detected stars)
pub const CATALOG_STAR_MULTIPLIER: usize = 100000;

// =============================================================================
// HTTP Client Constants (VizieR Connection)
// =============================================================================

/// Maximum idle connections per host for VizieR HTTP client
pub const VIZIER_POOL_MAX_IDLE_PER_HOST: usize = 16;

/// Idle timeout for VizieR HTTP connections in seconds
pub const VIZIER_POOL_IDLE_TIMEOUT_SECS: u64 = 60;

/// Request timeout for VizieR HTTP requests in seconds
pub const VIZIER_REQUEST_TIMEOUT_SECS: u64 = 30;

/// TCP keepalive interval for VizieR connections in seconds
pub const VIZIER_TCP_KEEPALIVE_SECS: u64 = 30;

/// VizieR TAP service URL
pub const VIZIER_TAP_URL: &str = "http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync";

/// VizieR catalog table name (UCAC5)
pub const VIZIER_CATALOG_TABLE: &str = "I/322A/out";

// =============================================================================
// Image Rendering Constants
// =============================================================================

/// Default low percentile for histogram stretch
pub const HISTOGRAM_STRETCH_LOW_PERCENTILE: f64 = 0.01;

/// Default high percentile for histogram stretch
pub const HISTOGRAM_STRETCH_HIGH_PERCENTILE: f64 = 0.995;

/// Pixel value stretch factor for PNG export
pub const PIXEL_STRETCH_FACTOR: u32 = 4;

/// Size of crosshair marker in pixels (for center marking)
pub const CROSSHAIR_MARKER_SIZE: u32 = 20;

/// Radius of circle marker at image center
pub const CENTER_CIRCLE_RADIUS: u32 = 10;

// =============================================================================
// Numerical Constants
// =============================================================================

/// SVD solver tolerance for least squares fitting
pub const SVD_TOLERANCE: f64 = 1e-14;
