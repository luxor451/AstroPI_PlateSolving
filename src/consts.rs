//! Constants used throughout the plate solving library.

// =============================================================================
// Quad Matching Constants
// =============================================================================

/// Tolerance for matching star quads (normalized distance difference)
/// Higher = more lenient matching (more matches but more false positives)
/// Lower = stricter matching (fewer matches but more accurate)
/// 0.001 is very strict, 0.005 is lenient
pub const MATCHED_TOLERANCE: f64 = 0.002;

/// Minimum number of matched quads required for a valid solution
/// Higher = more confident but harder to solve
/// Need enough to reject false positives through consensus
pub const MIN_MATCHED_QUADS: usize = 40;

/// Maximum graph hops for catalog quad generation (prevents memory explosion)
pub const MAX_CATALOG_QUAD_HOPS: usize = 2;

// =============================================================================
// Spiral Search Constants
// =============================================================================

/// Maximum number of spiral search iterations before giving up
pub const MAX_SPIRAL_ITERATIONS: usize = 20000;

/// Number of positions to process in parallel during spiral search
pub const SPIRAL_BATCH_SIZE: usize = 16;

// =============================================================================
// Star Matching Refinement Constants
// =============================================================================

/// Maximum iterations for star matching refinement
pub const MAX_STAR_MATCH_ITERATIONS: usize = 15;

/// Initial tolerance for star matching in arcseconds
/// Higher gives more initial matches to refine from
pub const INITIAL_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 120.0;

/// Minimum tolerance for star matching in arcseconds
pub const MIN_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 3.0;

/// Factor by which tolerance is reduced each iteration (tolerance /= factor^iteration)
/// Lower = slower reduction, more gradual refinement
pub const STAR_MATCH_TOLERANCE_REDUCTION_FACTOR: f64 = 1.3;

/// Convergence threshold for scale/rotation changes (radians)
/// Smaller = tighter convergence
pub const CONVERGENCE_THRESHOLD: f64 = 0.00001;

// =============================================================================
// Re-centering Constants
// =============================================================================

/// Maximum iterations for re-centering refinement
pub const MAX_RE_CENTER_ITERATIONS: usize = 15;

/// Convergence threshold for position during re-centering (degrees)
pub const RE_CENTER_CONVERGENCE_DEG: f64 = 0.00001;

// =============================================================================
// Star Detection Constants
// =============================================================================

/// Star detection threshold (fraction of max pixel value)
pub const STAR_THRESHOLD: f32 = 0.15;

/// Minimum number of connected pixels to consider a detection as a star
pub const MIN_STAR_PIXELS: u32 = 5;

/// Number of nearest neighbors for star graph construction
/// Higher = more quads generated (better matching but more memory/time)
pub const STAR_GRAPH_K_NEIGHBORS: usize = 6;

// =============================================================================
// Camera/Telescope Constants
// =============================================================================

/// Runtime-configurable camera/telescope parameters.
///
/// These values can be set via the `/update_camera_settings` endpoint
/// so that different optical setups work without recompilation.
#[derive(Debug, Clone)]
pub struct CameraConfig {
    /// Camera pixel size in microns
    pub pixel_size_micron: f64,
    /// Telescope focal length in millimeters
    pub focal_length_mm: f64,
    /// Camera bit depth
    pub bitdepth: u16,
}

impl Default for CameraConfig {
    fn default() -> Self {
        Self {
            pixel_size_micron: 6.0,
            focal_length_mm: 714.0,
            bitdepth: 14,
        }
    }
}

impl CameraConfig {
    /// Compute the expected pixel scale in arcsec/pixel.
    ///
    /// Formula: 206.265 * pixel_size_μm / focal_length_mm
    pub fn expected_pixel_scale(&self) -> f64 {
        ARCSEC_PER_RADIAN * self.pixel_size_micron / self.focal_length_mm
    }
}

/// Default camera pixel size in microns (used when no config is provided)
pub const PIXEL_SIZE_MICRON: f64 = 6.0;

/// Default telescope focal length in millimeters
pub const TELESCOPE_FOCAL_LENGTH: f64 = 714.0;

/// Default camera bit depth
pub const BITDEPTH: u16 = 14;

/// Conversion factor from focal length to pixel scale (arcsec/pixel = 206.265 * pixel_size_um / focal_length_mm)
pub const ARCSEC_PER_RADIAN: f64 = 206.265;

/// Expected pixel scale (arcsec/pixel) = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGTH
/// Used to validate solutions - reject if scale is too far from expected
/// Note: Actual scale from solving may differ slightly due to optical effects
/// For 714mm f/l + 6μm pixels: 206.265 * 6.0 / 714 = 1.733 arcsec/pixel
/// For 1200mm f/l + 6μm pixels: 206.265 * 6.0 / 1200 = 1.031 arcsec/pixel
/// Using mid-point to support multiple optical setups
pub const EXPECTED_PIXEL_SCALE: f64 = 1.4;

/// Tolerance for pixel scale validation (fraction of expected scale)
/// e.g., 0.5 means accept scales within ±50% of expected
/// 1.4 * 0.5 = 0.7, so range is 0.7 to 2.1 arcsec/pixel (covers 714mm to 1200mm setups)
pub const SCALE_TOLERANCE_FRACTION: f64 = 0.50;

/// Multiplier for maximum stars to fetch from catalog (relative to detected stars)
pub const CATALOG_STAR_MULTIPLIER: usize = 2;

// =============================================================================
// HTTP Client Constants (VizieR Connection)
// =============================================================================
#[allow(dead_code)]
/// Maximum idle connections per host for VizieR HTTP client
pub const VIZIER_POOL_MAX_IDLE_PER_HOST: usize = 16;
#[allow(dead_code)]
/// Idle timeout for VizieR HTTP connections in seconds
pub const VIZIER_POOL_IDLE_TIMEOUT_SECS: u64 = 60;
#[allow(dead_code)]
/// Request timeout for VizieR HTTP requests in seconds
pub const VIZIER_REQUEST_TIMEOUT_SECS: u64 = 30;
#[allow(dead_code)]
/// TCP keepalive interval for VizieR connections in seconds
pub const VIZIER_TCP_KEEPALIVE_SECS: u64 = 30;
#[allow(dead_code)]
/// VizieR TAP service URL
pub const VIZIER_TAP_URL: &str = "http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync";
#[allow(dead_code)]
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
pub const SVD_TOLERANCE: f64 = 1e-16;
