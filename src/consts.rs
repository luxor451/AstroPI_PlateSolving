//! Constants used throughout the plate solving library.

/// Tolerance for matching star quads (normalized distance difference)
pub const MATCHED_TOLERANCE: f64 = 0.001;

/// Maximum number of spiral search iterations before giving up
pub const MAX_SPIRAL_ITERATIONS: usize = 250;

/// Minimum number of matched quads required for a valid solution
pub const MIN_MATCHED_QUADS: usize = 3;

/// Maximum iterations for star matching refinement
pub const MAX_STAR_MATCH_ITERATIONS: usize = 10;

/// Initial tolerance for star matching in arcseconds
pub const INITIAL_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 60.0;

/// Minimum tolerance for star matching in arcseconds
pub const MIN_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 5.0;

/// Star detection threshold (fraction of max pixel value)
pub const STAR_THRESHOLD: f32 = 0.15;

/// Camera pixel size in microns
pub const PIXEL_SIZE_MICRON: f64 = 6.0;

/// Telescope focal length in millimeters
pub const TELESCOPE_FOCAL_LENGTH: f64 = 1200.0;

/// Camera bit depth
pub const BITDEPTH: u16 = 14;
