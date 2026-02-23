mod consts;
mod coordinate;
mod parse_catalog;
mod printing;
mod solver;
mod star_quads;
mod platesolve;
mod star_finding;

#[cfg(test)]
mod tests;

// Re-export public API
pub use consts::*;
pub use consts::CameraConfig;
pub use coordinate::{Arcdegrees, CoordinateEquatorial, RaHoursMinutesSeconds, Star};
pub use platesolve::{convert_cr3_to_dng, solve_plate, solve_plate_with_options, PlateSolvingResult, ImageAnalysisResult, TransformCoefficients};
pub use printing::{dng_to_png, cr3_to_png};
