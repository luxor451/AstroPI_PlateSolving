use std::path::Path;

mod coordinate;
mod parse_catalog;
mod printing;
mod solver;
mod star_quads;
mod platesolve;

#[cfg(test)]
mod tests;

use coordinate::*;
use platesolve::*;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize logger - can be controlled via RUST_LOG env var
    // e.g., RUST_LOG=debug cargo run, RUST_LOG=trace cargo run
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    // Define the initial coordinate estimate
    let initial_coord = CoordinateEquatorial::new(
        RaHoursMinutesSeconds::new(14, 15, 21.4),
        Arcdegrees::new(54, 1, 14.40),
    );

    // Path to the DNG file to analyze
    let file_path = Path::new("test.dng");

    // Solve the plate
    let result = solve_plate(file_path, &initial_coord)?;

    // Display results
    if result.coeffs_x.is_some() {
        let ra_rad = result.optical_axis_ra;
        let dec_rad = result.optical_axis_dec;
        let solved_coord = CoordinateEquatorial::from_radians(ra_rad, dec_rad);

        println!("RA:  {}", solved_coord.ra);
        println!("Dec: {}", solved_coord.dec);
        println!("Position angle: {:.2}°", result.rotation_deg);
    } else {
        eprintln!("Plate solving failed: not enough matches.");
    }

    Ok(())
}