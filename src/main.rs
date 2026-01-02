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

    // Input file - can be CR3 or DNG
    let input_file = Path::new("test_img/IMG_8993.CR3");
    let dng_file = Path::new("test.dng");

    // Convert CR3 to DNG if needed
    let created_dng = if input_file.extension().map(|e| e.to_ascii_lowercase()) == Some("cr3".into()) {
        convert_cr3_to_dng(input_file, dng_file)?;
        true
    } else {
        false
    };

    // Define the initial coordinate estimate
    let initial_coord = CoordinateEquatorial::new(
        RaHoursMinutesSeconds::new(14, 15, 21.4),
        Arcdegrees::new(54, 1, 14.40),
    );

    // Solve the plate
    let result = solve_plate(dng_file, &initial_coord)?;

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

    // Clean up temporary DNG file
    if created_dng {
        if let Err(e) = std::fs::remove_file(dng_file) {
            eprintln!("Warning: could not delete temporary DNG file: {}", e);
        }
    }

    Ok(())
}