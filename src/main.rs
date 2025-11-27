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
    // Define the initial coordinate estimate
    let initial_coord = CoordinateEquatorial::new(
        RaHoursMinutesSeconds::new(14, 01, 12.5),
        Arcdegrees::new(54, 20, 56.0),
    );

    // Path to the DNG file to analyze
    let file_path = Path::new("test.dng");

    println!("Starting plate solving for: {}", file_path.display());
    println!("Initial coordinates: {}", initial_coord);
    println!();

    // Solve the plate
    let result = solve_plate(file_path, &initial_coord, true)?;

    // Display results
    println!();
    println!("=== Plate Solving Results ===");
    println!("Matched quads: {}", result.matched_quads_count);

    if let (Some(coeffs_x), Some(coeffs_y)) = (result.coeffs_x, result.coeffs_y) {
        println!("Transformation coefficients:");
        println!("  X: {:?}", coeffs_x);
        println!("  Y: {:?}", coeffs_y);
        println!();
        println!("Plate solving successful!");
    } else {
        println!();
        println!("Not enough matches to determine transformation coefficients.");
        println!("Plate solving incomplete.");
    }

    Ok(())
}
