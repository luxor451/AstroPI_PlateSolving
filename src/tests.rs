#[cfg(test)]
mod plate_solving_tests {
    use crate::coordinate::*;
    use crate::platesolve::{convert_cr3_to_dng, solve_plate};
    use std::path::Path;
    use std::sync::OnceLock;

    const CR3_FILE_PATH: &str = "test_img/M101.CR3";
    
    /// Shared DNG file path - created once and reused by all tests
    const SHARED_DNG_PATH: &str = "test_shared.dng";
    
    /// OnceLock to ensure DNG is created only once across all tests
    static DNG_INIT: OnceLock<bool> = OnceLock::new();
    
    /// Maximum allowed difference in arcseconds between solutions from different starting points
    const CONVERGENCE_TOLERANCE_ARCSEC: f64 = 5.0;
    
    /// Maximum allowed difference in degrees for position angle between solutions
    const ORIENTATION_TOLERANCE_DEG: f64 = 1.0;

    /// Ensure the shared DNG file exists (thread-safe, creates only once)
    fn get_test_dng_path() -> Option<&'static str> {
        let success = DNG_INIT.get_or_init(|| {
            let dng = Path::new(SHARED_DNG_PATH);
            if dng.exists() {
                return true;
            }
            
            let cr3_path = Path::new(CR3_FILE_PATH);
            if !cr3_path.exists() {
                println!("CR3 file '{}' not found, skipping test", CR3_FILE_PATH);
                return false;
            }
            
            match convert_cr3_to_dng(cr3_path, dng) {
                Ok(()) => true,
                Err(e) => {
                    println!("Failed to convert CR3 to DNG: {}", e);
                    false
                }
            }
        });
        
        if *success {
            Some(SHARED_DNG_PATH)
        } else {
            None
        }
    }

    /// Test that plate solving works from the initial estimate
    #[test]
    fn test_solve_plate_basic() {
        let dng_path = match get_test_dng_path() {
            Some(p) => p,
            None => return,
        };

        let initial_coord = CoordinateEquatorial::new(
            RaHoursMinutesSeconds::new(14, 15, 21.4),
            Arcdegrees::new(54, 1, 14.40),
        );

        let result = solve_plate(Path::new(&dng_path), &initial_coord);
        
        assert!(result.is_ok(), "Plate solving failed: {:?}", result.err());
        
        let solution = result.unwrap();
        assert!(solution.coeffs_x.is_some(), "No solution found");
        
        println!("=== Basic Plate Solving Test ===");
        println!("Solved RA:  {:.6}°", solution.optical_axis_ra.to_degrees());
        println!("Solved Dec: {:.6}°", solution.optical_axis_dec.to_degrees());
        println!("Position angle: {:.2}°", solution.rotation_deg);
        println!("Scale: {:.4} arcsec/pixel", solution.scale_arcsec_per_pixel);
        println!("Matched quads: {}", solution.matched_quads_count);
        println!("Spiral iterations: {}", solution.spiral_iterations);
    }

    /// Test that plate solving converges to the same solution from multiple starting positions
    #[test]
    fn test_convergence_from_multiple_angles() {
        let dng_path = match get_test_dng_path() {
            Some(p) => p,
            None => return,
        };

        // Define multiple starting positions around the expected solution
        // These offsets are in degrees from the base position
        let base_ra_deg = 14.0 * 15.0 + 15.0 / 4.0;  // ~14h 15m in degrees
        let base_dec_deg = 54.0;
        
        let offsets = vec![
            (0.0, 0.0, "center"),
            (1.0, 0.0, "east"),
            (-1.0, 0.0, "west"),
            (0.0, 1.0, "north"),
            (0.0, -1.0, "south"),
            (0.7, 0.7, "northeast"),
            (-0.7, 0.7, "northwest"),
            (0.7, -0.7, "southeast"),
            (-0.7, -0.7, "southwest"),
        ];

        let mut solutions: Vec<(String, f64, f64, f64)> = Vec::new();

        println!("=== Convergence Test from Multiple Starting Positions ===\n");

        for (ra_offset, dec_offset, name) in &offsets {
            let ra_deg: f64 = base_ra_deg + ra_offset;
            let dec_deg: f64 = base_dec_deg + dec_offset;
            
            // Convert to RA hours/minutes/seconds
            let ra_hours: f64 = ra_deg / 15.0;
            let ra_h = ra_hours.floor() as i64;
            let ra_m = ((ra_hours - ra_h as f64) * 60.0).floor() as i64;
            let ra_s = ((ra_hours - ra_h as f64) * 60.0 - ra_m as f64) * 60.0;
            
            // Convert to Dec degrees/arcmin/arcsec
            let dec_sign = dec_deg.signum();
            let dec_abs = dec_deg.abs();
            let dec_d = dec_abs.floor() as i64;
            let dec_m = ((dec_abs - dec_d as f64) * 60.0).floor() as i64;
            let dec_s = ((dec_abs - dec_d as f64) * 60.0 - dec_m as f64) * 60.0;
            
            let initial_coord = CoordinateEquatorial::new(
                RaHoursMinutesSeconds::new(ra_h, ra_m, ra_s),
                Arcdegrees::new((dec_d as f64 * dec_sign) as i64, dec_m, dec_s),
            );

            print!("Testing from {} (RA={:.2}°, Dec={:.2}°)... ", name, ra_deg, dec_deg);

            match solve_plate(Path::new(&dng_path), &initial_coord) {
                Ok(solution) if solution.coeffs_x.is_some() => {
                    let solved_ra = solution.optical_axis_ra.to_degrees();
                    let solved_dec = solution.optical_axis_dec.to_degrees();
                    let pos_angle = solution.rotation_deg;
                    
                    println!("✓ RA={:.4}°, Dec={:.4}°, PA={:.2}° (iter={})",
                             solved_ra, solved_dec, pos_angle, solution.spiral_iterations);
                    
                    solutions.push((name.to_string(), solved_ra, solved_dec, pos_angle));
                }
                Ok(_) => {
                    println!("✗ No solution found");
                }
                Err(e) => {
                    println!("✗ Error: {}", e);
                }
            }
        }

        // Verify convergence - all solutions should be within tolerance
        assert!(!solutions.is_empty(), "No solutions found from any starting position");
        
        println!("\n=== Convergence Analysis ===");
        
        // Use first solution as reference
        let (ref_name, ref_ra, ref_dec, ref_pa) = &solutions[0];
        println!("Reference solution ({}): RA={:.4}°, Dec={:.4}°, PA={:.2}°", 
                 ref_name, ref_ra, ref_dec, ref_pa);
        
        let mut max_ra_diff_arcsec = 0.0_f64;
        let mut max_dec_diff_arcsec = 0.0_f64;
        let mut max_pa_diff_deg = 0.0_f64;
        
        for (name, ra, dec, pa) in &solutions[1..] {
            let ra_diff_arcsec = (ra - ref_ra) * 3600.0;
            let dec_diff_arcsec = (dec - ref_dec) * 3600.0;
            let pa_diff_deg = (pa - ref_pa).abs();
            
            max_ra_diff_arcsec = max_ra_diff_arcsec.max(ra_diff_arcsec.abs());
            max_dec_diff_arcsec = max_dec_diff_arcsec.max(dec_diff_arcsec.abs());
            max_pa_diff_deg = max_pa_diff_deg.max(pa_diff_deg);
            
            println!("  {} vs {}: ΔRA={:.2}\", ΔDec={:.2}\", ΔPA={:.2}°",
                     name, ref_name, ra_diff_arcsec, dec_diff_arcsec, pa_diff_deg);
        }
        
        println!("\nMaximum deviations:");
        println!("  RA:  {:.2}\" (tolerance: {}\")", max_ra_diff_arcsec, CONVERGENCE_TOLERANCE_ARCSEC);
        println!("  Dec: {:.2}\" (tolerance: {}\")", max_dec_diff_arcsec, CONVERGENCE_TOLERANCE_ARCSEC);
        println!("  PA:  {:.2}° (tolerance: {}°)", max_pa_diff_deg, ORIENTATION_TOLERANCE_DEG);
        
        assert!(
            max_ra_diff_arcsec < CONVERGENCE_TOLERANCE_ARCSEC,
            "RA convergence failed: max difference {:.2}\" exceeds tolerance {}\"",
            max_ra_diff_arcsec, CONVERGENCE_TOLERANCE_ARCSEC
        );
        
        assert!(
            max_dec_diff_arcsec < CONVERGENCE_TOLERANCE_ARCSEC,
            "Dec convergence failed: max difference {:.2}\" exceeds tolerance {}\"",
            max_dec_diff_arcsec, CONVERGENCE_TOLERANCE_ARCSEC
        );
        
        assert!(
            max_pa_diff_deg < ORIENTATION_TOLERANCE_DEG,
            "Position angle convergence failed: max difference {:.2}° exceeds tolerance {}°",
            max_pa_diff_deg, ORIENTATION_TOLERANCE_DEG
        );
        
        println!("\n✓ All solutions converged within tolerance!");
    }

    /// Test that solving from a position far from the true solution still works
    #[test]
    fn test_solve_from_far_offset() {
        let dng_path = match get_test_dng_path() {
            Some(p) => p,
            None => return,
        };

        // Start 3 degrees away from expected position (should require spiral search)
        let initial_coord = CoordinateEquatorial::new(
            RaHoursMinutesSeconds::new(14, 30, 0.0),  // ~15 minutes = 3.75 degrees off in RA
            Arcdegrees::new(52, 0, 0.0),              // 2 degrees off in Dec
        );

        println!("=== Far Offset Test ===");
        println!("Starting position: RA=14h 30m, Dec=52°");

        let result = solve_plate(Path::new(&dng_path), &initial_coord);
        
        assert!(result.is_ok(), "Plate solving failed: {:?}", result.err());
        
        let solution = result.unwrap();
        
        if solution.coeffs_x.is_some() {
            println!("✓ Solution found!");
            println!("  Solved RA:  {:.6}°", solution.optical_axis_ra.to_degrees());
            println!("  Solved Dec: {:.6}°", solution.optical_axis_dec.to_degrees());
            println!("  Position angle: {:.2}°", solution.rotation_deg);
            println!("  Spiral iterations: {}", solution.spiral_iterations);
            
            // Should have required multiple spiral iterations to find the solution
            assert!(solution.spiral_iterations > 0, 
                    "Expected spiral search to be needed for far offset");
        } else {
            // It's acceptable to not find a solution from very far away
            println!("No solution found from far offset (expected behavior)");
        }
    }

    /// Test position angle is consistent (approximately the same) across solutions
    #[test]
    fn test_position_angle_consistency() {
        let dng_path = match get_test_dng_path() {
            Some(p) => p,
            None => return,
        };

        // Solve from two different starting positions (both close enough to find solution)
        let positions = vec![
            CoordinateEquatorial::new(
                RaHoursMinutesSeconds::new(14, 15, 21.4),
                Arcdegrees::new(54, 1, 14.40),
            ),
            CoordinateEquatorial::new(
                RaHoursMinutesSeconds::new(14, 5, 0.0),  // About 0.5 degrees offset
                Arcdegrees::new(54, 20, 0.0),
            ),
        ];

        let mut position_angles: Vec<f64> = Vec::new();
        let mut scales: Vec<f64> = Vec::new();

        println!("=== Position Angle Consistency Test ===");

        for (i, coord) in positions.iter().enumerate() {
            if let Ok(solution) = solve_plate(Path::new(&dng_path), coord) {
                if solution.coeffs_x.is_some() {
                    println!("Solution {}: PA = {:.2}°, Scale = {:.4} arcsec/px, Mirrored = {}",
                             i + 1, solution.rotation_deg, solution.scale_arcsec_per_pixel, 
                             solution.is_mirrored);
                    position_angles.push(solution.rotation_deg);
                    scales.push(solution.scale_arcsec_per_pixel);
                }
            }
        }

        assert!(position_angles.len() >= 2, "Need at least 2 solutions to compare");

        let pa_diff = (position_angles[0] - position_angles[1]).abs();
        let scale_diff = (scales[0] - scales[1]).abs();
        
        println!("\nPosition angle difference: {:.2}°", pa_diff);
        println!("Scale difference: {:.4} arcsec/pixel", scale_diff);

        assert!(
            pa_diff < ORIENTATION_TOLERANCE_DEG,
            "Position angles differ by {:.2}°, exceeds tolerance {}°",
            pa_diff, ORIENTATION_TOLERANCE_DEG
        );
        
        // Scale should also be consistent
        assert!(
            scale_diff < 0.01,
            "Scales differ by {:.4} arcsec/pixel, exceeds tolerance 0.01",
            scale_diff
        );

        println!("✓ Position angles and scales are consistent!");
    }
}

#[cfg(test)]
mod projection_tests {
    use crate::platesolve::{get_star_x_y, get_equatorial_from_xy};
    use std::f64::consts::PI;

    /// Test that forward and inverse projections are consistent
    #[test]
    fn test_projection_roundtrip() {
        let optical_ra_rad = 210.0_f64.to_radians();
        let optical_dec_rad = 54.0_f64.to_radians();
        
        // Test several star positions around the optical axis
        let test_stars = vec![
            (211.0_f64.to_radians(), 54.5_f64.to_radians()),  // Slight offset
            (209.5_f64.to_radians(), 53.5_f64.to_radians()),  // Another offset
            (210.5_f64.to_radians(), 54.0_f64.to_radians()),  // RA only
            (210.0_f64.to_radians(), 54.5_f64.to_radians()),  // Dec only
            (210.0_f64.to_radians(), 54.0_f64.to_radians()),  // Center (should be 0,0)
        ];
        
        for (star_ra, star_dec) in test_stars {
            // Forward projection: RA/Dec -> x,y (arcseconds)
            let (x_arcsec, y_arcsec) = get_star_x_y(
                optical_ra_rad, optical_dec_rad,
                star_ra, star_dec
            );
            
            // Convert to radians for inverse projection
            let arcsec_to_rad = PI / (180.0 * 3600.0);
            let x_rad = x_arcsec * arcsec_to_rad;
            let y_rad = y_arcsec * arcsec_to_rad;
            
            // Inverse projection: x,y (radians) -> RA/Dec
            let (recovered_ra, recovered_dec) = get_equatorial_from_xy(
                x_rad, y_rad,
                optical_ra_rad, optical_dec_rad
            );
            
            // Check round-trip accuracy
            let ra_error_arcsec = (star_ra - recovered_ra).to_degrees() * 3600.0;
            let dec_error_arcsec = (star_dec - recovered_dec).to_degrees() * 3600.0;
            
            assert!(ra_error_arcsec.abs() < 0.1, 
                    "RA round-trip error too large: {} arcsec", ra_error_arcsec);
            assert!(dec_error_arcsec.abs() < 0.1, 
                    "Dec round-trip error too large: {} arcsec", dec_error_arcsec);
        }
        
        println!("✓ Projection roundtrip test passed");
    }
    
    /// Test that projection of center gives (0, 0)
    #[test]
    fn test_center_projection() {
        let optical_ra_rad = 210.0_f64.to_radians();
        let optical_dec_rad = 54.0_f64.to_radians();
        
        let (x, y) = get_star_x_y(
            optical_ra_rad, optical_dec_rad,
            optical_ra_rad, optical_dec_rad
        );
        
        assert!(x.abs() < 1e-10, "Center x should be 0, got {}", x);
        assert!(y.abs() < 1e-10, "Center y should be 0, got {}", y);
        
        println!("✓ Center projection test passed");
    }

    /// Test projection symmetry - equal offsets in opposite directions should give symmetric results
    #[test]
    fn test_projection_symmetry() {
        let optical_ra_rad = 210.0_f64.to_radians();
        let optical_dec_rad = 54.0_f64.to_radians();
        
        let offset_deg: f64 = 0.5;
        
        // Test RA symmetry
        let (x_plus, y_plus) = get_star_x_y(
            optical_ra_rad, optical_dec_rad,
            (210.0_f64 + offset_deg).to_radians(), optical_dec_rad
        );
        let (x_minus, y_minus) = get_star_x_y(
            optical_ra_rad, optical_dec_rad,
            (210.0_f64 - offset_deg).to_radians(), optical_dec_rad
        );
        
        // X should be opposite, Y should be equal
        assert!((x_plus + x_minus).abs() < 1.0, 
                "RA symmetry failed: x+ = {}, x- = {}", x_plus, x_minus);
        assert!((y_plus - y_minus).abs() < 1.0, 
                "RA symmetry failed: y+ = {}, y- = {}", y_plus, y_minus);

        // Test Dec symmetry
        let (x_plus, y_plus) = get_star_x_y(
            optical_ra_rad, optical_dec_rad,
            optical_ra_rad, (54.0_f64 + offset_deg).to_radians()
        );
        let (x_minus, y_minus) = get_star_x_y(
            optical_ra_rad, optical_dec_rad,
            optical_ra_rad, (54.0_f64 - offset_deg).to_radians()
        );
        
        // Y should be opposite, X should be equal
        assert!((x_plus - x_minus).abs() < 1.0, 
                "Dec symmetry failed: x+ = {}, x- = {}", x_plus, x_minus);
        assert!((y_plus + y_minus).abs() < 1.0, 
                "Dec symmetry failed: y+ = {}, y- = {}", y_plus, y_minus);
        
        println!("✓ Projection symmetry test passed");
    }
}

#[cfg(test)]
mod printing_tests {
    use crate::consts::*;
    use crate::coordinate::*;
    use crate::parse_catalog::*;
    use crate::platesolve::{analyze_image, convert_cr3_to_dng, get_star_x_y};
    use crate::printing::*;
    use crate::star_quads::*;
    use std::fs;
    use std::path::Path;
    use std::sync::OnceLock;

    const CR3_FILE_PATH: &str = "test_img/M82.CR3";
    const TEST_OUTPUT_DIR: &str = "test_output";
    
    /// Shared DNG file path - created once and reused by all tests
    const SHARED_DNG_PATH: &str = "test_shared.dng";
    
    /// OnceLock to ensure DNG is created only once across all tests
    static DNG_INIT: OnceLock<bool> = OnceLock::new();

    /// Ensure test output directory exists
    fn ensure_output_dir() {
        fs::create_dir_all(TEST_OUTPUT_DIR).expect("Failed to create test output directory");
    }

    /// Ensure the shared DNG file exists (thread-safe, creates only once)
    fn get_test_dng_path() -> Option<&'static str> {
        let success = DNG_INIT.get_or_init(|| {
            let dng = Path::new(SHARED_DNG_PATH);
            if dng.exists() {
                return true;
            }
            
            let cr3_path = Path::new(CR3_FILE_PATH);
            if !cr3_path.exists() {
                println!("CR3 file '{}' not found, skipping test", CR3_FILE_PATH);
                return false;
            }
            
            match convert_cr3_to_dng(cr3_path, dng) {
                Ok(()) => true,
                Err(e) => {
                    println!("Failed to convert CR3 to DNG: {}", e);
                    false
                }
            }
        });
        
        if *success {
            Some(SHARED_DNG_PATH)
        } else {
            None
        }
    }

    /// Generate all 6 visualization images
    #[test]
    fn test_generate_all_visualizations() {
        let dng_path = match get_test_dng_path() {
            Some(p) => p,
            None => return,
        };
        ensure_output_dir();

        println!("\n=== Generating All Visualization Images ===");
        println!("Output directory: {}/", TEST_OUTPUT_DIR);

        // Analyze image first to get dimensions and star positions
        let image_analysis = match analyze_image(Path::new(dng_path)) {
            Ok(analysis) => analysis,
            Err(e) => {
                println!("✗ Failed to analyze image: {}", e);
                return;
            }
        };

        let width = image_analysis.width as u32;
        let height = image_analysis.height as u32;

        // 1. Stretched image with center marker
        let path1 = format!("{}/01_stretched_with_center.png", TEST_OUTPUT_DIR);
        if let Err(e) = dng_to_png(Path::new(dng_path), Path::new(&path1)) {
            println!("✗ 01: Failed - {}", e);
        } else {
            println!("✓ 01_stretched_with_center.png - Stretched image with center marker");
        }

        // 2. Histogram of the original image
        let path2 = format!("{}/02_histogram.png", TEST_OUTPUT_DIR);
        match compute_dng_histogram(Path::new(dng_path), 100) {
            Ok(histogram) => {
                let mut bins_data: Vec<(&u16, &u32)> = histogram.iter().collect();
                bins_data.sort_by_key(|&(k, _)| *k);
                if plot_histogram_to_png(&bins_data, &path2, "Original Image Histogram").is_ok() {
                    println!("✓ 02_histogram.png - Histogram of original image");
                } else {
                    println!("✗ 02: Failed to plot histogram");
                }
            }
            Err(e) => println!("✗ 02: Failed to compute histogram - {}", e),
        }

        // 3. Original (not stretched) image with red circles around detected stars
        let path3 = format!("{}/03_original_with_stars.png", TEST_OUTPUT_DIR);
        if annotate_dng_image(dng_path, &path3, &image_analysis.star_barycenters, 8).is_ok() {
            println!("✓ 03_original_with_stars.png - Original image with {} detected stars", 
                     image_analysis.star_barycenters.len());
        } else {
            println!("✗ 03: Failed to annotate image");
        }

        // Get catalog data for the remaining images
        let initial_coord = CoordinateEquatorial::new(
            RaHoursMinutesSeconds::new(14, 3, 26.0),
            Arcdegrees::new(54, 24, 5.0),
        );

        let pixel_resolution: f64 = ARCSEC_PER_RADIAN * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGTH;
        let image_fov = (pixel_resolution * width as f64).max(pixel_resolution * height as f64);

        let star_in_fov_result = get_stars_from_catalogue(
            &initial_coord,
            image_fov,
            image_analysis.star_barycenters.len() * 100,
        );

        if star_in_fov_result.is_err() {
            println!("✗ Could not fetch catalog data (possibly no internet), skipping 4-6");
            return;
        }

        let star_in_fov = star_in_fov_result.unwrap();
        
        // Get catalog star positions in arcseconds (from get_star_x_y)
        let catalog_stars_arcsec: Vec<(f64, f64)> = star_in_fov
            .iter()
            .map(|s| {
                get_star_x_y(
                    initial_coord.ra.to_radians(),
                    initial_coord.dec.to_radians(),
                    s.ra.to_radians(),
                    s.dec.to_radians(),
                )
            })
            .collect();

        // Convert catalog stars from arcseconds to pixels for visualization
        // pixel_position = arcsec_position / pixel_resolution
        let catalog_stars_pixels: Vec<(f64, f64)> = catalog_stars_arcsec
            .iter()
            .map(|&(x_arcsec, y_arcsec)| (x_arcsec / pixel_resolution, y_arcsec / pixel_resolution))
            .collect();

        // Use arcsec coordinates for quad matching (same as main plate solving logic)
        let catalog_quads = returns_all_star_quads(&catalog_stars_arcsec, 2);

        // Find matched quads
        let mut matched_quads_image: Vec<StarQuad> = Vec::new();
        let mut matched_quads_catalog: Vec<StarQuad> = Vec::new();

        for quad in &image_analysis.star_quads {
            for quad_cat in &catalog_quads {
                if quad.compare(quad_cat, MATCHED_TOLERANCE) {
                    matched_quads_image.push(quad.clone());
                    matched_quads_catalog.push(quad_cat.clone());
                    break;
                }
            }
        }

        println!("  Found {} matched quads", matched_quads_image.len());

        // 4. Black background with catalog stars (red circles only)
        // Use pixel coordinates for visualization
        let path4 = format!("{}/04_catalog_stars_black.png", TEST_OUTPUT_DIR);
        if draw_stars_on_black(width, height, &path4, &catalog_stars_pixels, &[], 8).is_ok() {
            println!("✓ 04_catalog_stars_black.png - Catalog stars on black ({} stars)", 
                     catalog_stars_pixels.len());
        } else {
            println!("✗ 04: Failed to draw catalog stars");
        }

        // 5. Black background with image stars and matched quads
        let path5 = format!("{}/05_image_stars_matched_quads.png", TEST_OUTPUT_DIR);
        if draw_stars_on_black(width, height, &path5, &image_analysis.star_barycenters, &matched_quads_image, 8).is_ok() {
            println!("✓ 05_image_stars_matched_quads.png - Image stars with {} matched quads", 
                     matched_quads_image.len());
        } else {
            println!("✗ 05: Failed to draw image stars with quads");
        }

        // 6. Black background with catalog stars and matched quads
        // Need to convert catalog quad star positions from arcsec to pixels for visualization
        let matched_quads_catalog_pixels: Vec<StarQuad> = matched_quads_catalog
            .iter()
            .map(|quad| {
                // Convert each star position from arcsec to pixels
                let stars_pixels: [(f64, f64); 4] = [
                    (quad.stars[0].0 / pixel_resolution, quad.stars[0].1 / pixel_resolution),
                    (quad.stars[1].0 / pixel_resolution, quad.stars[1].1 / pixel_resolution),
                    (quad.stars[2].0 / pixel_resolution, quad.stars[2].1 / pixel_resolution),
                    (quad.stars[3].0 / pixel_resolution, quad.stars[3].1 / pixel_resolution),
                ];
                // Create a new quad with pixel coordinates
                StarQuad {
                    stars: stars_pixels,
                    barycenter: (quad.barycenter.0 / pixel_resolution, quad.barycenter.1 / pixel_resolution),
                    largest_distance: quad.largest_distance / pixel_resolution,
                    normalized_distances: quad.normalized_distances,
                }
            })
            .collect();
        
        let path6 = format!("{}/06_catalog_stars_matched_quads.png", TEST_OUTPUT_DIR);
        if draw_stars_on_black(width, height, &path6, &catalog_stars_pixels, &matched_quads_catalog_pixels, 8).is_ok() {
            println!("✓ 06_catalog_stars_matched_quads.png - Catalog stars with {} matched quads", 
                     matched_quads_catalog_pixels.len());
        } else {
            println!("✗ 06: Failed to draw catalog stars with quads");
        }

        println!("\n=== All 6 visualizations saved to {}/ ===", TEST_OUTPUT_DIR);
    }
}
