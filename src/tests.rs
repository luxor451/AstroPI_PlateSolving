#[cfg(test)]
mod png_generation_tests {
    use crate::coordinate::*;
    use crate::parse_catalog::*;
    use crate::printing::*;
    use crate::star_quads::*;
    use crate::platesolve::{analyze_image, get_star_x_y};
    use std::collections::HashMap;
    use std::path::Path;

    const TEST_FILE_PATH: &str = "test.dng";
    const MATCHED_TOLERANCE: f64 = 0.0009;

    #[test]
    fn test_generate_histogram_png() {
        let file_path = Path::new(TEST_FILE_PATH);
        if !file_path.exists() {
            println!("Test file '{}' not found, skipping test", TEST_FILE_PATH);
            return;
        }

        let result = analyze_image(file_path, false);
        assert!(
            result.is_ok(),
            "Failed to process DNG file: {:?}",
            result.err()
        );

        let image_analysis = result.unwrap();

        // Create histogram
        let mut histogram: HashMap<u16, u32> = HashMap::new();
        for &(_, _, pixel_value) in &image_analysis.centered_pixels {
            *histogram.entry(pixel_value).or_insert(0) += 1;
        }

        let mut bins_data: Vec<(&u16, &u32)> = histogram.iter().collect();
        bins_data.sort_by_key(|&(bin_start, _)| *bin_start);

        let result = plot_histogram_to_png(&bins_data, "test_histogram.png", "Pixel Value Histogram");
        assert!(result.is_ok(), "Failed to generate histogram PNG");
        
        // Verify file was created
        assert!(
            Path::new("test_histogram.png").exists(),
            "Histogram PNG file was not created"
        );
        println!("✓ Generated test_histogram.png");
    }

    #[test]
    fn test_generate_pixel_matrix_png() {
        let file_path = Path::new(TEST_FILE_PATH);
        if !file_path.exists() {
            println!("Test file '{}' not found, skipping test", TEST_FILE_PATH);
            return;
        }

        let result = analyze_image(file_path, false);
        assert!(
            result.is_ok(),
            "Failed to process DNG file: {:?}",
            result.err()
        );

        let image_analysis = result.unwrap();

        // Rebuild 2D matrix
        let mut pixel_matrix = vec![vec![0u16; image_analysis.width]; image_analysis.height];
        let cx = (image_analysis.width / 2) as i32;
        let cy = (image_analysis.height / 2) as i32;
        for &(x, y, val) in &image_analysis.centered_pixels {
            let c = (x + cx) as usize;
            let r = (y + cy) as usize;
            if r < image_analysis.height && c < image_analysis.width {
                pixel_matrix[r][c] = val;
            }
        }

        let result = save_pixel_matrix_to_png(&pixel_matrix, "test_pixel_matrix.png");
        assert!(result.is_ok(), "Failed to generate pixel matrix PNG");
        
        assert!(
            Path::new("test_pixel_matrix.png").exists(),
            "Pixel matrix PNG file was not created"
        );
        println!("✓ Generated test_pixel_matrix.png");
    }

    #[test]
    fn test_generate_annotated_stars_png() {
        let file_path = Path::new(TEST_FILE_PATH);
        if !file_path.exists() {
            println!("Test file '{}' not found, skipping test", TEST_FILE_PATH);
            return;
        }

        let result = analyze_image(file_path, false);
        assert!(
            result.is_ok(),
            "Failed to process DNG file: {:?}",
            result.err()
        );

        let image_analysis = result.unwrap();

        let result = annotate_stars_on_image(
            &image_analysis.centered_pixels,
            image_analysis.width,
            image_analysis.height,
            &image_analysis.star_barycenters,
            &image_analysis.star_quads,
            "test_annotated_stars.png",
            15,
        );
        assert!(result.is_ok(), "Failed to generate annotated stars PNG");
        
        assert!(
            Path::new("test_annotated_stars.png").exists(),
            "Annotated stars PNG file was not created"
        );
        println!("✓ Generated test_annotated_stars.png");
    }

    #[test]
    fn test_generate_catalog_comparison_png() {
        let file_path = Path::new(TEST_FILE_PATH);
        if !file_path.exists() {
            println!("Test file '{}' not found, skipping test", TEST_FILE_PATH);
            return;
        }

        let initial_coord = CoordinateEquatorial::new(
            RaHoursMinutesSeconds::new(14, 01, 12.5),
            Arcdegrees::new(54, 20, 56.0),
        );

        let result = analyze_image(file_path, false);
        if result.is_err() {
            println!("Could not process DNG file, skipping test");
            return;
        }

        let image_analysis = result.unwrap();

        let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT;
        let image_fov_x = pixel_resolution * image_analysis.width as f64;
        let image_fov_y = pixel_resolution * image_analysis.height as f64;
        let image_fov = image_fov_x.max(image_fov_y);

        // Try to get catalog stars (may fail if no internet connection)
        let star_in_fov_result = get_stars_from_catalogue(
            &initial_coord,
            image_fov * 1.0,
            image_analysis.star_barycenters.len() * 100,
        );

        if star_in_fov_result.is_err() {
            println!("Could not fetch catalog data (possibly no internet), skipping catalog comparison test");
            return;
        }

        let star_in_fov = star_in_fov_result.unwrap();
        let vec_star = star_in_fov
            .iter()
            .map(|s| {
                get_star_x_y(
                    initial_coord.ra.to_radians(),
                    initial_coord.dec.to_radians(),
                    s.ra.to_radians(),
                    s.dec.to_radians(),
                )
            })
            .collect::<Vec<(f64, f64)>>();

        let star_quad_from_cat = returns_all_star_quads(&vec_star, 2);

        // Find matched quads
        let mut matched_quads: Vec<(&StarQuad, &StarQuad)> = Vec::new();
        for quad in &image_analysis.star_quads {
            for quad_cat in &star_quad_from_cat {
                if quad.compare(quad_cat, MATCHED_TOLERANCE) {
                    matched_quads.push((quad, quad_cat));
                }
            }
        }

        let mut matched_quad_image: Vec<StarQuad> = Vec::new();
        let mut matched_quad_cat: Vec<StarQuad> = Vec::new();

        for elem in matched_quads.iter() {
            matched_quad_image.push(elem.0.clone());
            matched_quad_cat.push(elem.1.clone());
        }

        // Generate catalog visualization
        let result = annotate_stars_on_image(
            &[],
            (image_analysis.width as f64 * 1.0) as usize,
            (image_analysis.width as f64 * 1.0) as usize,
            &vec_star,
            &matched_quad_cat,
            "test_catalog_stars.png",
            15,
        );
        assert!(result.is_ok(), "Failed to generate catalog stars PNG");
        assert!(
            Path::new("test_catalog_stars.png").exists(),
            "Catalog stars PNG file was not created"
        );
        println!("✓ Generated test_catalog_stars.png");

        // Generate image visualization with matched quads
        let result = annotate_stars_on_image(
            &[],
            (image_analysis.width as f64 * 1.0) as usize,
            (image_analysis.width as f64 * 1.0) as usize,
            &image_analysis.star_barycenters,
            &matched_quad_image,
            "test_image_matched_quads.png",
            15,
        );
        assert!(result.is_ok(), "Failed to generate image matched quads PNG");
        assert!(
            Path::new("test_image_matched_quads.png").exists(),
            "Image matched quads PNG file was not created"
        );
        println!("✓ Generated test_image_matched_quads.png");
    }
}

