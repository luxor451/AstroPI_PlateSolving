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

/// Converts a DNG file to PNG with histogram stretching for better visibility.
///
/// The function applies a percentile-based stretch to enhance contrast,
/// mapping the lower percentile to black and the upper percentile to white.
///
/// # Arguments
///
/// * `input_path` - Path to the input DNG file
/// * `output_path` - Path for the output PNG file
/// * `low_percentile` - Lower percentile for black point (e.g., 0.01 for 1%)
/// * `high_percentile` - Upper percentile for white point (e.g., 0.99 for 99%)
///
/// # Returns
///
/// Result indicating success or an error
fn dng_to_png_stretched(
    input_path: &Path,
    output_path: &Path,
    low_percentile: f64,
    high_percentile: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    use image::{Rgb, RgbImage};

    // Load the raw DNG file
    let raw_image = rawloader::decode_file(input_path)?;

    // Get the raw data as u16
    let data = match raw_image.data {
        rawloader::RawImageData::Integer(ref data) => data,
        rawloader::RawImageData::Float(_) => {
            return Err("Float raw data not supported".into());
        }
    };

    let width = raw_image.width;
    let height = raw_image.height;

    // Convert to grayscale by simple debayering (averaging 2x2 blocks)
    // This is a simple approach - for better quality, use proper demosaicing
    let gray_width = width / 2;
    let gray_height = height / 2;

    let mut gray_values: Vec<u16> = Vec::with_capacity(gray_width * gray_height);

    for y in 0..gray_height {
        for x in 0..gray_width {
            let bayer_x = x * 2;
            let bayer_y = y * 2;

            // Average the 2x2 Bayer block
            let p00 = data[bayer_y * width + bayer_x] as u32;
            let p01 = data[bayer_y * width + bayer_x + 1] as u32;
            let p10 = data[(bayer_y + 1) * width + bayer_x] as u32;
            let p11 = data[(bayer_y + 1) * width + bayer_x + 1] as u32;

            let avg = ((p00 + p01 + p10 + p11) / 4) as u16;
            gray_values.push(avg);
        }
    }

    // Calculate histogram for percentile stretching
    let mut sorted_values = gray_values.clone();
    sorted_values.sort_unstable();

    let low_idx = ((sorted_values.len() as f64) * low_percentile) as usize;
    let high_idx = ((sorted_values.len() as f64) * high_percentile) as usize;
    let high_idx = high_idx.min(sorted_values.len() - 1);

    let black_point = sorted_values[low_idx] as f64;
    let white_point = sorted_values[high_idx] as f64;

    // Apply histogram stretch and convert to 8-bit
    let range = (white_point - black_point).max(1.0);

    let mut img = RgbImage::new(gray_width as u32, gray_height as u32);

    for (i, &value) in gray_values.iter().enumerate() {
        let x = (i % gray_width) as u32;
        let y = (i / gray_width) as u32;

        // Apply stretch: map [black_point, white_point] to [0, 255]
        let stretched = ((value as f64 - black_point) / range * 255.0)
            .clamp(0.0, 255.0) as u8;

        img.put_pixel(x, y, Rgb([stretched, stretched, stretched]));
    }

    // Draw a crosshair marker at the center of the image
    let center_x = gray_width as u32 / 2;
    let center_y = gray_height as u32 / 2;
    let marker_size = 20u32;
    let marker_color = Rgb([255u8, 0u8, 0u8]); // Bright red marker

    // Draw horizontal line of the crosshair
    for dx in 0..=marker_size {
        if center_x >= dx {
            img.put_pixel(center_x - dx, center_y, marker_color);
        }
        if center_x + dx < gray_width as u32 {
            img.put_pixel(center_x + dx, center_y, marker_color);
        }
    }

    // Draw vertical line of the crosshair
    for dy in 0..=marker_size {
        if center_y >= dy {
            img.put_pixel(center_x, center_y - dy, marker_color);
        }
        if center_y + dy < gray_height as u32 {
            img.put_pixel(center_x, center_y + dy, marker_color);
        }
    }

    // Draw a small circle around the center for better visibility
    let circle_radius = 10u32;
    for angle in 0..360 {
        let rad = (angle as f64).to_radians();
        let cx = (center_x as f64 + circle_radius as f64 * rad.cos()).round() as u32;
        let cy = (center_y as f64 + circle_radius as f64 * rad.sin()).round() as u32;
        if cx < gray_width as u32 && cy < gray_height as u32 {
            img.put_pixel(cx, cy, marker_color);
        }
    }

    // Save as PNG
    img.save(output_path)?;

    Ok(())
}

/// Converts a DNG file to PNG with default histogram stretching parameters.
///
/// Uses 1% and 99.5% percentiles for the stretch.
#[allow(dead_code)]
fn dng_to_png(input_path: &Path, output_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    dng_to_png_stretched(input_path, output_path, 0.01, 0.995)
}

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

    // Optionally export the DNG as a stretched PNG for visualization
    let png_output = Path::new("test_stretched.png");
    let _ = dng_to_png(file_path, png_output);

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