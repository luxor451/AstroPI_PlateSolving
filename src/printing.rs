
use image::{Rgb, RgbImage};
use imageproc::drawing::{draw_hollow_circle_mut, draw_line_segment_mut};
use plotters::prelude::*;
use std::error::Error;
use std::path::Path;

use crate::star_quads::*;

/// Options for rendering DNG images
#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub enum DngRenderMode {
    /// Histogram stretch with percentile clipping (default: 1% low, 99.5% high)
    Stretched { low_percentile: f64, high_percentile: f64 },
    /// Linear scaling (no stretch, just map raw range to 0-255)
    Linear,
    /// Gamma correction with specified gamma value
    Gamma(f64),
}

impl Default for DngRenderMode {
    fn default() -> Self {
        DngRenderMode::Stretched { low_percentile: 0.01, high_percentile: 0.995 }
    }
}

/// Loads a DNG file and converts it to an RGB image with the specified render mode.
/// This is the base function that all annotation functions should use.
fn load_dng_as_rgb(dng_path: &Path, mode: DngRenderMode) -> Result<RgbImage, Box<dyn Error>> {
    // Load the raw DNG file
    let raw_image = rawloader::decode_file(dng_path)?;

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

    // Calculate min/max for different modes
    let (black_point, white_point) = match mode {
        DngRenderMode::Stretched { low_percentile, high_percentile } => {
            let mut sorted_values = gray_values.clone();
            sorted_values.sort_unstable();
            
            let low_idx = ((sorted_values.len() as f64) * low_percentile) as usize;
            let high_idx = ((sorted_values.len() as f64) * high_percentile) as usize;
            let high_idx = high_idx.min(sorted_values.len() - 1);
            
            (sorted_values[low_idx] as f64, sorted_values[high_idx] as f64)
        }
        DngRenderMode::Linear | DngRenderMode::Gamma(_) => {
            let min_val = *gray_values.iter().min().unwrap_or(&0) as f64;
            let max_val = *gray_values.iter().max().unwrap_or(&65535) as f64;
            (min_val, max_val)
        }
    };

    let range = (white_point - black_point).max(1.0);

    let mut img = RgbImage::new(gray_width as u32, gray_height as u32);

    for (i, &value) in gray_values.iter().enumerate() {
        let x = (i % gray_width) as u32;
        let y = (i / gray_width) as u32;

        // Normalize to 0-1 range
        let normalized = ((value as f64 - black_point) / range).clamp(0.0, 1.0);

        // Apply gamma if needed
        let gamma_corrected = match mode {
            DngRenderMode::Gamma(gamma) => normalized.powf(1.0 / gamma),
            _ => normalized,
        };

        let pixel_value = (gamma_corrected * 255.0).clamp(0.0, 255.0) as u8;
        img.put_pixel(x, y, Rgb([pixel_value, pixel_value, pixel_value]));
    }

    Ok(img)
}
#[allow(dead_code)]
/// Returns the dimensions (width, height) of a DNG image after debayering.
pub fn get_dng_dimensions(dng_path: &Path) -> Result<(usize, usize), Box<dyn Error>> {
    let raw_image = rawloader::decode_file(dng_path)?;
    Ok((raw_image.width / 2, raw_image.height / 2))
}
#[allow(dead_code)]
pub fn plot_histogram_to_png(
    bins_data: &Vec<(&u16, &u32)>,
    output_path: &str,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if bins_data.is_empty() {
        eprintln!("No data to plot for histogram.");
        return Ok(());
    }

    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_count = bins_data
        .iter()
        .map(|&(_, count)| *count)
        .max()
        .unwrap_or(1) as f64;
    // Determine the max bin start value for x-axis range
    let max_bin_start = bins_data.last().map_or(16300, |&(bin_val, _)| *bin_val);

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(80)
        .build_cartesian_2d(0u32..(max_bin_start as u32 + 100u32), 0f64..max_count * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Pixel Value Bins")
        .y_desc("Count")
        .axis_desc_style(("sans-serif", 15))
        .draw()?;

    chart.draw_series(bins_data.iter().map(|&(bin_start, count)| {
        let x0 = *bin_start as u32;
        let x1 = (bin_start + 99) as u32; // Bin width is 100
        let y = *count as f64;
        let mut rect = Rectangle::new([(x0, 0.0), (x1, y)], BLUE.filled());
        rect.set_margin(0, 0, 1, 1); // Add small margin between bars
        rect
    }))?;

    root.present()?;
    // println!("Histogram saved to {}", output_path);
    Ok(())
}
#[allow(dead_code)]
pub fn save_pixel_matrix_to_png(
    matrix: &Vec<Vec<u16>>,
    output_path: &str,
) -> Result<(), image::ImageError> {
    if matrix.is_empty() {
        return Err(image::ImageError::Parameter(
            image::error::ParameterError::from_kind(
                image::error::ParameterErrorKind::DimensionMismatch,
            ),
        ));
    }

    let height = matrix.len();
    let width = matrix[0].len();

    if width == 0 {
        return Err(image::ImageError::Parameter(
            image::error::ParameterError::from_kind(
                image::error::ParameterErrorKind::DimensionMismatch,
            ),
        ));
    }

    // Flatten the Vec<Vec<u16>> into Vec<u16> for ImageBuffer.
    let mut flat_data: Vec<u16> = Vec::with_capacity(width * height);
    for row in matrix {
        if row.len() != width {
            return Err(image::ImageError::Parameter(
                image::error::ParameterError::from_kind(
                    image::error::ParameterErrorKind::DimensionMismatch,
                ),
            ));
        }
        flat_data.extend_from_slice(row);
    }

    // Stretch pixel values by a factor of 4 and clamp to u16::MAX.
    for pixel_value in flat_data.iter_mut() {
        *pixel_value = (*pixel_value as u32 * 4).min(u16::MAX as u32) as u16;
    }

    // Create an ImageBuffer for Luma<u16> (16-bit grayscale).
    // `ImageBuffer::from_raw` takes ownership of `flat_data`.
    // It returns Some(image_buffer) if successful, None if data length doesn't match dimensions.
    let image_buffer = match image::ImageBuffer::<image::Luma<u16>, Vec<u16>>::from_raw(
        width as u32,
        height as u32,
        flat_data,
    ) {
        Some(buffer) => buffer,
        None => {
            // This case implies flat_data.len() != width * height * num_channels,
            // which should be caught by earlier checks or indicate an internal logic error.
            return Err(image::ImageError::Parameter(
                image::error::ParameterError::from_kind(
                    image::error::ParameterErrorKind::DimensionMismatch, // Data length issue
                ),
            ));
        }
    };

    // Save the image buffer to a PNG file.
    // The image crate infers the format from the file extension.
    image_buffer.save(output_path)?;

    // println!("PNG image saved to {}", output_path);
    Ok(())
}
#[allow(dead_code)]
/// Annotates an existing image file with star positions and quad connections.
/// 
/// Opens an existing PNG file and overlays red circles at each star position
/// and white lines connecting stars in quads. The star positions should be
/// in centered coordinates (relative to image center).
///
/// # Arguments
///
/// * `input_image_path` - Path to the existing image file to annotate
/// * `output_path` - Path where the annotated image will be saved
/// * `star_positions` - Star positions in centered coordinates (relative to image center)
/// * `star_quads` - Star quads with connection information
/// * `circle_radius` - Radius of circles to draw around stars
pub fn annotate_existing_image<P: AsRef<Path>, Q: AsRef<Path>>(
    input_image_path: P,
    output_path: Q,
    star_positions: &[StarPosXY],
    star_quads: &[StarQuad],
    circle_radius: i32,
) -> Result<(), Box<dyn Error>> {
    let mut img_rgb = image::open(&input_image_path)?.to_rgb8();
    let width = img_rgb.width() as i32;
    let height = img_rgb.height() as i32;
    let cx = width / 2;
    let cy = height / 2;

    // Draw circles at star positions
    let red = Rgb([255, 0, 0]);
    for &(bx, by) in star_positions {
        let center_x = bx.round() as i32 + cx;
        let center_y = by.round() as i32 + cy;

        if center_x >= 0 && center_x < width && center_y >= 0 && center_y < height {
            draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red);
        }
    }

    // Draw lines for star quads
    let white_color = Rgb([255u8, 255u8, 255u8]);
    for quad in star_quads {
        if quad.stars.len() > 1 {
            let central_star = quad.stars[0];
            let central_point = (
                (central_star.0 as i32 + cx) as f32,
                (central_star.1 as i32 + cy) as f32,
            );

            for neighbor_star in quad.stars.iter().skip(1) {
                let neighbor_point = (
                    (neighbor_star.0 as i32 + cx) as f32,
                    (neighbor_star.1 as i32 + cy) as f32,
                );
                draw_line_segment_mut(&mut img_rgb, central_point, neighbor_point, white_color);
            }
        }
    }

    img_rgb.save(&output_path)?;
    Ok(())
}
#[allow(dead_code)]
/// Reconstructs a grayscale PNG from your centered pixels, then overlays
/// red circles at each `(x,y)` barycenter, and re­saves it.
pub fn annotate_stars_on_image<P: AsRef<Path>>(
    centered_pixels: &[(i32, i32, u16)],
    width: usize,
    height: usize,
    star_positions: &[StarPosXY],
    star_quads: &[StarQuad],
    output_path: P,
    circle_radius: i32,
) -> Result<(), Box<dyn Error>> {
    // 1) rebuild a 2D u16 matrix
    let mut pixel_matrix = vec![vec![0u16; width]; height];
    let cx = (width / 2) as i32;
    let cy = (height / 2) as i32;
    for &(x, y, val) in centered_pixels {
        let c = (x + cx) as usize;
        let r = (y + cy) as usize;
        if r < height && c < width {
            pixel_matrix[r][c] = val;
        }
    }

    // 2) save as grayscale
    save_pixel_matrix_to_png(
        &pixel_matrix,
        output_path.as_ref().to_str().ok_or("Invalid output path")?,
    )?;

    // 3) reopen, draw circles, re-save
    if !star_positions.is_empty() {
        let mut img_rgb = image::open(&output_path)?.to_rgb8();
        let red = Rgb([255, 0, 0]);
        for &(bx, by) in star_positions {
            let center_x = bx.round() as i32 + cx;
            let center_y = by.round() as i32 + cy;

            // println!("Drawing circle at ({}, {}) with radius {}", center_x, center_y, circle_radius);

            draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red);
        }

        // Draw lines for star quads
        let white_color = Rgb([255u8, 255u8, 255u8]);
        // println!("Drawing lines for {} star quad(s).", star_quads.len());
        for quad in star_quads {
            if quad.stars.len() > 1 {
                let central_star = quad.stars[0];
                // Convert centered star coordinates to image coordinates
                let central_point = (
                    (central_star.0 as i32 + cx) as f32,
                    (central_star.1 as i32 + cy) as f32,
                );

                for neighbor_star in quad.stars.iter().skip(1) {
                    let neighbor_point = (
                        (neighbor_star.0 as i32 + cx) as f32,
                        (neighbor_star.1 as i32 + cy) as f32,
                    );
                    draw_line_segment_mut(&mut img_rgb, central_point, neighbor_point, white_color);
                }
            }
        }

        img_rgb.save(&output_path)?;
        // println!("Successfully drew circles and quad lines on {} and re-saved.", output_path.as_ref().display());
    }

    Ok(())
}

#[allow(dead_code)]
/// Annotates a DNG image directly with star positions (red circles only, no quads).
/// 
/// Loads the DNG file with linear scaling and overlays red circles at star positions.
/// Note: The DNG is debayered (2x2 averaging) so the output image is half the raw dimensions.
/// Star positions should be in raw image centered coordinates - they will be scaled by 0.5.
///
/// # Arguments
///
/// * `dng_path` - Path to the DNG file to annotate
/// * `output_path` - Path where the annotated image will be saved
/// * `star_positions` - Star positions in centered coordinates (relative to raw image center)
/// * `circle_radius` - Radius of circles to draw around stars
pub fn annotate_dng_image<P: AsRef<Path>, Q: AsRef<Path>>(
    dng_path: P,
    output_path: Q,
    star_positions: &[StarPosXY],
    circle_radius: i32,
) -> Result<(), Box<dyn Error>> {
    // Use Linear mode to show the original image without histogram stretching
    let mut img_rgb = load_dng_as_rgb(dng_path.as_ref(), DngRenderMode::Linear)?;
    let width = img_rgb.width() as i32;
    let height = img_rgb.height() as i32;
    let cx = width / 2;
    let cy = height / 2;

    // Draw circles at star positions
    // Note: Star positions are in raw image coordinates (before debayering).
    // The debayered image is half the size, so we scale positions by 0.5.
    let red = Rgb([255, 0, 0]);
    for &(bx, by) in star_positions {
        // Scale from raw coordinates to debayered coordinates
        let scaled_x = bx * 0.5;
        let scaled_y = by * 0.5;
        let center_x = scaled_x.round() as i32 + cx;
        let center_y = scaled_y.round() as i32 + cy;

        if center_x >= 0 && center_x < width && center_y >= 0 && center_y < height {
            draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red);
        }
    }

    img_rgb.save(output_path.as_ref())?;
    Ok(())
}
#[allow(dead_code)]
/// Draws stars and quads on a black background.
///
/// # Arguments
///
/// * `width` - Image width in pixels
/// * `height` - Image height in pixels
/// * `output_path` - Path where the image will be saved
/// * `star_positions` - Star positions in centered coordinates (relative to image center)
/// * `star_quads` - Star quads with connection information (can be empty)
/// * `circle_radius` - Radius of circles to draw around stars
pub fn draw_stars_on_black<Q: AsRef<Path>>(
    width: u32,
    height: u32,
    output_path: Q,
    star_positions: &[StarPosXY],
    star_quads: &[StarQuad],
    circle_radius: i32,
) -> Result<(), Box<dyn Error>> {
    let mut img_rgb = RgbImage::new(width, height);
    // Image is already black (all zeros)
    
    let cx = (width / 2) as i32;
    let cy = (height / 2) as i32;

    // Draw circles at star positions
    let red = Rgb([255, 0, 0]);
    for &(bx, by) in star_positions {
        let center_x = bx.round() as i32 + cx;
        let center_y = by.round() as i32 + cy;

        if center_x >= 0 && center_x < width as i32 && center_y >= 0 && center_y < height as i32 {
            draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red);
        }
    }

    // Draw lines for star quads
    let white_color = Rgb([255u8, 255u8, 255u8]);
    for quad in star_quads {
        if quad.stars.len() > 1 {
            let central_star = quad.stars[0];
            let central_point = (
                (central_star.0 as i32 + cx) as f32,
                (central_star.1 as i32 + cy) as f32,
            );

            for neighbor_star in quad.stars.iter().skip(1) {
                let neighbor_point = (
                    (neighbor_star.0 as i32 + cx) as f32,
                    (neighbor_star.1 as i32 + cy) as f32,
                );
                draw_line_segment_mut(&mut img_rgb, central_point, neighbor_point, white_color);
            }
        }
    }

    img_rgb.save(output_path.as_ref())?;
    Ok(())
}
#[allow(dead_code)]
/// Computes histogram bins from a DNG file.
///
/// Returns a HashMap where keys are bin starts (multiples of bin_size) and values are counts.
pub fn compute_dng_histogram(dng_path: &Path, bin_size: u16) -> Result<std::collections::HashMap<u16, u32>, Box<dyn Error>> {
    use std::collections::HashMap;
    
    let raw_image = rawloader::decode_file(dng_path)?;
    let data = match raw_image.data {
        rawloader::RawImageData::Integer(ref data) => data,
        rawloader::RawImageData::Float(_) => {
            return Err("Float raw data not supported".into());
        }
    };

    let mut histogram: HashMap<u16, u32> = HashMap::new();
    for &value in data.iter() {
        let bin = (value / bin_size) * bin_size;
        *histogram.entry(bin).or_insert(0) += 1;
    }

    Ok(histogram)
}

/// Converts a DNG file to PNG with the specified render mode.
///
/// # Arguments
///
/// * `input_path` - Path to the input DNG file
/// * `output_path` - Path for the output PNG file
/// * `mode` - Render mode (stretched, linear, gamma)
///
/// # Returns
///
/// Result indicating success or an error
pub fn dng_to_png_with_mode(
    input_path: &Path,
    output_path: &Path,
    mode: DngRenderMode,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut img = load_dng_as_rgb(input_path, mode)?;
    
    // Draw a crosshair marker at the center of the image
    let center_x = img.width() / 2;
    let center_y = img.height() / 2;
    let marker_size = 20u32;
    let marker_color = Rgb([255u8, 0u8, 0u8]); // Bright red marker

    // Draw horizontal line of the crosshair
    for dx in 0..=marker_size {
        if center_x >= dx {
            img.put_pixel(center_x - dx, center_y, marker_color);
        }
        if center_x + dx < img.width() {
            img.put_pixel(center_x + dx, center_y, marker_color);
        }
    }

    // Draw vertical line of the crosshair
    for dy in 0..=marker_size {
        if center_y >= dy {
            img.put_pixel(center_x, center_y - dy, marker_color);
        }
        if center_y + dy < img.height() {
            img.put_pixel(center_x, center_y + dy, marker_color);
        }
    }

    // Draw a small circle around the center for better visibility
    let circle_radius = 10u32;
    for angle in 0..360 {
        let rad = (angle as f64).to_radians();
        let cx = (center_x as f64 + circle_radius as f64 * rad.cos()).round() as u32;
        let cy = (center_y as f64 + circle_radius as f64 * rad.sin()).round() as u32;
        if cx < img.width() && cy < img.height() {
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
pub fn dng_to_png(input_path: &Path, output_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    dng_to_png_with_mode(input_path, output_path, DngRenderMode::default())
}