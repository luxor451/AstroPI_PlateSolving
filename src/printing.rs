use image::Rgb; 
use imageproc::drawing::{draw_hollow_circle_mut, draw_line_segment_mut};
use plotters::prelude::*;
use std::error::Error;
use std::path::Path;

use crate::star_quads::*;

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
