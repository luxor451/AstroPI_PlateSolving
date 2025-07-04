use std::path::Path;
use std::collections::HashMap; // Added VecDeque

use image::Rgb; // For defining the red color
use imageproc::drawing::{draw_hollow_circle_mut, draw_line_segment_mut}; // For drawing circles


mod star_quads;
use star_quads::*;


fn get_quad_from_file(file_path: &Path, verbose : bool) -> Result<Vec<StarQuad>, Box<dyn std::error::Error>> {
    // Get the pixel matrix from the DNG file
    let pixel_matrix = get_pixel_matrix_from_dng(file_path)?;

    println!("Nb of pixel above 16000: {}", pixel_matrix.iter().flatten().filter(|&&x| x > 16000).count());

    println!("Calculating star barycenters...");
    let star_barycenters = calculate_star_barycenters(&pixel_matrix);
    println!("Found {} star barycenters.", star_barycenters.len());


    println!("calculating star quads...");
    let star_quads = returns_all_star_quads(&star_barycenters, 3);
    println!("Found {} star quads.", star_quads.len());

    

    if verbose {
        let output_image_path = "output_image.png";

        // First, save the grayscale image from the pixel matrix.
        // This ensures "output_image.png" exists and has the up-to-date content from pixel_matrix.
        save_pixel_matrix_to_png(&pixel_matrix, output_image_path)?;

        if !star_barycenters.is_empty() {
            // Now, open the saved image, draw circles, and save it again.
            let dynamic_image = image::open(output_image_path)?;
            // Convert to Rgb8. The original image is Luma16. 
            // to_rgb8() will convert it to an 8-bit RGB image, suitable for drawing colored shapes.
            let mut img_rgb = dynamic_image.to_rgb8(); 
            
            let red_color = Rgb([255u8, 0u8, 0u8]);
            let circle_radius = 20; // Radius of the circle in pixels; adjust as needed for visibility.

            println!("Drawing circles for {} detected star(s) on {}.", star_barycenters.len(), output_image_path);
            for (bary_x, bary_y) in &star_barycenters {
                // Barycenter coordinates (x, y) from calculation are (column_index, row_index).
                // imageproc's draw_hollow_circle_mut expects (x, y) where x is horizontal (column) and y is vertical (row).
                let center_x = bary_x.round() as i32;
                let center_y = bary_y.round() as i32;
                
                draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red_color);
            }

            // Draw lines for star quads
            let white_color = Rgb([255u8, 255u8, 255u8]);
            println!("Drawing lines for {} star quad(s).", star_quads.len());
            for quad in &star_quads {
                if quad.stars.len() > 1 {
                    let central_star = quad.stars[0];
                    let central_point = (central_star.0 as f32, central_star.1 as f32);

                    for neighbor_star in quad.stars.iter().skip(1) {
                        let neighbor_point = (neighbor_star.0 as f32, neighbor_star.1 as f32);
                        draw_line_segment_mut(
                            &mut img_rgb,
                            central_point,
                            neighbor_point,
                            white_color,
                        );
                    }
                }
            }

            // Save the modified image, overwriting the previous grayscale version.
            img_rgb.save(output_image_path)?;
            println!("Successfully drew circles and quad lines on {} and re-saved.", output_image_path);
        } else {
            println!("No star barycenters to draw, so {} remains unchanged after initial save.", output_image_path);
        }



        // Create a histogram of pixel values
        let mut histogram: HashMap<u16, u32> = HashMap::new();
        for row in &pixel_matrix {
            for &pixel_value in row {
                *histogram.entry(pixel_value).or_insert(0) += 1;
            }
        }

        // Convert histogram to a sorted vector of tuples (bin_start, count)
        let mut bins_data: Vec<(&u16, &u32)> = histogram.iter().collect();
        bins_data.sort_by_key(|&(bin_start, _)| *bin_start);

        // Plot the histogram to a PNG file
        plot_histogram_to_png(&bins_data, "histogram.png", "Pixel Value Histogram")?;
    };

    return Ok(star_quads);
}


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let file_path = Path::new("test.dng"); // Replace with your DNG file path

    let (x, y) = get_image_size(file_path)?;
    
    let star_quads = get_quad_from_file(file_path, false)?;

    let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT; 

    println!("Image size: {}x{} pixels", x, y);
    println!("Pixel resolution: {:.6} arcsec", pixel_resolution);

    println!("Total FOV {:.6} arcmin by {:.6} arcmin", 
        (x as f64 * pixel_resolution / 60.0), 
        (y as f64 * pixel_resolution / 60.0));
    
    Ok(())
}