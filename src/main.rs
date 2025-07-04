use std::path::Path;
use std::collections::HashMap; // Added VecDeque

use image::Rgb; // For defining the red color
use imageproc::drawing::{draw_hollow_circle_mut, draw_line_segment_mut}; // For drawing circles


mod star_quads;
mod parse_catalog;
use star_quads::*;
use parse_catalog::*; 



struct arcdegrees {
    pub degrees: i64, // in degrees
    pub arcminutes: i64, // in minutes
    pub arcseconds: f64, // in seconds
}


struct ra_hours_minutes_seconds {
    pub hours: i64, // in hours
    pub minutes: i64, // in minutes
    pub seconds: f64, // in seconds
}

struct coordinate_equatorial {
    pub ra: ra_hours_minutes_seconds, // in degrees
    pub dec: arcdegrees, // in degrees
}


struct coordinate_x_y {
    pub x: f64, 
    pub y: f64, 
}

impl arcdegrees {
    pub fn new(degrees: i64, arcminutes: i64, arcseconds: f64) -> Self {
        arcdegrees {
            degrees,
            arcminutes,
            arcseconds,
        }
    }

    pub fn to_degrees(&self) -> f64 {
        self.degrees as f64 + self.arcminutes as f64 / 60.0 + self.arcseconds / 3600.0
    }
}

impl ra_hours_minutes_seconds {
    pub fn new(hours: i64, minutes: i64, seconds: f64) -> Self {
        ra_hours_minutes_seconds {
            hours,
            minutes,
            seconds,
        }
    }

    pub fn to_arcdegrees(&self) -> arcdegrees {
        let total_seconds = (self.hours * 3600 + self.minutes * 60 + self.seconds as i64) as f64;
        let degrees = total_seconds / 240.0; // 1 hour = 15 degrees, so 1 second = 15/3600 degrees
        let arcminutes = (degrees.fract() * 60.0) as i64;
        let arcseconds = (degrees.fract() * 3600.0).fract() * 60.0;

        arcdegrees::new(degrees as i64, arcminutes, arcseconds)
    }

    pub fn to_degrees(&self) -> f64 {
        let total_seconds = (self.hours * 3600 + self.minutes * 60) as f64 + self.seconds * 60.0;
        total_seconds / 240.0 // 1 hour = 15 degrees, so 1 second = 15/3600 degrees
    }
}

impl coordinate_equatorial {
    pub fn new(ra: ra_hours_minutes_seconds, dec: arcdegrees) -> Self {
        coordinate_equatorial { ra, dec }
    }
}




fn get_quad_from_file(file_path: &Path, verbose : bool) -> Result<Vec<StarQuad>, Box<dyn std::error::Error>> {
    // Get image dimensions first
    let (width, height) = get_image_size(file_path)?;

    // Get the pixel matrix from the DNG file with centered coordinates
    let centered_pixels = get_pixel_matrix_from_dng(file_path)?;

    if verbose { 
        println!("Nb of pixel above 16000: {}", centered_pixels.iter().filter(|&&(_, _, val)| val > 16000).count());
    }
    
    if verbose { println!("Calculating star barycenters...")};
    let star_barycenters = calculate_star_barycenters(&centered_pixels, width, height);
    if verbose { println!("Found {} star barycenters.", star_barycenters.len())};


    if verbose { println!("calculating star quads...") };
    let star_quads = returns_all_star_quads(&star_barycenters, 3);
    if verbose {println!("Found {} star quads.", star_quads.len())}

    

    if verbose {
        let output_image_path = "output_image.png";

        // Reconstruct the 2D pixel matrix for saving and histogram
        let mut pixel_matrix = vec![vec![0u16; width]; height];
        let cx = (width / 2) as i32;
        let cy = (height / 2) as i32;
        for &(x, y, val) in &centered_pixels {
            let c = (x + cx) as usize;
            let r = (y + cy) as usize;
            if r < height && c < width {
                pixel_matrix[r][c] = val;
            }
        }

        // First, save the grayscale image from the pixel matrix.
        save_pixel_matrix_to_png(&pixel_matrix, output_image_path)?;

        if !star_barycenters.is_empty() {
            // Now, open the saved image, draw circles, and save it again.
            let dynamic_image = image::open(output_image_path)?;
            let mut img_rgb = dynamic_image.to_rgb8(); 
            
            let red_color = Rgb([255u8, 0u8, 0u8]);
            let circle_radius = 20; // Radius of the circle in pixels; adjust as needed for visibility.

            println!("Drawing circles for {} detected star(s) on {}.", star_barycenters.len(), output_image_path);
            for (bary_x, bary_y) in &star_barycenters {
                // Convert centered barycenter coordinates to image coordinates (top-left origin)
                let center_x = (bary_x.round() as i32) + cx;
                let center_y = (bary_y.round() as i32) + cy;
                
                draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red_color);
            }

            // Draw lines for star quads
            let white_color = Rgb([255u8, 255u8, 255u8]);
            println!("Drawing lines for {} star quad(s).", star_quads.len());
            for quad in &star_quads {
                if quad.stars.len() > 1 {
                    let central_star = quad.stars[0];
                    // Convert centered star coordinates to image coordinates
                    let central_point = ((central_star.0 as i32 + cx) as f32, (central_star.1 as i32 + cy) as f32);

                    for neighbor_star in quad.stars.iter().skip(1) {
                        let neighbor_point = ((neighbor_star.0 as i32 + cx) as f32, (neighbor_star.1 as i32 + cy) as f32);
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
        for &(_, _, pixel_value) in &centered_pixels {
            *histogram.entry(pixel_value).or_insert(0) += 1;
        }

        // Convert histogram to a sorted vector of tuples (bin_start, count)
        let mut bins_data: Vec<(&u16, &u32)> = histogram.iter().collect();
        bins_data.sort_by_key(|&(bin_start, _)| *bin_start);

        // Plot the histogram to a PNG file
        plot_histogram_to_png(&bins_data, "histogram.png", "Pixel Value Histogram")?;
    };

    return Ok(star_quads);
}

fn get_star_x_y(optical_axis_right_ascension_radians: f64, 
                optical_axis_declination_radians : f64, 
                star_right_ascension_radians: f64,
                star_declination_radians: f64 
                ) -> coordinate_x_y {

    let sin_optical_axis_dec = optical_axis_declination_radians.sin();
    let cos_optical_axis_dec = optical_axis_declination_radians.cos();

    let sin_star_dec = star_declination_radians.sin();
    let cos_star_dec = star_declination_radians.cos();

    let ra_delta = star_right_ascension_radians - optical_axis_right_ascension_radians;
    let sin_delta_ra = ra_delta.sin();
    let cos_delta_ra = ra_delta.cos();

    let delta_value = (cos_optical_axis_dec * cos_star_dec * cos_delta_ra 
                            + sin_optical_axis_dec * sin_star_dec) / (3600.0 * (180.0 / std::f64::consts::PI));

    let star_x = - cos_star_dec * sin_delta_ra / delta_value;
    let star_y = - (sin_optical_axis_dec * cos_star_dec  * cos_delta_ra - cos_optical_axis_dec * sin_star_dec) / delta_value;

    coordinate_x_y { x : star_x, y : star_y }
}


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let initial_coord = coordinate_equatorial::new( 
        ra_hours_minutes_seconds::new(14, 01, 12.5), 
        arcdegrees::new(54, 20, 56.0)
    );

    let file_path = Path::new("test.dng"); // Replace with your DNG file path

    let (images_x_size, image_y_size) = get_image_size(file_path)?;
    
    let star_quads = get_quad_from_file(file_path, true)?;

    let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT; 

    println!("Image size: {}x{} pixels", images_x_size, image_y_size);
    println!("Pixel resolution: {:.6} arcsec", pixel_resolution);

    let image_resolution_x = pixel_resolution * images_x_size as f64; // Convert to arcsec
    let image_resolution_y = pixel_resolution * image_y_size as f64; // Convert to arcsec

    println!("Total FOV {:.6} arcmin by {:.6} arcmin", 
        image_resolution_x, 
        image_resolution_y);



    let initial_coord_ra_deg = initial_coord.ra.to_degrees();
    let initial_coord_dec_deg = initial_coord.dec.to_degrees();

    let star_in_fov = find_stars_in_fov(
        "catalog.dat", // Replace with your catalog file path
        initial_coord_ra_deg,
        initial_coord_dec_deg,
        image_resolution_x, // Convert arcmin to arcsec
    )?;

    println!("Found {} stars in the FOV.", star_in_fov.len());
    for star in &star_in_fov {
        let star_x_y = get_star_x_y(
            initial_coord_ra_deg.to_radians(),
            initial_coord_dec_deg.to_radians(),
            star.ra.to_radians(),
            star.dec.to_radians(),
        );
        println!("star x : {:.6}, star y: {:.6}", star_x_y.x, star_x_y.y);
    }
    
    Ok(())
}