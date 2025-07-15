use std::path::Path;
use std::collections::HashMap; 

mod star_quads;
mod parse_catalog;
mod coordinate;
mod printing;
use star_quads::*;
use parse_catalog::*; 
use coordinate::*;
use printing::*;

const MATCHED_TOLERANCE: f64 = 0.0009; // Tolerance for matching star quads


fn get_quad_from_file(file_path: &Path, verbose : bool) -> Result<(Vec<StarQuad>, Vec<(i32, i32, u16)>, Vec<(f64, f64)>), Box<dyn std::error::Error>> {
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
    let star_quads = returns_all_star_quads(&star_barycenters, 1);
    if verbose {println!("Found {} star quads.", star_quads.len())}

    

    if verbose {
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
    }

    return Ok((star_quads, centered_pixels, star_barycenters));
}

fn get_equatorial_from_xy(
    x: f64, 
    y: f64, 
    optical_axis_ra_rad: f64, 
    optical_axis_dec_rad: f64
) -> (f64, f64) {
    let sin_dec0 = optical_axis_dec_rad.sin();
    let cos_dec0 = optical_axis_dec_rad.cos();

    let rho = (x.powi(2) + y.powi(2)).sqrt();
    let c = rho.atan();
    let sin_c = c.sin();
    let cos_c = c.cos();

    let dec = if rho == 0.0 {
        optical_axis_dec_rad
    } else {
        (cos_c * sin_dec0 + y * sin_c * cos_dec0 / rho).asin()
    };

    let ra = if rho == 0.0 {
        optical_axis_ra_rad
    } else {
        optical_axis_ra_rad + (x * sin_c).atan2(rho * cos_dec0 * cos_c - y * sin_dec0 * sin_c)
    };

    (ra, dec)
}

fn get_star_x_y(optical_axis_right_ascension_radians: f64, 
                optical_axis_declination_radians : f64, 
                star_right_ascension_radians: f64,
                star_declination_radians: f64 
                ) -> StarPosXY {

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

    return (star_x, star_y)
}



fn solve_transformation(
    matched_quads: &[(&StarQuad, &StarQuad)],
    pixel_resolution: f64
) -> Result<(f64, f64, f64, f64, f64, f64), Box<dyn std::error::Error>> {
    let mut image_points = Vec::new();
    let mut catalog_points = Vec::new();

    for (image_quad, catalog_quad) in matched_quads.iter() {
        for i in 0..image_quad.stars.len() {
            // The image stars are already in pixel coordinates relative to center
            let image_star = image_quad.stars[i];
            // The catalog stars are in standard coordinates, need to convert to pixels
            let catalog_star = catalog_quad.stars[i];

            image_points.push(image_star);
            catalog_points.push(catalog_star);
        }
    }

    if image_points.len() < 3 {
        return Err("Not enough matching points to solve transformation (need at least 3).".into());
    }

    // Create the M matrix for the least-squares problem
    let m = nalgebra::DMatrix::from_fn(image_points.len(), 3, |r, c| {
        match c {
            0 => image_points[r].0, // cameraX
            1 => image_points[r].1, // cameraY
            2 => 1.0,               // constant for C1/C2
            _ => 0.0,
        }
    });

    // Create the b_x and b_y vectors
    let b_x = nalgebra::DVector::from_iterator(catalog_points.len(), catalog_points.iter().map(|p| p.0));
    let b_y = nalgebra::DVector::from_iterator(catalog_points.len(), catalog_points.iter().map(|p| p.1));

    // Solve for the transformation coefficients using SVD-based least squares
    let svd = nalgebra::SVD::new(m, true, true);
    let x_coeffs = svd.solve(&b_x, 1e-18).map_err(|_| "Failed to solve for X coefficients")?;
    let y_coeffs = svd.solve(&b_y, 1e-18).map_err(|_| "Failed to solve for Y coefficients")?;

    let (a1, b1, c1) = (x_coeffs[0], x_coeffs[1], x_coeffs[2]);
    let (a2, b2, c2) = (y_coeffs[0], y_coeffs[1], y_coeffs[2]);

    Ok((a1, b1, c1, a2, b2, c2))
}


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let initial_coord = CoordinateEquatorial::new( 
        RaHoursMinutesSeconds::new(14, 01, 12.5), 
        Arcdegrees::new(54, 20, 56.0)
    );

    let file_path = Path::new("test.dng"); // Replace with your DNG file path

    let (image_x_size, image_y_size) = get_image_size(file_path)?;
    
    let (star_quads, centered_pixels, star_barycenters) = get_quad_from_file(file_path, false)?;

    println!("Found {} star quads in the image.", star_quads.len());

    let nb_of_star = star_barycenters.len();

    println!("Found {} star barycenters in the image.", nb_of_star);

    // annotate_stars_on_image(
    //     &centered_pixels,
    //     image_x_size,
    //     image_y_size,
    //     &star_barycenters,
    //     &star_quads,
    //     "annotated_image.png",
    //     15,
    // )?;

    let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT; 



    let image_fov_x = pixel_resolution * image_x_size as f64; // Convert to arcsec
    let image_fov_y = pixel_resolution * image_y_size as f64; // Convert to arcsec

    let image_fov = image_fov_x.max(image_fov_y);

    let star_in_fov = get_stars_from_catalogue(&initial_coord, image_fov * 1.0, nb_of_star * 100000)?;


    let vec_star = star_in_fov.iter()
        .map(|s| (get_star_x_y(
            initial_coord.ra.to_radians(),
            initial_coord.dec.to_radians(),
            s.ra.to_radians(),
            s.dec.to_radians(),
        )))
        .collect::<Vec<(f64, f64)>>();


    let star_quad_from_cat = returns_all_star_quads(&vec_star, 2);

    println!("Found {} star quads from catalogue.", star_quad_from_cat.len());

    let mut matches = 0;

    let mut matched_quads: Vec<(&StarQuad, &StarQuad)> = Vec::new();

    //TODO : implement while loop to compare each quad from the image with each quad from the catalogue

    for quad in &star_quads {
        for quad_cat in &star_quad_from_cat {
            if quad.compare(quad_cat, MATCHED_TOLERANCE) {
                matched_quads.push((quad, quad_cat));
                matches += 1;
            }
        }
    }

    let mut matched_quad_image: Vec<StarQuad> = Vec::new();
    let mut matched_quad_cat: Vec<StarQuad> = Vec::new();

    for elem in matched_quads.iter() {
        matched_quad_image.push(elem.0.clone());
        matched_quad_cat.push(elem.1.clone());
    }

    println!("Found {} matches between image and catalogue.", matches);


    annotate_stars_on_image(
        &[],
        (image_x_size as f64 * 1.0) as usize,
        (image_x_size as f64 * 1.0) as usize,
        &vec_star,
        &matched_quad_cat,
        "cat.png",
        15,
    )?;

    annotate_stars_on_image(
        &[],
        (image_x_size as f64 * 1.0) as usize,
        (image_x_size as f64 * 1.0) as usize,
        &star_barycenters,
        &matched_quad_image,
        "img.png",
        15,
    )?;


    // if matches > 2 {
    //     let (a1, b1, c1, a2, b2, c2) = solve_transformation(&matched_quads, pixel_resolution)?;
    //     println!("Solved transformation matrix:");
    //     println!("starDatabaseX = {:.4} * cameraX + {:.4} * cameraY + {:.4}", a1, b1, c1);
    //     println!("starDatabaseY = {:.4} * cameraX + {:.4} * cameraY + {:.4}", a2, b2, c2);

    //     // The center of the image in camera coordinates is (0, 0)
    //     let center_x_proj = c1; // a1*0 + b1*0 + c1
    //     let center_y_proj = c2; // a2*0 + b2*0 + c2

    //     let (center_ra_rad, center_dec_rad) = get_equatorial_from_xy(
    //         center_x_proj, 
    //         center_y_proj, 
    //         initial_coord.ra.to_radians(), 
    //         initial_coord.dec.to_radians()
    //     );

    //     let solved_coord = CoordinateEquatorial::from_radians(center_ra_rad, center_dec_rad);

    //     println!("\nSolved Image Center Coordinates:");
    //     println!("Right Ascension: {}", solved_coord.ra);
    //     println!("Declination: {}", solved_coord.dec);

    // } else {
    //     println!("Not enough matches found to determine the transformation matrix.");
    // }
    
    Ok(())
}
    
