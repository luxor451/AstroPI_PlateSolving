use std::path::Path;
use std::collections::HashMap; // Added VecDeque



mod star_quads;
mod parse_catalog;
mod coordinate;
mod printing;
use star_quads::*;
use parse_catalog::*; 
use coordinate::*;
use printing::*;


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
    let star_quads = returns_all_star_quads(&star_barycenters, 3);
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


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let initial_coord = CoordinateEquatorial::new( 
        RaHoursMinutesSeconds::new(14, 01, 12.5), 
        Arcdegrees::new(54, 20, 56.0)
    );

    let file_path = Path::new("test.dng"); // Replace with your DNG file path

    let (image_x_size, image_y_size) = get_image_size(file_path)?;
    
    let (star_quads, centered_pixels, star_barycenters) = get_quad_from_file(file_path, false)?;

    let nb_of_star = star_barycenters.len();

    annotate_stars_on_image(
        &centered_pixels,
        image_x_size,
        image_y_size,
        &star_barycenters,
        &star_quads,
        "annotated_image.png",
        15,
    )?;

    let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT; 



    let image_fov_x = pixel_resolution * image_x_size as f64; // Convert to arcsec
    let image_fov_y = pixel_resolution * image_y_size as f64; // Convert to arcsec

    let image_fov = image_fov_x.max(image_fov_y);

    let star_in_fov = get_stars_from_catalogue(&initial_coord, image_fov * 1.5, nb_of_star * 3)?;


    let vec_star = star_in_fov.iter()
        .map(|s| (get_star_x_y(
            initial_coord.ra.to_radians(),
            initial_coord.dec.to_radians(),
            s.ra.to_radians(),
            s.dec.to_radians(),
        )))
        .collect::<Vec<(f64, f64)>>();


    let star_quad_from_cat = returns_all_star_quads(&vec_star, 6);

    annotate_stars_on_image(
        &[],
        (image_x_size as f64 * 1.5) as usize,
        (image_x_size as f64 * 1.5) as usize,
        &vec_star,
        &star_quad_from_cat,
        "annotated_image_with_catalog.png",
        15,
    )?;
    
    Ok(())
}