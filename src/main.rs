use std::collections::HashMap;
use std::path::Path;

mod coordinate;
mod parse_catalog;
mod printing;
mod solver;
mod star_quads;
use coordinate::*;
use parse_catalog::*;
use printing::*;
use solver::*;
use star_quads::*;

const MATCHED_TOLERANCE: f64 = 0.0009; // Tolerance for matching star quads

fn get_quad_from_file(
    file_path: &Path,
    verbose: bool,
) -> Result<(Vec<StarQuad>, Vec<(i32, i32, u16)>, Vec<(f64, f64)>), Box<dyn std::error::Error>> {
    // Get image dimensions first
    let (width, height) = get_image_size(file_path)?;

    // Get the pixel matrix from the DNG file with centered coordinates
    let centered_pixels = get_pixel_matrix_from_dng(file_path)?;

    if verbose {
        println!(
            "Nb of pixel above 16000: {}",
            centered_pixels
                .iter()
                .filter(|&&(_, _, val)| val > 16000)
                .count()
        );
    }

    if verbose {
        println!("Calculating star barycenters...")
    };
    let star_barycenters = calculate_star_barycenters(&centered_pixels, width, height);
    if verbose {
        println!("Found {} star barycenters.", star_barycenters.len())
    };

    if verbose {
        println!("calculating star quads...")
    };
    let star_quads = returns_all_star_quads(&star_barycenters, 1);
    if verbose {
        println!("Found {} star quads.", star_quads.len())
    }

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
    optical_axis_dec_rad: f64,
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

fn get_star_x_y(
    optical_axis_right_ascension_radians: f64,
    optical_axis_declination_radians: f64,
    star_right_ascension_radians: f64,
    star_declination_radians: f64,
) -> StarPosXY {
    let sin_optical_axis_dec = optical_axis_declination_radians.sin();
    let cos_optical_axis_dec = optical_axis_declination_radians.cos();

    let sin_star_dec = star_declination_radians.sin();
    let cos_star_dec = star_declination_radians.cos();

    let ra_delta = star_right_ascension_radians - optical_axis_right_ascension_radians;
    let sin_delta_ra = ra_delta.sin();
    let cos_delta_ra = ra_delta.cos();

    let delta_value = (cos_optical_axis_dec * cos_star_dec * cos_delta_ra
        + sin_optical_axis_dec * sin_star_dec)
        / (3600.0 * (180.0 / std::f64::consts::PI));

    let star_x = -cos_star_dec * sin_delta_ra / delta_value;
    let star_y = -(sin_optical_axis_dec * cos_star_dec * cos_delta_ra
        - cos_optical_axis_dec * sin_star_dec)
        / delta_value;

    return (star_x, star_y);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let initial_coord = CoordinateEquatorial::new(
        RaHoursMinutesSeconds::new(14, 01, 12.5),
        Arcdegrees::new(54, 20, 56.0),
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

    let star_in_fov =
        get_stars_from_catalogue(&initial_coord, image_fov * 1.0, nb_of_star * 100000)?;

    let vec_star = star_in_fov
        .iter()
        .map(|s| {
            (get_star_x_y(
                initial_coord.ra.to_radians(),
                initial_coord.dec.to_radians(),
                s.ra.to_radians(),
                s.dec.to_radians(),
            ))
        })
        .collect::<Vec<(f64, f64)>>();

    let star_quad_from_cat = returns_all_star_quads(&vec_star, 2);

    println!(
        "Found {} star quads from catalogue.",
        star_quad_from_cat.len()
    );

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

    // annotate_stars_on_image(
    //     &[],
    //     (image_x_size as f64 * 1.0) as usize,
    //     (image_x_size as f64 * 1.0) as usize,
    //     &vec_star,
    //     &matched_quad_cat,
    //     "cat.png",
    //     15,
    // )?;

    // annotate_stars_on_image(
    //     &[],
    //     (image_x_size as f64 * 1.0) as usize,
    //     (image_x_size as f64 * 1.0) as usize,
    //     &star_barycenters,
    //     &matched_quad_image,
    //     "img.png",
    //     15,
    // )?;

    if matches >= 3 {
        let img_quad_barycenters_x = matched_quad_image
            .iter()
            .map(|quad| quad.barycenter.0)
            .collect::<Vec<f64>>();

        let img_quad_barycenters_y = matched_quad_image
            .iter()
            .map(|quad| quad.barycenter.1)
            .collect::<Vec<f64>>();

        let cat_quad_barycenters_x = matched_quad_cat
            .iter()
            .map(|quad| quad.barycenter.0)
            .collect::<Vec<f64>>();

        let cat_quad_barycenters_y = matched_quad_cat
            .iter()
            .map(|quad| quad.barycenter.1)
            .collect::<Vec<f64>>();

        let coeffs_x = solve_projection(
            &cat_quad_barycenters_x,
            &img_quad_barycenters_x,
            &img_quad_barycenters_y,
        );

        let coeffs_y = solve_projection(
            &cat_quad_barycenters_y,
            &img_quad_barycenters_x,
            &img_quad_barycenters_y,
        );

        println!("Coefficients for X: {:?}", coeffs_x);
        println!("Coefficients for Y: {:?}", coeffs_y);
    } else {
        println!("Not enough matches found to determine the transformation matrix.");
    }

    Ok(())
}
