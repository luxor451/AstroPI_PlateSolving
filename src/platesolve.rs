use crate::coordinate::*;
use crate::parse_catalog::*;
use crate::solver::*;
use crate::star_quads::*;
use std::path::Path;

const MATCHED_TOLERANCE: f64 = 0.0009;

/// Result of extracting star quads from an image
pub struct ImageAnalysisResult {
    pub star_quads: Vec<StarQuad>,
    pub centered_pixels: Vec<(i32, i32, u16)>,
    pub star_barycenters: Vec<(f64, f64)>,
    pub width: usize,
    pub height: usize,
}

/// Result of plate solving
pub struct PlateSolvingResult {
    pub matched_quads_count: usize,
    pub coeffs_x: Option<(f64, f64, f64)>,
    pub coeffs_y: Option<(f64, f64, f64)>,
    pub optical_axis_ra: f64,
    pub optical_axis_dec: f64,
}

/// Extracts star quads and pixel data from a DNG image file.
///
/// # Arguments
///
/// * `file_path` - Path to the DNG file to analyze.
/// * `verbose` - If true, prints progress information to stdout.
///
/// # Returns
///
/// An `ImageAnalysisResult` containing star quads, pixel data, barycenters, and image dimensions.
pub fn analyze_image(
    file_path: &Path,
    verbose: bool,
) -> Result<ImageAnalysisResult, Box<dyn std::error::Error>> {
    let (width, height) = get_image_size(file_path)?;
    let centered_pixels = get_pixel_matrix_from_dng(file_path)?;

    if verbose {
        println!(
            "Nb of pixels above 16000: {}",
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
        println!("Calculating star quads...")
    };
    let star_quads = returns_all_star_quads(&star_barycenters, 1);
    if verbose {
        println!("Found {} star quads.", star_quads.len())
    }

    Ok(ImageAnalysisResult {
        star_quads,
        centered_pixels,
        star_barycenters,
        width,
        height,
    })
}

/// Converts equatorial coordinates (RA/Dec) to image plane coordinates (x/y).
///
/// Uses gnomonic projection with the optical axis as the projection center.
///
/// # Arguments
///
/// * `optical_axis_right_ascension_radians` - RA of the optical axis in radians.
/// * `optical_axis_declination_radians` - Dec of the optical axis in radians.
/// * `star_right_ascension_radians` - RA of the star in radians.
/// * `star_declination_radians` - Dec of the star in radians.
///
/// # Returns
///
/// A tuple `(x, y)` representing the star's position on the image plane.
pub fn get_star_x_y(
    optical_axis_right_ascension_radians: f64,
    optical_axis_declination_radians: f64,
    star_right_ascension_radians: f64,
    star_declination_radians: f64,
) -> (f64, f64) {
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

    (star_x, star_y)
}

/// Converts image plane coordinates (x/y) to equatorial coordinates (RA/Dec).
///
/// Inverse of gnomonic projection using the optical axis as reference.
///
/// # Arguments
///
/// * `x` - X coordinate on the image plane.
/// * `y` - Y coordinate on the image plane.
/// * `optical_axis_ra_rad` - RA of the optical axis in radians.
/// * `optical_axis_dec_rad` - Dec of the optical axis in radians.
///
/// # Returns
///
/// A tuple `(ra, dec)` in radians representing the equatorial coordinates.
pub fn get_equatorial_from_xy(
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

/// Performs plate solving on a DNG image to determine astrometric solution.
///
/// Matches star quads from the image against a catalog to compute transformation
/// coefficients between image and sky coordinates.
///
/// # Arguments
///
/// * `file_path` - Path to the DNG file.
/// * `initial_coord` - Initial estimate of the image center coordinates.
/// * `verbose` - If true, prints detailed progress information.
///
/// # Returns
///
/// A `PlateSolvingResult` containing matched quad count, transformation coefficients,
/// and optical axis coordinates. Returns error if image processing or catalog query fails.
pub fn solve_plate(
    file_path: &Path,
    initial_coord: &CoordinateEquatorial,
    verbose: bool,
) -> Result<PlateSolvingResult, Box<dyn std::error::Error>> {
    // Analyze the image
    let image_analysis = analyze_image(file_path, verbose)?;

    if verbose {
        println!("Found {} star quads in the image.", image_analysis.star_quads.len());
        println!("Found {} star barycenters in the image.", image_analysis.star_barycenters.len());
    }

    // Calculate field of view
    let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT;
    let image_fov_x = pixel_resolution * image_analysis.width as f64;
    let image_fov_y = pixel_resolution * image_analysis.height as f64;
    let image_fov = image_fov_x.max(image_fov_y);

    // Get stars from catalog
    let nb_of_stars = image_analysis.star_barycenters.len();
    let star_in_fov = get_stars_from_catalogue(initial_coord, image_fov * 1.0, nb_of_stars * 100000)?;

    if verbose {
        println!("Retrieved {} stars from catalog.", star_in_fov.len());
    }

    // Convert catalog stars to x,y coordinates
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

    // Generate star quads from catalog
    let star_quad_from_cat = returns_all_star_quads(&vec_star, 2);

    if verbose {
        println!("Found {} star quads from catalog.", star_quad_from_cat.len());
    }

    // Match quads between image and catalog
    let mut matched_quads: Vec<(&StarQuad, &StarQuad)> = Vec::new();
    for quad in &image_analysis.star_quads {
        for quad_cat in &star_quad_from_cat {
            if quad.compare(quad_cat, MATCHED_TOLERANCE) {
                matched_quads.push((quad, quad_cat));
            }
        }
    }

    let matched_count = matched_quads.len();
    if verbose {
        println!("Found {} matches between image and catalog.", matched_count);
    }

    // Calculate transformation coefficients if enough matches
    let (coeffs_x, coeffs_y) = if matched_count >= 3 {
        let mut matched_quad_image: Vec<StarQuad> = Vec::new();
        let mut matched_quad_cat: Vec<StarQuad> = Vec::new();

        for elem in matched_quads.iter() {
            matched_quad_image.push(elem.0.clone());
            matched_quad_cat.push(elem.1.clone());
        }

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

        if verbose {
            println!("Coefficients for X: {:?}", coeffs_x);
            println!("Coefficients for Y: {:?}", coeffs_y);
        }

        (Some(coeffs_x), Some(coeffs_y))
    } else {
        if verbose {
            println!("Not enough matches found to determine the transformation matrix.");
        }
        (None, None)
    };

    Ok(PlateSolvingResult {
        matched_quads_count: matched_count,
        coeffs_x,
        coeffs_y,
        optical_axis_ra: initial_coord.ra.to_radians(),
        optical_axis_dec: initial_coord.dec.to_radians(),
    })
}
