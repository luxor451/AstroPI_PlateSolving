use crate::coordinate::*;
use crate::parse_catalog::*;
use crate::solver::*;
use crate::star_quads::*;
use std::path::Path;

const MATCHED_TOLERANCE: f64 = 0.00085;

/// Maximum number of spiral search iterations before giving up
const MAX_SPIRAL_ITERATIONS: usize = 100;

/// Minimum number of matched quads required for a valid solution
const MIN_MATCHED_QUADS: usize = 3;

/// Maximum iterations for star matching refinement
const MAX_STAR_MATCH_ITERATIONS: usize = 10;

/// Initial tolerance for star matching in arcseconds
const INITIAL_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 60.0;

/// Minimum tolerance for star matching in arcseconds
const MIN_STAR_MATCH_TOLERANCE_ARCSEC: f64 = 5.0;

/// Result of extracting star quads from an image
pub struct ImageAnalysisResult {
    pub star_quads: Vec<StarQuad>,
    pub centered_pixels: Vec<(i32, i32, u16)>,
    pub star_barycenters: Vec<(f64, f64)>,
    pub width: usize,
    pub height: usize,
}

/// Transformation coefficients for coordinate system projection
/// 
/// These coefficients define the affine transformation from image coordinates
/// to catalog projected coordinates:
/// - X_catalog = A1 * x_img + B1 * y_img + C1
/// - Y_catalog = A2 * x_img + B2 * y_img + C2
#[derive(Debug, Clone, Copy)]
pub struct TransformCoefficients {
    pub a1: f64, // X-axis: coefficient for image x
    pub b1: f64, // X-axis: coefficient for image y
    pub c1: f64, // X-axis: translation offset
    pub a2: f64, // Y-axis: coefficient for image x
    pub b2: f64, // Y-axis: coefficient for image y
    pub c2: f64, // Y-axis: translation offset
}

impl TransformCoefficients {
    /// Creates new transformation coefficients from two sets of (A, B, C) values.
    pub fn new(coeffs_x: (f64, f64, f64), coeffs_y: (f64, f64, f64)) -> Self {
        Self {
            a1: coeffs_x.0,
            b1: coeffs_x.1,
            c1: coeffs_x.2,
            a2: coeffs_y.0,
            b2: coeffs_y.1,
            c2: coeffs_y.2,
        }
    }

    /// Transforms image coordinates to catalog projected plane coordinates.
    /// 
    /// # Arguments
    /// 
    /// * `img_x` - X coordinate in image space (centered pixel coordinates).
    /// * `img_y` - Y coordinate in image space (centered pixel coordinates).
    /// 
    /// # Returns
    /// 
    /// A tuple `(cat_x, cat_y)` in the catalog projected coordinate system.
    pub fn image_to_catalog(&self, img_x: f64, img_y: f64) -> (f64, f64) {
        let cat_x = self.a1 * img_x + self.b1 * img_y + self.c1;
        let cat_y = self.a2 * img_x + self.b2 * img_y + self.c2;
        (cat_x, cat_y)
    }

    /// Transforms catalog projected coordinates back to image coordinates.
    /// 
    /// Uses the inverse of the affine transformation matrix.
    /// 
    /// # Arguments
    /// 
    /// * `cat_x` - X coordinate in catalog projected space.
    /// * `cat_y` - Y coordinate in catalog projected space.
    /// 
    /// # Returns
    /// 
    /// A tuple `(img_x, img_y)` in image coordinate space.
    pub fn catalog_to_image(&self, cat_x: f64, cat_y: f64) -> (f64, f64) {
        // Inverse of the 2x2 transformation matrix [a1, b1; a2, b2]
        let det = self.a1 * self.b2 - self.b1 * self.a2;
        if det.abs() < 1e-14 {
            // Singular matrix, return origin
            return (0.0, 0.0);
        }
        
        // Subtract translation first
        let dx = cat_x - self.c1;
        let dy = cat_y - self.c2;
        
        // Apply inverse matrix
        let img_x = (self.b2 * dx - self.b1 * dy) / det;
        let img_y = (-self.a2 * dx + self.a1 * dy) / det;
        
        (img_x, img_y)
    }

    /// Returns the scale factor (pixels to arcseconds).
    pub fn scale(&self) -> f64 {
        ((self.a1.powi(2) + self.a2.powi(2)).sqrt() + (self.b1.powi(2) + self.b2.powi(2)).sqrt()) / 2.0
    }

    /// Returns the rotation angle in radians.
    pub fn rotation(&self) -> f64 {
        self.a2.atan2(self.a1)
    }
}

/// Result of plate solving
pub struct PlateSolvingResult {
    pub matched_quads_count: usize,
    pub coeffs_x: Option<(f64, f64, f64)>,
    pub coeffs_y: Option<(f64, f64, f64)>,
    pub transform: Option<TransformCoefficients>,
    pub optical_axis_ra: f64,
    pub optical_axis_dec: f64,
    /// Number of spiral search iterations performed (0 = found at initial position)
    pub spiral_iterations: usize,
    /// The search position where the solution was found (RA in degrees)
    pub solution_ra_deg: f64,
    /// The search position where the solution was found (Dec in degrees)
    pub solution_dec_deg: f64,
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

    // cos(c) = sin(δ₀)sin(δ) + cos(δ₀)cos(δ)cos(Δα) - the angular distance cosine
    let cos_c = sin_optical_axis_dec * sin_star_dec 
        + cos_optical_axis_dec * cos_star_dec * cos_delta_ra;
    
    // Standard gnomonic projection formulas:
    // x = cos(δ) * sin(Δα) / cos(c)
    // y = (cos(δ₀) * sin(δ) - sin(δ₀) * cos(δ) * cos(Δα)) / cos(c)
    // 
    // Scale factor to convert from radians to arcseconds: 3600 * (180/π) ≈ 206264.8
    let rad_to_arcsec = 3600.0 * (180.0 / std::f64::consts::PI);
    
    let star_x = cos_star_dec * sin_delta_ra / cos_c * rad_to_arcsec;
    let star_y = (cos_optical_axis_dec * sin_star_dec 
        - sin_optical_axis_dec * cos_star_dec * cos_delta_ra) / cos_c * rad_to_arcsec;

    (star_x, star_y)
}

/// Converts image plane coordinates (x/y) to equatorial coordinates (RA/Dec).
///
/// Inverse of gnomonic projection using the optical axis as reference.
/// Input coordinates should be tangent plane coordinates (dimensionless, in radians).
///
/// # Arguments
///
/// * `x` - X coordinate on the tangent plane (radians).
/// * `y` - Y coordinate on the tangent plane (radians).
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
    
    // For small angles, atan(rho) ≈ rho, but we need the proper formula
    // The gnomonic projection maps: tan(c) = rho, so c = atan(rho)
    let c = rho.atan();
    let sin_c = c.sin();
    let cos_c = c.cos();

    let dec = if rho < 1e-10 {
        optical_axis_dec_rad
    } else {
        (cos_c * sin_dec0 + y * sin_c * cos_dec0 / rho).asin()
    };

    let ra = if rho < 1e-10 {
        optical_axis_ra_rad
    } else {
        optical_axis_ra_rad + (x * sin_c).atan2(rho * cos_dec0 * cos_c - y * sin_dec0 * sin_c)
    };

    (ra, dec)
}


/// Converts image pixel coordinates to equatorial coordinates (RA/Dec).
///
/// This is the complete transformation pipeline:
/// 1. Transform pixel coordinates to catalog projected plane using the affine transform
/// 2. Convert projected plane coordinates to equatorial coordinates using inverse gnomonic projection
///
/// # Arguments
///
/// * `pixel_x` - X coordinate in centered pixel coordinates (0 = image center).
/// * `pixel_y` - Y coordinate in centered pixel coordinates (0 = image center).
/// * `transform` - The transformation coefficients from plate solving.
/// * `optical_axis_ra_rad` - RA of the optical axis in radians.
/// * `optical_axis_dec_rad` - Dec of the optical axis in radians.
///
/// # Returns
///
/// A tuple `(ra, dec)` in radians representing the equatorial coordinates.
pub fn pixel_to_equatorial(
    pixel_x: f64,
    pixel_y: f64,
    transform: &TransformCoefficients,
    optical_axis_ra_rad: f64,
    optical_axis_dec_rad: f64,
) -> (f64, f64) {
    // Step 1: Transform image coordinates to catalog projected plane
    let (cat_x, cat_y) = transform.image_to_catalog(pixel_x, pixel_y);
    
    // Step 2: Convert projected coordinates to equatorial
    // Note: The projected coordinates are in arcseconds, need to convert to radians
    let arcsec_to_rad = std::f64::consts::PI / (180.0 * 3600.0);
    let x_rad = cat_x * arcsec_to_rad;
    let y_rad = cat_y * arcsec_to_rad;
    
    get_equatorial_from_xy(x_rad, y_rad, optical_axis_ra_rad, optical_axis_dec_rad)
}

/// Calculates the equatorial coordinates of the image center.
///
/// The image center is at pixel coordinates (0, 0) in the centered coordinate system.
/// This function uses the transformation coefficients to find where the optical center
/// of the camera points in the sky.
///
/// # Arguments
///
/// * `transform` - The transformation coefficients from plate solving.
/// * `initial_optical_axis_ra_rad` - Initial RA estimate of the optical axis in radians.
/// * `initial_optical_axis_dec_rad` - Initial Dec estimate of the optical axis in radians.
///
/// # Returns
///
/// A tuple `(ra, dec)` in radians representing the true optical axis position.
pub fn calculate_image_center_equatorial(
    transform: &TransformCoefficients,
    initial_optical_axis_ra_rad: f64,
    initial_optical_axis_dec_rad: f64,
) -> (f64, f64) {
    // The image center is at pixel (0,0)
    // Transform gives us: cat_x = C1, cat_y = C2 when img_x = img_y = 0
    // These are the offsets in the catalog's projected coordinate system (arcseconds)
    
    let (cat_x, cat_y) = transform.image_to_catalog(0.0, 0.0);
    
    // Convert from arcseconds to radians for the inverse projection
    let arcsec_to_rad = std::f64::consts::PI / (180.0 * 3600.0);
    let x_rad = cat_x * arcsec_to_rad;
    let y_rad = cat_y * arcsec_to_rad;
    
    get_equatorial_from_xy(x_rad, y_rad, initial_optical_axis_ra_rad, initial_optical_axis_dec_rad)
}


/// Generates spiral search offsets around the center point.
/// 
/// The spiral moves outward in FOV-sized steps:
/// - First iteration: center (0, 0)
/// - Then expands in a spiral pattern: right, up, left, left, down, down, right, right, right, ...
/// 
/// # Arguments
/// 
/// * `max_iterations` - Maximum number of positions to generate.
/// 
/// # Returns
/// 
/// A vector of (delta_x, delta_y) offsets in FOV units.
fn generate_spiral_offsets(max_iterations: usize) -> Vec<(i32, i32)> {
    let mut offsets = Vec::with_capacity(max_iterations);
    
    // Start at center
    offsets.push((0, 0));
    
    if max_iterations == 1 {
        return offsets;
    }
    
    let mut x: i32 = 0;
    let mut y: i32 = 0;
    let mut dx: i32 = 1;
    let mut dy: i32 = 0;
    let mut steps_in_direction = 1;
    let mut steps_taken = 0;
    let mut direction_changes = 0;
    
    while offsets.len() < max_iterations {
        // Move one step
        x += dx;
        y += dy;
        offsets.push((x, y));
        steps_taken += 1;
        
        // Check if we need to turn
        if steps_taken == steps_in_direction {
            steps_taken = 0;
            direction_changes += 1;
            
            // Turn left (counterclockwise)
            let temp = dx;
            dx = -dy;
            dy = temp;
            
            // Increase step count every 2 turns
            if direction_changes % 2 == 0 {
                steps_in_direction += 1;
            }
        }
    }
    
    offsets
}

/// Result of attempting to match quads at a specific sky position.
struct MatchAttemptResult<'a> {
    matched_quads: Vec<(&'a StarQuad, StarQuad)>,
    search_coord_ra_deg: f64,
    search_coord_dec_deg: f64,
    catalog_stars_xy: Vec<(f64, f64)>,
}

/// Attempts to match image quads against catalog quads at a specific sky position.
/// 
/// # Arguments
/// 
/// * `image_quads` - Star quads extracted from the image.
/// * `search_ra_deg` - Right Ascension of the search center in degrees.
/// * `search_dec_deg` - Declination of the search center in degrees.
/// * `fov_arcsec` - Field of view in arcseconds.
/// * `max_stars` - Maximum number of stars to retrieve from catalog.
/// * `verbose` - If true, prints debug information.
/// 
/// # Returns
/// 
/// A `MatchAttemptResult` containing matched quad pairs and the search coordinates.
fn try_match_at_position<'a>(
    image_quads: &'a [StarQuad],
    search_ra_deg: f64,
    search_dec_deg: f64,
    fov_arcsec: f64,
    max_stars: usize,
    verbose: bool,
) -> Result<MatchAttemptResult<'a>, Box<dyn std::error::Error>> {
    // Create coordinate for this search position
    let search_coord = CoordinateEquatorial::from_degrees(search_ra_deg, search_dec_deg);
    
    // Get stars from catalog at this position
    let stars_in_fov = get_stars_from_catalogue(&search_coord, fov_arcsec * 1.0, max_stars)?;
    
    if verbose {
        println!("  Retrieved {} stars from catalog at RA={:.4}°, Dec={:.4}°", 
                 stars_in_fov.len(), search_ra_deg, search_dec_deg);
    }
    
    // Convert catalog stars to x,y coordinates
    let vec_star: Vec<(f64, f64)> = stars_in_fov
        .iter()
        .map(|s| {
            get_star_x_y(
                search_ra_deg.to_radians(),
                search_dec_deg.to_radians(),
                s.ra.to_radians(),
                s.dec.to_radians(),
            )
        })
        .collect();
    
    // Generate star quads from catalog
    let star_quad_from_cat = returns_all_star_quads(&vec_star, 2);
    
    if verbose {
        println!("  Generated {} catalog quads.", star_quad_from_cat.len());
    }
    
    // Match quads between image and catalog
    // Each image quad should only match one catalog quad (best match)
    // and each catalog quad should only be used once
    // Also avoid duplicate image quads (same barycenter)
    let mut matched_quads: Vec<(&StarQuad, StarQuad)> = Vec::new();
    let mut used_catalog_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut used_image_barycenters: std::collections::HashSet<(i64, i64)> = std::collections::HashSet::new();
    
    for quad in image_quads {
        // Create a key from the barycenter (rounded to avoid floating point issues)
        let barycenter_key = (
            (quad.barycenter.0 * 10.0).round() as i64,
            (quad.barycenter.1 * 10.0).round() as i64,
        );
        
        // Skip if we've already matched a quad with this barycenter
        if used_image_barycenters.contains(&barycenter_key) {
            continue;
        }
        
        let mut best_match: Option<(usize, &StarQuad, f64)> = None;
        
        for (cat_idx, quad_cat) in star_quad_from_cat.iter().enumerate() {
            // Skip if this catalog quad is already matched
            if used_catalog_indices.contains(&cat_idx) {
                continue;
            }
            
            if quad.compare(quad_cat, MATCHED_TOLERANCE) {
                // Calculate match quality (sum of squared differences in normalized distances)
                let match_quality: f64 = quad.normalized_distances
                    .iter()
                    .zip(quad_cat.normalized_distances.iter())
                    .map(|(d1, d2)| (d1 - d2).powi(2))
                    .sum();
                
                if best_match.is_none() || match_quality < best_match.unwrap().2 {
                    best_match = Some((cat_idx, quad_cat, match_quality));
                }
            }
        }
        
        // If we found a match for this image quad, record it
        if let Some((cat_idx, quad_cat, _)) = best_match {
            matched_quads.push((quad, quad_cat.clone()));
            used_catalog_indices.insert(cat_idx);
            used_image_barycenters.insert(barycenter_key);
        }
    }
    
    Ok(MatchAttemptResult {
        matched_quads,
        search_coord_ra_deg: search_ra_deg,
        search_coord_dec_deg: search_dec_deg,
        catalog_stars_xy: vec_star,
    })
}

/// Performs plate solving on a DNG image to determine astrometric solution.
///
/// Matches star quads from the image against a catalog to compute transformation
/// coefficients between image and sky coordinates. If the initial position doesn't
/// yield enough matches, performs a spiral search around the initial coordinates,
/// with each search position offset by the camera's field of view.
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
    solve_plate_with_options(file_path, initial_coord, verbose, MAX_SPIRAL_ITERATIONS)
}

/// Performs plate solving on a DNG image with configurable spiral search options.
///
/// Matches star quads from the image against a catalog to compute transformation
/// coefficients between image and sky coordinates. If the initial position doesn't
/// yield enough matches, performs a spiral search around the initial coordinates,
/// with each search position offset by the camera's field of view.
///
/// # Arguments
///
/// * `file_path` - Path to the DNG file.
/// * `initial_coord` - Initial estimate of the image center coordinates.
/// * `verbose` - If true, prints detailed progress information.
/// * `max_spiral_iterations` - Maximum number of spiral search positions to try (1 = no spiral search).
///
/// # Returns
///
/// A `PlateSolvingResult` containing matched quad count, transformation coefficients,
/// and optical axis coordinates. Returns error if image processing or catalog query fails.
pub fn solve_plate_with_options(
    file_path: &Path,
    initial_coord: &CoordinateEquatorial,
    verbose: bool,
    max_spiral_iterations: usize,
) -> Result<PlateSolvingResult, Box<dyn std::error::Error>> {
    // Analyze the image
    let image_analysis = analyze_image(file_path, verbose)?;

    if verbose {
        println!("Found {} star quads in the image.", image_analysis.star_quads.len());
        println!("Found {} star barycenters in the image.", image_analysis.star_barycenters.len());
    }

    // Calculate field of view in arcseconds
    let pixel_resolution: f64 = 206.265 * PIXEL_SIZE_MICRON / TELESCOPE_FOCAL_LENGHT;
    let image_fov_x_arcsec = pixel_resolution * image_analysis.width as f64;
    let image_fov_y_arcsec = pixel_resolution * image_analysis.height as f64;
    let image_fov = image_fov_x_arcsec.max(image_fov_y_arcsec);
    
    // Convert FOV to degrees for spiral search offsets
    let fov_x_deg = image_fov_x_arcsec / 3600.0;
    let fov_y_deg = image_fov_y_arcsec / 3600.0;
    
    if verbose {
        println!("Image FOV: {:.2}\" x {:.2}\" ({:.4}° x {:.4}°)", 
                 image_fov_x_arcsec, image_fov_y_arcsec, fov_x_deg, fov_y_deg);
    }

    let nb_of_stars = image_analysis.star_barycenters.len();
    let max_stars = nb_of_stars * 100000;
    
    // Get initial coordinates in degrees
    let (initial_ra_deg, initial_dec_deg) = initial_coord.to_degrees();
    
    // Generate spiral search pattern
    let spiral_offsets = generate_spiral_offsets(max_spiral_iterations.max(1));
    
    // Variables to store the best match result
    let mut best_matched_quads: Vec<(StarQuad, StarQuad)> = Vec::new();
    let mut best_search_ra_deg = initial_ra_deg;
    let mut best_search_dec_deg = initial_dec_deg;
    let mut best_catalog_stars_xy: Vec<(f64, f64)> = Vec::new();
    let mut solution_iteration: usize = 0;
    
    // Perform spiral search
    for (iteration, (offset_x, offset_y)) in spiral_offsets.iter().enumerate() {
        // Calculate the search position
        // Account for RA scaling with declination (cos(dec) factor)
        let dec_rad = initial_dec_deg.to_radians();
        let ra_scale = if dec_rad.cos().abs() > 0.001 { 
            1.0 / dec_rad.cos() 
        } else { 
            1.0 
        };
        
        let search_ra_deg = initial_ra_deg + (*offset_x as f64) * fov_x_deg * ra_scale;
        let search_dec_deg = initial_dec_deg + (*offset_y as f64) * fov_y_deg;
        
        // Clamp declination to valid range
        let search_dec_deg = search_dec_deg.clamp(-90.0, 90.0);
        // Normalize RA to 0-360 range
        let search_ra_deg = ((search_ra_deg % 360.0) + 360.0) % 360.0;
        
        if verbose {
            if iteration == 0 {
                println!("Attempting match at initial position...");
            } else {
                println!("Spiral search iteration {}: offset ({}, {}), RA={:.4}°, Dec={:.4}°", 
                         iteration, offset_x, offset_y, search_ra_deg, search_dec_deg);
            }
        }
        
        // Try to match at this position
        let match_result = try_match_at_position(
            &image_analysis.star_quads,
            search_ra_deg,
            search_dec_deg,
            image_fov,
            max_stars,
            verbose,
        )?;
        
        let matched_count = match_result.matched_quads.len();
        
        if verbose {
            println!("  Found {} matches at this position.", matched_count);
        }
        
        // Check if we found enough matches
        if matched_count >= MIN_MATCHED_QUADS {
            if verbose {
                println!("Found {} matches (>= {}), proceeding with solution.", 
                         matched_count, MIN_MATCHED_QUADS);
            }
            
            // Convert to owned quads for storage
            best_matched_quads = match_result.matched_quads
                .into_iter()
                .map(|(img, cat)| (img.clone(), cat))
                .collect();
            best_search_ra_deg = match_result.search_coord_ra_deg;
            best_search_dec_deg = match_result.search_coord_dec_deg;
            best_catalog_stars_xy = match_result.catalog_stars_xy;
            solution_iteration = iteration;
            break;
        }
        
        // Keep track of best result so far (even if below threshold)
        if matched_count > best_matched_quads.len() {
            best_matched_quads = match_result.matched_quads
                .into_iter()
                .map(|(img, cat)| (img.clone(), cat))
                .collect();
            best_search_ra_deg = match_result.search_coord_ra_deg;
            best_search_dec_deg = match_result.search_coord_dec_deg;
            best_catalog_stars_xy = match_result.catalog_stars_xy;
            solution_iteration = iteration;
        }
    }
    
    let matched_count = best_matched_quads.len();
    
    if verbose {
        if matched_count >= MIN_MATCHED_QUADS {
            println!("Solution found with {} matched quads at RA={:.4}°, Dec={:.4}°", 
                     matched_count, best_search_ra_deg, best_search_dec_deg);
        } else {
            println!("Spiral search completed. Best result: {} matches (need {} for solution).", 
                     matched_count, MIN_MATCHED_QUADS);
        }
    }

    // Calculate transformation coefficients if enough matches
    // We build an overdetermined system of linear equations:
    // For X-axis: X_catalog = A1 * x_img + B1 * y_img + C1
    // For Y-axis: Y_catalog = A2 * x_img + B2 * y_img + C2
    // With more than 3 matched quads, we use least squares (SVD) to find best fit
    // 
    // We use an iterative approach that matches individual stars:
    // 1. Get initial rough transformation from matched quad barycenters
    // 2. Use transformation to predict catalog positions for ALL image stars
    // 3. Match each image star to nearest catalog star within tolerance
    // 4. Re-fit transformation using individual star matches
    // 5. Iterate until convergence
    let (coeffs_x, coeffs_y, transform, optical_ra_deg, optical_dec_deg, final_matched_count) = if matched_count >= MIN_MATCHED_QUADS {
        
        // Step 1: Get initial transformation from quad barycenters
        let img_x: Vec<f64> = best_matched_quads.iter().map(|(q, _)| q.barycenter.0).collect();
        let img_y: Vec<f64> = best_matched_quads.iter().map(|(q, _)| q.barycenter.1).collect();
        let cat_x: Vec<f64> = best_matched_quads.iter().map(|(_, q)| q.barycenter.0).collect();
        let cat_y: Vec<f64> = best_matched_quads.iter().map(|(_, q)| q.barycenter.1).collect();
        
        let mut coeffs_x = solve_projection(&cat_x, &img_x, &img_y);
        let mut coeffs_y = solve_projection(&cat_y, &img_x, &img_y);
        let mut transform = TransformCoefficients::new(coeffs_x, coeffs_y);
        
        if verbose {
            println!("Initial transformation from {} quad barycenters:", best_matched_quads.len());
            println!("  Scale: {:.4} arcsec/pixel, Rotation: {:.4}°", 
                     transform.scale(), transform.rotation().to_degrees());
        }
        
        // Step 2-5: Iteratively refine by matching individual stars
        let mut matched_stars: Vec<((f64, f64), (f64, f64))> = Vec::new(); // (image_pos, catalog_pos)
        
        for iteration in 0..MAX_STAR_MATCH_ITERATIONS {
            // Use current transformation to predict catalog positions for all image stars
            let mut new_matched_stars: Vec<((f64, f64), (f64, f64))> = Vec::new();
            
            // Tolerance decreases more gradually: divide by 1.5 each iteration instead of 2
            let tolerance = (INITIAL_STAR_MATCH_TOLERANCE_ARCSEC / (1.5_f64).powi(iteration as i32))
                .max(MIN_STAR_MATCH_TOLERANCE_ARCSEC);
            
            for &(img_x, img_y) in &image_analysis.star_barycenters {
                // Predict where this image star should be in catalog coordinates
                let (pred_cat_x, pred_cat_y) = transform.image_to_catalog(img_x, img_y);
                
                // Find the nearest catalog star
                let mut best_dist = f64::MAX;
                let mut best_cat: Option<(f64, f64)> = None;
                
                for &(cat_x, cat_y) in &best_catalog_stars_xy {
                    let dx = pred_cat_x - cat_x;
                    let dy = pred_cat_y - cat_y;
                    let dist = (dx * dx + dy * dy).sqrt();
                    
                    if dist < best_dist {
                        best_dist = dist;
                        best_cat = Some((cat_x, cat_y));
                    }
                }
                
                // Accept match if within tolerance
                if best_dist < tolerance {
                    if let Some(cat_pos) = best_cat {
                        new_matched_stars.push(((img_x, img_y), cat_pos));
                    }
                }
            }
            
            if verbose {
                println!("Star matching iteration {}: {} matched stars (tolerance: {:.1}\")", 
                         iteration + 1, new_matched_stars.len(), tolerance);
            }
            
            // Need at least 3 matched stars to fit transformation
            if new_matched_stars.len() < 3 {
                if verbose {
                    println!("  Too few star matches, using previous result");
                }
                break;
            }
            
            // Re-fit transformation with matched stars
            let img_x: Vec<f64> = new_matched_stars.iter().map(|((x, _), _)| *x).collect();
            let img_y: Vec<f64> = new_matched_stars.iter().map(|((_, y), _)| *y).collect();
            let cat_x: Vec<f64> = new_matched_stars.iter().map(|(_, (x, _))| *x).collect();
            let cat_y: Vec<f64> = new_matched_stars.iter().map(|(_, (_, y))| *y).collect();
            
            let new_coeffs_x = solve_projection(&cat_x, &img_x, &img_y);
            let new_coeffs_y = solve_projection(&cat_y, &img_x, &img_y);
            let new_transform = TransformCoefficients::new(new_coeffs_x, new_coeffs_y);
            
            // Calculate RMS residual
            let rms: f64 = new_matched_stars.iter()
                .map(|((ix, iy), (cx, cy))| {
                    let (px, py) = new_transform.image_to_catalog(*ix, *iy);
                    (px - cx).powi(2) + (py - cy).powi(2)
                })
                .sum::<f64>() / new_matched_stars.len() as f64;
            let rms = rms.sqrt();
            
            if verbose {
                println!("  RMS residual: {:.2}\" Scale: {:.4}\"/pix Rotation: {:.2}°",
                         rms, new_transform.scale(), new_transform.rotation().to_degrees());
            }
            
            // Check for convergence (transformation doesn't change much)
            let scale_change = (new_transform.scale() - transform.scale()).abs();
            let rotation_change = (new_transform.rotation() - transform.rotation()).abs();
            
            // Update for next iteration
            coeffs_x = new_coeffs_x;
            coeffs_y = new_coeffs_y;
            transform = new_transform;
            matched_stars = new_matched_stars;
            
            if scale_change < 0.0001 && rotation_change < 0.0001 && iteration > 0 {
                if verbose {
                    println!("  Converged!");
                }
                break;
            }
        }
        
        let final_match_count = matched_stars.len();

        if verbose {
            println!();
            println!("Final transformation coefficients:");
            println!("  X-axis: A1={:.6}, B1={:.6}, C1={:.6}", 
                     coeffs_x.0, coeffs_x.1, coeffs_x.2);
            println!("  Y-axis: A2={:.6}, B2={:.6}, C2={:.6}", 
                     coeffs_y.0, coeffs_y.1, coeffs_y.2);
            println!("  Scale: {:.4} arcsec/pixel", transform.scale());
            println!("  Rotation: {:.4} degrees", transform.rotation().to_degrees());
            println!("  Matched stars: {}", final_match_count);
            
            // Calculate and display the true image center in equatorial coordinates
            let (center_ra, center_dec) = calculate_image_center_equatorial(
                &transform,
                best_search_ra_deg.to_radians(),
                best_search_dec_deg.to_radians(),
            );
            println!("Image center RA:  {:.6}° ({:.6} rad)", center_ra.to_degrees(), center_ra);
            println!("Image center Dec: {:.6}° ({:.6} rad)", center_dec.to_degrees(), center_dec);
        }

        (Some(coeffs_x), Some(coeffs_y), Some(transform), best_search_ra_deg, best_search_dec_deg, final_match_count)
    } else {
        if verbose {
            println!("Not enough matches found to determine the transformation matrix.");
            println!("Need at least {} matched quads, found {}.", MIN_MATCHED_QUADS, matched_count);
        }

        (None, None, None, initial_ra_deg, initial_dec_deg, 0)
    };

    // Calculate the refined optical axis position if we have a valid transform
    let (mut refined_ra, mut refined_dec) = if let Some(ref t) = transform {
        calculate_image_center_equatorial(
            t,
            optical_ra_deg.to_radians(),
            optical_dec_deg.to_radians(),
        )
    } else {
        (initial_coord.ra.to_radians(), initial_coord.dec.to_radians())
    };

    // Re-centering refinement: if the solved position is significantly different from the
    // search position, re-project catalog stars from the solved position and refine
    let mut final_coeffs_x = coeffs_x;
    let mut final_coeffs_y = coeffs_y;
    let mut final_transform = transform;
    let mut final_match_count = final_matched_count;

    if coeffs_x.is_some() {
        const RE_CENTER_THRESHOLD_DEG: f64 = 0.5; // Re-center if more than 0.5° off
        const MAX_RE_CENTER_ITERATIONS: usize = 3;
        
        let offset_deg = ((refined_ra.to_degrees() - optical_ra_deg).powi(2) 
                        + (refined_dec.to_degrees() - optical_dec_deg).powi(2)).sqrt();
        
        if offset_deg > RE_CENTER_THRESHOLD_DEG {
            if verbose {
                println!();
                println!("Re-centering: solved position is {:.2}° from search position", offset_deg);
            }
            
            for re_iter in 0..MAX_RE_CENTER_ITERATIONS {
                // Use the current solved position as the new projection center
                let new_center_ra_deg = refined_ra.to_degrees();
                let new_center_dec_deg = refined_dec.to_degrees();
                
                if verbose {
                    println!("Re-center iteration {}: projecting from RA={:.4}°, Dec={:.4}°",
                             re_iter + 1, new_center_ra_deg, new_center_dec_deg);
                }
                
                // Re-fetch catalog stars centered on the new position
                let new_match_result = try_match_at_position(
                    &image_analysis.star_quads,
                    new_center_ra_deg,
                    new_center_dec_deg,
                    image_fov,
                    max_stars,
                    false, // not verbose for re-centering
                )?;
                
                if new_match_result.matched_quads.len() >= MIN_MATCHED_QUADS {
                    // Get initial transformation from quad barycenters
                    let img_x: Vec<f64> = new_match_result.matched_quads.iter().map(|(q, _)| q.barycenter.0).collect();
                    let img_y: Vec<f64> = new_match_result.matched_quads.iter().map(|(q, _)| q.barycenter.1).collect();
                    let cat_x: Vec<f64> = new_match_result.matched_quads.iter().map(|(_, q)| q.barycenter.0).collect();
                    let cat_y: Vec<f64> = new_match_result.matched_quads.iter().map(|(_, q)| q.barycenter.1).collect();
                    
                    let mut iter_coeffs_x = solve_projection(&cat_x, &img_x, &img_y);
                    let mut iter_coeffs_y = solve_projection(&cat_y, &img_x, &img_y);
                    let mut iter_transform = TransformCoefficients::new(iter_coeffs_x, iter_coeffs_y);
                    
                    // Iteratively refine by matching individual stars
                    for iteration in 0..MAX_STAR_MATCH_ITERATIONS {
                        let mut new_matched_stars: Vec<((f64, f64), (f64, f64))> = Vec::new();
                        let tolerance = (INITIAL_STAR_MATCH_TOLERANCE_ARCSEC / (1.5_f64).powi(iteration as i32))
                            .max(MIN_STAR_MATCH_TOLERANCE_ARCSEC);
                        
                        for &(img_x, img_y) in &image_analysis.star_barycenters {
                            let (pred_cat_x, pred_cat_y) = iter_transform.image_to_catalog(img_x, img_y);
                            let mut best_dist = f64::MAX;
                            let mut best_cat: Option<(f64, f64)> = None;
                            
                            for &(cat_x, cat_y) in &new_match_result.catalog_stars_xy {
                                let dx = pred_cat_x - cat_x;
                                let dy = pred_cat_y - cat_y;
                                let dist = (dx * dx + dy * dy).sqrt();
                                if dist < best_dist {
                                    best_dist = dist;
                                    best_cat = Some((cat_x, cat_y));
                                }
                            }
                            
                            if best_dist < tolerance {
                                if let Some(cat_pos) = best_cat {
                                    new_matched_stars.push(((img_x, img_y), cat_pos));
                                }
                            }
                        }
                        
                        if new_matched_stars.len() < 3 {
                            break;
                        }
                        
                        let img_x: Vec<f64> = new_matched_stars.iter().map(|((x, _), _)| *x).collect();
                        let img_y: Vec<f64> = new_matched_stars.iter().map(|((_, y), _)| *y).collect();
                        let cat_x: Vec<f64> = new_matched_stars.iter().map(|(_, (x, _))| *x).collect();
                        let cat_y: Vec<f64> = new_matched_stars.iter().map(|(_, (_, y))| *y).collect();
                        
                        let new_coeffs_x = solve_projection(&cat_x, &img_x, &img_y);
                        let new_coeffs_y = solve_projection(&cat_y, &img_x, &img_y);
                        let new_transform = TransformCoefficients::new(new_coeffs_x, new_coeffs_y);
                        
                        let scale_change = (new_transform.scale() - iter_transform.scale()).abs();
                        let rotation_change = (new_transform.rotation() - iter_transform.rotation()).abs();
                        
                        iter_coeffs_x = new_coeffs_x;
                        iter_coeffs_y = new_coeffs_y;
                        iter_transform = new_transform;
                        final_match_count = new_matched_stars.len();
                        
                        if scale_change < 0.0001 && rotation_change < 0.0001 && iteration > 0 {
                            break;
                        }
                    }
                    
                    // Calculate new optical axis position
                    let (new_ra, new_dec) = calculate_image_center_equatorial(
                        &iter_transform,
                        new_center_ra_deg.to_radians(),
                        new_center_dec_deg.to_radians(),
                    );
                    
                    // Check for convergence
                    let delta_ra = (new_ra - refined_ra).to_degrees().abs();
                    let delta_dec = (new_dec - refined_dec).to_degrees().abs();
                    
                    if verbose {
                        println!("  New center: RA={:.6}°, Dec={:.6}° (Δ={:.4}°, {:.4}°)",
                                 new_ra.to_degrees(), new_dec.to_degrees(), delta_ra, delta_dec);
                    }
                    
                    // Update results
                    final_coeffs_x = Some(iter_coeffs_x);
                    final_coeffs_y = Some(iter_coeffs_y);
                    final_transform = Some(iter_transform);
                    refined_ra = new_ra;
                    refined_dec = new_dec;
                    
                    // Converged if position changed by less than 0.01°
                    if delta_ra < 0.01 && delta_dec < 0.01 {
                        if verbose {
                            println!("  Converged!");
                        }
                        break;
                    }
                }
            }
        }
    }

    Ok(PlateSolvingResult {
        matched_quads_count: final_match_count,
        coeffs_x: final_coeffs_x,
        coeffs_y: final_coeffs_y,
        transform: final_transform,
        optical_axis_ra: refined_ra,
        optical_axis_dec: refined_dec,
        spiral_iterations: solution_iteration,
        solution_ra_deg: best_search_ra_deg,
        solution_dec_deg: best_search_dec_deg,
    })
}
