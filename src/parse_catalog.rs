use std::error::Error;
use std::fs::File;
use std::path::Path;
use once_cell::sync::Lazy;
use crate::coordinate::*; // Assumes your Star struct is defined here

/// Represents a row in your local Gaia CSV catalogue.
/// Matches header: source_id,ra,dec,phot_g_mean_mag
#[derive(Debug, serde::Deserialize)]
struct LocalCatalogRow {
    ra: f64,
    dec: f64,
    phot_g_mean_mag: f64,
}

/// Queries a local CSV star catalog for stars within a field of view.
///
/// # Arguments
///
/// * `center_coordinate` - Center position (RA/Dec).
/// * `fov` - Field of view size in arcseconds.
/// * `nb_of_star` - Maximum number of stars to retrieve.
/// * `catalog_path` - Path to the local .csv file.
///
/// # Returns
///
/// A vector of `Star` objects sorted by brightness.
pub fn get_stars_from_catalogue(
    center_coordinate: &CoordinateEquatorial,
    fov: f64,
    nb_of_star: usize,
) -> Result<Vec<Star>, Box<dyn Error>> {
    let (center_ra, center_dec) = center_coordinate.to_degrees();
    
    // Convert FOV from arcseconds to degrees
    let fov_deg = fov / 3600.0;
    let half_fov = fov_deg / 2.0;

    log::debug!(
        "Searching local CSV for stars at RA: {:.6}, DEC: {:.6}, FOV: {:.6}°",
        center_ra, center_dec, fov_deg
    );

    // Load and cache the CSV into memory once (the file is sorted by ra, dec, mag).
    static CATALOG: Lazy<Vec<LocalCatalogRow>> = Lazy::new(|| {
        // Try several known locations, then fall back to a recursive search.
        let candidate_paths = [
            "AstroPI_PlateSolving/catalogue/gaia_local.csv",
            "catalogue/gaia_local.csv",
            "Gaia_dowloader/gaia_local.csv",
        ];

        let found_path = candidate_paths
            .iter()
            .map(Path::new)
            .find(|p| p.exists())
            .map(|p| p.to_path_buf())
            .or_else(|| {
                // Recursive search from cwd and common parent directories
                let search_roots = [
                    std::env::current_dir().ok(),
                    std::env::current_exe().ok().and_then(|p| p.parent().map(|d| d.to_path_buf())),
                ];
                for root in search_roots.iter().flatten() {
                    for entry in walkdir::WalkDir::new(root)
                        .max_depth(5)
                        .into_iter()
                        .filter_map(|e| e.ok())
                    {
                        if entry.file_name() == "gaia_local.csv" {
                            return Some(entry.into_path());
                        }
                    }
                }
                None
            });

        let path = match found_path {
            Some(p) => {
                println!("[catalog] Found catalog at: {}", p.display());
                p
            }
            None => {
                log::error!("Could not find gaia_local.csv in any known location");
                eprintln!("[catalog] ERROR: Could not find gaia_local.csv anywhere!");
                return Vec::new();
            }
        };

        let file = match File::open(&path) {
            Ok(f) => f,
            Err(e) => {
                log::error!("Failed to open catalog {}: {}", path.display(), e);
                return Vec::new();
            }
        };
        let mut rdr = csv::ReaderBuilder::new().has_headers(true).from_reader(file);
        let mut rows: Vec<LocalCatalogRow> = Vec::new();
        for res in rdr.deserialize() {
            match res {
                Ok(rec) => rows.push(rec),
                Err(e) => log::warn!("Skipping bad row in catalog: {}", e),
            }
        }
        println!("[catalog] Loaded {} stars from {}", rows.len(), path.display());
        rows
    });

    // We will collect candidates as tuples: (ra, dec, magnitude)
    let mut candidates: Vec<(f64, f64, f64)> = Vec::new();

    // Compute RA search interval(s), handling wraparound (0 - 360°)
    let ra_min = (center_ra - half_fov).rem_euclid(360.0);
    let ra_max = (center_ra + half_fov).rem_euclid(360.0);

    // Helper: lower_bound (first index with ra >= value)
    let lower_bound = |arr: &Vec<LocalCatalogRow>, value: f64| -> usize {
        match arr.binary_search_by(|r| r.ra.partial_cmp(&value).unwrap_or(std::cmp::Ordering::Equal)) {
            Ok(i) => i,
            Err(i) => i,
        }
    };

    // Helper: upper_bound (first index with ra > value)
    let upper_bound = |arr: &Vec<LocalCatalogRow>, value: f64| -> usize {
        match arr.binary_search_by(|r| r.ra.partial_cmp(&value).unwrap_or(std::cmp::Ordering::Equal)) {
            Ok(mut i) => {
                // advance past equal elements
                while i < arr.len() && (arr[i].ra - value).abs() < std::f64::EPSILON {
                    i += 1;
                }
                i
            }
            Err(i) => i,
        }
    };

    let arr = &*CATALOG;

    // If catalog failed to load, fall back to empty set
    if arr.is_empty() {
        log::warn!("Catalog is empty, returning no stars.");
        return Ok(Vec::new());
    }

    // Determine contiguous RA ranges to scan
    let ranges: Vec<(usize, usize)> = if ra_min <= ra_max {
        // single interval [ra_min, ra_max]
        let start = lower_bound(arr, ra_min);
        let end = upper_bound(arr, ra_max);
        vec![(start, end)]
    } else {
        // wrap-around: [ra_min, 360) and [0, ra_max]
        let start1 = lower_bound(arr, ra_min);
        let end1 = arr.len();
        let start2 = 0usize;
        let end2 = upper_bound(arr, ra_max);
        vec![(start1, end1), (start2, end2)]
    };

    // Scan the small subset(s) and filter by exact RA/Dec box
    for (start, end) in ranges {
        for rec in &arr[start..end] {
            // Dec filter
            if (rec.dec - center_dec).abs() > half_fov {
                continue;
            }

            // RA delta (circular)
            let mut delta_ra = (rec.ra - center_ra).abs();
            if delta_ra > 180.0 {
                delta_ra = 360.0 - delta_ra;
            }
            if delta_ra > half_fov {
                continue;
            }

            candidates.push((rec.ra, rec.dec, rec.phot_g_mean_mag));
        }
    }

    // 3. Sort by Magnitude (Brightest stars have lower magnitude values)
    // f64 doesn't implement Ord, so we use partial_cmp
    candidates.sort_by(|a, b| {
        a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal)
    });

    // 4. Transform to your internal Star struct and take the top N
    let result_stars: Vec<Star> = candidates
        .into_iter()
        .take(nb_of_star)
        .map(|(ra, dec, _mag)| Star { ra, dec })
        .collect();

    log::info!(
        "Found {} stars in local CSV catalog within FOV.", result_stars.len()
    );
    Ok(result_stars)
}