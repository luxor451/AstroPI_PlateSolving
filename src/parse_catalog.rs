use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Represents a star with its coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Star {
    /// Right Ascension in degrees.
    pub ra: f64,
    /// Declination in degrees.
    pub dec: f64,
}

/// Finds stars from a catalog file within a specified square field of view.
///
/// # Arguments
///
/// * `catalog_path` - The path to the catalog.dat file.
/// * `center_ra` - The Right Ascension of the center of the FOV in degrees.
/// * `center_dec` - The Declination of the center of the FOV in degrees.
/// * `fov_arcsec` - The width and height of the square FOV in arcseconds.
///
/// # Returns
///
/// A `Result` containing a vector of `Star` structs that are within the FOV,
/// or an `io::Error` if the file cannot be read.
pub fn find_stars_in_fov(
    catalog_path: &str,
    center_ra: f64,
    center_dec: f64,
    fov_arcsec: f64,
) -> io::Result<Vec<Star>> {
    let path = Path::new(catalog_path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let fov_deg = fov_arcsec / 3600.0;
    let half_fov_deg = fov_deg / 2.0;

    // Define the search area
    let min_ra = center_ra - half_fov_deg;
    let max_ra = center_ra + half_fov_deg;
    let min_dec = center_dec - half_fov_deg;
    let max_dec = center_dec + half_fov_deg;

    let mut stars_in_fov = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('|').collect();

        // The format seems to have RA at index 2 and Dec at index 3
        if parts.len() > 3 {
            if let (Ok(ra), Ok(dec)) = (parts[2].trim().parse::<f64>(), parts[3].trim().parse::<f64>()) {
                // Check if the star is within the square FOV
                if ra >= min_ra && ra <= max_ra && dec >= min_dec && dec <= max_dec {
                    stars_in_fov.push(Star { ra, dec });
                }
            }
        }
    }

    Ok(stars_in_fov)
}




#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_find_stars_in_fov() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "9537 00375 1| |10.001|-20.001|...").unwrap(); // Inside
        writeln!(file, "9537 00379 1| |10.002|-19.999|...").unwrap(); // Inside
        writeln!(file, "9537 00380 1| |10.1|-20.0|...").unwrap();     // Outside (RA too high)
        writeln!(file, "9537 00387 1| |9.9|-20.0|...").unwrap();      // Outside (RA too low)
        writeln!(file, "9537 00388 1| |10.0|-19.9|...").unwrap();     // Outside (Dec too high)
        writeln!(file, "9537 00389 1| |10.0|-20.1|...").unwrap();     // Outside (Dec too low)
        writeln!(file, "bad line").unwrap();                         // Invalid line

        let catalog_path = file.path().to_str().unwrap();
        let center_ra = 10.0;
        let center_dec = -20.0;
        // 1 degree = 3600 arcseconds. 0.01 deg = 36 arcsec.
        let fov_arcsec = 36.0; // +/- 0.005 degrees from center

        let stars = find_stars_in_fov(catalog_path, center_ra, center_dec, fov_arcsec).unwrap();

        assert_eq!(stars.len(), 2);
        assert_eq!(stars[0], Star { ra: 10.001, dec: -20.001 });
        assert_eq!(stars[1], Star { ra: 10.002, dec: -19.999 });
    }
}