use reqwest::blocking::Client;
use std::collections::HashMap;
use std::sync::OnceLock;
use std::time::Duration;

use crate::consts::*;
use crate::coordinate::*;

/// Global HTTP client for VizieR queries - reuses connections across requests
static VIZIER_CLIENT: OnceLock<Client> = OnceLock::new();

/// Returns a shared HTTP client configured for VizieR queries.
/// The client maintains a connection pool and reuses connections for better performance.
fn get_vizier_client() -> &'static Client {
    VIZIER_CLIENT.get_or_init(|| {
        Client::builder()
            .pool_max_idle_per_host(VIZIER_POOL_MAX_IDLE_PER_HOST)
            .pool_idle_timeout(Duration::from_secs(VIZIER_POOL_IDLE_TIMEOUT_SECS))
            .timeout(Duration::from_secs(VIZIER_REQUEST_TIMEOUT_SECS))
            .tcp_keepalive(Duration::from_secs(VIZIER_TCP_KEEPALIVE_SECS))
            .build()
            .expect("Failed to create HTTP client")
    })
}

/// Queries the UCAC5 star catalog via VizieR TAP service for stars within a field of view.
///
/// # Arguments
///
/// * `center_coordinate` - Center position (RA/Dec) of the search region.
/// * `fov` - Field of view size in arcseconds.
/// * `nb_of_star` - Maximum number of stars to retrieve, ordered by brightness (Vmag).
///
/// # Returns
///
/// A vector of `Star` objects containing RA/Dec coordinates in degrees.
pub fn get_stars_from_catalogue(
    center_coordinate: &CoordinateEquatorial,
    fov: f64,
    nb_of_star: usize,
) -> Result<Vec<Star>, Box<dyn std::error::Error>> {
    let client = get_vizier_client();

    let (ra, dec) = center_coordinate.to_degrees();
    let fov = fov / 3600.0; // Convert arcseconds to degrees

    log::debug!(
        "Searching for stars in FOV centered at RA: {:.6}, DEC: {:.6} with FOV: {:.6} degrees",
        ra, dec, fov
    );

    let query = format!(
        r#"
        SELECT TOP {nb_of_star} *
        FROM "{catalog}"
        WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), BOX('ICRS', {ra}, {dec}, {fov}, {fov}))
        ORDER BY Vmag
        "#,
        nb_of_star = nb_of_star,
        catalog = VIZIER_CATALOG_TABLE,
        ra = ra,
        dec = dec,
        fov = fov
    );

    let mut form = HashMap::new();
    form.insert("REQUEST", "doQuery");
    form.insert("LANG", "ADQL");
    form.insert("FORMAT", "json");
    form.insert("QUERY", &query);

    let json: serde_json::Value = client
        .post(VIZIER_TAP_URL)
        .form(&form)
        .send()?
        .json()?;

    let metadata = json["metadata"].as_array().unwrap();
    let ra_idx = metadata
        .iter()
        .position(|m| m["name"] == "RAJ2000")
        .unwrap();
    let dec_idx = metadata
        .iter()
        .position(|m| m["name"] == "DEJ2000")
        .unwrap();

    let mut res = Vec::with_capacity(nb_of_star);

    if let Some(data) = json["data"].as_array() {
        for entry in data {
            let ra = entry[ra_idx].as_f64().unwrap();
            let dec = entry[dec_idx].as_f64().unwrap();
            res.push(Star { ra, dec });
        }
    }

    Ok(res)
}
