use reqwest::blocking::Client;
use std::collections::HashMap;

use crate::coordinate::*;

pub fn get_stars_from_catalogue(
    center_coordinate: &CoordinateEquatorial,
    fov: f64,
    nb_of_star: usize,
) -> Result<Vec<Star>, Box<dyn std::error::Error>> {
    let client = Client::new();

    let (ra, dec) = center_coordinate.to_degrees();
    let fov = fov / 3600.0; // Convert arcseconds to degrees

    println!(
        "Searching for stars in FOV centered at RA: {:.6}, DEC: {:.6} with FOV: {:.6} degrees",
        ra, dec, fov
    );

    let query = format!(
        r#"
        SELECT TOP {nb_of_star} *
        FROM "I/322A/out"
        WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), BOX('ICRS', {ra}, {dec}, {fov}, {fov}))
        ORDER BY Vmag
        "#,
        nb_of_star = nb_of_star,
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
        .post("http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync")
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
