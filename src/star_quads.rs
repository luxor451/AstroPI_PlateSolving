use rawloader::{decode_file, RawImageData};
use std::path::Path;
use std::collections::VecDeque; 
use plotters::prelude::*;


use std::fmt;
pub const STAR_THREESHOLD: f32 = 0.2;
pub const PIXEL_SIZE_MICRON: f64 = 6.0;
pub const TELESCOPE_FOCAL_LENGHT: f64 = 1200.0;




const BITDEPTH: u16 = 14; 


type StarBarycenter = (f64, f64);

#[derive(Debug, Clone)]
pub struct StarQuad {
    pub stars : Vec<StarBarycenter>,
    pub barycenter: (f64, f64),
    pub largest_distance: f64,
    pub normalized_distances: Vec<f64>,
}

pub fn max_from_vec(vec: &Vec<f64>) -> f64 {
    let mut res = 0.0;
    for &value in vec {
        if value > res {
            res = value;
        }
    };

    return res
}

impl fmt::Display for StarQuad {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "  Barycenter: ({:.2}, {:.2})", self.barycenter.0, self.barycenter.1)?;
        writeln!(f, "  Largest Distance: {:.2}", self.largest_distance)?;
        write!(f, "  Stars: [")?;
        for (i, star) in self.stars.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "({:.2}, {:.2})", star.0, star.1)?;
        }
        writeln!(f, "]")?;
        write!(f, "  Normalized Distances: [")?;
        for (i, dist) in self.normalized_distances.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:.4}", dist)?;
        }
        write!(f, "]")
    }
}



impl StarQuad {
    pub fn new(stars : Vec<StarBarycenter>) -> Self {

        let n = stars.len();

        let x = stars.iter().map(|(x, _)| *x).sum::<f64>() / n as f64;
        let y = stars.iter().map(|(_, y)| *y).sum::<f64>() / n as f64;

        let mut distances = Vec::new();
        
        for i in 0..n {
            for j in (i + 1)..n {
                distances.push(distance_between_stars(&stars[i], &stars[j]));
            }
        }

        let biggest_distance = max_from_vec(&distances);

        let mut normalized_distances: Vec<f64> = distances.iter().map(|&d| d / biggest_distance).collect();

        normalized_distances.sort_by(|a, b| b.partial_cmp(a).unwrap());

        StarQuad {
            barycenter: (x, y),
            largest_distance: biggest_distance,
            normalized_distances: normalized_distances,
            stars : stars,
        }
        
    }

}





pub fn plot_histogram_to_png(
    bins_data: &Vec<(&u16, &u32)>,
    output_path: &str,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if bins_data.is_empty() {
        eprintln!("No data to plot for histogram.");
        return Ok(());
    }

    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_count = bins_data.iter().map(|&(_, count)| *count).max().unwrap_or(1) as f64;
    // Determine the max bin start value for x-axis range
    let max_bin_start = bins_data.last().map_or(16300, |&(bin_val, _)| *bin_val);


    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(80)
        .build_cartesian_2d(0u32..(max_bin_start as u32 + 100u32), 0f64..max_count * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Pixel Value Bins")
        .y_desc("Count")
        .axis_desc_style(("sans-serif", 15))
        .draw()?;

    chart.draw_series(
        bins_data.iter().map(|&(bin_start, count)| {
            let x0 = *bin_start as u32;
            let x1 = (bin_start + 99) as u32; // Bin width is 100
            let y = *count as f64;
            let mut rect = Rectangle::new([(x0, 0.0), (x1, y)], BLUE.filled());
            rect.set_margin(0,0,1,1); // Add small margin between bars
            rect
        }),
    )?;

    root.present()?;
    println!("Histogram saved to {}", output_path);
    Ok(())
}

pub fn get_image_size(file_path: &Path) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    // Attempt to decode the RAW file (e.g., DNG)
    let raw_image = decode_file(file_path)?;

    // Return the width and height of the image
    Ok((raw_image.width, raw_image.height))
}

/// Reads a DNG (or other RAW) into a flat list of pixels whose
/// (x, y) coordinates are measured from the center of the frame.
/// 
/// Returned Vec has entries `(xc, yc, value)` where
///   xc = column_index – (width/2)
///   yc = row_index    – (height/2)
/// so that `(0,0)` is the exact center of the image.
pub fn get_pixel_matrix_from_dng(
    file_path: &Path,
) -> Result<Vec<(i32, i32, u16)>, Box<dyn std::error::Error>> {
    // decode DNG into raw image
    let raw = decode_file(file_path)?;
    let w = raw.width as usize;
    let h = raw.height as usize;
    let tot = w.checked_mul(h).ok_or("image too large")?;
    if tot == 0 {
        return Ok(Vec::new());
    }

    // pull out u16 data (integer or float→scaled)
    let flat: Vec<u16> = match raw.data {
        RawImageData::Integer(v) => v,
        RawImageData::Float(v) => {
            let maxval = (2u32.pow(BITDEPTH as u32) - 1) as f32;
            v.into_iter()
             .map(|f| (f.clamp(0.0, 1.0) * maxval) as u16)
             .collect()
        }
    };

    if flat.len() < tot {
        return Err("Incomplete image data".into());
    }

    // compute center offsets
    let cx = (w / 2) as i32;
    let cy = (h / 2) as i32;

    // build Vec<(x_centered, y_centered, value)>
    let mut out = Vec::with_capacity(tot);
    for r in 0..h {
        for c in 0..w {
            let idx = r * w + c;
            let val = flat[idx];
            out.push(((c as i32) - cx, (r as i32) - cy, val));
        }
    }

    Ok(out)
}


pub fn save_pixel_matrix_to_png(matrix: &Vec<Vec<u16>>,output_path: &str) -> Result<(), image::ImageError> {
    if matrix.is_empty() {
        return Err(image::ImageError::Parameter(
            image::error::ParameterError::from_kind(
                image::error::ParameterErrorKind::DimensionMismatch,
            ),
        ));
    }

    let height = matrix.len();
    let width = matrix[0].len();

    if width == 0 { 
        return Err(image::ImageError::Parameter(
            image::error::ParameterError::from_kind(
                image::error::ParameterErrorKind::DimensionMismatch,
            ),
        ));
    }

    // Flatten the Vec<Vec<u16>> into Vec<u16> for ImageBuffer.
    let mut flat_data: Vec<u16> = Vec::with_capacity(width * height);
    for row in matrix {
        if row.len() != width {
            return Err(image::ImageError::Parameter(
                image::error::ParameterError::from_kind(
                    image::error::ParameterErrorKind::DimensionMismatch, 
                ),
            ));
        }
        flat_data.extend_from_slice(row);
    }

    // Stretch pixel values by a factor of 4 and clamp to u16::MAX.
    for pixel_value in flat_data.iter_mut() {
        *pixel_value = (*pixel_value as u32 * 4).min(u16::MAX as u32) as u16;
    }

    // Create an ImageBuffer for Luma<u16> (16-bit grayscale).
    // `ImageBuffer::from_raw` takes ownership of `flat_data`.
    // It returns Some(image_buffer) if successful, None if data length doesn't match dimensions.
    let image_buffer = match image::ImageBuffer::<image::Luma<u16>, Vec<u16>>::from_raw(width as u32, height as u32, flat_data) {
        Some(buffer) => buffer,
        None => {
            // This case implies flat_data.len() != width * height * num_channels,
            // which should be caught by earlier checks or indicate an internal logic error.
            return Err(image::ImageError::Parameter(
                image::error::ParameterError::from_kind(
                    image::error::ParameterErrorKind::DimensionMismatch, // Data length issue
                ),
            ));
        }
    };

    // Save the image buffer to a PNG file.
    // The image crate infers the format from the file extension.
    image_buffer.save(output_path)?;

    println!("PNG image saved to {}", output_path);
    Ok(())
}


pub fn calculate_star_barycenters(
    pixels: &[(i32, i32, u16)],
    width: usize,
    height: usize,
) -> Vec<StarBarycenter> {
    if pixels.is_empty() {
        return Vec::new();
    }

    let star_threshold = ((2.0_f32.powi(BITDEPTH as i32) - 1.0) * STAR_THREESHOLD) as u16;

    // Create a map for quick lookup of pixel values by centered coordinates
    let pixel_map: std::collections::HashMap<(i32, i32), u16> =
        pixels.iter().map(|&(x, y, val)| ((x, y), val)).collect();

    let mut visited: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();
    let mut barycenters = Vec::new();

    let cx = (width / 2) as i32;
    let cy = (height / 2) as i32;

    for r_idx in 0..height {
        for c_idx in 0..width {
            let x = c_idx as i32 - cx;
            let y = r_idx as i32 - cy;
            let coord = (x, y);

            if !visited.contains(&coord) {
                if let Some(&pixel_value) = pixel_map.get(&coord) {
                    if pixel_value > star_threshold {
                        // Found a new potential star, start BFS
                        let mut weighted_x_sum_star = 0.0;
                        let mut weighted_y_sum_star = 0.0;
                        let mut total_mass_star = 0.0;
                        let mut star_pixels = 0;

                        let mut queue = VecDeque::new();
                        queue.push_back(coord);
                        visited.insert(coord);

                        while let Some((curr_x, curr_y)) = queue.pop_front() {
                            if let Some(&val) = pixel_map.get(&(curr_x, curr_y)) {
                                let mass = val as f64;
                                weighted_x_sum_star += curr_x as f64 * mass;
                                weighted_y_sum_star += curr_y as f64 * mass;
                                total_mass_star += mass;
                                star_pixels += 1;

                                // Check 8-connectivity neighbors
                                for dy in -1..=1 {
                                    for dx in -1..=1 {
                                        if dx == 0 && dy == 0 {
                                            continue;
                                        }
                                        let neighbor_coord = (curr_x + dx, curr_y + dy);
                                        if !visited.contains(&neighbor_coord) {
                                            if let Some(&neighbor_val) = pixel_map.get(&neighbor_coord) {
                                                if neighbor_val > star_threshold {
                                                    visited.insert(neighbor_coord);
                                                    queue.push_back(neighbor_coord);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if total_mass_star > 0.0 && star_pixels > 2 { // Ensure star is at least 3 pixels
                            barycenters.push((
                                weighted_x_sum_star / total_mass_star,
                                weighted_y_sum_star / total_mass_star,
                            ));
                        }
                    }
                }
            }
        }
    }
    barycenters
}

pub fn distance_between_stars(star1: &StarBarycenter, star2: &StarBarycenter,) -> f64 {
    let dx = star1.0 - star2.0;
    let dy = star1.1 - star2.1;
    (dx * dx + dy * dy).sqrt()
}


pub fn calculate_distance_matrix(
    barycenters: &[(f64, f64)],
) -> Vec<Vec<f64>> {
    let n = barycenters.len();
    let mut distance_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            if i != j {
                distance_matrix[i][j] = distance_between_stars(&barycenters[i], &barycenters[j]);
            } else {
                distance_matrix[i][j] = std::f64::MAX;
            }
        }
    }
    distance_matrix
}


pub fn returns_all_star_quads(stars: &[StarBarycenter], number_of_star_by_quads: usize) -> Vec<StarQuad> {
    if stars.len() < number_of_star_by_quads {
        return Vec::new();
    }
    let distance_mat: Vec<Vec<f64>> = calculate_distance_matrix(stars);
    let n: usize = stars.len();

    let mut star_quads_vec: Vec<StarQuad> = Vec::new();

    for i in 0..n {
        let mut star_with_distances: Vec<(usize, f64)> = distance_mat[i]
            .iter()
            .enumerate()
            .filter(|(idx, _)| *idx != i)
            .map(|(idx, &dist)| (idx, dist))
            .collect();

        // Sort other stars by distance (ascending)
        star_with_distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Take the current star and the `number_of_star_by_quads - 1` closest ones
        let mut closest_stars: Vec<StarBarycenter> = Vec::with_capacity(number_of_star_by_quads);
        closest_stars.push(stars[i]); // Always include the current star

        for (star_idx, _) in star_with_distances.iter().take(number_of_star_by_quads - 1) {
            closest_stars.push(stars[*star_idx]);
        }
        
        if closest_stars.len() == number_of_star_by_quads {
            star_quads_vec.push(StarQuad::new(closest_stars));
        }
    }

    star_quads_vec
}


