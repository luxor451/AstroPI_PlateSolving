use rawloader::{decode_file, RawImageData};
use std::path::Path;
use std::collections::VecDeque; 
use plotters::prelude::*;


use std::fmt;

pub const STAR_THREESHOLD: f32 = 0.2;
pub const PIXEL_SIZE_MICRON: f64 = 6.0;
pub const TELESCOPE_FOCAL_LENGHT: f64 = 1200.0;
pub const IMAGE_BINNING_FACTOR: u8 = 1;




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

pub fn get_pixel_matrix_from_dng(file_path: &Path) -> Result<Vec<Vec<u16>>, Box<dyn std::error::Error>> {
    // Attempt to decode the RAW file (e.g., DNG)
    let raw_image = decode_file(file_path)?;

    let width = raw_image.width;
    let height = raw_image.height;
    let total_pixels = width * height;

    // Handle images with zero area (e.g., width or height is 0)
    if total_pixels == 0 {
        // Return an empty matrix if the image has no pixels
        return Ok(Vec::new());
    }

    // Extract pixel data into a 1D Vec<u16>.
    // Integer data is used directly.
    // Float data is assumed to be normalized (0.0-1.0) and is scaled to u16 range [0, 16383].
    // This scaling (to 2^14 - 1) matches the histogram processing in the provided main function.
    let pixel_data_1d: Vec<u16> = match raw_image.data {
        RawImageData::Integer(data) => {
            if data.len() < total_pixels {
                return Err(format!(
                    "Integer pixel data is insufficient. Expected at least {} elements, but found {}.",
                    total_pixels,
                    data.len()
                )
                .into());
            }
            // Take the first 'total_pixels' elements, in case the buffer is larger.
            data.into_iter().take(total_pixels).collect()
        }
        RawImageData::Float(data) => {
            if data.len() < total_pixels {
                return Err(format!(
                    "Float pixel data is insufficient. Expected at least {} elements, but found {}.",
                    total_pixels,
                    data.len()
                )
                .into());
            }
            // Scale float values (assumed 0.0-1.0) to u16 range (0-16383) and collect.
            data.into_iter()
                .take(total_pixels)
                .map(|val| (val.clamp(0.0, 1.0) * (2.0_f32.powi(BITDEPTH as i32) - 1.0)) as u16)
                .collect()
        }
    };

    // Reshape the 1D pixel data (Vec<u16>) into a 2D matrix (Vec<Vec<u16>>)
    // The matrix will have 'height' rows, each row containing 'width' pixel values.
    let mut matrix: Vec<Vec<u16>> = Vec::with_capacity(height);
    for i in 0..height {
        let row_start_index = i * width;
        let row_end_index = row_start_index + width;
        // Create a new vector for the row from a slice of pixel_data_1d.
        // This is safe because pixel_data_1d is confirmed to have total_pixels elements.
        matrix.push(pixel_data_1d[row_start_index..row_end_index].to_vec());
    }

    Ok(matrix)


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

    // Stretch pixel values by a factor of 10 and clamp to u16::MAX.
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


pub fn calculate_star_barycenters(pixel_matrix: &Vec<Vec<u16>>) -> Vec<StarBarycenter> {
    let mut barycenters = Vec::new();
    if pixel_matrix.is_empty() || pixel_matrix[0].is_empty() {
        return barycenters;
    }

    let height = pixel_matrix.len();
    let width = pixel_matrix[0].len();
    let star_threshold = ((2.0_f32.powi(BITDEPTH as i32) - 1.0) * STAR_THREESHOLD) as u16; // Adjusted threshold for star detection
    let mut visited = vec![vec![false; width]; height];

    for r in 0..height {
        for c in 0..width {
            if pixel_matrix[r][c] > star_threshold && !visited[r][c] {
                // Found a new potential star, start BFS
                let mut weighted_x_sum_star = 0.0;
                let mut weighted_y_sum_star = 0.0;
                let mut total_mass_star = 0.0;
                
                let mut queue = VecDeque::new();
                queue.push_back((r, c));
                visited[r][c] = true;

                while let Some((curr_r, curr_c)) = queue.pop_front() {
                    let pixel_value = pixel_matrix[curr_r][curr_c];
                    // No need to check threshold again if only adding valid pixels to queue,
                    // but good for safety if logic changes.
                    // if pixel_value > star_threshold { // Already ensured by initial check and neighbor check
                    let mass = pixel_value as f64;
                    weighted_x_sum_star += curr_c as f64 * mass; // x is column index
                    weighted_y_sum_star += curr_r as f64 * mass; // y is row index
                    total_mass_star += mass;

                    // Check 8-connectivity neighbors
                    for dr in -1..=1 {
                        for dc in -1..=1 {
                            if dr == 0 && dc == 0 {
                                continue;
                            }

                            let nr = curr_r as isize + dr;
                            let nc = curr_c as isize + dc;

                            if nr >= 0 && nr < height as isize && nc >= 0 && nc < width as isize {
                                let nr_usize = nr as usize;
                                let nc_usize = nc as usize;
                                if pixel_matrix[nr_usize][nc_usize] > star_threshold && !visited[nr_usize][nc_usize] {
                                    visited[nr_usize][nc_usize] = true;
                                    queue.push_back((nr_usize, nc_usize));
                                }
                            }
                        }
                    }
                    // }
                }

                if total_mass_star > 0.0 {
                    barycenters.push((
                        weighted_x_sum_star / total_mass_star,
                        weighted_y_sum_star / total_mass_star,
                    ));
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

pub fn returns_all_star_quads(stars: &[StarBarycenter], number_of_star_by_quads : usize) -> Vec<StarQuad> {
    let distance_mat: Vec<Vec<f64>> = calculate_distance_matrix(stars);
    let n: usize = stars.len();

    let mut star_quads_vec: Vec<StarQuad> = Vec::new();

    for i in 0..n {
        let distance_to_current_star = &distance_mat[i];

        let mut star_with_distances: Vec<(StarBarycenter, f64)>  = Vec::new();

        for j in 0..n {
            star_with_distances.push((stars[j], distance_to_current_star[j]));
        }

        // Sort star_with_distances by distance (ascending)
        star_with_distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Take the first `number_of_star_by_quads` stars (excluding the current star itself)
        let mut closest_stars: Vec<StarBarycenter> = Vec::new();

        closest_stars.push(stars[i]); // Always include the current star

        for (star, _) in star_with_distances.iter().take(number_of_star_by_quads) {
                closest_stars.push(*star);
        }
        star_quads_vec.push(StarQuad::new(closest_stars));
    }

    star_quads_vec

}


