use rawloader::{decode_file, RawImageData};
use std::f64::MAX;
use std::path::Path;
use std::collections::{HashMap, VecDeque}; // Added VecDeque
use plotters::prelude::*;
use std::fs::File;
use std::io::{self, Write, BufWriter};
use image::Rgb; // For defining the red color
use imageproc::drawing::{draw_hollow_circle_mut, draw_line_segment_mut}; // For drawing circles

use std::fmt;

const BITDEPTH: u16 = 14; // Assuming 14-bit depth for DNG files, adjust as 

const STAR_THREESHOLD: f32 = 0.2;

type StarBarycenter = (f64, f64);

#[derive(Debug, Clone)]
struct StarQuad {
    stars : Vec<StarBarycenter>,
    barycenter: (f64, f64),
    largest_distance: f64,
    normalized_distances: Vec<f64>,
}

fn max_from_vec(vec: &Vec<f64>) -> f64 {
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





fn plot_histogram_to_png(
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


fn get_pixel_matrix_from_dng(file_path: &Path) -> Result<Vec<Vec<u16>>, Box<dyn std::error::Error>> {
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


fn save_pixel_matrix_to_png(matrix: &Vec<Vec<u16>>,output_path: &str) -> Result<(), image::ImageError> {
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


fn calculate_star_barycenters(pixel_matrix: &Vec<Vec<u16>>) -> Vec<StarBarycenter> {
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

fn distance_between_stars(star1: &StarBarycenter, star2: &StarBarycenter,) -> f64 {
    let dx = star1.0 - star2.0;
    let dy = star1.1 - star2.1;
    (dx * dx + dy * dy).sqrt()
}


fn calculate_distance_matrix(
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

fn print_matrix(matrix: &Vec<Vec<f64>>) {
    for row in matrix {
        for value in row {
            print!("{:.2} ", value);
        }
        println!();
    }
}

fn returns_all_star_quads(stars: &[StarBarycenter], number_of_star_by_quads : usize) -> Vec<StarQuad> {
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




fn main() -> Result<(), Box<dyn std::error::Error>> {






    let file_path = Path::new("test.dng"); // Replace with your DNG file path

    // Get the pixel matrix from the DNG file
    let pixel_matrix = get_pixel_matrix_from_dng(file_path)?;

    println!("Nb of pixel above 16000: {}", pixel_matrix.iter().flatten().filter(|&&x| x > 16000).count());

    println!("Calculating star barycenters...");
    let star_barycenters = calculate_star_barycenters(&pixel_matrix);
    println!("Found {} star barycenters.", star_barycenters.len());


    println!("calculating star quads...");
    let star_quads = returns_all_star_quads(&star_barycenters, 3);
    println!("Found {} star quads.", star_quads.len());

    let print = false;
    if print {
        let output_image_path = "output_image.png";

        // First, save the grayscale image from the pixel matrix.
        // This ensures "output_image.png" exists and has the up-to-date content from pixel_matrix.
        save_pixel_matrix_to_png(&pixel_matrix, output_image_path)?;

        if !star_barycenters.is_empty() {
            // Now, open the saved image, draw circles, and save it again.
            let dynamic_image = image::open(output_image_path)?;
            // Convert to Rgb8. The original image is Luma16. 
            // to_rgb8() will convert it to an 8-bit RGB image, suitable for drawing colored shapes.
            let mut img_rgb = dynamic_image.to_rgb8(); 
            
            let red_color = Rgb([255u8, 0u8, 0u8]);
            let circle_radius = 20; // Radius of the circle in pixels; adjust as needed for visibility.

            println!("Drawing circles for {} detected star(s) on {}.", star_barycenters.len(), output_image_path);
            for (bary_x, bary_y) in &star_barycenters {
                // Barycenter coordinates (x, y) from calculation are (column_index, row_index).
                // imageproc's draw_hollow_circle_mut expects (x, y) where x is horizontal (column) and y is vertical (row).
                let center_x = bary_x.round() as i32;
                let center_y = bary_y.round() as i32;
                
                draw_hollow_circle_mut(&mut img_rgb, (center_x, center_y), circle_radius, red_color);
            }

            // Draw lines for star quads
            let white_color = Rgb([255u8, 255u8, 255u8]);
            println!("Drawing lines for {} star quad(s).", star_quads.len());
            for quad in &star_quads {
                if quad.stars.len() > 1 {
                    let central_star = quad.stars[0];
                    let central_point = (central_star.0 as f32, central_star.1 as f32);

                    for neighbor_star in quad.stars.iter().skip(1) {
                        let neighbor_point = (neighbor_star.0 as f32, neighbor_star.1 as f32);
                        draw_line_segment_mut(
                            &mut img_rgb,
                            central_point,
                            neighbor_point,
                            white_color,
                        );
                    }
                }
            }

            // Save the modified image, overwriting the previous grayscale version.
            img_rgb.save(output_image_path)?;
            println!("Successfully drew circles and quad lines on {} and re-saved.", output_image_path);
        } else {
            println!("No star barycenters to draw, so {} remains unchanged after initial save.", output_image_path);
        }



        // // Create a histogram of pixel values
        // let mut histogram: HashMap<u16, u32> = HashMap::new();
        // for row in &pixel_matrix {
        //     for &pixel_value in row {
        //         *histogram.entry(pixel_value).or_insert(0) += 1;
        //     }
        // }

        // // Convert histogram to a sorted vector of tuples (bin_start, count)
        // let mut bins_data: Vec<(&u16, &u32)> = histogram.iter().collect();
        // bins_data.sort_by_key(|&(bin_start, _)| *bin_start);

        // // Plot the histogram to a PNG file
        // plot_histogram_to_png(&bins_data, "histogram.png", "Pixel Value Histogram")?;
    }
    Ok(())
}