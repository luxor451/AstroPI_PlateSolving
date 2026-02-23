use crate::consts::*;
use crate::star_quads::StarPosXY;
use std::collections::VecDeque;
use std::time::Instant;

/// Calculates the barycenters (centroids) of star regions in an image.
///
/// Uses flood-fill (BFS) on bright pixels to identify connected star regions
/// and computes the intensity-weighted centroid for each.
///
/// # Arguments
///
/// * `pixels` - A slice of tuples `(x, y, value)` representing pixel coordinates
///              (centered at image center) and their brightness values.
/// * `width` - The width of the image in pixels.
/// * `height` - The height of the image in pixels.
///
/// # Returns
///
/// A `Vec<StarPosXY>` containing the (x, y) coordinates of each detected star's barycenter,
/// with coordinates centered at the image center. Returns an empty vector if no stars are found
/// or if the input is empty.
pub fn calculate_star_barycenters(
    pixels: &[(i32, i32, u16)],
    width: usize,
    height: usize,
    cam: &CameraConfig,
) -> Vec<StarPosXY> {
    let solve_start = Instant::now();
    if pixels.is_empty() || width == 0 || height == 0 {
        return Vec::new();
    }

    let star_threshold = ((2.0_f32.powi(cam.bitdepth as i32) - 1.0) * STAR_THRESHOLD) as u16;
    
    let cx = (width / 2) as i32;
    let cy = (height / 2) as i32;

    // Use a flat grid for O(1) pixel lookups instead of HashMap
    // Grid stores pixel values, 0 means no data or below threshold
    let mut grid: Vec<u16> = vec![0; width * height];
    
    // Pre-filter: collect only bright pixels and their positions
    let mut bright_pixels: Vec<(usize, usize)> = Vec::new();
    
    for &(x, y, val) in pixels {
        // Convert from centered coordinates back to grid indices
        let col = (x + cx) as usize;
        let row = (y + cy) as usize;
        
        if col < width && row < height {
            let idx = row * width + col;
            grid[idx] = val;
            
            if val > star_threshold {
                bright_pixels.push((col, row));
            }
        }
    }

    // Track visited pixels using a bitset for memory efficiency
    let mut visited: Vec<bool> = vec![false; width * height];
    let mut barycenters = Vec::new();
    
    // Reusable queue to avoid allocations
    let mut queue: VecDeque<(usize, usize)> = VecDeque::with_capacity(256);

    // Only iterate over bright pixels, not all pixels
    for (start_col, start_row) in bright_pixels {
        let start_idx = start_row * width + start_col;
        
        if visited[start_idx] {
            continue;
        }

        // Found a new potential star, start BFS
        let mut weighted_x_sum = 0.0;
        let mut weighted_y_sum = 0.0;
        let mut total_mass = 0.0;
        let mut star_pixels = 0u32;

        queue.clear();
        queue.push_back((start_col, start_row));
        visited[start_idx] = true;

        while let Some((col, row)) = queue.pop_front() {
            let idx = row * width + col;
            let val = grid[idx];
            
            // Convert back to centered coordinates for barycenter calculation
            let x = col as f64 - cx as f64;
            let y = row as f64 - cy as f64;
            
            let mass = val as f64;
            weighted_x_sum += x * mass;
            weighted_y_sum += y * mass;
            total_mass += mass;
            star_pixels += 1;

            // Check 8-connectivity neighbors using bounds checking
            let row_min = row.saturating_sub(1);
            let row_max = (row + 1).min(height - 1);
            let col_min = col.saturating_sub(1);
            let col_max = (col + 1).min(width - 1);

            for nr in row_min..=row_max {
                for nc in col_min..=col_max {
                    if nr == row && nc == col {
                        continue;
                    }
                    
                    let neighbor_idx = nr * width + nc;
                    
                    if !visited[neighbor_idx] && grid[neighbor_idx] > star_threshold {
                        visited[neighbor_idx] = true;
                        queue.push_back((nc, nr));
                    }
                }
            }
        }

        if total_mass > 0.0 && star_pixels > MIN_STAR_PIXELS {
            barycenters.push((
                weighted_x_sum / total_mass,
                weighted_y_sum / total_mass,
            ));
        }
    }
    let solve_elapsed = solve_start.elapsed();
    println!("Stars found in: {:.3} s", solve_elapsed.as_secs_f64());
    barycenters
}
