use itertools::Itertools;
use rawloader::{RawImageData, decode_file};
use std::collections::VecDeque;
use std::fmt;
use std::path::Path;

pub const STAR_THREESHOLD: f32 = 0.15;
pub const PIXEL_SIZE_MICRON: f64 = 6.0;
pub const TELESCOPE_FOCAL_LENGHT: f64 = 1200.0;

const BITDEPTH: u16 = 14;

pub type StarPosXY = (f64, f64);

#[derive(Debug, Clone)]
pub struct StarQuad {
    pub stars: [StarPosXY; 4],
    pub barycenter: (f64, f64),
    pub largest_distance: f64,
    pub normalized_distances: [f64; 6],
}

impl fmt::Display for StarQuad {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "  Barycenter: ({:.2}, {:.2})",
            self.barycenter.0, self.barycenter.1
        )?;
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
        write!(f, "]\n")
    }
}

impl StarQuad {
    /// Creates a new `StarQuad` from four stars.
    ///
    /// Computes the barycenter, all pairwise distances, and normalizes distances
    /// by the largest distance. Normalized distances are sorted in descending order
    /// for use in quad matching.
    ///
    /// # Arguments
    ///
    /// * `stars` - An array of 4 `StarPosXY` tuples representing star coordinates.
    ///
    /// # Returns
    ///
    /// A new `StarQuad` instance with computed properties.
    pub fn new(stars: [StarPosXY; 4]) -> Self {
        let n = stars.len();

        let x = stars.iter().map(|(x, _)| *x).sum::<f64>() / n as f64;
        let y = stars.iter().map(|(_, y)| *y).sum::<f64>() / n as f64;

        let distances: [f64; 6] = [
            distance_between_stars(&stars[0], &stars[1]),
            distance_between_stars(&stars[0], &stars[2]),
            distance_between_stars(&stars[0], &stars[3]),
            distance_between_stars(&stars[1], &stars[2]),
            distance_between_stars(&stars[1], &stars[3]),
            distance_between_stars(&stars[2], &stars[3]),
        ];

        let biggest_distance = distances.into_iter().fold(0.0, f64::max);

        let mut normalized_distances: [f64; 6] = distances.map(|d| d / biggest_distance);

        normalized_distances.sort_by(|a, b| b.partial_cmp(a).unwrap());

        StarQuad {
            barycenter: (x, y),
            largest_distance: biggest_distance,
            normalized_distances: normalized_distances,
            stars: stars,
        }
    }

    /// Compares two star quads for similarity based on normalized distances.
    ///
    /// # Arguments
    ///
    /// * `other` - The star quad to compare against.
    /// * `tolerance` - Maximum allowed difference between normalized distances.
    ///
    /// # Returns
    ///
    /// `true` if at least 5 out of 6 normalized distances match within tolerance.
    pub fn compare(&self, other: &StarQuad, tolerance: f64) -> bool {
        let matches = self
            .normalized_distances
            .iter()
            .zip(other.normalized_distances.iter())
            .filter(|(d1, d2)| (*d1 - *d2).abs() <= tolerance)
            .count();

        matches >= 5
    }
}

pub struct StarGraph {
    pub stars: Vec<StarPosXY>,
    pub adj: Vec<Vec<usize>>,
}

impl StarGraph {
    /// Creates a new StarGraph connecting each star to its 3 nearest neighbors.
    ///
    /// # Arguments
    ///
    /// * `stars` - Slice of star positions to build the graph from.
    ///
    /// # Returns
    ///
    /// A `StarGraph` with adjacency lists for nearest neighbor connections.
    pub fn new(stars: &[StarPosXY]) -> Self {
        let n = stars.len();
        if n == 0 {
            return StarGraph {
                stars: Vec::new(),
                adj: Vec::new(),
            };
        }

        let distance_matrix = calculate_distance_matrix(stars);
        let mut adj = vec![Vec::new(); n];

        for i in 0..n {
            let neighbors = distance_matrix[i]
                .iter()
                .enumerate()
                .sorted_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(idx, _)| idx)
                .take(3)
                .collect::<Vec<usize>>();
            adj[i] = neighbors;
        }

        StarGraph {
            stars: stars.to_vec(),
            adj,
        }
    }

    /// Finds all stars within n hops from a starting star using BFS.
    ///
    /// # Arguments
    ///
    /// * `start_node_idx` - Index of the starting star.
    /// * `n` - Maximum number of hops (graph edges) to traverse.
    ///
    /// # Returns
    ///
    /// `Some(Vec<StarPosXY>)` containing stars within distance n, or `None` if index is invalid.
    pub fn find_stars_within_distance(
        &self,
        start_node_idx: usize,
        n: usize,
    ) -> Option<Vec<StarPosXY>> {
        if start_node_idx >= self.stars.len() {
            return None;
        }

        let mut queue: VecDeque<(usize, usize)> = VecDeque::new();
        queue.push_back((start_node_idx, 0));

        let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
        visited.insert(start_node_idx);

        let mut result_indices = Vec::new();
        result_indices.push(start_node_idx); // Add the starting node (distance 0)

        while let Some((current_node, distance)) = queue.pop_front() {
            if distance >= n {
                continue; // Don't explore further from this node
            }

            for &neighbor in &self.adj[current_node] {
                if !visited.contains(&neighbor) {
                    visited.insert(neighbor);
                    result_indices.push(neighbor);
                    queue.push_back((neighbor, distance + 1));
                }
            }
        }

        Some(
            result_indices
                .into_iter()
                .map(|idx| self.stars[idx])
                .collect(),
        )
    }
}

/// Extracts image dimensions from a DNG file.
///
/// # Arguments
///
/// * `file_path` - Path to the DNG file.
///
/// # Returns
///
/// A tuple `(width, height)` representing the image dimensions in pixels.
pub fn get_image_size(file_path: &Path) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let raw_image = decode_file(file_path)?;
    Ok((raw_image.width, raw_image.height))
}

/// Extracts pixel data from a DNG file with coordinates centered at the image center.
///
/// # Arguments
///
/// * `file_path` - Path to the DNG file.
///
/// # Returns
///
/// A vector of `(x, y, value)` tuples where coordinates are centered at the image center.
pub fn get_pixel_matrix_from_dng(
    file_path: &Path,
) -> Result<Vec<(i32, i32, u16)>, Box<dyn std::error::Error>> {
    let raw = decode_file(file_path)?;
    let w = raw.width as usize;
    let h = raw.height as usize;
    let tot = w.checked_mul(h).ok_or("image too large")?;
    if tot == 0 {
        return Ok(Vec::new());
    }

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

    let cx = (w / 2) as i32;
    let cy = (h / 2) as i32;

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


/// Calculates the barycenters (centroids) of star regions in an image.
///
///
/// **Note:** This is a naive approach that should be refined for better performance.
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
) -> Vec<StarPosXY> {
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
                                            if let Some(&neighbor_val) =
                                                pixel_map.get(&neighbor_coord)
                                            {
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

                        if total_mass_star > 0.0 && star_pixels > 2 {
                            // Ensure star is at least 3 pixels
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

/// Calculates the Euclidean distance between two stars.
///
/// # Arguments
///
/// * `star1` - First star position (x, y).
/// * `star2` - Second star position (x, y).
///
/// # Returns
///
/// The distance between the two stars.
pub fn distance_between_stars(star1: &StarPosXY, star2: &StarPosXY) -> f64 {
    let dx = star1.0 - star2.0;
    let dy = star1.1 - star2.1;
    (dx * dx + dy * dy).sqrt()
}

/// Computes the pairwise distance matrix for a set of star positions.
///
/// # Arguments
///
/// * `barycenters` - Slice of star positions.
///
/// # Returns
///
/// A 2D matrix where `matrix[i][j]` is the distance from star i to star j.
/// Diagonal elements are set to `f64::MAX`.
pub fn calculate_distance_matrix(barycenters: &[(f64, f64)]) -> Vec<Vec<f64>> {
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

/// Generates all possible star quads from a set of stars using graph-based neighbor search.
///
/// # Arguments
///
/// * `stars` - Slice of star positions.
/// * `n` - Maximum number of graph hops to search for neighboring stars.
///
/// # Returns
///
/// A vector of all `StarQuad` instances formed from combinations of nearby stars.
pub fn returns_all_star_quads(stars: &[StarPosXY], n: usize) -> Vec<StarQuad> {
    let mut res = Vec::new();

    let star_graph = StarGraph::new(stars);

    for i in 0..star_graph.stars.len() {
        if let Some(nearby_stars) = star_graph.find_stars_within_distance(i, n) {
            let combination = nearby_stars.iter().combinations(4);
            for quad in combination {
                if quad.len() == 4 {
                    let star_quad = StarQuad::new([*quad[0], *quad[1], *quad[2], *quad[3]]);
                    res.push(star_quad);
                }
            }
        }
    }
    res
}
