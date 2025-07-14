use rawloader::{decode_file, RawImageData};
use std::path::Path;
use std::collections::VecDeque; 
use std::fmt;
use itertools::Itertools;


pub const STAR_THREESHOLD: f32 = 0.2;
pub const PIXEL_SIZE_MICRON: f64 = 6.0;
pub const TELESCOPE_FOCAL_LENGHT: f64 = 1200.0;

const BITDEPTH: u16 = 14; 

pub type StarPosXY = (f64, f64);

#[derive(Debug, Clone)]
pub struct StarQuad {
    pub stars : [StarPosXY; 4],
    pub barycenter: (f64, f64),
    pub largest_distance: f64,
    pub normalized_distances: [f64; 6],
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
        write!(f, "]\n")
    }
}




impl StarQuad {
    /**
     * Creates a new `StarQuad` from a given vector of stars.
     *
     * This function calculates the properties of the star quad. It first computes the
     * barycenter (geometric center) of the provided stars. Then, it calculates all
     * pairwise distances between the stars, identifies the largest distance, and uses it
     * to normalize all other distances. The normalized distances are then sorted in
     * descending order.
     *
     * # Arguments
     *
     * * `stars` - A vector of `StarPosXY` tuples, where each tuple represents the (x, y) coordinates of a star.
     *
     * # Returns
     *
     * A new `StarQuad` instance containing the original stars, their barycenter, the largest
     * pairwise distance, and a sorted vector of normalized distances.
     */
    pub fn new(stars : [StarPosXY; 4]) -> Self {

        let n = stars.len();

        let x = stars.iter().map(|(x, _)| *x).sum::<f64>() / n as f64;
        let y = stars.iter().map(|(_, y)| *y).sum::<f64>() / n as f64;

        let distances : [f64; 6] = [
            distance_between_stars(&stars[0], &stars[1]),
            distance_between_stars(&stars[0], &stars[2]),
            distance_between_stars(&stars[0], &stars[3]),

            distance_between_stars(&stars[1], &stars[2]),
            distance_between_stars(&stars[1], &stars[3]),

            distance_between_stars(&stars[2], &stars[3])
        ];
        
        

        let biggest_distance = distances.into_iter().fold(0.0, f64::max);

        let mut normalized_distances: [f64; 6] = distances.map(|d| d / biggest_distance);

        normalized_distances.sort_by(|a, b| b.partial_cmp(a).unwrap());

        StarQuad {
            barycenter: (x, y),
            largest_distance: biggest_distance,
            normalized_distances: normalized_distances,
            stars : stars,
        }
        
    }

    pub fn compare(&self, other: &StarQuad, tolerance: f64) -> bool {
        let matches = self.normalized_distances.iter()
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
    /// Creates a new StarGraph connecting each star to its `num_neighbors` nearest neighbors.
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
                .take(4)
                .collect::<Vec<usize>>();
            adj[i] = neighbors;
        }

        StarGraph {
            stars: stars.to_vec(),
            adj,
        }
    }

    /// Finds all stars at a distance of n hops from a given start_node_idx.
    /// Uses Breadth-First Search (BFS).
    ///
    /// # Arguments
    /// * `start_node_idx` - The index of the starting star in the graph's `stars` vector.
    /// * `n` - The number of hops (edges) away from the start node.
    ///
    /// # Returns
    /// A `Vec<StarPosXY>` containing the stars at exactly distance `n`.
    /// Returns `None` if `start_node_idx` is out of bounds.
    pub fn find_stars_at_distance_n(&self, start_node_idx: usize, n: usize) -> Option<Vec<StarPosXY>> {
        if start_node_idx >= self.stars.len() {
            return None;
        }

        if n == 0 {
            return Some(vec![self.stars[start_node_idx]]);
        }

        let mut queue: VecDeque<(usize, usize)> = VecDeque::new();
        queue.push_back((start_node_idx, 0));

        let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
        visited.insert(start_node_idx);

        let mut result_indices = Vec::new();

        while let Some((current_node, distance)) = queue.pop_front() {
            if distance == n {
                result_indices.push(current_node);
                continue; // Don't explore further from this node
            }

            if distance > n {
                break; // Optimization: no need to check deeper nodes
            }

            for &neighbor in &self.adj[current_node] {
                if !visited.contains(&neighbor) {
                    visited.insert(neighbor);
                    queue.push_back((neighbor, distance + 1));
                }
            }
        }

        Some(result_indices.into_iter().map(|idx| self.stars[idx]).collect())
    }
}




pub fn get_image_size(file_path: &Path) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let raw_image = decode_file(file_path)?;
    Ok((raw_image.width, raw_image.height))
}

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

pub fn distance_between_stars(star1: &StarPosXY, star2: &StarPosXY,) -> f64 {
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






pub fn returns_all_star_quads(stars: &[StarPosXY], n: usize) -> Vec<StarQuad> {
    let mut res = Vec::new();

    let distance_matrix = calculate_distance_matrix(stars);

    for row in distance_matrix.iter() {
        let clossest_stars = row.iter().enumerate().sorted_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());

        let nth_clossest = clossest_stars.map(|(idx, _)| stars[idx]).take(4).collect::<Vec<_>>().try_into().unwrap();

        res.push(StarQuad::new(nth_clossest));
    }

    res
}

