use nalgebra::{DMatrix, DVector};

use crate::consts::*;

/// Solves for transformation coefficients mapping image coordinates to reference coordinates.
///
/// Uses SVD to solve the problem: `A*x = b`
///
/// # Arguments
///
/// * `ref_quad` - Reference coordinate values (catalog positions).
/// * `img_quad_x` - X coordinates from the image.
/// * `img_quad_y` - Y coordinates from the image.
///
/// # Returns
///
/// A tuple `(x_1, x_2, x_3)` representing the transformation coefficients.
pub fn solve_projection(
    ref_quad: &Vec<f64>,
    img_quad_x: &Vec<f64>,
    img_quad_y: &Vec<f64>,
) -> (f64, f64, f64) {
    let nb_of_row = img_quad_x.len();
    // Create the matrix A and vector b for the linear system Ax = b
    let mut data = Vec::with_capacity(nb_of_row * 3);
    for i in 0..nb_of_row {
        data.push(img_quad_x[i]);
        data.push(img_quad_y[i]);
        data.push(1.0);
    }
    let a = DMatrix::from_row_slice(nb_of_row, 3, &data);

    let b = DVector::from_column_slice(&ref_quad);

    // Solve Ax = b for x using SVD
    let svd = a.svd(true, true);
    let x = svd.solve(&b, SVD_TOLERANCE).expect("Failed to solve for x");

    return (x[0], x[1], x[2]);
}
