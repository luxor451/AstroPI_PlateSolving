fn print_matrix(matrix: &[[f64; 3]; 3]) {
    for row in matrix.iter() {
        for &value in row.iter() {
            print!("{:8.2} ", value);
        }
        println!();
    }
}