# Plate Solving

A Rust implementation of astronomical plate solving for determining the precise celestial coordinates of stars in astrophotography images.

## Overview

Plate solving is the process of matching star patterns in an image to known star catalogs, allowing you to determine exactly where your telescope is pointing. This tool analyzes DNG (Digital Negative) image files, detects stars, and matches their geometric patterns against the UCAC5 star catalog to compute astrometric solutions.

## How It Works

This implementation is based on the **star quad matching algorithm** described in [How Plate Solving Works](https://olegignat.com/how-plate-solving-works/) by Oleg Ignat.

### Algorithm Overview

1. **Star Detection**: Extract bright pixels from the DNG image and calculate star barycenters (centroids) using a flood-fill approach on connected regions above a brightness threshold.

2. **Quad Formation**: Group stars into "quads" (sets of 4 stars). For each star, find nearby neighbors using a graph-based search and form all possible 4-star combinations.

3. **Geometric Fingerprinting**: For each quad, compute:
   - The barycenter (center point)
   - All 6 pairwise distances between the 4 stars
   - Normalized distances (divided by the largest distance)
   
   These normalized distances create a scale-invariant "fingerprint" of the star pattern.

4. **Catalog Query**: Query the UCAC5 catalog via VizieR TAP service to retrieve stars in the estimated field of view.

5. **Catalog Quads**: Generate quads from the catalog stars using gnomonic projection to convert celestial coordinates (RA/Dec) to image plane coordinates.

6. **Pattern Matching**: Compare image quads against catalog quads by matching their normalized distance fingerprints within a tolerance threshold.

7. **Transformation Solving**: Use matched quad pairs to solve for the affine transformation coefficients between image coordinates and celestial coordinates using Singular Value Decomposition (SVD).

## Features

- **DNG Support**: Direct processing of raw Digital Negative files
- **UCAC5 Catalog**: Queries the comprehensive UCAC5 star catalog
- **Robust Matching**: Scale-invariant geometric pattern matching
- **Gnomonic Projection**: Accurate coordinate transformations
- **Graph-based Neighbor Search**: Efficient quad generation using nearest-neighbor graphs

## Installation

### Prerequisites

- Rust (edition 2024 or later)
- Internet connection (for catalog queries)

### Build

```bash
git clone https://github.com/luxor451/AstroPI
cargo build --release
```

## Usage

### Basic Usage

Edit `src/main.rs` to configure your target image and initial coordinate estimate:

```rust
let initial_coord = CoordinateEquatorial::new(
    RaHoursMinutesSeconds::new(14, 01, 12.5),  // RA: 14h 01m 12.5s
    Arcdegrees::new(54, 20, 56.0),              // Dec: +54° 20' 56"
);

let file_path = Path::new("test.dng");
```

Run the plate solver:

```bash
cargo run --release
```


Test the plate solver:

```bash
cargo test
```

This will create some debug png to see the quads and histogram of the image

### Output

The program will output:
- Number of detected stars and quads
- Number of matched quads between image and catalog
- Transformation coefficients (if enough matches are found)
- Success/failure status

Example output:
```
Starting plate solving for: test.dng
Initial coordinates: RA: 14h 1m 12.5s, Dec: +54° 20' 56.0"

Found 150 star barycenters in the image.
Found 2500 star quads in the image.
Found 3200 star quads from catalogue.
Found 45 matches between image and catalogue.

=== Plate Solving Results ===
Matched quads: 45
Transformation coefficients:
  X: (0.0012, -0.0003, 150.2)
  Y: (0.0003, 0.0012, -75.8)

Plate solving successful!
```

## Project Structure

```
src/
├── main.rs           # Main entry point
├── coordinate.rs     # Coordinate types (RA/Dec, projection functions)
├── star_quads.rs     # Star detection, quad formation, matching
├── parse_catalog.rs  # UCAC5 catalog queries via VizieR
├── solver.rs         # SVD-based transformation solver
├── platesolve.rs     # High-level plate solving workflow
├── printing.rs       # Visualization utilities
└── tests.rs          # Test suite
```

## Configuration

Key constants can be adjusted in the source files:

- `STAR_THRESHOLD` in `star_quads.rs`: Brightness threshold for star detection
- `MATCHED_TOLERANCE` in `main.rs`: Tolerance for quad matching (default: 0.0009)
- `PIXEL_SIZE_MICRON` and `TELESCOPE_FOCAL_LENGTH` in `star_quads.rs`: Camera/telescope specifications

## Technical Details

### Dependencies

- **rawloader**: DNG/RAW image decoding
- **nalgebra**: Linear algebra (SVD solver)
- **reqwest**: HTTP client for catalog queries
- **serde/serde_json**: JSON parsing
- **itertools**: Iterator utilities
- **image/imageproc**: Image processing utilities
- **plotters**: Visualization (optional)

### Coordinate Systems

- **Equatorial Coordinates**: Right Ascension (RA) and Declination (Dec) in the ICRS reference frame
- **Image Coordinates**: Pixel positions centered at the image center
- **Gnomonic Projection**: Maps the celestial sphere onto a flat tangent plane

### Performance Considerations

The current star detection algorithm uses a naive flood-fill approach and should be optimized for production use. Future improvements may include:
- More efficient connected component labeling
- Local maxima detection
- Hierarchical region growing
- GPU acceleration

## References

- [How Plate Solving Works](https://olegignat.com/how-plate-solving-works/) - Oleg Ignat's excellent explanation of the quad matching algorithm
- [UCAC5 Catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=I/322A) - USNO CCD Astrograph Catalog
- [VizieR TAP Service](http://tapvizier.u-strasbg.fr/) - Astronomical catalog query service
