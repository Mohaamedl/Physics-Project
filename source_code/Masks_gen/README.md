# Phase Mask Generation for Optical Beams

This folder contains MATLAB functions to generate phase masks for various types of optical beams, including Laguerre-Gaussian, Hermite-Gaussian, and Bessel-Gaussian beams. Additionally, auxiliary functions are provided to support these main functions. 

## Overview

### Main Functions

1. **Laguerre-Gaussian Phase Masks**
   - **File**: `laguerre_gauss_phase_mask.m`
   - **Description**: Generates phase masks for Laguerre-Gaussian beams.

2. **Hermite-Gaussian Phase Masks**
   - **File**: `hermite_gauss_phase_mask.m`
   - **Description**: Generates phase masks for Hermite-Gaussian beams.

3. **Bessel-Gaussian Phase Masks**
   - **File**: `bessel_gauss_phase_mask.m`
   - **Description**: Generates phase masks for Bessel-Gaussian beams.

### Auxiliary Functions

1. **Generalized Laguerre Polynomial**
   - **File**: `fastLaguerre.m`
   - **Description**: Computes generalized Laguerre polynomials.

2. **2D Grid Creation**
   - **File**: `grid2D.m`
   - **Description**: Creates a 2D grid matrix to support the phase mask.

3. **Parameter Management**
   - **File**: `setParameters.m`
   - **Description**: Manages parameters for the functions.

### Testing

- **File**: `mainTest.m`
- **Description**: Contains code to test the functions and verify their functionality.

## How to Use

1. **Setup**: Ensure all required files are in the same directory or properly referenced in your MATLAB path.

2. **Generate Phase Masks**: Use the corresponding main function (`laguerre_gauss_phase_mask.m`, `hermite_gauss_phase_mask.m`, or `bessel_gauss_phase_mask.m`) to generate the desired phase masks. Adjust parameters as needed using `setParameters.m`.

3. **Test**: Run `mainTest.m` to see examples of how the functions can be used and verify their outputs.

## Contributing

If you have any improvements or suggestions, feel free to create a pull request or open an issue.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

Happy experimenting with optical beams!
