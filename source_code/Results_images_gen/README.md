# Experimental Results Image Generation

This folder contains scripts used to generate the images of experimental results shown in the report of this work. The scripts are organized into different categories based on the type of beams and masks. Each category includes scripts for generating the phase masks, expected results, and experimental results.

## Folder Structure

### 1. Laguerre-Gaussian Beams (LG)

- **Script**: `LG.m`
  - **Description**: Generates images showing the phase masks sent to the SLM, the expected results, and the experimental results for Laguerre-Gaussian beams.
- **Images**: Contains experimental result images related to Laguerre-Gaussian beams.

### 2. Hermite-Gaussian Beams (HG)

- **Script**: `HG.m`
  - **Description**: Generates images showing the phase masks sent to the SLM, the expected results, and the experimental results for Hermite-Gaussian beams.
- **Images**: Contains experimental result images related to Hermite-Gaussian beams.

### 3. Bessel-Gaussian Beams (BG)

- **Script**: `BG.m`
  - **Description**: Generates images showing the phase masks sent to the SLM, the expected results, and the experimental results for Bessel-Gaussian beams.
- **Images**: Contains experimental result images related to Bessel-Gaussian beams.

### 4. Forked Masks and Orbital Angular Momentum Beams (FM)

- **Scripts**:
  - **`FM_progression.m`**: Shows the results of the progression of the diffraction pattern separation into diffraction orders 0 and 1.
  - **`OAMs.m`**: Displays results obtained for beams of different orders of Orbital Angular Momentum (OAM).
- **Images**: Contains experimental result images related to forked masks and OAM beams.

### 5. Alignment (ALIGN)

- **Script**: `align.m`
  - **Description**: Generates images showing the results of linearization as the polarization of the beam varies.
- **Images**: Contains experimental result images related to polarization variation.

## How to Use

1. **Run the Scripts**: Execute the relevant `.m` file for the beam type or mask category of interest to generate the images. The scripts will produce images of the phase masks, expected results, and experimental results.
   
2. **View Results**: Check the corresponding folder for the generated images of the experimental results.

## Contributing

For any improvements or suggestions, please feel free to create a pull request or open an issue.

## License

This project is licensed under the MIT License. See the [LICENSE](../../LICENSE) file in the parent directory for more details.

