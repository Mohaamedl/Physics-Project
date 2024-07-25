# Physics-Project

## Structured Light Generation Using Spatial Light Modulators

**Final Year Project in Physics**

**Author**: MOHAMED HADDADI

This study explores the generation and analysis of structured light beams, with a specific focus on Hermite-Gaussian, Laguerre-Gaussian, and Bessel-Gaussian beams, as well as their solutions carrying Orbital Angular Momentum (OAM). Using Spatial Light Modulators (SLMs), these beams are experimentally produced, and their spatial intensity is captured by a Laser Beam Profile (LBP). Experimental results are then compared with their theoretical counterparts. The primary objective of this project includes improving the purity of the generated beams to enhance their applicability in optical communication systems. Additionally, a graphical user interface has been developed, allowing for the creation of these beams, real-time characterization, and comparison with theoretical intensity profiles.

## Structure

- **Report**: The project report is available in the root directory of this repository:
  - [report.docx](report.docx)
  - [report.pdf](report.pdf)

- **Presentation**: The PowerPoint presentation used for the final project presentation is available:
  - [presentation.pptx](presentation.pptx)

- **Resources**:
  - **Books**: Refer to the `books` folder for the books used in this project.
  - **Papers**: The `papers` folder contains relevant research papers.
  - **Images**: Diagrams and images used in the report can be found in the `images` folder.

- **Data**: The `Data` folder includes intensity profile data collected during the project:
  - Data is organized into subfolders named by the collection dates.
  - Images collected are named using abbreviations: L-G, H-G, and B-G.

- **Source Code**: Scripts used to generate phase masks, send them to the SLM, simulate results, and obtain experimental results are located in the `source_code` folder. This folder is further divided into:
  - **`masks_gen`**: Scripts for generating phase masks for the mentioned beams.
  - **`GUI`**: The graphical user interface that integrates all scripts into an easy-to-use and intuitive interface.
  - **`results_images_gen`**: Scripts used to generate images of the results obtained, which are included in the report.

Each subfolder within `source_code` contains its own `README.md` file for detailed information on its contents and usage. Please refer to these files for more specific instructions.

## Installation and Usage

1. **Installation**: Follow the instructions in the respective `README.md` files within the `source_code` folder for installation and usage details of the scripts and GUI.

2. **Usage**: For detailed instructions on how to use the provided scripts and GUI, refer to the `README.md` files in the `masks_gen`, `GUI`, and `results_images_gen` folders.

## Contributing

For improvements or suggestions, please create a pull request or open an issue.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file in the parent directory for more details.

---

Enjoy exploring the structured light generation and analysis!
