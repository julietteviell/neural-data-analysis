neural-data-analysis

This repository contains all MATLAB code and processed data required to reproduce the figures of the manuscript. Each figure has its own folder with the scripts and the .mat data files needed to run it. Example data from a single mouse is provided so every script runs out of the box.

Repository structure

Each figN/ folder is self-contained:





Scripts (*.m): MATLAB scripts that load the processed data and generate the panels of the corresponding figure.



Functions (*.m): MATLAB functions that are needed for the corresponding scripts.



Data (*.mat): pre-processed data files. For each figure, the data from one example mouse is provided either in the main branch (m31), or in the figure folder when needed, so the scripts can be run without the full dataset.

Requirements





MATLAB R2024b or later

Usage





Clone the repository:

 git clone https://github.com/<user>/<repo>.git
 cd <repo>



Open MATLAB (R2024b+) and cd into the repository folder.



Run the script of the figure you want to reproduce, e.g.:

 run('fig1/fig1_script.m')

The script loads the example .mat files provided in the same folder and generates the figure panels. No additional setup is required for the example data.



Note: The provided data are from a single example mouse, meant to demonstrate the pipeline. Figures in the manuscript are based on the full dataset across animals; contact the corresponding author for access.

Data

The .mat files contain processed data (not raw recordings). Raw data are available from the corresponding author upon reasonable request.

License

This project is licensed under the MIT License — see LICENSE for details.

Contact

For questions about the code or data, please open an issue or contact the corresponding author.
