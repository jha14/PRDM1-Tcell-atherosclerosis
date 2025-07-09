# Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation

This repository contains the code and resources necessary to reproduce the main findings of the paper:  
**"Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation."**

## 🗂 Repository structure
├── data/ # Processed datasets (access via Zenodo)
├── codes/ # Numbered R scripts used in the analysis
└── figures/ # Output plots (optional, if applicable)


- **data/**: Processed data used in the analyses can be accessed via Zenodo:  
  🔗 [https://doi.org/10.5281/zenodo.15709332](https://doi.org/10.5281/zenodo.15709332)
- **codes/**: Contains ordered R scripts (e.g., `01_preprocessing.R`, `02_analysis.R`, etc.) for reproducing figures and results presented in the manuscript.

## ⚙️ Environment

The scripts were developed and tested in the following environment:

- R version: `>= 4.2.0`
- Required R packages:
  - `Seurat`
  - `dplyr`
  - `ggplot2`
  - `clusterProfiler`
  - `WGCNA`
  - `ComplexHeatmap`
  - *(add any other packages used)*


## ▶️ How to run

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/PRDM1-Tcell-atherosclerosis.git
   cd PRDM1-Tcell-atherosclerosis

2. Download the processed data from Zenodo and place it in the data/ directory.
3. Run the scripts sequentially from the codes/ folder:

  ```bash
  source("codes/01_preprocessing.R")
  source("codes/02_analysis.R")
  ...

4. Output figures and tables will be saved to the working directory or subfolders as specified.

📄 Citation
If you use this code or data in your research, please cite the original paper:

Han J, et al. (2025). Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation. Genome Medicine. [DOI pending]

BibTeX:
@article{han2025prdm1,
  title={Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation},
  author={Han, Jin and others},
  journal={Genome Medicine},
  year={2025},
  doi={10.1186/s13073-025-00xxx}
}
