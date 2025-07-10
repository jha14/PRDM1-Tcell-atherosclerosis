# Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation

This repository contains the code and resources necessary to reproduce the main findings of the paper:  
**"Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation"**

## 🗂 Repository structure
```bash
.
├── data/        # Processed datasets (access via Zenodo)
├── codes/       # Numbered R scripts used in the analysis
├── results/     # Output plots
├── LICENSE
└── README.md
```
- **data/**: Processed data used in the analyses can be accessed via Zenodo:  
  🔗 [https://doi.org/10.5281/zenodo.15852128](https://doi.org/10.5281/zenodo.15852128)
  
- **codes/**: Contains ordered R scripts (e.g., `00_functions.R`, `01_wgcna.R`, etc.) for reproducing figures and results presented in the paper.

## ⚙️ Environment

The scripts were developed and tested in the following environment:

```
| Software / Package    | Version   |
|-----------------------|-----------|
| R                     | ≥ 3.6.0   |
| WGCNA                 | 1.73      |
| bnlearn               | 4.5       |
| clusterProfiler       | 3.12.0    |
| GOSemSim              | 2.20.0    |
| GENIE3                | 1.6.0     |
| minet                 | 3.42.0    |
| Seurat                | 5.1.0     |
```

## ▶️ How to run

1. Clone this repository:

   ```bash
   git clone https://github.com/jha14/PRDM1_Tcell_atherosclerosis.git
   cd PRDM1_Tcell_atherosclerosis

3. Download the processed data from Zenodo and place it in the data/ directory.
4. Run the scripts sequentially from the codes/ folder:

   ```bash
   source("codes/00_functions.R")
   source("codes/01_wgcna.R")
   ...

6. Output figures and tables will be saved to the results/ folder as specified.

## 📄 Citation
If you use this code or data in your research, please cite the original paper:

Han Jin, Sanne L. Maas, et al. (2025). Identification of a PRDM1-regulated T cell network to regulate atherosclerotic plaque inflammation. Genome Medicine. [DOI pending]

(*Han Jin and Sanne L. Maas contributed equally to this work*)
