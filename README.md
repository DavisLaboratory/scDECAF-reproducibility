# scDECAF Reproducibility Repository

This repository contains notebooks and supporting materials to reproduce analyses from the manuscript:

**[Identification of cell types, states and programs by learning gene set representations](https://www.biorxiv.org/content/10.1101/2023.09.08.556842v1)**

It is designed to complement the main software package [scDECAF](https://github.com/DavisLaboratory/scDECAF) by providing end-to-end, executable examples.

---

## Software Requirements

For installation instructions and environment setup, please refer to the main [scDECAF repository](https://github.com/DavisLaboratory/scDECAF).

---

## Example Notebooks

The following Jupyter notebooks illustrate how to reproduce key results and analyses:

* **Pathway/gene signature screening and scoring**
  [Kang et al. (2018): 25K PBMC single cells](https://github.com/DavisLaboratory/scDECAF-reproducibility/blob/master/kang_pbmc/kang_pbmc.ipynb)
* **Optimization of sparsity operator** 
  [Experimentation with the shrinkage penalty and gene set screening results in Kang et al. (2018)](https://github.com/DavisLaboratory/scDECAF-reproducibility/blob/master/kang_pbmc/sparse_mode_effect.ipynb)
* **PMBC COVID-19 analysis**
  [Combining reference atlas mapping and Milo analysis with scDECAF gene set screening](https://github.com/DavisLaboratory/scDECAF-reproducibility/blob/master/cite_pbmc/TotalVI_scDECAF_analysis-addMilo.ipynb)
* **Drug2cell analysis with scDECAF**
  [Running scDECAF with pre-computed Drug2cell scores in HECOA Organoid Atlas]

More examples will be added as the manuscript evolves.

---

## Usage

1. Clone this repository:

   ```bash
   git clone https://github.com/DavisLaboratory/scDECAF-reproducibility.git
   cd scDECAF-reproducibility
   ```
2. Ensure that you have installed `scDECAF` and its dependencies (see [installation guide](https://github.com/DavisLaboratory/scDECAF)).
3. Open the notebooks with Jupyter Lab or Jupyter Notebook:

   ```bash
   jupyter lab
   ```
4. Run through the cells to reproduce figures and results.

---

## Data Availability

Data required for running the notebooks can be accessed from the sources cited in the paper (e.g., Kang et al. 2018). Where possible, direct download links are provided inside the notebooks.

---

## Citation

If you use these materials, please cite:

> Hediyehzadeh, Whitfield et al., *Identification of cell types, states and programs by learning gene set representations*, bioRxiv (2023). [DOI: 10.1101/2023.09.08.556842](https://doi.org/10.1101/2023.09.08.556842)

---

## License

This repository follows the same license as the main [scDECAF repository](https://github.com/DavisLaboratory/scDECAF).
