## scgenokit: Single-Cell Genomic Data Analysis Toolkit

# 🔬 scgenokit: Single-Cell Genomic Data Analysis Toolkit

**scgenokit** is a Python-based toolkit designed to simplify and streamline the analysis of single-cell RNA sequencing (scRNA-seq) data. It provides a collection of modular and user-friendly functions for performing key steps in the scRNA-seq analysis pipeline, including:

- ✅ Quality control  
- ✅ Clustering  
- ✅ Dimensionality reduction  
- ✅ Differential expression analysis  
- ✅ Gene set enrichment analysis  

The goal of this toolkit is to make single-cell data analysis more **accessible** and **reproducible** for both biologists and bioinformaticians. Whether you're working with raw count data or pre-processed datasets, **scgenokit** helps you efficiently explore cellular heterogeneity and extract meaningful biological insights.

---

## 🚀 Quick Start

### Installation

```bash
pip install -e .
```

### Example Usage
```
from scgenokit.qc import filter_cells
from scgenokit.plot import plot_umap

# Example: filter cells with low gene counts
filtered_data = filter_cells(raw_data)

# Plot UMAP
plot_umap(filtered_data)
```


## 🛠️ Features

📦 Modular functions for each analysis step
📈 Built-in visualizations (e.g., UMAP, PCA)
🔎 Differential expression and enrichment tools
🔄 Compatible with AnnData and Scanpy workflows
🧪 Automated tests and high code quality


## 📂 Project Structure
scgenokit/
├── __init__.py
├── qc.py
├── clustering.py
├── plot.py
├── de.py

tests/
├── test_qc.py
├── test_clustering.py
├── test_plot.py
└── ...

## 👨‍💻 Developers

This project is developed by Yuejun Xu and Jingyun Liu as part of a data science software development course, with a focus on:

📁 Modular code design
🔄 Version control with Git
✅ Automated testing with pytest
📚 Clear documentation and reproducibility

## 🙌 Acknowledgments

Special thanks to the teaching team of the course for their guidance on software best practices and reproducible scientific computing.
