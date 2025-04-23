# 🔬 scgenokit: Single-Cell Genomic Data Analysis Toolkit

**scgenokit** is a Python-based toolkit designed to simplify and streamline the analysis of single-cell RNA sequencing (scRNA-seq) data. It provides a collection of modular and user-friendly functions for performing key steps in the scRNA-seq analysis pipeline, including:

- ✅ Quality control: Remove low-quality cells and genes to ensure reliable downstream analysis.  
- ✅ Clustering: Group similar cells based on expression profiles to reveal underlying cell types or states.  
- ✅ Dimensionality reduction: Reduce high-dimensional gene expression data into interpretable 2D or 3D space using PCA, UMAP, or t-SNE.  
- ✅ Differential expression analysis: Identify genes that are significantly up- or down-regulated between different clusters or conditions.  
- ✅ Gene set enrichment analysis: Highlight biological processes or pathways enriched in specific gene sets to interpret cell functions.  

The goal of this toolkit is to make single-cell data analysis more **accessible** and **reproducible** for both biologists and bioinformaticians. Whether you're working with raw count data or pre-processed datasets, **scgenokit** helps you efficiently explore cellular heterogeneity and extract meaningful biological insights.

---

## 📦 Installation

You can install the toolkit using `pip`:
```bash
pip install scgenokit
```
Or clone the repository directly:
```
git clone https://github.com/jingyun20/821-final-project.git
cd 821-final-project
pip install -e .
```

## 🚀 Usage
Here is a quick example:
```
from scgenokit.qc import filter_cells
from scgenokit.plot import plot_umap

adata = qc.load_and_filter("example_data.h5ad")
adata = clustering.run_leiden(adata)
qc.plot_qc_metrics(adata)
reduce_dim.plot_umap(adata)
```

## 🛠️ Features

- 📦 Modular functions for each analysis step
- 📈 Built-in visualizations (e.g., UMAP, PCA)
- 🔎 Differential expression and enrichment tools
- 🔄 Compatible with AnnData and Scanpy workflows
- 🧪 Automated tests and high code quality


## 🧪 Testing
To run tests:
```
pytest
```

## 👨‍💻 Developers

This project is developed by **Yuejun Xu** and **Jingyun Liu** as part of a software tools for data science (BIOSTAT 821) course, with a focus on:

- Modular code design
- Version control with Git
- Automated testing with pytest
- Clear documentation and reproducibility

## 🙌 Acknowledgments

Special thanks to the teaching team of the course for their guidance on software best practices and reproducible scientific computing. We are especially grateful to Dr. Patrick Wang (Instructor), Jason (Teaching Assistant), and Jonathan (Teaching Assistant) for their thoughtful feedback and technical mentorship, which have been invaluable in shaping this project.
