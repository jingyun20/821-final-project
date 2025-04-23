import scanpy as sc
import gseapy as gp


# ---------- Quality Control ----------
def qc_filtering(adata, min_genes=200, max_genes=2500, max_mt=5.0):
    """Filter cells based on gene count and mitochondrial gene percentage."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")  # O(m)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, percent_top=[]
    )  # O(nm)

    adata = adata[adata.obs["n_genes_by_counts"] > min_genes, :]  # O(n)
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes, :]  # O(n)
    adata = adata[adata.obs["pct_counts_mt"] < max_mt, :]  # O(n)

    return adata  # Total: O(nm)


# ---------- Normalization ----------
def normalize_log(adata):
    """Normalize total counts per cell and log-transform the data."""
    sc.pp.normalize_total(adata, target_sum=1e4)  # O(nm)
    sc.pp.log1p(adata)  # O(nm)
    return adata  # Total: O(nm)


# ---------- Dimensionality Reduction ----------
def run_pca(adata, n_comps=50):
    """Run PCA on the normalized data."""
    max_comps = min(n_comps, adata.shape[0], adata.shape[1])  # O(1)

    # arpack solver (default in scanpy) fails when n_components == min(n_samples, n_features)
    if max_comps == min(adata.shape[0], adata.shape[1]):  # O(1)
        max_comps -= 1

    sc.pp.pca(adata, n_comps=max_comps)  # O(nm * k), where k = n_comps
    return adata


def run_umap(adata):
    """Run UMAP on the PCA-reduced data."""
    sc.pp.neighbors(adata, n_neighbors=15)  # O(n^2)
    sc.tl.umap(adata, min_dist=0.1)  # O(n log n)
    return adata


# ---------- Clustering ----------
def leiden_clustering(adata, resolution=1.0):
    """Cluster the cells using the Leiden algorithm."""
    sc.tl.leiden(adata, resolution=resolution)  # O(n log n)
    return adata


# ---------- Differential Expression ----------
def differential_expression(adata, groupby="leiden", method="t-test"):
    """Run DE analysis across groups."""
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)  # O(nm)
    return adata


# ---------- Gene Set Enrichment Analysis ----------
def run_gsea(gene_list, gene_sets="KEGG_2016", organism="Human"):
    """Run GSEA using gseapy."""
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism=organism,
        outdir=None,
        no_plot=True,
    )  # O(g * s)
    return enr.results


# ---------- Optional: Quick pipeline ----------
def run_full_pipeline(adata):
    """Run full pipeline from QC to clustering (example only)."""
    adata = qc_filtering(adata)  # O(nm)
    adata = normalize_log(adata)  # O(nm)
    adata = run_pca(adata)  # O(nm * k)
    adata = run_umap(adata)  # O(n^2 + n log n)
    adata = leiden_clustering(adata)  # O(n log n)
    return adata


# ---------- Overall Summary ----------
# This module provides a complete scRNA-seq preprocessing and analysis workflow,
# including quality control, normalization, dimensionality reduction (PCA/UMAP),
# clustering (Leiden), differential expression analysis, and gene set enrichment (GSEA).
#
# Time complexity of the full pipeline (excluding GSEA):
#   - Dominated by UMAP’s graph construction: O(n^2)
#   - Total: O(nm + n^2 + n log n)
#
# Where:
#   n = number of cells
#   m = number of genes
#   k = number of PCA components (typically << m)
#   g = number of genes in gene list (GSEA)
#   s = number of gene sets (GSEA)
