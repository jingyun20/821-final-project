import scanpy as sc
import gseapy as gp


# ---------- Quality Control ----------
def qc_filtering(adata, min_genes=200, max_genes=2500, max_mt=5.0):
    """Filter cells based on gene count and mitochondrial gene percentage."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=[])

    adata = adata[adata.obs["n_genes_by_counts"] > min_genes, :]
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes, :]
    adata = adata[adata.obs["pct_counts_mt"] < max_mt, :]

    return adata


# ---------- Normalization ----------
def normalize_log(adata):
    """Normalize total counts per cell and log-transform the data."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


# ---------- Dimensionality Reduction ----------
def run_pca(adata, n_comps=50):
    """Run PCA on the normalized data."""
    max_comps = min(n_comps, adata.shape[0], adata.shape[1])

    # arpack solver (default in scanpy) fails when n_components == min(n_samples, n_features)
    if max_comps == min(adata.shape[0], adata.shape[1]):
        max_comps -= 1

    sc.pp.pca(adata, n_comps=max_comps)
    return adata


def run_umap(adata):
    """Run UMAP on the PCA-reduced data."""
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.1)
    return adata


# ---------- Clustering ----------
def leiden_clustering(adata, resolution=1.0):
    """Cluster the cells using the Leiden algorithm."""
    sc.tl.leiden(adata, resolution=resolution)
    return adata


# ---------- Differential Expression ----------
def differential_expression(adata, groupby="leiden", method="t-test"):
    """Run DE analysis across groups."""
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)
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
    )
    return enr.results


# ---------- Optional: Quick pipeline ----------
def run_full_pipeline(adata):
    """Run full pipeline from QC to clustering (example only)."""
    adata = qc_filtering(adata)
    adata = normalize_log(adata)
    adata = run_pca(adata)
    adata = run_umap(adata)
    adata = leiden_clustering(adata)
    return adata
