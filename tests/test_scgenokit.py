"""Unit test package for scgenokit."""

import pytest
import numpy as np
import pandas as pd
from anndata import AnnData

from scgenokit.scgenokit import (
    qc_filtering,
    normalize_log,
    run_pca,
    run_umap,
    leiden_clustering,
)


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for testing."""
    X = np.random.poisson(1, (20, 100)).astype(float)
    var_names = [f"Gene{i}" for i in range(100)]
    obs_names = [f"Cell{i}" for i in range(20)]

    var = pd.DataFrame(index=var_names)
    var.index.values[:5] = [f"MT-Gene{i}" for i in range(5)]

    adata = AnnData(X=X, var=var, obs=pd.DataFrame(index=obs_names))
    return adata


def test_qc_filtering(mock_adata):
    """Test the quality control filtering function."""
    adata_filtered = qc_filtering(mock_adata)
    assert adata_filtered.shape[0] <= mock_adata.shape[0]
    assert "pct_counts_mt" in adata_filtered.obs


def test_normalize_log(mock_adata):
    """Test the normalization and log transformation function."""
    adata = normalize_log(mock_adata.copy())
    assert np.all(adata.X >= 0)


def test_run_pca(mock_adata):
    """Test the PCA function."""
    adata = normalize_log(mock_adata.copy())
    adata = run_pca(adata)
    assert "X_pca" in adata.obsm


def test_run_umap(mock_adata):
    """Test the UMAP function."""
    adata = normalize_log(mock_adata.copy())
    adata = run_pca(adata)
    adata = run_umap(adata)
    assert "X_umap" in adata.obsm


def test_leiden_clustering(mock_adata):
    """Test the Leiden clustering function."""
    adata = normalize_log(mock_adata.copy())
    adata = run_pca(adata)
    adata = run_umap(adata)
    adata = leiden_clustering(adata)
    assert "leiden" in adata.obs
