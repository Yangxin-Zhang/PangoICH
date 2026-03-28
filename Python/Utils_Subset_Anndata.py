
# Python/Utils_Subset_Anndata.py

import h5py
import numpy as np
import scanpy as sc
import pandas as pd
import anndata
from scipy.sparse import csr_matrix

def Subset_Anndata_From_H5ad_Py(H5ad_path,subset_var,var_value,subset_path):

    adata = sc.read_h5ad(H5ad_path, backed='r')

    cell_indices = np.where(adata.obs[subset_var] == var_value)[0]

    f = h5py.File(H5ad_path, "r")
    f_X = f["X"]

    n_genes = f_X.attrs.get('shape', [0, 0])[1]
    indptr = f_X["indptr"][:]

    new_data = []
    new_indices = []
    new_indptr = [0]

    for old_idx in cell_indices:
        start = indptr[old_idx]
        end = indptr[old_idx + 1]

        if end > start:
            cell_data = f_X["data"][start:end]
            cell_cols = f_X["indices"][start:end]

            new_data.extend(cell_data)
            new_indices.extend(cell_cols)

        new_indptr.append(len(new_data))

    X_subset = csr_matrix((new_data, new_indices, new_indptr),
                          shape=(len(cell_indices), n_genes))

    adata_subset = adata[cell_indices, :].to_memory()
    adata_subset.X = X_subset
    
    anndata.settings.allow_write_nullable_strings = True
    adata_subset.write(subset_path,compression='gzip')

  
