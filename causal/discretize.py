import pandas as pd
import scanpy as sc
from tqdm import tqdm
import json

NUM_LEVELS = 4

def process(name, source_path):
    adata = sc.read(source_path)
    data = adata.to_df("lognorm")
    
    # filter cells and genes with at least 5% non-zero values
    row_zeros = (data.values == 0).mean(axis=1)
    col_zeros = (data.values == 0).mean(axis=0)
    adata = adata[row_zeros < 0.95, col_zeros < 0.95].copy()
    data = adata.to_df("lognorm")

    # discretize expression levels by percentile
    levels = pd.CategoricalDtype(range(1 + NUM_LEVELS), ordered=True)
    binmaps = dict()
    for col_name in tqdm(data.columns, desc="Discretize"):
        nonzero_values = data.loc[data[col_name] > 0, col_name]
        nonzero_values, bins = pd.cut(nonzero_values, bins=NUM_LEVELS, labels=range(1, 1 + NUM_LEVELS), retbins=True, duplicates='drop')
        binmaps[col_name] = bins.tolist()
        data.loc[nonzero_values.index, col_name] = nonzero_values.astype(int)
        data[col_name] = data[col_name].astype(levels)

    # add demographic data
    for col_name in ["gravity", "age", "clustertype"]:
        col = adata.obs[col_name].astype('category')
        data[col_name] = col

    # save data
    json.dump(binmaps, open(f"data/binmaps_{name}.json", "w"), indent=4)
    data = data.reset_index(drop=True)
    data.to_feather(f"data/discretized_{name}.feather")

process("402", "../data/402-2.h5ad")
process("403", "../data/403-2.h5ad")
