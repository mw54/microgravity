import scanpy as sc
from tqdm import tqdm
import numpy as np

def process(name, source_path, target_path):
    adata = sc.read(source_path)

    # construct metadata
    celltypes = sorted(adata.obs["clustertype"].unique().tolist()) + ["Overall"]
    ages = sorted(adata.obs["age"].unique().tolist())
    metrics = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'mean', 'std', 'pct']
    indices = {
        "celltype": dict(zip(celltypes, range(len(celltypes)))),
        "age": dict(zip(ages, range(len(ages)))),
        "metric": dict(zip(metrics, range(len(metrics))))
    }

    dfs = list() # dimension [celltype, age, gene, metric]
    for celltype in tqdm(celltypes, desc=name):
        dfs.append(list())
        for age in ages:
            # stratify dataset by experimental group
            if celltype == "Overall":
                subset = adata[(adata.obs["age"] == age), :].copy()
            else:
                subset = adata[(adata.obs["clustertype"] == celltype) & (adata.obs["age"] == age), :].copy()
            
            # compute DEG w.r.t. experimental conditions
            sc.tl.rank_genes_groups(subset, groupby="gravity", reference="1G", method="wilcoxon", key_added=f"{celltype}-{age}", layer="lognorm", use_raw=False)
            df = sc.get.rank_genes_groups_df(subset, group=None, key=f"{celltype}-{age}").set_index("names", drop=True)
            df = df.loc[adata.var_names.values]

            # compute expression metrics
            data = subset.layers["lognorm"].toarray()
            df["mean"] = np.mean(data, axis=0)
            df["std"] = np.std(data, axis=0)
            df["pct"] = np.where(data > 0, 1, 0).mean(axis=0) # percentage of cells with expression > 0
            dfs[indices["celltype"][celltype]].append(df.loc[:, metrics].to_numpy())

    # permute dfs to satisfy adata.varm format
    dfs = np.array(dfs).transpose(2, 0, 1, 3) # dimension [gene, celltype, age, metrics]
    
    # store results and indices
    adata.varm["DEG"] = dfs
    adata.uns["DEG_indices"] = indices
    adata.write(target_path)

process("402", "../data/402-1.h5ad", "../data/402-2.h5ad")
process("403", "../data/403-1.h5ad", "../data/403-2.h5ad")
process("404", "../data/404-1.h5ad", "../data/404-2.h5ad")
process("405", "../data/405-1.h5ad", "../data/405-2.h5ad")
