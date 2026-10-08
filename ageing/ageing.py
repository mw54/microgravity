import scanpy as sc
import pandas as pd
import numpy as np

def process(name, source_path):
    adata = sc.read(source_path)

    # derive instrinsic ageing signature
    subset = adata[(adata.obs["gravity"] == "1g"), :].copy()
    subset.obs.value_counts()
    sc.tl.rank_genes_groups(subset, groupby="age", reference="young", method="wilcoxon", key_added="ageing", layer="lognorm", use_raw=False)
    intrinsic = sc.get.rank_genes_groups_df(subset, group=None, key="ageing").set_index("names", drop=True)
    intrinsic = intrinsic.loc[adata.var_names.values]
    lognorm = subset.layers["lognorm"].toarray()
    intrinsic["mean"] = np.mean(lognorm, axis=0)
    intrinsic["std"] = np.std(lognorm, axis=0)
    intrinsic["pct"] = np.where(lognorm > 0, 1, 0).mean(axis=0)

    # load gravity DEGs
    celltypes = adata.uns["DEG_indices"]["celltype"]
    ages = adata.uns["DEG_indices"]["age"]
    metrics = adata.uns["DEG_indices"]["metric"]
    degs = adata.varm["DEG"]

    # determine ageing signatures
    ageing = pd.DataFrame({
        "scores": intrinsic["scores"],
        "pvals_adj": intrinsic["pvals_adj"],
        "intrinsic": intrinsic["logfoldchanges"],
        "old": degs[:,celltypes["Overall"],ages["old"],metrics["logfoldchanges"]],
        "young": degs[:,celltypes["Overall"],ages["young"],metrics["logfoldchanges"]],
    }, index=intrinsic.index)

    # save data
    intrinsic.to_csv(f"data/intrinsic_{name}.csv", index=True)
    ageing.to_csv(f"data/ageing_{name}.csv", index=True)


process("402", "../data/402-2.h5ad")
process("403", "../data/403-2.h5ad")
process("404", "/groups/xlu/zgu4/osd404/data/404-2.h5ad")
process("405", "/groups/xlu/zgu4/osd405/data/405-2.h5ad")
