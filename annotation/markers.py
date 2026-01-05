import scanpy as sc
import constants
import json

def process(name, source_path, target_path):
    adata = sc.read(source_path)
    clustertypes = {str(vv): k for k, v in constants.clustertypes[name].items() for vv in v}
    adata.obs["clustertype"] = adata.obs["leiden"].map(clustertypes)

    sc.tl.rank_genes_groups(adata, groupby="clustertype", method="wilcoxon", layer="lognorm", use_raw=False, key_added="celltype_markers")
    df = sc.get.rank_genes_groups_df(adata, None, pval_cutoff=0.05, log2fc_min=0.5, key="celltype_markers")
    df.to_csv(f"data/markers_{name}.csv", index=False)
    summ = df.groupby('group').apply(lambda x: x.head(50)["names"].tolist()).to_dict()
    json.dump(summ, open(f"data/markers_{name}.json", "w"))
    adata.write(target_path)

process("402", "../data/402-0.h5ad", "../data/402-1.h5ad")
process("403", "../data/403-0.h5ad", "../data/403-1.h5ad")
process("404", "../data/404-0.h5ad", "../data/404-1.h5ad")
process("405", "../data/405-0.h5ad", "../data/405-1.h5ad")
