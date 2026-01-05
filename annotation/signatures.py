import scanpy as sc
import json

def process(name, source_path):
    adata = sc.read(source_path)
    df = sc.get.rank_genes_groups_df(adata, None, pval_cutoff=0.05, log2fc_min=0.5)
    summ = df.groupby('group').apply(lambda x: x.head(100)["names"].tolist()).to_dict()
    df.to_csv(f"data/signatures_{name}.csv", index=False)
    json.dump(summ, open(f"data/signatures_{name}.json", "w"))

process("402", "../data/402-0.h5ad")
process("403", "../data/403-0.h5ad")
process("404", "../data/404-0.h5ad")
process("405", "../data/405-0.h5ad")
