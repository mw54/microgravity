import pandas as pd
import scanpy as sc
import gseapy as gp
from tqdm import tqdm
import json

gene_sets = ["data/m2.cp.v2025.1.Mm.symbols.gmt"]

def process(name, source_path):
    adata = sc.read(source_path)

    # get metadata
    celltypes = sorted(adata.obs["clustertype"].unique().tolist()) + ["Overall"]
    ages = sorted(adata.obs["age"].unique().tolist())

    dfs = list()
    for celltype in tqdm(celltypes, desc=name):
        for age in ages:
            # stratify dataset by experimental group
            if celltype == "Overall":
                subset = adata[adata.obs["age"] == age, :].copy()
            else:
                subset = adata[(adata.obs["clustertype"] == celltype) & (adata.obs["age"] == age), :].copy()

            # ensure the reference is 1G
            subset.obs['gravity'] = pd.Categorical(subset.obs['gravity'], categories=["uG", "1G"], ordered=True)

            # compute GSEA
            results = gp.gsea(data=subset.to_df("lognorm").transpose(), # row -> genes, column-> samples
                gene_sets=gene_sets,
                cls=subset.obs.gravity,
                permutation_num=1000,
                permutation_type='phenotype',
                outdir=f"data/{name}/{celltype}-{age}",
                method='s2n', # signal_to_noise
                threads=32,
                no_plot=True
            )

            # format results
            df = results.res2d.copy()
            df["Celltype"] = celltype
            df["Age"] = age
            dfs.append(df)

    # format and save concatenated results
    dfs = pd.concat(dfs)
    dfs["Gene Set"] = dfs.apply(lambda row: row["Term"].split('__')[0], axis=1)
    dfs["Database"] = dfs.apply(lambda row: row["Term"].split("__")[1].split("_")[0], axis=1)
    dfs["Term"] = dfs.apply(lambda row: " ".join(row["Term"].split("__")[1].split("_")[1:]), axis=1)
    dfs["|NES|"] = dfs.apply(lambda row: abs(row["NES"]), axis=1)
    dfs["Direction"] = "ns"
    dfs.loc[(dfs["FDR q-val"] < 0.05) & (dfs["NES"] > 1.5), "Direction"] = "up"
    dfs.loc[(dfs["FDR q-val"] < 0.05) & (dfs["NES"] < -1.5), "Direction"] = "down"
    dfs = dfs.sort_values(by=["Celltype", "Age", "Direction", "|NES|"], ascending=[True, True, False, False])
    dfs.to_csv(f"data/gsea_{name}.csv", index=False)

process("402", "../data/402-2.h5ad")
process("403", "../data/403-2.h5ad")
process("404", "../data/404-2.h5ad")
process("405", "../data/405-2.h5ad")
