import scanpy as sc
import anndata as ad
import pandas as pd
import constants

def process(name, target_path):
    adatas = dict()
    counts = dict()

    # load data
    for sample, path in constants.samples[name].items():
        adata = sc.read_10x_mtx(path)
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
        adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
        adatas[sample] = adata

    # concatenate all data
    adata = ad.concat(adatas, label="sample")
    adata.raw = adata.copy()
    adata.obs['group'] = adata.obs['sample'].map(constants.groups)
    adata.obs["gravity"] = adata.obs["group"].map(constants.gravity)
    adata.obs["age"] = adata.obs["group"].map(constants.age)
    adata.obs_names_make_unique()
    counts["raw"] = adata.obs["sample"].value_counts()

    # filter droplets
    adata.obs["conditional"] = False
    adata.obs.loc[
        (adata.obs["total_counts"] < constants.thresholds["upper_total_counts"]) &
        (adata.obs["pct_counts_mt"] < constants.thresholds["upper_pct_mt"]) &
        (adata.obs["pct_counts_ribo"] < constants.thresholds["upper_pct_ribo"]) &
        (adata.obs["pct_counts_hb"] < constants.thresholds["upper_pct_hb"]),
        "conditional"
    ] = True

    adata.obs["droplet"] = "other"
    adata.obs.loc[adata.obs["n_genes_by_counts"] <= constants.thresholds["upper_n_genes_empty"], "droplet"] = "empty"
    adata.obs.loc[
        (adata.obs["log1p_n_genes_by_counts"] > constants.thresholds["lower_log_n_genes"]) &
        (adata.obs["pct_counts_in_top_100_genes"] < constants.thresholds["upper_pct_top_100"]) &
        (adata.obs["total_counts"] < constants.thresholds["upper_total_counts"]) &
        (adata.obs["pct_counts_mt"] < constants.thresholds["upper_pct_mt"]) &
        (adata.obs["pct_counts_ribo"] < constants.thresholds["upper_pct_ribo"]) &
        (adata.obs["pct_counts_hb"] < constants.thresholds["upper_pct_hb"]),
        "droplet"
    ] = "cell"
    

    # cell counts
    counts["cell"] = adata.obs.loc[adata.obs["droplet"] == "cell", "sample"].value_counts()
    counts["empty"] = adata.obs.loc[adata.obs["droplet"] == "empty", "sample"].value_counts()
    counts["other"] = adata.obs.loc[adata.obs["droplet"] == "other", "sample"].value_counts()

    counts["cell|condition"] = adata.obs.loc[(adata.obs["droplet"] == "cell") & (adata.obs["conditional"] == True), "sample"].value_counts()
    counts["empty|condition"] = adata.obs.loc[(adata.obs["droplet"] == "empty") & (adata.obs["conditional"] == True), "sample"].value_counts()
    counts["other|condition"] = adata.obs.loc[(adata.obs["droplet"] == "other") & (adata.obs["conditional"] == True), "sample"].value_counts()

    pd.DataFrame(counts).to_csv(f"data/counts_{name}.csv", index=True)

    # normalization and transformation
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    adata.write(target_path)

process("402", "contamination/data/402-0.h5ad")
process("403", "contamination/data/403-0.h5ad")
process("404", "contamination/data/404-0.h5ad")
process("405", "contamination/data/405-0.h5ad")
