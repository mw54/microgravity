import scanpy as sc
import anndata as ad
import pandas as pd
import constants

sc.settings.n_jobs = 16

def process(name, target_path):
    adatas = dict()
    cell_counts = dict()

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
    cell_counts["raw"] = adata.obs["sample"].value_counts()

    # filter cells
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[
        (adata.obs["log1p_n_genes_by_counts"] > constants.thresholds["lower_log_n_genes"]) &
        (adata.obs["total_counts"] < constants.thresholds["upper_total_counts"]) &
        (adata.obs["pct_counts_in_top_100_genes"] < constants.thresholds["upper_pct_top_100"]) &
        (adata.obs["pct_counts_mt"] < constants.thresholds["upper_pct_mt"]) &
        (adata.obs["pct_counts_ribo"] < constants.thresholds["upper_pct_ribo"]) &
        (adata.obs["pct_counts_hb"] < constants.thresholds["upper_pct_hb"])
        ,:
    ].copy()
    sc.pp.scrublet(adata, batch_key="sample")
    adata = adata[~adata.obs['predicted_doublet'], :].copy()
    cell_counts["filtered"] = adata.obs["sample"].value_counts()

    # save cell counts
    df = pd.DataFrame(cell_counts)
    df.to_csv(f"data/counts_{name}.csv", index=True)

    # normalization and transformation
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    # feature selection
    sc.pp.highly_variable_genes(adata, batch_key="sample")
    sc.pl.highly_variable_genes(adata, save=f"_{name}.png")

    # dimension reduction
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=64)
    sc.pl.pca_variance_ratio(adata, n_pcs=64, log=True, save=f"_{name}.png")
    sc.pl.pca(
        adata,
        color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        save=f"_{name}.png"
    )

    # batch effect check
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="sample", size=2, save=f"_samples_{name}.png")
    sc.pl.umap(adata, color="group", size=2, save=f"_groups_{name}.png")
    sc.pl.umap(adata, color="gravity", size=2, save=f"_gravity_{name}.png")
    sc.pl.umap(adata, color="age", size=2, save=f"_age_{name}.png")

    # clustering
    sc.tl.leiden(adata, resolution=1.0, flavor="igraph", random_state=42)
    sc.pl.umap(adata, color=["leiden"], legend_loc='on data', save=f"_cluster_{name}.png")

    # reassessment
    sc.pl.umap(
        adata,
        color=[
            'doublet_score', 
            'n_genes_by_counts',
            'log1p_n_genes_by_counts',
            'total_counts',
            'log1p_total_counts',
            'pct_counts_in_top_50_genes',
            'pct_counts_in_top_100_genes',
            'pct_counts_in_top_200_genes',
            'pct_counts_in_top_500_genes',
            'total_counts_mt',
            'log1p_total_counts_mt',
            'pct_counts_mt',
            'total_counts_ribo',
            'log1p_total_counts_ribo',
            'pct_counts_ribo',
            'total_counts_hb',
            'log1p_total_counts_hb',
            'pct_counts_hb'
        ],
        wspace=0.1,
        hspace=0.2,
        size=3,
        ncols=3,
        save=f"_reassessment_{name}.png"
    )

    # differentially expressed genes
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", layer="lognorm", use_raw=False)
    sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=10, save=f"top10_{name}.png")

    # save data
    adata.write(target_path)

process("402", "../data/402-0.h5ad")
process("403", "../data/403-0.h5ad")
process("404", "../data/404-0.h5ad")
process("405", "../data/405-0.h5ad")
