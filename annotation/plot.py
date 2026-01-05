import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import constants

organ_map = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def plot_stackedbar(title, data, path):
    counts = dict()
    for item in data["clustertype"].unique():
        counts[item] = data.loc[data["clustertype"] == item, "group"].value_counts()

    counts = pd.DataFrame(counts)
    counts.index.name = "clustertype"
    counts.columns.name = "group"

    counts = counts[sorted(data["clustertype"].unique().tolist())]
    counts = counts.apply(lambda row: row / sum(row), axis=1)

    ax = counts.plot(
        kind='bar',
        stacked=True,
        figsize=(5, 5),
        xlabel="group",
        ylabel='proportion'
    )
    ax.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()

def plot_dotplot(title, data, markers, path):
    nrows = 2
    ncols = int(len(markers) / nrows + 1 / nrows)
    fig, axes = plt.subplots(
        ncols=ncols,
        nrows=nrows,
        sharey=True,  # All subplots share the Y-axis scale
        figsize=(13, 9),
        dpi=300
    )
    fig.subplots_adjust(wspace=0.1)
    fig.subplots_adjust(hspace=0.3)

    data = data.copy()
    data.loc[data["logfoldchanges"] < 0, "logfoldchanges"] = np.nan
    data['logfoldchanges'] = data['logfoldchanges'].clip(-5, 5)
    data['nlog10padj'] = data['nlog10padj'].clip(0, 300)

    for i, clustertype in enumerate(markers):
        subset = data[data['gene'].isin(markers[clustertype])]
        ax = axes[int(i / ncols), int(i % ncols)]
        sns.scatterplot(
            data=subset,
            x='gene',
            y='clustertype',
            hue='logfoldchanges',
            palette='Reds',
            hue_norm=(-5, 5),
            size_norm=(0, 300),
            size='nlog10padj',
            sizes=(0, 100),
            ax=ax,
            legend=False if i + 1 != ncols else True
        )
        ax.set_title(clustertype, size=8)
        ax.set_xlim(-0.5, len(markers[clustertype]) - 0.5)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', rotation=90)
        if i % ncols != 0:
            ax.tick_params(axis='y', left=False)
        if i + 1 == ncols:
            sns.move_legend(ax, "center left", bbox_to_anchor=(1.02, 0.5))
            legend = ax.get_legend()
            for text in legend.texts:
                if text.get_text() == "logfoldchanges":
                    text.set_text(r"$\log_2$(FE) this/other")
                if text.get_text() == "nlog10padj":
                    text.set_text(r"$-\log_{10}$(p-adj)")

    # remove unused panels
    for i in range(len(markers), nrows * ncols):
        ax = axes[int(i / ncols), int(i % ncols)]
        ax.remove()

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot(name, source_path):
    adata = sc.read(source_path)

    # select relevant genes for expression metrics
    markers = list()
    for v in constants.markers[name].values():
        for vv in v:
            if vv not in markers:
                markers.append(vv)
    expr = adata[:,markers].to_df(layer="lognorm")
    expr["clustertype"] = expr.index.map(adata.obs["clustertype"])

    # compute expression mean and percentage
    mean_expr = expr.groupby("clustertype").apply(lambda x: x.mean(axis=0))
    pct_expr = expr.groupby("clustertype").apply(lambda x: (x > 0).mean(axis=0))

    # create annotation markers dotplot CSV
    dotplot = sc.get.rank_genes_groups_df(adata, None, key="celltype_markers")
    dotplot = dotplot[["group", "names", "logfoldchanges", "pvals_adj"]]
    dotplot.rename(columns={"group": "clustertype", "names": "gene"}, inplace=True)
    dotplot = dotplot.loc[dotplot["gene"].isin(markers),]
    min_padj = dotplot.loc[dotplot["pvals_adj"] > 0, "pvals_adj"].min()
    dotplot["nlog10padj"] = dotplot["pvals_adj"].apply(lambda x: -np.log10(x) if x > 0 else -np.log10(min_padj))
    dotplot["mean_expr"] = dotplot.apply(lambda row: mean_expr.loc[row["clustertype"], row["gene"]], axis=1)
    dotplot["pct_expr"] = dotplot.apply(lambda row: pct_expr.loc[row["clustertype"], row["gene"]], axis=1)
    dotplot["marker"] = dotplot.apply(lambda row: row["gene"] in constants.markers[name][row["clustertype"]] if row["clustertype"] in constants.markers[name] else False, axis=1)
    dotplot.to_csv(f"data/dotplot_{name}.csv", index=False)

    # create celltype UMAP CSV
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index)
    umap["leiden"] = adata.obs["leiden"]
    umap["clustertype"] = adata.obs["clustertype"]
    umap["sample"] = adata.obs["sample"]
    umap["group"] = adata.obs["group"]
    umap.to_csv(f"data/umap_{name}.csv", index=False)

    # celltype UMAP preview
    sc.pl.umap(adata, color="clustertype", size=2, title=organ_map[name], show=False)
    plt.savefig(f"figures/umap_celltype_{name}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # marker dotplot preview
    plot_dotplot(organ_map[name], dotplot, constants.markers[name], f"figures/dotplot_significance_{name}.png")

    # cell type proportion preview
    plot_stackedbar(organ_map[name], umap, f"figures/stackedbar_celltype_proportion_{name}.png")
    

plot("402", "../data/402-1.h5ad")
plot("403", "../data/403-1.h5ad")
plot("404", "../data/404-1.h5ad")
plot("405", "../data/405-1.h5ad")
