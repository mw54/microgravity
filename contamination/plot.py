import scanpy as sc
import pandas as pd
import constants
import matplotlib.pyplot as plt
import seaborn as sns

organ_map = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def dotplot_droplet(title, data, path):
    droplets = data["droplet"].unique().tolist()
    genes = data["gene"].unique().tolist()
    fig, axes = plt.subplots(
        ncols=len(droplets),
        nrows=1,
        sharey=True,  # All subplots share the Y-axis scale
        figsize=(5, 5),
        dpi=300
    )
    fig.subplots_adjust(wspace=0.1)
    for i, droplet in enumerate(droplets):
        subset = data[data["droplet"] == droplet].copy()
        subset["mean"] = subset["mean"].clip(0, 8)
        ax = axes[i]
        sns.scatterplot(
            data=subset,
            x='gene',
            y='sample',
            hue='mean',
            palette='Blues',
            hue_norm=(0, 8),
            size_norm=(0, 1),
            size='pct',
            sizes=(0, 100),
            ax=ax,
            legend=False if i < len(droplets) - 1 else True
        )
        ax.set_title(droplet, size=8)
        ax.set_xlim(-0.5, len(genes) - 0.5)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', rotation=90)
        if i > 0:
            ax.tick_params(axis='y', left=False)
        if i == len(droplets) - 1:
            sns.move_legend(ax, "center left", bbox_to_anchor=(1.02, 0.5))

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def dotplot_sample(title, data, path):
    samples = data["sample"].unique().tolist()
    genes = data["gene"].unique().tolist()
    fig, axes = plt.subplots(
        ncols=len(samples),
        nrows=1,
        sharex=True,
        sharey=True,
        figsize=(20, 5),
        dpi=300
    )
    for i, sample in enumerate(samples):
        subset = data[data["sample"] == sample].copy()
        subset["mean"] = subset["mean"].clip(0, 8)
        ax = axes[i]
        sns.scatterplot(
            data=subset,
            x='gene',
            y='clustertype',
            hue='mean',
            palette='Blues',
            hue_norm=(0, 8),
            size_norm=(0, 1),
            size='pct',
            sizes=(0, 100),
            ax=ax,
            legend=False if i < len(samples) - 1 else True
        )
        ax.set_title(sample, size=8)
        ax.set_xlim(-0.5, len(genes) - 0.5)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', rotation=90)
        if i > 0:
            ax.tick_params(axis='y', left=False)
        if i == len(samples) - 1:
            sns.move_legend(ax, "center left", bbox_to_anchor=(1.02, 0.5))

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def plot_droplet(name, source_path):
    adata = sc.read(source_path)

    # get expression dataframe
    expr = adata[:,constants.genes].to_df(layer="lognorm")
    expr[["sample", "droplet", "conditional"]] = adata.obs[["sample", "droplet", "conditional"]]

    # unconditional expression
    mean_expr = expr.groupby(["sample", "droplet"]).apply(lambda x: x.mean(axis=0)).melt(None, constants.genes, "gene", "mean", ignore_index=False).set_index("gene", append=True)
    pct_expr = expr.groupby(["sample", "droplet"]).apply(lambda x: (x > 0).mean(axis=0)).melt(None, constants.genes, "gene", "pct", ignore_index=False).set_index("gene", append=True)
    df = pd.concat([mean_expr, pct_expr], axis=1).reset_index()
    df.to_csv(f"data/unconditional_{name}.csv", index=False)
    dotplot_droplet(f"Unfiltered, {organ_map[name]}", df, f"figures/dotplot_unfiltered_{name}.png")

    # conditional expression
    mean_expr = expr.loc[expr["conditional"]].groupby(["sample", "droplet"]).apply(lambda x: x.mean(axis=0)).melt(None, constants.genes, "gene", "mean", ignore_index=False).set_index("gene", append=True)
    pct_expr = expr.loc[expr["conditional"]].groupby(["sample", "droplet"]).apply(lambda x: (x > 0).mean(axis=0)).melt(None, constants.genes, "gene", "pct", ignore_index=False).set_index("gene", append=True)
    df = pd.concat([mean_expr, pct_expr], axis=1).reset_index()
    df.to_csv(f"data/conditional_{name}.csv", index=False)
    dotplot_droplet(f"Filtered, {organ_map[name]}", df, f"figures/dotplot_filtered_{name}.png")

def plot_sample(name, source_path):
    adata = sc.read(source_path)

    expr = adata[:,constants.genes].to_df(layer="lognorm")
    expr[["sample", "clustertype"]] = adata.obs[["sample", "clustertype"]]
    mean_expr = expr.groupby(["sample", "clustertype"]).apply(lambda x: x.mean(axis=0)).melt(None, constants.genes, "gene", "mean", ignore_index=False).set_index("gene", append=True)
    pct_expr = expr.groupby(["sample", "clustertype"]).apply(lambda x: (x > 0).mean(axis=0)).melt(None, constants.genes, "gene", "pct", ignore_index=False).set_index("gene", append=True)
    df = pd.concat([mean_expr, pct_expr], axis=1).reset_index()
    df.to_csv(f"data/sample_{name}.csv", index=False)

    dotplot_sample(f"Filtered, {organ_map[name]}", df, f"figures/dotplot_celltype_{name}.png")
    

plot_droplet("402", "data/402-0.h5ad")
plot_droplet("403", "data/403-0.h5ad")
plot_droplet("404", "data/404-0.h5ad")
plot_droplet("405", "data/405-0.h5ad")

plot_sample("402", "data/402-2.h5ad")
plot_sample("403", "data/403-2.h5ad")
plot_sample("404", "data/404-2.h5ad")
plot_sample("405", "data/405-2.h5ad")