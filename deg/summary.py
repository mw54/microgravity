import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

genes = ["S100a8", "S100a9", "Hba-a2", "Hbb-bs"]
ages = ["old", "young"]
organs = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def plot_dotplot(data, path):
    fig, axes = plt.subplots(
        ncols=5,
        nrows=3,
        figsize=(20, 4),
        dpi=300,
        gridspec_kw={'height_ratios': [1, 6, 6], 'width_ratios': [6, 6, 6, 6, 1]}
    )
    for i, age in enumerate(ages):
        for j, organ in enumerate(organs):
            subset = data.loc[(data['age'] == age) & (data["organ"] == organ)].copy()
            gene_cate = pd.CategoricalDtype(categories=genes, ordered=True)
            cell_cate = pd.CategoricalDtype(categories=sorted(subset["celltype"].unique()), ordered=True)
            subset["gene"] = subset["gene"].astype(gene_cate)
            subset["celltype"] = subset["celltype"].astype(cell_cate)
            ax = axes[i+1, j]
            sns.scatterplot(
                data=subset,
                x='celltype',
                y='gene',
                hue='logfoldchanges',
                palette='RdBu_r',
                hue_norm=(-3.5, 3.5),
                size_norm=(0, 5),
                size='nlog10padj',
                sizes=(0, 100),
                ax=ax,
                legend=False
            )
            ax.set_ylim(-0.5, len(genes) - 0.5)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='x',  rotation=90)
            if i+1 < 2:
                ax.set_xticks([])
                ax.tick_params(axis='x', bottom=False)
            if j > 0:
                ax.set_yticks([])
                ax.tick_params(axis='y', left=False)

    for i, age in enumerate(ages):
        ax = axes[i+1, 4]
        ax.axis('off')
        ax.text(
            0.0, 0.5, # Position the text at the center of the axes (0.5, 0.5)
            age.capitalize(),
            transform=ax.transAxes,
            fontsize=12,
            ha='center',
            va='center',
            rotation=270
        )

    for j, organ in enumerate(organs):
        ax = axes[0, j]
        ax.axis('off')
        ax.text(
            0.5, 0.0,
            organs[organ],
            transform=ax.transAxes,
            fontsize=12,
            ha='center',
            va='center',
            rotation=0
        )

    axes[0, 4].axis('off')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot(names, source_paths):
    dfs = list()
    for name, path in zip(names, source_paths):
        df = pd.read_csv(path)
        df["organ"] = name
        dfs.append(df)
    dfs = pd.concat(dfs, axis=0, ignore_index=True)
    dfs = dfs.loc[dfs["gene"].isin(genes)]
    plot_dotplot(dfs, "figures/summary.png")

plot(["402", "403", "404", "405"], ["data/volcano_402.csv", "data/volcano_403.csv", "data/volcano_404.csv", "data/volcano_405.csv"])