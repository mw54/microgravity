import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

terms = {
    "Mitochondria": ["AEROBIC RESPIRATION AND RESPIRATORY ELECTRON TRANSPORT", "ELECTRON TRANSPORT CHAIN", "OXIDATIVE PHOSPHORYLATION", "MITOCHONDRIAL TRANSLATION"],
    "Proliferation": ["G1 S TRANSITION", "APC C MEDIATED DEGRADATION OF CELL CYCLE PROTEINS", "MITOTIC G2 G2 M PHASES", "DNA REPLICATION", "CDK MEDIATED PHOSPHORYLATION AND REMOVAL OF CDC6"],
    "Motility": ["RHO GTPASE CYCLE", "RAC1 GTPASE CYCLE", "CDC42 GTPASE CYCLE", "RHOA GTPASE CYCLE", "CELL EXTRACELLULAR MATRIX INTERACTIONS"],
    "Immunity": ["SIGNALING BY THE B CELL RECEPTOR BCR", "DOWNSTREAM TCR SIGNALING", "NEUTROPHIL DEGRANULATION", "TOLLLIKE RECEPTOR SIGNALING", "FCGAMMA RECEPTOR FCGR DEPENDENT PHAGOCYTOSIS"]
}

organs = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def plot_dotplot(title, data, path):
    fig, axes = plt.subplots(
        ncols=5,
        nrows=5,
        figsize=(20, 7),
        dpi=300,
        gridspec_kw={'height_ratios': [1, 6, 6, 6, 6], 'width_ratios': [6, 6, 6, 6, 1]}
    )
    for i, aspect in enumerate(terms):
        for j, organ in enumerate(organs):
            subset = data.loc[(data['Term'].isin(terms[aspect])) & (data["Organ"] == organ)].copy()
            term_cate = pd.CategoricalDtype(categories=sorted(subset["Term"].unique()), ordered=True)
            cell_cate = pd.CategoricalDtype(categories=sorted(subset["Celltype"].unique()), ordered=True)
            subset["Term"] = subset["Term"].astype(term_cate)
            subset["Celltype"] = subset["Celltype"].astype(cell_cate)
            ax = axes[i+1, j]
            sns.scatterplot(
                data=subset,
                x='Celltype',
                y='Term',
                hue='NES',
                palette='RdBu_r',
                hue_norm=(-3.5, 3.5),
                size_norm=(0, 5),
                size='nlog10padj',
                sizes=(0, 100),
                ax=ax,
                legend=False
            )
            ax.set_ylim(-0.5, len(terms[aspect]) - 0.5)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='x',  rotation=90)
            if i+1 < 4:
                ax.set_xticks([])
                ax.tick_params(axis='x', bottom=False)
            if j > 0:
                ax.set_yticks([])
                ax.tick_params(axis='y', left=False)

    for i, aspect in enumerate(terms):
        ax = axes[i+1, 4]
        ax.axis('off')
        ax.text(
            0.0, 0.5, # Position the text at the center of the axes (0.5, 0.5)
            aspect,
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

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def plot(names, source_paths):
    dfs = list()
    for name, path in zip(names, source_paths):
        df = pd.read_csv(path)
        df["Organ"] = name
        dfs.append(df)
    dfs = pd.concat(dfs, axis=0, ignore_index=True)
    min_padj = dfs.loc[dfs["FDR q-val"] > 0, "FDR q-val"].min()
    dfs["nlog10padj"] = dfs["FDR q-val"].apply(lambda x: -np.log10(x) if x > 0 else -np.log10(min_padj))
    dfs = dfs.loc[dfs["Term"].isin([t for aspect in terms.values() for t in aspect])]
    for age in dfs[ "Age"].unique():
        plot_dotplot(age.capitalize(), dfs.loc[dfs["Age"] == age], f"figures/summary_{age}.png")

plot(["402", "403", "404", "405"], ["data/gsea_402.csv", "data/gsea_403.csv", "data/gsea_404.csv", "data/gsea_405.csv"])