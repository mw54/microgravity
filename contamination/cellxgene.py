import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def dotplot_droplet(title, data, path):
    organs = data["organ"].unique().tolist()
    genes = data["gene"].unique().tolist()
    fig, axes = plt.subplots(
        ncols=len(organs),
        nrows=1,
        sharey=True,  # All subplots share the Y-axis scale
        figsize=(16, 5),
        dpi=300,
        gridspec_kw={'width_ratios': [len(data.loc[data["organ"] == organ, "celltype"].unique().tolist()) for organ in organs]}
    )
    for i, organ in enumerate(organs):
        subset = data[data["organ"] == organ].copy()
        subset["mean"] = subset["mean"].clip(0, 8)
        ax = axes[i]
        sns.scatterplot(
            data=subset,
            x='celltype',
            y='gene',
            hue='mean',
            palette='Blues',
            hue_norm=(0, 8),
            size_norm=(0, 1),
            size='pct',
            sizes=(0, 100),
            ax=ax,
            legend=False if i < len(organs) - 1 else True
        )
        ax.set_title(organ, size=8)
        ax.set_xlim(-0.5, len(data.loc[data["organ"] == organ, "celltype"].unique().tolist()) - 0.5)
        ax.set_ylim(-0.5, len(genes) - 0.5)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', rotation=90)
        if i > 0:
            ax.tick_params(axis='y', left=False)
        if i == len(organs) - 1:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1.0))

    fig.suptitle(title)
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.02)
    plt.savefig(path)
    plt.close()

def plot(name, source_path):
    data = pd.read_csv(source_path)
    data["pct"] = data.apply(lambda row: row["Number of Cells Expressing Genes"] / row["Cell Count"], axis=1)
    df = data[["Tissue", "Cell Type", "Gene Symbol", "Expression", "pct"]]
    df = df.rename(columns={"Tissue": "organ", "Cell Type": "celltype", "Gene Symbol": "gene", "Expression": "mean", "pct": "pct"})
    dotplot_droplet(name, df, f"dotplot_{name.lower()}.png")

plot("CELLxGENE", "data/cellxgene.csv")
