import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

organ_map = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def plot_dotplot(title, data, path):
    plt.figure(figsize=(11, 6), dpi=300)
    ax = sns.scatterplot(
        x='Celltype', y='Term', data=data,
        hue='NES',
        palette='RdBu_r',
        hue_norm=(-3.5, 3.5),
        size_norm=(0, 5),
        size='nlog10padj',
        sizes=(0, 100),
    )

    # customize legend labels
    legend = ax.get_legend()
    for text in legend.texts:
        if text.get_text() == "nlog10padj":
            text.set_text(r"$-\log_{10}$(p-adj)")
    sns.move_legend(ax, "center left", bbox_to_anchor=(1.02, 0.5))
    
    plt.title(title)
    plt.xticks(rotation=90, ha="center")
    plt.xlabel('Cell Type')
    plt.ylabel('Top Terms')

    original_labels = [label.get_text() for label in ax.get_yticklabels()]
    wrapped_labels = [(label[:57] + '...') if len(label) > 60 else label for label in original_labels]
    ax.set_yticklabels(wrapped_labels)

    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot_barplot(title, data, path):
    plt.figure(figsize=(10, 5), dpi=300)
    ax = sns.barplot(
        data=data,
        x='NES',
        y='Term',
        palette="Blues", # Color palette for the bars based on hue_norm
        hue="nlog10padj", # Use -log10(fdr) for coloring
        hue_norm=(0, 5),
        legend=True,
    )
    ax.legend_.set_title(r"$-\log_{10}$(p-adj)")

    # truncate overflowing y tick labels
    original_labels = [label.get_text() for label in ax.get_yticklabels()]
    wrapped_labels = [(label[:77] + '...') if len(label) > 80 else label for label in original_labels]
    ax.set_yticklabels(wrapped_labels)

    plt.title(title)
    plt.xlabel("Normalized Enrichment Score")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot(name, source_path):
    dfs = pd.read_csv(source_path)
    min_padj = dfs.loc[dfs["FDR q-val"] > 0, "FDR q-val"].min()
    dfs["nlog10padj"] = dfs["FDR q-val"].apply(lambda x: -np.log10(x) if x > 0 else -np.log10(min_padj))
    # dfs.to_csv(source_path)

    # top GSEA terms
    terms = dfs[dfs["Direction"] != "ns"].groupby("Age").apply(
        lambda x: x.groupby("Celltype").apply(
            lambda x: x.groupby("Direction").apply(
                lambda x: x.sort_values("|NES|", ascending=False)["Term"].head(30).to_list()
            ).to_dict()
        ).to_dict()
    ).to_dict()
    json.dump(terms, open(f"data/terms_{name}.json", "w"))

    # tally significant terms across cell types
    tally = dfs.groupby("Age").apply(lambda x: x.loc[x["Direction"] != "ns", "Term"].value_counts()).reset_index()
    tally.to_csv(f"data/tally_{name}.csv", index=False)

    # create barplot for individual celltype-age
    for celltype in dfs["Celltype"].unique():
        for age in dfs["Age"].unique():
            # take the corresponding subset, sort by NES and take the two extreme 15 terms
            subset = dfs[(dfs["Celltype"] == celltype) & (dfs["Age"] == age)].copy()
            subset = pd.concat([subset.sort_values("NES", ascending=True).head(15), subset.sort_values("NES", ascending=True).tail(15)], axis=0)
            subset = subset.loc[~subset.duplicated("Term")]
            plot_barplot(f"{organ_map[name]}, {age.capitalize()}, {celltype}", subset, f"figures/bar_{celltype.replace(" ", "_").lower()}_{age}_{name}.png")

    # create dotplot for top tallied terms
    for age in dfs["Age"].unique():
        plot_dotplot(f"{name} {age}", dfs.loc[(dfs["Age"] == age) & (dfs["Term"].isin(tally.loc[tally["Age"] == age, "Term"].head(20).tolist()))], f"figures/dotplot_{age.lower().replace(" ", "_")}_{name}.png")
    
plot("402", "data/gsea_402.csv")
plot("403", "data/gsea_403.csv")
plot("404", "data/gsea_404.csv")
plot("405", "data/gsea_405.csv")
