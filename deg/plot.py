import scanpy as sc
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib_venn import venn2
from adjustText import adjust_text

organ_map = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def plot_volcano(title, data, path, top_labels=10):
    plt.figure(figsize=(5, 5), dpi=300)
    ax = sns.scatterplot(
        x='logfoldchanges', y='nlog10padj', data=data,
        hue='direction',
        palette={'up': 'red', 'down': 'blue', 'ns': 'grey'},
        alpha=0.6, s=20,
        hue_order=['up', 'down', 'ns']
    )
    plt.title(title)
    plt.xlabel(r'$\log_2$(Fold Change) uG/1G')
    plt.ylabel(r'$-\log_{10}$(Adjusted p-Value)')

    # add text label for top significant DEGs
    data["|scores|"] = data["scores"].abs()
    labels = data.loc[data["direction"] != "ns"].groupby("direction").apply(lambda x: x.sort_values("|scores|", ascending=False).head(top_labels))
    texts = []
    for _, row in labels.iterrows():
        texts.append(ax.text(row['logfoldchanges'], row['nlog10padj'], row["gene"], fontsize=9))
    adjust_text(
        texts,
        arrowprops=dict(arrowstyle='-', color='gray', lw=0.5)
    )
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot_dotplot(title, data, path):
    plt.figure(figsize=(7, 6), dpi=300)
    data = data.copy()
    data['logfoldchanges'] = data['logfoldchanges'].clip(-3.5, 3.5)
    data['nlog10padj'] = data['nlog10padj'].clip(0, 5)
    ax = sns.scatterplot(
        x='celltype', y='gene', data=data,
        hue='logfoldchanges',
        palette='RdBu_r',
        hue_norm=(-3.5, 3.5),
        size_norm=(0, 5),
        size='nlog10padj',
        sizes=(0, 100),
        legend=True
    )
    sns.move_legend(ax, "center left", bbox_to_anchor=(1.02, 0.5))
    legend = ax.get_legend()
    legend.texts[0].set_text(r"$\log_2$(FC) uG/1G")
    legend.texts[6].set_text(r"$-\log_{10}$(p-adj)")
    plt.title(title)
    plt.xticks(rotation=90, ha="center")
    plt.xlabel('Cell Type')
    plt.ylabel(f'Top DEGs in {title}')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot(name, source_path):
    adata = sc.read(source_path)

    # get metadata
    celltypes = adata.uns["DEG_indices"]["celltype"]
    ages = adata.uns["DEG_indices"]["age"]
    metrics = adata.uns["DEG_indices"]["metric"]
    dfs = adata.varm["DEG"]

    # get DEG data
    degs = list()
    for celltype, celltype_index in celltypes.items():
        for age, age_index in ages.items():
            deg = pd.DataFrame(dfs[:, celltype_index, age_index, list(metrics.values())], columns=list(metrics.keys()))
            deg["gene"] = adata.var_names
            deg["celltype"] = celltype
            deg["age"] = age
            degs.append(deg)
    
    # create dataframe
    degs = pd.concat(degs, axis=0)

    # determine significance
    degs["direction"] = "ns"
    degs.loc[(degs["pvals_adj"] < 0.05) & (degs["logfoldchanges"] > 0.5), "direction"] = "up"
    degs.loc[(degs["pvals_adj"] < 0.05) & (degs["logfoldchanges"] < -0.5), "direction"] = "down"
    
    # calculate -log10(padj)
    min_padj = degs.loc[degs["pvals_adj"] > 0, "pvals_adj"].min()
    degs["nlog10padj"] = degs["pvals_adj"].apply(lambda x: -np.log10(x) if x > 0 else -np.log10(min_padj))
    degs.to_csv(f"data/volcano_{name}.csv", index=False)

    # select significant DEGs
    sigs = degs.groupby("celltype").apply(lambda x: x.groupby('age').apply(lambda x:x.loc[x["direction"] != "ns", "gene"].tolist()).to_dict()).to_dict()
    json.dump(sigs, open(f"data/venn_{name}.json", "w"))

    # tally significant genes
    tally = degs.groupby("age").apply(lambda x: x.loc[x["direction"] != "ns", "gene"].value_counts()).reset_index()
    tally.to_csv(f"data/tally_{name}.csv", index=False)

    # create preview volcano plot
    for celltype in degs["celltype"].unique():
        for age in degs["age"].unique():
            plot_volcano(f"{organ_map[name]}, {age.capitalize()}, {celltype}", degs.loc[(degs["celltype"] == celltype) & (degs["age"] == age)], f"figures/volcano_{celltype.lower().replace(" ", "_")}_{age.lower().replace(" ", "_")}_{name}.png")

    # create dotplot for top tallied DEGs
    for age in degs["age"].unique():
        plot_dotplot(f"{organ_map[name]}, {age.capitalize()}", degs.loc[(degs["age"] == age) & (degs["gene"].isin(tally.loc[tally["age"] == age, "gene"].head(20).tolist()))], f"figures/dotplot_{age.lower().replace(" ", "_")}_{name}.png")

plot("402", "../data/402-2.h5ad")
plot("403", "../data/403-2.h5ad")
plot("404", "../data/404-2.h5ad")
plot("405", "../data/405-2.h5ad")


    
