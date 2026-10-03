import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

organ_map = {
    "402": "Femur",
    "403": "Humerus",
    "404": "Blood",
    "405": "Spleen"
}

def plot_scatter(title, data, var, path):
    data = data.copy()
    data["|scores|"] = data["scores"].abs()
    data = data.sort_values("|scores|", ascending=True)
    data['nlog10padj'] = data['nlog10padj'].clip(0, 5)
    mask = (data['pvals_adj'] < 0.05)
    slope, intercept = np.polyfit(data.loc[mask, 'intrinsic'], data.loc[mask, var], 1)
    r2 = np.corrcoef(data.loc[mask, 'intrinsic'], data.loc[mask, var])[0, 1]**2

    plt.figure(figsize=(5, 5), dpi=300)
    sns.scatterplot(x='intrinsic', y=var, data=data, hue='nlog10padj', palette='Blues', s=20, legend=True)
    plt.axline((0, intercept), (-intercept / slope, 0), color="black", linestyle="--", label=f"$R^2={r2:.2f}$")
    plt.legend(title=r"$-\log_{10}$(p-adj)")

    plt.title(title)
    plt.xlim(-6, 6)
    plt.xlabel(r'$\log_2$(Fold Change) Old/Young')
    plt.ylabel(r'$-\log_2$(Fold Change) uG/1G')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot_volcano(title, data, path, top_labels=10):
    plt.figure(figsize=(5, 5), dpi=300)
    ax = sns.scatterplot(
        x='intrinsic', y='nlog10padj', data=data, hue='direction',
        palette={'up': 'red', 'down': 'blue', 'ns': 'grey'},
        alpha=0.6, s=20,
        hue_order=['up', 'down', 'ns']
    )
    plt.xlim(-7.5, 7.5)
    plt.title(title)
    plt.xlabel(r'$\log_2$(Fold Change) Old/Young')
    plt.ylabel(r'$-\log_{10}$(Adjusted p-Value)')

    # add text label for top significant DEGs
    data["|scores|"] = data["scores"].abs()
    labels = data.loc[data["direction"] != "ns"].groupby("direction").apply(lambda x: x.sort_values("|scores|", ascending=False).head(top_labels))
    texts = [ax.text(row['intrinsic'], row['nlog10padj'], row["names"], fontsize=9) for _, row in labels.iterrows()]
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def plot(name, source_path):
    ageing = pd.read_csv(source_path)
    min_padj = ageing.loc[ageing["pvals_adj"] > 0, "pvals_adj"].min()
    ageing["nlog10padj"] = ageing["pvals_adj"].apply(lambda x: -np.log10(x) if x > 0 else -np.log10(min_padj))
    ageing["direction"] = "ns"
    ageing.loc[(ageing["pvals_adj"] < 0.05) & (ageing["intrinsic"] > 0.5), "direction"] = "up"
    ageing.loc[(ageing["pvals_adj"] < 0.05) & (ageing["intrinsic"] < -0.5), "direction"] = "down"

    plot_volcano(f"{organ_map[name]}, Ageing", ageing, f"figures/volcano_{name}.png")
    for age in ["old", "young"]:
        plot_scatter(f"{organ_map[name]}, {age.capitalize()}", ageing, age, f"figures/scatter_{age}_{name}.png")


plot("402", "data/ageing_402.csv")
plot("403", "data/ageing_403.csv")
plot("404", "data/ageing_404.csv")
plot("405", "data/ageing_405.csv")
    
