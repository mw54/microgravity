import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import multiprocessing as mp
import constants

def plot(args):
    name, sample, path = args

    adata = sc.read_10x_mtx(path)
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    sc.pl.violin(
        adata,
        [
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
        jitter=0.4,
        multi_panel=True,
        show=False
    )
    plt.savefig(f"figures/violin_{name}_{sample}.png", dpi=300, bbox_inches='tight')
    plt.close()
    return sample, adata

def process(name):
    adatas = dict()
    args = [(name, sample, path) for sample, path in constants.samples[name].items()]
    with mp.Pool(processes=4) as pool:
        for sample, adata in pool.imap_unordered(plot, args):
            adatas[sample] = adata

    adata = ad.concat(adatas, label="sample")
    sc.pl.violin(
        adata,
        [
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
        stripplot=False,
        multi_panel=True,
        show=False
    )
    plt.savefig(f"figures/violin_{name}_overall.png", dpi=300, bbox_inches='tight')
    plt.close()

process("402")
process("403")
process("404")
process("405")
