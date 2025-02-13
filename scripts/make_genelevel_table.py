import anndata as ad
import pandas as pd
import click

import scgenome


@click.command()
@click.argument('input_bin_cn_h5ad')
@click.argument('genes_gtf')
@click.argument('output_gene_cn_csv')
def make_genelevel_table(input_bin_cn_h5ad, genes_gtf, output_gene_cn_csv):
    adata = ad.read_h5ad(input_bin_cn_h5ad)

    genes = scgenome.tl.read_ensemble_genes_gtf(genes_gtf)

    cols = [
        'hscr.a1', 'hscr.a2',
        'major_cn', 'minor_cn',
        'modal.a1', 'modal.a2',
    ]

    adata_genes = scgenome.tl.aggregate_genes(adata, genes, agg_layers=cols)

    data = []
    for col in cols:
        data.append(adata_genes.to_df(layer=col).melt(ignore_index=False, value_name=col).set_index('gene_id', append=True))
    data = pd.concat(data, axis=1)

    data['gene_name'] = adata_genes.var.loc[data.index.get_level_values(1), 'gene_name'].values
    data['ploidy'] = adata_genes.obs.loc[data.index.get_level_values(0), 'ploidy'].values
    data['n_wgd'] = adata_genes.obs.loc[data.index.get_level_values(0), 'n_wgd'].values

    total_cn = (data['minor_cn'] + data['major_cn']).round()

    # Default Neutral state
    data['gistic_value'] = 0

    # Shallow deletion as anything less than 2**(n_wgd + 1)
    data.loc[total_cn < 2**(data['n_wgd'] + 1), 'gistic_value'] = -1

    # Deep deletion as homozygous deletion, overrides shallow deletion
    data.loc[total_cn == 0, 'gistic_value'] = -2

    # Gain as anything greater than 2**(n_wgd + 1)
    data.loc[total_cn > 2**(data['n_wgd'] + 1), 'gistic_value'] = 1

    # Amplification as anything greater than 2**(n_wgd + 2), overrides gain
    data.loc[total_cn > 2**(data['n_wgd'] + 2), 'gistic_value'] = 2

    # Make pivot table
    data = data[['gistic_value']].reset_index().pivot(index='gene_id', columns='sample_id', values='gistic_value')
    assert data.isna().sum().sum() == 0
    
    # Save
    data.to_csv(output_gene_cn_csv, index=True, header=True, sep='\t')

if __name__ == "__main__":
    make_genelevel_table()
