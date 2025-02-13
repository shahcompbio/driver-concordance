import pandas as pd
import click

import scgenome


@click.command()
@click.argument('input_table')
@click.argument('input_events')
@click.argument('genes_gtf')
@click.argument('output_table')
def make_genelevel_table(input_table, input_events, genes_gtf, output_table):
    event_map = {
        'hlamp': 2,
        'homdel': -2,
    }

    genes = scgenome.tl.read_ensemble_genes_gtf(genes_gtf)
    genes_df = genes.as_df()
    # id2name = dict(zip(genes_df['gene_id'], genes_df['gene_name']))
    name2id = dict(zip(genes_df['gene_name'], genes_df['gene_id']))

    data = pd.read_table(input_table, index_col=0)
    assert data.isna().sum().sum() == 0

    events = pd.read_csv(input_events)
    for _, row in events.iterrows():
        sample_id = row['sample']
        gene_name = row['gene']
        gene_id = name2id[gene_name]
        event = row['event']
        value = event_map[event]
        data.loc[gene_id, sample_id] = value

    data.to_csv(output_table, index=True, header=True, sep='\t')

if __name__ == "__main__":
    make_genelevel_table()
