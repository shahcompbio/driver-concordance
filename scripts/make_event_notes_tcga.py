import pickle

import click
import anndata as ad
import pandas as pd


def get_cohort_table(gdata):
    obs = gdata.obs.copy().reset_index()
    obs_ixs = 'sample_id,case,cancer_type'.split(',')
    obs = obs[obs_ixs].rename(columns={'case':'case_id'})
    return obs

def make_stats_summary(events):
    tts = events['cancer_type'].unique()
    data = []
    for tt in tts:
        tt = tt.strip()
        tt_events = events[events['cancer_type']==tt].copy()
        for gene, eventdf in tt_events.groupby('gene'):
            cnt = eventdf.shape[0]
            field = [gene, cnt, tt]
            data.append(field)
        # tt_stats[tt].set_index('gene', inplace=True)
        # tt_stats[tt].sort_values('count', ascending=False, inplace=True)
    df = pd.DataFrame(data, columns=['gene', 'count', 'cancer_type'])
    df.sort_values(['cancer_type', 'count'], ascending=[True, False], inplace=True)
    return df

@click.command()
@click.option('--cn_anndata', type=str, help="input gene-level cnvs anndata")
@click.option('--mutation_anndata', type=str, help="input annotated anndata for mutations")
@click.option('--events_cn', type=str, help="input events table for cnvs")
@click.option('--events_mutation', type=str, help="events table for mutations")
@click.option('--cohort', type=str, help="output cohort samples table")
@click.option('--summary', type=str, help="output tumor-model comparison summary")
@click.option('--counts', type=str, help="output tumor type counts from available data")
@click.option('--genes_to_exclude', type=str, nargs=2, help="genes to exclude")
def make_samples_table(cn_anndata, mutation_anndata, events_cn, events_mutation,
                       cohort, summary, counts, genes_to_exclude):
    gdata_cn = ad.read_h5ad(cn_anndata)
    gdata_mutation = ad.read_h5ad(mutation_anndata)
    shared_samples = sorted(set(gdata_cn.obs.index.tolist()) & set(gdata_mutation.obs.index.tolist()))
    gdata_cn = gdata_cn[shared_samples]
    gdata_mutation = gdata_mutation[shared_samples]
    obs = get_cohort_table(gdata_cn)
    obs.to_csv(cohort, index=False)

    cnts = obs['cancer_type'].value_counts().reset_index()
    cnts.columns = ['cancer_type', 'count']
    cnts.set_index('cancer_type', inplace=True)
    cnts.to_csv(counts, index=True)

    event_cols = ['sample_id', 'gene', 'event', 'cancer_type']
    events_cn = pd.read_csv(events_cn)[event_cols]
    events_mutation = pd.read_csv(events_mutation)[event_cols]
    events = pd.concat([events_cn, events_mutation]).reset_index(drop=True).drop_duplicates()
    events = events[events['sample_id'].isin(shared_samples)]
    pivot_summary = make_stats_summary(events)
    pivot_summary.to_csv(summary, index=False)

if __name__ == "__main__":
    make_samples_table()
