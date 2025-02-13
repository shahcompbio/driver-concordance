import pickle

import click
import anndata as ad
import pandas as pd

from gdan_hcmi_tools.process import VEP_meta


def make_notes_from_events(events):
    events = events.copy()
    event_cols = ['sample', 'sample_type', 'gene', 'event', 'case', 'cancer_type']
    for event_col in event_cols:
        assert event_col in events.columns, f'event_col {event_col} not in \nevents.columns: {events.columns}'
    events['sample_type'] = events['sample_type'].replace({'model_expanded':'model'})
    events = events[event_cols].reset_index(drop=True).drop_duplicates()
    events['event_rank'] = events['event'].map(VEP_meta.event_rank)
    events = events.sort_values(['sample', 'gene', 'event_rank'], ascending=False).reset_index(drop=True)
    for (case_id, gene, event), edf in events.groupby(['case', 'gene', 'event']):
        sample_types = set(edf['sample_type'])
        note = 'NA'
        if sample_types == {'tumor', 'model'} or sample_types == {'tumor', 'model', 'model_expanded'}:
            note = 'both'
        elif sample_types == {'tumor'} or sample_types == {'tumor', 'model_expanded'}:
            note = 'tumor-only'
        elif sample_types == {'model'} or sample_types == {'model', 'model_expanded'}:
            note = 'model-only'
        elif sample_types == {'model_expanded'}:
            continue
        else:
            raise ValueError(sample_types)
        events.loc[edf.index, 'note'] = note
    case_events = events[['case', 'cancer_type', 'gene', 'event', 'note']].drop_duplicates()
    ixs = []
    for (case_id, gene), gdf in case_events.groupby(['case', 'gene']):
        ixs.append(gdf.index[0])
    case_events = case_events.loc[ixs].reset_index(drop=True)
    notes = case_events[['cancer_type', 'gene', 'event', 'note']].value_counts()
    notes = notes.reset_index().rename(columns={0:'count'})
    return notes

def make_summary_from_notes(notes):
    notes = notes.sort_values('count', ascending=False).reset_index(drop=True)
    summary = (notes[['cancer_type', 'gene', 'note', 'count']]
        .groupby(['cancer_type', 'gene', 'note'])
        .agg({'count':'sum'})
        .sort_values('count', ascending=False)
        .reset_index())
    return summary

def get_summary_per_tumor_type(notes, count_both_as=1):
    df = pd.DataFrame()
    for cancer_type, cnotes in notes.groupby('cancer_type'):
        stats = cnotes.pivot(index=['gene'], columns=['note']).fillna(0)['count']
        if 'both' in stats.columns:
            stats['both'] *= count_both_as
        else:
            stats['both'] = 0
        for note_type in ['model-only', 'tumor-only']:
            if note_type not in stats.columns:
                stats[note_type] = 0
        stats['total'] = stats.sum(axis=1)
        stats.reset_index(inplace=True)
        stats.sort_values('total', ascending=False, inplace=True)
        stats['cancer_type'] = cancer_type
        df = pd.concat([df, stats])
    return df

def get_tumor_type_counts(adata):
    tts = sorted(adata.obs['cancer_type'].unique())
    data = []
    for _, tt in enumerate(tts):
        if tt == 'Others': continue
        tt_subset = adata.obs[adata.obs['cancer_type'] == tt]
        tt_count = tt_subset.index.str.slice(0, 17).unique().shape[0]
        data.append([tt, tt_count])
    df = pd.DataFrame(data, columns=['cancer_type', 'count'])
    df.sort_values('cancer_type', inplace=True)
    df.sort_values('count', ascending=False, inplace=True)
    return df


@click.command()
@click.option('--cohort', type=str, help="input cohort pickle path")
@click.option('--mutation_anndata', type=str, help="output annotated anndata for mutations")
@click.option('--events_path', type=str, help="input mutation events table")
@click.option('--notes_path', type=str, help="output tumor-model event comparison notes table")
@click.option('--cohort_table', type=str, help="output cohort obs data")
@click.option('--summary', type=str, help="output tumor-model comparison summary")
@click.option('--counts', type=str, help="output tumor type counts from available data")
@click.option('--genes_to_exclude', type=str, nargs=2, help="genes to exclude")
def make_event_notes(cohort, mutation_anndata, events_path,  
        notes_path, cohort_table, summary, counts, genes_to_exclude):
    with open(cohort, 'rb') as handle:
        cohort = pickle.load(handle)
    adata = ad.read_h5ad(mutation_anndata)
    adata = adata[adata.obs.index.isin(cohort.sample_ids)]
    adata.obs.to_csv(cohort_table, index=True)

    # Concordance
    event_cols = ['sample', 'sample_type', 'gene', 'event', 'case', 'cancer_type']
    events = pd.read_csv(events_path)
    events = events[events['sample'].isin(cohort.sample_ids)]
    # events['sample_type'] = events['sample_type'].replace({'model_expanded':'model'})
    events = events[event_cols].reset_index(drop=True).drop_duplicates()
    events = events[~events['gene'].isin(genes_to_exclude)]

    notes = make_notes_from_events(events)
    notes.to_csv(notes_path, index=False)
    _summary = make_summary_from_notes(notes)
    pivot_summary = get_summary_per_tumor_type(_summary)
    tumor_type_counts = get_tumor_type_counts(adata)
    pivot_summary.to_csv(summary, index=False)
    tumor_type_counts.to_csv(counts, index=False)


if __name__ == "__main__":
    make_event_notes()