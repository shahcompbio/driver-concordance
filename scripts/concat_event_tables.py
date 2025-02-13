import click
import numpy as np
import pandas as pd

def normalize_variant_table(input_path, add_change=False):
    """
    Normalizes a variant table by ensuring required columns are present and adds missing optional columns.

    Parameters
    ----------
    input_path : str or Path
        Path to the CSV file containing variant data with at least the required columns.

    Returns
    -------
    DataFrame
        A DataFrame containing the normalized variant data with the following columns:
        - 'sample': Unique identifier for the sample.
        - 'case': Identifier for the case or patient.
        - 'sample_type': Type of sample (e.g., tumor, normal).
        - 'cancer_type': Type of cancer associated with the sample.
        - 'gene': Gene involved in the variant.
        - 'event': Type of event (e.g., mutation type).
        - 'variant': Formatted as 'gene:event', generated if not already present.
        - 'impact': Predicted impact of the variant (default NaN if not provided).
        - 'vaf': Variant allele frequency (default NaN if not provided).
    """
    df = pd.read_csv(input_path)
    if add_change:
        cols = ['sample', 'case', 'sample_type', 'cancer_type', 'gene', 'event', 'variant', 'change', 'impact', 'vaf']
    else:
        cols = ['sample', 'case', 'sample_type', 'cancer_type', 'gene', 'event', 'variant', 'impact', 'vaf']
    cols_req = ['sample', 'case', 'gene', 'event', 'cancer_type', 'sample_type']
    for col in cols_req:
        assert col in df.columns, f'{col} not in {df.columns}'

    if 'variant' not in df.columns:
        df['variant'] = df['gene'] + ':' + df['event']
    if add_change:
        if 'change' not in df.columns:
            df['change'] = df['variant'].replace(':', ' ')
    if 'impact' not in df.columns:
        df['impact'] = np.nan
    if 'vaf' not in df.columns:
        df['vaf'] = np.nan
    return df[cols]

@click.command()
@click.option('-i', '--input', nargs=3, type=click.Path(exists=True), required=True, help="input event file(s)")
@click.option('-o', '--output', type=str, required=True, help="output concatenated table")
def process_vep(input, output):
    events = []
    for input_path in input:
        data = normalize_variant_table(input_path, add_change=False)
        events.append(data)
    events = pd.concat(events)
    events.to_csv(output, index=False)

if __name__ == "__main__":
    process_vep()