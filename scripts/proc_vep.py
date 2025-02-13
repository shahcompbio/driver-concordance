import re

import click
import numpy as np
import pandas as pd

from gdan_hcmi_tools.process import VEP_meta, read_vcf, transform_aa_change


def make_variant_data(vcf_path, driver_genes, sample_id):
    """
    Generate a DataFrame of variant data for a given sample by processing a VCF file.

    This function reads variant data for a specified sample from a VCF file, filters it for driver genes, and maps each variant
    to a specific event type based on VEP annotations. Each unique variant is recorded with metadata, including its variant allele frequency (VAF).

    Parameters
    ----------
    vcf_path_df : pd.DataFrame
        A DataFrame containing sample IDs as indices and paths to their respective VCF files in a column named 'path'.
    driver_genes : list
        A list of driver gene symbols to filter for, ensuring only variants in driver genes are recorded.
    sample_id : str
        The sample identifier for which to process variant data.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing variant information for the sample, with columns:
        - 'case': a truncated identifier of the sample (first 17 characters of `sample_id`).
        - 'sample': the full sample ID.
        - 'variant': the variant details in "chrom:pos:ref:alt" format.
        - 'gene': the gene symbol where the variant was found.
        - 'event': the type of event, derived from VEP consequence mapping.
        - 'impact': the variant impact, based on VEP annotations.
        - 'vaf': the variant allele frequency.
    """
    vep_config = VEP_meta()
    variant_data = []
    variant_tuple_saved = {}
    case_id = sample_id[:17]
    df = read_vcf(vcf_path)
    for _, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        var = f'{chrom}:{pos}:{ref}:{alt}'
        res = re.search('CSQ=([^;]+)', row['info'])
        flag_promoter_mis = 'Promoter_Mis' in row['info']
        text = res.groups()[0]
        vep_isoforms = text.split(',')
        fmt = row['FORMAT'].split(':')
        gt = row[sample_id].split(':')
        gt = dict(zip(fmt, gt))
        vaf = np.nan
        if 'AF' in gt:
            vaf = sum([float(f) for f in gt['AF'].split(',')])
        elif 'Mutect2Multi_AF' in gt:
            vaf = sum([float(f) for f in gt['Mutect2Multi_AF'].split(',')])
        vep_results = [s.split('|') for s in vep_isoforms]
        for _, vep_result in enumerate(vep_results):
            vep_data = {vep_config.vep_fmts[i]:vep_result[i] for i in range(vep_config.fmt_size)}
            gene_symbol = vep_data['SYMBOL']
            cdna_change = vep_data['HGVSc'].split(':')[-1]
            ptn_change = vep_data['HGVSp'].split(':')[-1].replace('%3D', '=')
            if gene_symbol not in set(driver_genes): # filter out non drivers
                continue
            impact = vep_data['IMPACT']
            consequences = vep_data['Consequence'].split('&')
            for consequence in consequences:
                change = transform_aa_change(ptn_change)
                if 'splice' in consequence:
                    change = cdna_change
                change = f'{gene_symbol} {change}'
                flag_non_promoter = (consequence != 'upstream_gene_variant' and consequence in vep_config.consequences_map)
                flag_promoter = (consequence == 'upstream_gene_variant' and flag_promoter_mis)
                flag_save = flag_non_promoter or flag_promoter
                if flag_save:
                    vkey = vep_config.consequences_map[consequence]
                    event = vep_config.layer_map[vkey]
                    variant_tuple = (case_id, sample_id, var, gene_symbol, event, impact, vaf)
                    if variant_tuple not in variant_tuple_saved:
                        variant_tuple_saved[variant_tuple] = True
                        variant_data.append(variant_tuple)
    variant_data = pd.DataFrame(variant_data, columns=['case', 'sample', 'variant', 'gene', 'event', 'impact', 'vaf'])
    return variant_data


@click.command()
@click.option('--vcf', type=str, help="input table of sample - path map")
@click.option('--drivers_mutation_path', type=str, help="input driver genes table for simple somatic mutations")
@click.option('--sample_id', type=str, help="sample id")
@click.option('--variant_data_path', type=str, help="output variant table")
def process_vep(vcf, drivers_mutation_path, sample_id, variant_data_path):
    drivers_mutation_df = pd.read_table(drivers_mutation_path, comment='#')
    drivers_mutation = drivers_mutation_df['Gene'].unique().tolist()

    variant_data = make_variant_data(vcf, drivers_mutation, sample_id)
    variant_data.to_csv(variant_data_path, index=False)

if __name__ == "__main__":
    process_vep()