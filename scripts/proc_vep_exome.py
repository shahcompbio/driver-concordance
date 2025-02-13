import click
import re
import os

import pandas as pd

from gdan_hcmi_tools.process import VEP_meta, read_vcf, extract_vaf, transform_aa_change

def make_variant_data(vcf_path, driver_genes, sample_id):
    """
    Processes a VCF file to extract relevant variant information for driver genes.
    
    Parameters:
    - vcf_dir (str): Directory path where VCF files are stored.
    - driver_genes (list): List of driver genes for filtering.
    - sample_id (str): Unique sample identifier.
    
    Returns:
    - DataFrame: Filtered variant data with columns for case, sample, variant, gene, event, impact, and VAF.
    """
    vep_config = VEP_meta()
    case_id = sample_id[:17]
    assert os.path.exists(vcf_path), f'ERROR: {vcf_path} does not exist for sample ID {sample_id}'

    df = read_vcf(vcf_path)
    variant_data = []
    variant_seen = set()  # For tracking duplicates
    driver_gene_set = set(driver_genes)  # Convert once for efficient lookup

    for _, row in df.iterrows():
        variant = parse_variant(row)
        vaf = extract_vaf(row, sample_id)
        vep_isoforms = extract_vep_data(row['info'])
        
        for vep_isoform in vep_isoforms:
            vep_data = parse_vep_fields(vep_isoform, vep_config)
            gene_symbol = vep_data['SYMBOL']
            cdna_change = vep_data['HGVSc'].split(':')[-1]
            ptn_change = vep_data['HGVSp'].split(':')[-1].replace('%3D', '=')
            if not vep_data or vep_data['Allele'] == '__UNKNOWN__':
                continue
            
            if not is_driver_gene(gene_symbol, driver_gene_set):
                continue
            
            impact = vep_data['IMPACT']
            if impact not in vep_config.impact_selected:
                continue
            
            for consequence in vep_data['Consequence'].split('&'):
                change = transform_aa_change(ptn_change)
                if 'splice' in consequence:
                    change = cdna_change
                change = f'{gene_symbol} {change}'

                if should_save_variant(gene_symbol, consequence, impact, vep_config):
                    event = get_event_type(consequence, vep_config)
                    variant_tuple = (case_id, sample_id, variant, gene_symbol, event, impact, vaf)
                    
                    if variant_tuple not in variant_seen:
                        variant_seen.add(variant_tuple)
                        variant_data.append(variant_tuple)

    return pd.DataFrame(variant_data, columns=['case', 'sample', 'variant', 'gene', 'event', 'impact', 'vaf'])

def parse_variant(row):
    """Extracts and formats variant details."""
    chrom, pos, ref, alt = row['chrom'], row['pos'], row['ref'], row['alt'].replace(',__UNKNOWN__', '')
    return f'{chrom}:{pos}:{ref}:{alt}'

def extract_vep_data(info_field):
    """Extracts VEP annotations from the info field of a VCF row."""
    match = re.search('CSQ=([^;]+)', info_field)
    return match.groups()[0].split(',') if match else []

def parse_vep_fields(vep_isoform, vep_config):
    """Parses VEP fields and maps them to a dictionary based on VEP configuration."""
    fields = vep_isoform.split('|')
    if len(fields) < vep_config.fmt_size:
        return None
    return {vep_config.vep_fmts[i]: fields[i] for i in range(vep_config.fmt_size)}

def is_driver_gene(gene_symbol, driver_gene_set):
    """Checks if a gene symbol is in the set of driver genes."""
    return gene_symbol in driver_gene_set

def should_save_variant(gene_symbol, consequence, impact, vep_config):
    """Determines if a variant should be saved based on gene, consequence, and impact."""
    is_non_promoter = (consequence != 'upstream_gene_variant' and consequence in vep_config.consequences_map)
    is_promoter = (gene_symbol == 'TERT' and consequence == 'upstream_gene_variant' and impact in vep_config.impact_selected)
    return is_non_promoter or is_promoter

def get_event_type(consequence, vep_config):
    """Maps a consequence to an event type based on the VEP configuration."""
    return vep_config.layer_map[vep_config.consequences_map[consequence]]


@click.command()
@click.option('--vcf', type=str, help="input vcf dir")
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