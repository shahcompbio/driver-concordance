import glob
import os

import click
import pandas as pd

def get_vcf_path_table(vcf_dir):
    """
    Generates a table of VCF file paths and sample identifiers from a specified directory of VCF files.

    Parameters
    ----------
    vcf_dir : str or Path
        Path to the directory containing VCF files with a ".vep.vcf" extension. 
        These files should represent paired tumor and normal samples.

    Returns
    -------
    DataFrame
        A DataFrame indexed by 'tumor_short' with columns:
        - 'tumor_id': Full identifier for the tumor sample.
        - 'normal_id': Full identifier for the normal sample.
        - 'tumor_short': Shortened identifier for the tumor sample.
        - 'normal_short': Shortened identifier for the normal sample.
        - 'path': Absolute path to the VCF file.
    """
    vcf_paths = glob.glob(f'{vcf_dir}/*.vep.vcf')
    vcf_path_df = pd.DataFrame(columns=['tumor_id', 'path'])
    for vcf_path in vcf_paths:
        _, filename = os.path.split(vcf_path)
        tumor_id = filename.replace('.vep.vcf', '')
        vcf_path = os.path.abspath(vcf_path)
        field = [tumor_id, vcf_path]
        vcf_path_df.loc[vcf_path_df.shape[0]] = field
    vcf_path_df = vcf_path_df.sort_values('tumor_id').reset_index(drop=True)
    vcf_path_df.set_index('tumor_id', inplace=True)
    return vcf_path_df

@click.command()
@click.option('--vcf_dir', type=str, help="input directory with simple somatic mutation vcf files after VEP annotation")
@click.option('--vcf_df_path', type=str, help="output table of sample - path map")
def extract_vep_sample_ids(vcf_dir, vcf_df_path):
    vcf_path_df = get_vcf_path_table(vcf_dir)
    vcf_path_df.to_csv(vcf_df_path, index=True)

if __name__ == "__main__":
    extract_vep_sample_ids()