import pickle

import click

from gdan_hcmi_tools.process import Cohort, read_tm_link

@click.command()
@click.option('--links_path', type=str, help="input tumor-model linkage table")
@click.option('--ct_path', type=str, help="input cancer type table")
@click.option('--pp_path', type=str, help="input purity ploidy table path")
@click.option('--cohort_data', type=str, help="output pickle path for the processed cohort data")
def make_cohort_data(links_path, pp_path, ct_path, cohort_data):
    links = read_tm_link(links_path)
    cohort = Cohort(links)
    cohort.add_purity_ploidy(pp_path)
    cohort.add_cancer_type(ct_path)
    with open(cohort_data, 'wb') as handle:
        pickle.dump(cohort, handle)
        
if __name__ == "__main__":
    make_cohort_data()