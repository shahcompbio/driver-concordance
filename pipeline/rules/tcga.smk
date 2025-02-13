# Create plots related to HCMI-TCGA comparison

rule make_genelevel_anndata_tcga:
    input:
        h5ad = os.path.join(base, config['TCGA']['cn']['anndata']['annotated']),
        curated_events = os.path.join(base, config['HCMI']['cn']['curated_events']), # reuse HCMI curated genes
        refflat = config['refflat'],
    output:
        h5ad = os.path.join(base, config['TCGA']['cn']['anndata']['genelevel']),
        events_cn = os.path.join(base, config['TCGA']['cn']['events']),
    shell:
        '{python_bin} ../scripts/make_genelevel_anndata_tcga.py -i {input.h5ad} '
        '-d {input.curated_events} -r {input.refflat} -og {output.h5ad} -oe {output.events_cn}'

rule make_snvindel_anndata_tcga:
    input:
        maf_path_table = os.path.join(base, config['TCGA']['mutation']['maf_table']),
        drivers_mutation = os.path.join(base, config['HCMI']['mutation']['driver_genes']), # reuse HCMI driver genes
    output:
        adata_mutation = os.path.join(base, config['TCGA']['mutation']['anndata']['genelevel']),
        events_mutation = os.path.join(base, config['TCGA']['mutation']['events']),
    shell:
        '{python_bin} ../scripts/make_snvindel_anndata_tcga.py '
        '--maf_path_table {input.maf_path_table} '
        '--drivers_mutation_path {input.drivers_mutation} '
        '--mutation_anndata {output.adata_mutation} '
        '--mutation_events {output.events_mutation}'

rule make_sample_list_tcga:
    input:
        adata_cn = os.path.join(base, config['TCGA']['cn']['anndata']['genelevel']),
        adata_mutation = os.path.join(base, config['TCGA']['mutation']['anndata']['genelevel']),
        events_cn = os.path.join(base, config['TCGA']['cn']['events']),
        events_mutation = os.path.join(base, config['TCGA']['mutation']['events']),
    output:
        samples_table = os.path.join(base, config['TCGA']['cohort_table']), # also act as supplementary table for fig3A
        summary = os.path.join(base, config['TCGA']['summary']),
        counts = os.path.join(base, config['TCGA']['counts']),
    params:
        genes_to_exclude = ' '.join(config['HCMI']['genes_to_exclude']), # reuse HCMI genes_to_exclude
    shell:
        '{python_bin} ../scripts/make_event_notes_tcga.py '
        '--cn_anndata {input.adata_cn} --mutation_anndata {input.adata_mutation} '
        '--events_cn {input.events_cn} --events_mutation {input.events_mutation} '
        '--summary {output.summary} --counts {output.counts} '
        '--cohort {output.samples_table} --genes_to_exclude {params.genes_to_exclude}'
        
rule make_data_for_fig3A:
    input:
        summary = os.path.join(base, config['TCGA']['summary']),
        counts = os.path.join(base, config['TCGA']['counts']),
        notes_hcmi = os.path.join(base, config['HCMI']['notes']),
        cohort_hcmi = os.path.join(base, config['HCMI']['cohort_table']),
    output:
        freq = os.path.join(base, config['HCMI']['freq']),
    shell:
        '{python_bin} ../scripts/make_data_for_fig3A.py --summary {input.summary} --counts {input.counts} '
        '--notes_hcmi {input.notes_hcmi} --cohort_hcmi {input.cohort_hcmi} '
        '--freq {output.freq}'

rule plot_fig3A:
    input:
        freq = os.path.join(base, config['HCMI']['freq']),
    output:
        pdf = os.path.join(base, config['plots']['fig3A']),
    params:
        top_n_genes = 5,
        min_samples = 10,
    shell:
        '{python_bin} ../scripts/plot_fig3A.py --freq {input.freq} --fig_path {output.pdf} '
        '--min_samples_per_cancer_type {params.min_samples} --n_top_genes_to_show {params.top_n_genes}'

rule plot_fig3A_supplementary:
    input:
        freq = os.path.join(base, config['HCMI']['freq']),
    output:
        pdf = os.path.join(base, config['plots']['fig3A_supplementary']),
    params:
        top_n_genes = 5,
        min_samples = 1, # show all cohort as long as min 1 sample exists
        figwidth = 12,
        figheight = 8,
        nrow = 3,
    shell:
        '{python_bin} ../scripts/plot_fig3A.py --freq {input.freq} --fig_path {output.pdf} '
        '--min_samples_per_cancer_type {params.min_samples} --n_top_genes_to_show {params.top_n_genes} '
        '--nrow {params.nrow} --figwidth {params.figwidth} --figheight {params.figheight}'
