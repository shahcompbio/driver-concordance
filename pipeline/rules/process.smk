# Create plots related to tumor-model comparison

def get_vep_vcf_path(wildcards):
    paths = VEP_PATHS.copy()
    assert 'tumor_short' in paths.columns
    assert wildcards.sample in paths['tumor_short'].values, wildcards.sample
    paths = paths[paths['tumor_short']==wildcards.sample]
    assert paths.shape[0] == 1, paths
    path = paths['path'].squeeze()
    return path

def get_vep_vcf_path_exome(wildcards):
    paths = VEP_EXOME_PATHS.copy()
    assert wildcards.sample in paths['tumor_id'].values, wildcards.sample
    paths = paths[paths['tumor_id']==wildcards.sample]
    assert paths.shape[0] == 1, paths
    path = paths['path'].squeeze()
    return path

rule make_genelevel_anndata_hcmi:
    input:
        h5ad_consensus = os.path.join(base, config['HCMI']['cn']['anndata']['annotated']),
        h5ad_remixt = os.path.join(base, config['HCMI']['cn']['anndata_remixt']['annotated']),
        curated_events = os.path.join(base, config['HCMI']['cn']['curated_events']),
        refflat = config['refflat'],
        cohort_data = os.path.join(base, config['HCMI']['cohort']),
    output:
        h5ad_consensus = os.path.join(base, config['HCMI']['cn']['anndata']['genelevel']),
        h5ad_remixt = os.path.join(base, config['HCMI']['cn']['anndata_remixt']['genelevel']),
        events = os.path.join(base, config['HCMI']['cn']['events']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/make_genelevel_anndata_hcmi.py -ic {input.h5ad_consensus} -ir {input.h5ad_remixt} '
        '-d {input.curated_events} -r {input.refflat} -c {input.cohort_data} -og {output.h5ad_consensus} -or {output.h5ad_remixt} -oe {output.events}'

rule extract_ssm_sample_ids:
    input:
        vcf_dir = os.path.join(base, config['HCMI']['mutation']['vcf_dir']),
    output:
        vcf_path_df = os.path.join(base, config['HCMI']['mutation']['vcf_table']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/extract_vep_sample_ids.py --vcf_dir {input.vcf_dir} --vcf_df_path {output.vcf_path_df}'

rule make_cohort_data:
    input:
        pp = os.path.join(base, config['HCMI']['purity_ploidy_table']),
        links = os.path.join(base, config['HCMI']['pairs_map']),
        cancer_types = os.path.join(base, config['HCMI']['metadata']),
    output:
        cohort_data = os.path.join(base, config['HCMI']['cohort']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/make_cohort_data.py --links_path {input.links} '
        '--ct_path {input.cancer_types} --pp_path {input.pp} '
        '--cohort_data {output.cohort_data}'

rule process_vep:
    input:
        vcf = get_vep_vcf_path,
        drivers_mutation = os.path.join(base, config['HCMI']['mutation']['driver_genes']),
    output:
        variant_data_path = os.path.join(base, config['HCMI']['mutation']['variant_data_dir'], '{sample}.csv'),
    resources:
        mem_mb=4000,
    shell:
        '{python_bin} ../scripts/proc_vep.py --vcf {input.vcf} --drivers_mutation_path {input.drivers_mutation} '
        '--sample_id {wildcards.sample} --variant_data_path {output.variant_data_path}'

rule process_vep_exome:
    input:
        vcf = get_vep_vcf_path_exome,
        drivers_mutation = os.path.join(base, config['HCMI']['mutation']['driver_genes']),
    output:
        variant_data_path = os.path.join(base, config['HCMI']['mutation']['variant_data_dir_exome'], '{sample}.csv'),
    shell:
        '{python_bin} ../scripts/proc_vep_exome.py --vcf {input.vcf} --drivers_mutation_path {input.drivers_mutation} '
        '--sample_id {wildcards.sample} --variant_data_path {output.variant_data_path}'

rule aggregate_variant_table:
    input:
        cohort_data = os.path.join(base, config['HCMI']['cohort']),
        variant_data_paths = expand(
            os.path.join(base, config['HCMI']['mutation']['variant_data_dir'], '{sample}.csv'), 
            sample=VEP_SAMPLES
        ),
    output:
        aggregated_variants = os.path.join(base, config['HCMI']['mutation']['events']),
    params:
        variant_data_dir = os.path.join(base, config['HCMI']['mutation']['variant_data_dir']),
    shell:
        '{python_bin} ../scripts/aggregate_variant_tables.py '
        '--cohort_data {input.cohort_data} --variant_data_dir {params.variant_data_dir} '
        '--aggregated_output {output.aggregated_variants}'

rule aggregate_variant_table_exome:
    input:
        cohort_data = os.path.join(base, config['HCMI']['cohort']),
        variant_data_paths = expand(
            os.path.join(base, config['HCMI']['mutation']['variant_data_dir_exome'], '{sample}.csv'), 
            sample=VEP_EXOME_SAMPLES,
        ),
    output:
        aggregated_variants = os.path.join(base, config['HCMI']['mutation']['events_exome']),
    params:
        variant_data_dir = os.path.join(base, config['HCMI']['mutation']['variant_data_dir_exome']),
    shell:
        '{python_bin} ../scripts/aggregate_variant_tables.py '
        '--cohort_data {input.cohort_data} --variant_data_dir {params.variant_data_dir} '
        '--aggregated_output {output.aggregated_variants}'

rule make_snvindel_anndata:
    input:
        cohort_data = os.path.join(base, config['HCMI']['cohort']),
        drivers_mutation = os.path.join(base, config['HCMI']['mutation']['driver_genes']),
        variant_data_paths = expand(os.path.join(base, config['HCMI']['mutation']['variant_data_dir'], '{sample}.csv'), sample=VEP_SAMPLES),
    output:
        adata_mutation = os.path.join(base, config['HCMI']['mutation']['anndata']['genelevel']),
    params:
        variant_data_dir = os.path.join(base, config['HCMI']['mutation']['variant_data_dir']),
    shell:
        '{python_bin} ../scripts/make_snvindel_anndata.py '
        '--cohort_data {input.cohort_data} '
        '--drivers_mutation_path {input.drivers_mutation} '
        '--variant_data_dir {params.variant_data_dir} '
        '--mutation_anndata {output.adata_mutation}'

rule concatenate_event_tables:
    input:
        cn = os.path.join(base, config['HCMI']['cn']['events']),
        snvindel = os.path.join(base, config['HCMI']['mutation']['events']),
        snvindel_exome = os.path.join(base, config['HCMI']['mutation']['events_exome']),
    output:
        events = os.path.join(base, config['HCMI']['events']),
    shell:
        '{python_bin} ../scripts/concat_event_tables.py -i {input.cn} {input.snvindel} {input.snvindel_exome} -o {output.events}'

rule make_samples_list:
    input:
        cohort = os.path.join(base, config['HCMI']['cohort']),
        adata_cn = os.path.join(base, config['HCMI']['cn']['anndata']['genelevel']),
        adata_mutation = os.path.join(base, config['HCMI']['mutation']['anndata']['genelevel']),
        qc = os.path.join(base, config['HCMI']['samples_to_exclude']),
    output:
        cohort = os.path.join(base, config['HCMI']['cohort_final']),
    shell:
        '{python_bin} ../scripts/make_final_samples_list.py --cohort {input.cohort} '
        "--cn_anndata {input.adata_cn} --mutation_anndata {input.adata_mutation} --qc '{input.qc}' "
        '-o {output.cohort}'

rule make_event_notes_summary:
    input:
        cohort = os.path.join(base, config['HCMI']['cohort_final']),
        adata_mutation = os.path.join(base, config['HCMI']['mutation']['anndata']['genelevel']),
        events = os.path.join(base, config['HCMI']['events']),
    output:
        notes = os.path.join(base, config['HCMI']['notes']),
        cohort_table = os.path.join(base, config['HCMI']['cohort_table']),
        summary = os.path.join(base, config['HCMI']['summary']),
        counts = os.path.join(base, config['HCMI']['counts']),
    params:
        genes_to_exclude = ' '.join(config['HCMI']['genes_to_exclude'])
    shell:
        '{python_bin} ../scripts/make_event_notes.py --cohort {input.cohort} '
        '--mutation_anndata {input.adata_mutation} --events_path {input.events} '
        '--notes_path {output.notes} --summary {output.summary} --counts {output.counts} '
        '--cohort_table {output.cohort_table} '
        '--genes_to_exclude {params.genes_to_exclude}'

rule plot_fig2E:
    input:
        summary = os.path.join(base, config['HCMI']['summary']),
        counts = os.path.join(base, config['HCMI']['counts']),
        drivers_cn = os.path.join(base, config['HCMI']['cn']['driver_genes']),
    output:
        pdf = os.path.join(base, config['plots']['fig2E']),
    shell:
        '{python_bin} ../scripts/plot_fig2E.py --summary {input.summary} --counts {input.counts} '
        '--cn_drivers_path {input.drivers_cn} '
        '--fig_path {output.pdf}'


rule cn_comparison_consensus:
    input:
        pairs_map = os.path.join(base, config['HCMI']['pairs_map']),
        anndata_cn = os.path.join(base, config['HCMI']['cn']['anndata']['annotated']),
    output:
        pairwise = os.path.join(base, config['HCMI']['cn']['pairwise']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/make_data_for_fig2C.py '
        '--input_tm_link {input.pairs_map} '
        '--input_cn_adata {input.anndata_cn} '
        '--output_pairwise_csv {output.pairwise}'

rule plot_fig2C:
    input:
        pairwise = os.path.join(base, config['HCMI']['cn']['pairwise']),
    output:
        pdf1 = os.path.join(base, config['plots']['fig2C_1']),
        pdf2 = os.path.join(base, config['plots']['fig2C_2']),
    shell:
        '{python_bin} ../scripts/plot_fig2C.py --pairwise {input.pairwise} -o1 {output.pdf1} -o2 {output.pdf2}'
    
rule plot_fig2D:
    input:
        pairwise = os.path.join(base, config['HCMI']['cn']['pairwise']),
    output:
        pdf = os.path.join(base, config['plots']['fig2D']),
    shell:
        '{python_bin} ../scripts/plot_fig2D.py --pairwise {input.pairwise} -o {output.pdf}'


rule make_genelevel_table:
    input:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata']['annotated']),
        gtf = config['gtf'],
    output:
        table = os.path.join(base, config['HCMI']['cn']['genelevel']['cnv']),
    resources:
        mem_mb=40000,
    shell:
        'python ../scripts/make_genelevel_table.py {input.h5ad} {input.gtf} {output.table}'

rule update_genelevel_table:
    input:
        table = os.path.join(base, config['HCMI']['cn']['genelevel']['cnv']),
        events = os.path.join(base, config['HCMI']['cn']['events']),
        gtf = config['gtf'],
    output:
        table = os.path.join(base, config['HCMI']['cn']['genelevel']['cnv_merged']),
    resources:
        mem_mb=40000,
    shell:
        'python ../scripts/update_genelevel_table.py {input.table} {input.events} {input.gtf} {output.table}'
