# Make anndata h5ad files

def get_hcmi_consensus_path(wildcards):
    paths = PATHS_CONSENSUS[PATHS_CONSENSUS['sample_id']==wildcards.sample]
    assert paths.shape[0] == 1, paths
    path = paths['path'].squeeze()
    return path

def get_hcmi_remixt_path(wildcards):
    paths = PATHS_REMIXT[PATHS_REMIXT['sample_id']==wildcards.sample]
    assert paths.shape[0] == 1, paths
    path = paths['path'].squeeze()
    return path

rule make_cntable_hcmi_consensus:
    input:
        cn = get_hcmi_consensus_path,
    output:
        cn = os.path.join(base, config['HCMI']['cn']['consensus']['intermediate_dir'], '{sample}.cn.csv'),
    shell:
        '{python_bin} ../scripts/make_cntable_hcmi_consensus.py -i {input.cn} -b 50000 -o {output.cn} -s {wildcards.sample}'

rule make_cntable_hcmi_remixt:
    input:
        cn = get_hcmi_remixt_path,
    output:
        cn = os.path.join(base, config['HCMI']['cn']['remixt']['intermediate_dir'], '{sample}.cn.csv'),
    shell:
        '{python_bin} ../scripts/make_cntable_hcmi_remixt.py -i {input.cn} -b 50000 -o {output.cn} -s {wildcards.sample}'

rule make_anndata_hcmi_consensus:
    input:
        expand(
            os.path.join(base, config['HCMI']['cn']['consensus']['intermediate_dir'], '{sample}.cn.csv'),
            sample=SAMPLES_CONSENSUS,
        )
    output:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata']['segment']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/make_anndata_hcmi_consensus.py {input} -o {output.h5ad}'

rule make_anndata_hcmi_remixt:
    input:
        expand(
            os.path.join(base, config['HCMI']['cn']['remixt']['intermediate_dir'], '{sample}.cn.csv'),
            sample=SAMPLES_REMIXT,
        )
    output:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata_remixt']['segment']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/make_anndata_hcmi_remixt.py {input} -o {output.h5ad}'

rule make_anndata_tcga:
    input:
        cn_dir = os.path.join(base, config['TCGA']['cn']['paths']),
    output:
        h5ad = os.path.join(base, config['TCGA']['cn']['anndata']['segment']),
    params:
        bin_size = 50000,
    resources:
        mem_mb=240000,
    shell:
        '{python_bin} ../scripts/make_anndata_tcga.py -c {input.cn_dir} -o {output.h5ad} '
        '-b {params.bin_size}'

rule add_obs_to_anndata_tcga:
    input:
        h5ad = os.path.join(base, config['TCGA']['cn']['anndata']['segment']),
        metadata = os.path.join(base, config['TCGA']['metadata']),
    output:
        h5ad = os.path.join(base, config['TCGA']['cn']['anndata']['annotated']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/add_obs_to_anndata.py -i {input.h5ad} '
        '-c {input.metadata} -s TCGA '
        '-o {output.h5ad}'

rule add_obs_to_anndata_hcmi_consensus:
    input:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata']['segment']),
        metadata = os.path.join(base, config['HCMI']['metadata']),
        pairs_map = os.path.join(base, config['HCMI']['pairs_map']),
        analysis_tracker = os.path.join(base, config['HCMI']['analysis_tracker']),
        purity_ploidy_table = os.path.join(base, config['HCMI']['purity_ploidy_table'])
    output:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata']['annotated']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/add_obs_to_anndata.py -i {input.h5ad} '
        '-s HCMI-consensus '
        '-c {input.metadata} '
        '-p {input.pairs_map} '
        '-pp {input.purity_ploidy_table} '
        '-m {input.analysis_tracker} '
        '-o {output.h5ad}'

rule add_obs_to_anndata_hcmi_remixt:
    input:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata_remixt']['segment']),
        metadata = os.path.join(base, config['HCMI']['metadata']),
        pairs_map = os.path.join(base, config['HCMI']['pairs_map']),
        analysis_tracker = os.path.join(base, config['HCMI']['analysis_tracker']),
        purity_ploidy_table = os.path.join(base, config['HCMI']['purity_ploidy_table'])
    output:
        h5ad = os.path.join(base, config['HCMI']['cn']['anndata_remixt']['annotated']),
    resources:
        mem_mb=80000,
    shell:
        '{python_bin} ../scripts/add_obs_to_anndata.py -i {input.h5ad} '
        '-s HCMI-remixt '
        '-c {input.metadata} '
        '-p {input.pairs_map} '
        '-pp {input.purity_ploidy_table} '
        '-m {input.analysis_tracker} '
        '-o {output.h5ad}'
