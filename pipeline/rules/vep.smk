rule maf_to_vcf_broad:
    input:
        maf=os.path.join(data_dir, 'broad/FF_exome/renamed/{sample}.maf')
    output:
        vcf=os.path.join(base, 'maf2vcf_exome/{sample}.vcf')
    params:
        ref=config['fasta'],
        outdir=os.path.join(base, 'maf2vcf_exome'),
    singularity:
        'docker://cgrlab/vcf2maf:latest',
    shell:
        'perl /opt/vcf2maf/maf2vcf.pl --input-maf {input.maf} --output-dir {params.outdir} --ref-fasta {params.ref}'

rule run_vep_broad:
    input:
        vcf=os.path.join(base, 'data/maf2vcf_exome/{sample}.vcf'),
    output:
        vcf=os.path.join(base, 'data/vep_vcf_exome/{sample}.vep.vcf'),
    params:
        vep_dir=config['vep']['data'],
        vep_fasta=config['vep']['ref_fasta'],
    resources:
        mem_mb = 24000,
    singularity:
        'docker://ensemblorg/ensembl-vep:release_105.0',
    shell:
        'bash ../scripts/run_vep.sh {input.vcf} {output.vcf} {params.vep_dir} {params.vep_fasta}'

rule run_vep:
    input:
        vcf=os.path.join(data_dir, 'nygc/plusReadSupportAndM2MultiPass/{pair}.plusReadSupportAndM2MultiPass.intersection.vcf'),
    output:
        vcf=os.path.join(base, 'data/vep_vcf/{pair}.vep.vcf'),
    params:
        vep_dir=config['vep']['data'],
        vep_fasta=config['vep']['ref_fasta'],
    resources:
        mem_mb = 24000,
    singularity:
        'docker://ensemblorg/ensembl-vep:release_105.0',
    shell:
        'bash ../scripts/run_vep.sh {input.vcf} {output.vcf} {params.vep_dir} {params.vep_fasta}'