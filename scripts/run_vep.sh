[ $# -ne 4 ] && { echo -e "\nUsage: $0 <in.vcf> <out.vcf> <vep.dir> <vep.fasta>\n"; exit 1; }
in_vcf=$1
out_vcf=$2
vep_dir=$3
vep_fasta=$4

vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats \
    --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein \
    --biotype --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number \
    --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir ${vep_dir} \
    --fasta ${vep_fasta} \
    --input_file $in_vcf \
    --output_file $out_vcf \
    --polyphen b --af --af_1kg --regulatory \
    --buffer_size 500 --force_overwrite
