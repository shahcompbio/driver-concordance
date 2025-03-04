cluster_fmt="sbatch --partition=componc_cpu --cpus-per-task={threads} --mem={resources.mem_mb} --job-name={rule}.{wildcards} --error=logs/{rule}/{rule}.{wildcards}.%j.err --output=logs/{rule}/{rule}.{wildcards}.%j.out --time=24:00:00"
cmd="snakemake --executor cluster-generic"
cmd="$cmd --cluster-generic-submit-cmd \"$cluster_fmt\""
cmd="$cmd --profile profile/"
cmd="$cmd --singularity-args \"--cleanenv --bind /data1 --bind /home\""
cmd="$cmd -p"
cmd="$cmd --jobs 600"
cmd="$cmd --configfile config.yaml"
cmd="$cmd --use-singularity"
cmd="$cmd --dry-run"
cmd="$cmd --forceall"
cmd="$cmd --rulegraph"

# echo $cmd
eval $cmd
