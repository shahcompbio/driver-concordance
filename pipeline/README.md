# GDAN HCMI Project: Driver event concordance pipeline

## Overview
This project is part of the **Genomic Data Analysis Network (GDAN)** initiative within the **Human Cancer Models Initiative (HCMI)**. The primary goal is to assess the concordance of driver mutations between paired **tumor** and **tumor-model** whole genome sequencing (WGS) samples.

## Usage
```bash
cluster_fmt="sbatch --partition=componc_cpu --cpus-per-task={threads} --mem={resources.mem_mb} --job-name={rule}.{wildcards} --error=logs/{rule}/{rule}.{wildcards}.%j.err --output=logs/{rule}/{rule}.{wildcards}.%j.out --time=24:00:00"
cmd="snakemake --executor cluster-generic"
cmd="$cmd --cluster-generic-submit-cmd \"$cluster_fmt\""
cmd="$cmd --profile profile/"
cmd="$cmd --singularity-args \"--cleanenv --bind /path/to/data\""
cmd="$cmd -p"
cmd="$cmd --jobs 600"
cmd="$cmd --configfile config.yaml"
cmd="$cmd --use-singularity"
cmd="$cmd --forceall"
eval $cmd
```

## Contact
For questions or collaboration, feel free to reach out to Seongmin Choi (chois7@mskcc.org).