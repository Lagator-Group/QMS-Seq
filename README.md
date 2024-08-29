# QMS-Seq
A pipeline for Quantitative Mutational Scan sequencing to identify resistance mutations across a whole bacterial genome

# Usage 
```
sh sweep_script.sh manifest.tsv reference.fasta
```
## manifest.tsv:
```
sample_1_name      /path/to/forward/read/sample_1_read_1.fastq.gz      /path/to/reverse/read/sample_1_read_2.fastq.gz
sample_2_name      /path/to/forward/read/sample_2_read_1.fastq.gz      /path/to/reverse/read/sample_2_read_2.fastq.gz
```
