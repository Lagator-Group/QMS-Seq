# QMS-Seq: Quantitative Microbial Sequencing Pipeline

**QMS-Seq** is a bioinformatics pipeline designed to identify low-frequency single nucleotide variants (SNVs) and indels from sequencing data. It specifically utilizes GATK's Base Quality Score Recalibration (BQSR) guided by initial high-confidence calls to improve the accuracy of subsequent **LoFreq** variant calling.

See the publication in [Nature Communications](https://www.nature.com/articles/s41467-025-56050-2) for more details.

## Workflow Overview

1. **Read QC & Trimming:** `TrimGalore` filters low-quality bases (Q < 30) and removes adapters.
2. **Alignment:** `Bowtie2` aligns reads to a provided reference.
3. **Bootstrap Recalibration:**
* Initial SNVs are called via `bcftools` (> 1% prevalence).
* These SNVs serve as "known sites" for GATK `BaseRecalibrator` to model platform-specific errors.
4. **Variant Calling:** `LoFreq` performs the final sensitive variant calling (including indels) on the recalibrated data.


## Installation & Dependencies

Ensure the following tools are in your `PATH`:

* **Alignment/Processing:** `bowtie2`, `samtools`, `picard`, `bcftools`
* **QC:** `trim_galore`, `fastqc`
* **Variant Calling:** `GATK` (v4+), `lofreq`


## Usage

### 1. Prepare the Manifest File

The pipeline requires a tab-separated `manifest.tsv` with three columns (no header):
`Sample_ID` | `Path_to_R1` | `Path_to_R2`

**Example `manifest.tsv`:**

```text
SampleA    data/SampleA_L001_R1.fq.gz    data/SampleA_L001_R2.fq.gz
SampleB    data/SampleB_L001_R1.fq.gz    data/SampleB_L001_R2.fq.gz

```

### 2. Run the Script

Execute the script by providing the manifest and the reference FASTA:

```bash
bash sweep_script.sh manifest.tsv reference.fasta

```

---

## Output Structure

All results are organized in the `QMS/` directory:

* `TRIMMING/`: Cleaned FastQ files and FastQC reports.
* `ALIGNING/`: Reference index and initial BAM files.
* `EXPECTED/`: VCFs used for BQSR bootstrapping.
* `RECALIBRATING/`: Recalibrated BAMs and GATK tables.
* `LOFREQ/`: **Final output files** (`*_variants.vcf`).

## Note on Computational Resources

The pipeline uses multiple threads for alignment (`-p 4`) and sorting (`-@ 8`). Ensure your environment has sufficient CPU and RAM for GATK operations, which are memory-intensive.

