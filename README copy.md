# HPV Amplicon Typing Pipeline (Demonstration)

This repository demonstrates a complete pipeline for genotyping Human Papillomavirus (HPV) from NGS amplicon sequencing data using publicly available SRA samples. It performs quality control, read alignment, quantification, summarization, and visualization of HPV genotypes.

## 🧪 Key Features

* Uses [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads to all HPV genomes
* Quantifies mapped reads per HPV reference genome
* Identifies dominant HPV type per sample
* Generates summary report, barplot, and heatmap for easy visualization

## 📁 Project Structure

```
.
├── References/                  # Contains HPV reference genome and Bowtie2 index
├── SRR11860xx/                 # Per-sample directories containing FASTQ files
├── hpv_typing_results/         # Output directory for BAMs, plots, and result tables
├── scripts/
│   └── HPV_genotyping_pipeline_demo.py   # Main Python pipeline script
├── run_hpv_genotyping_pipeline.sh        # SLURM submission script
├── download_hpv_amplicons.sh             # Download script (optional)
└── README.md
```

## ⚙️ Major Pipeline Steps

1. **Index HPV reference genome** using Bowtie2
2. **Trim reads** with `fastp`
3. **Align reads** using `bowtie2`
4. **Convert and sort BAMs**, then index using `samtools`
5. **Count mapped reads** per HPV reference using `samtools idxstats`
6. **Generate reports**:

   * `HPV_genotype_report.tsv`: All genotypes per sample
   * `HPV_genotype_report_condensed.tsv`: Genotypes with non-zero support only
   * `HPV_dominant_summary.tsv`: Dominant HPV type and read count per sample
7. **Create plots**:

   * `HPV_genotype_barplot.png`: Stacked barplot
   * `HPV_genotype_heatmap.png`: Annotated heatmap

## 🧬 Input Data

This demo uses the following SRA samples:

```
SRR1186008 - SRR1186016
```

You can customize the pipeline to analyze your own HPV amplicon data.

## 📦 Reference

Reference genome downloaded from [PAVE](https://pave.niaid.nih.gov/search/search_database) on June 28, 2025, by selecting **all HPV genomes**, saved as:

```
References/HPV_pave_Downloaded_06_28_2025.fa
```

## 🚀 Running the Pipeline

Submit to Biowulf:

```bash
sbatch run_hpv_genotyping_pipeline.sh
```

Or run locally:

```bash
python scripts/HPV_genotyping_pipeline_demo.py
```

## 📊 Example Output

**Dominant HPV Summary (HPV\_dominant\_summary.tsv):**

```
	Dominant_HPV_Type	Mapped_Reads
SRR1186008	HPV16REF|lcl|Human	1029
SRR1186009	HPV16REF|lcl|Human	1014
...
```

**Genotype Report (condensed):**

```
	HPV16REF|lcl|Human	HPV18REF|lcl|Human
SRR1186008	1029	0
SRR1186009	1014	0
...
```

---

Developed by: **Patrick Xuechun Zhao, Ph.D.**

---
