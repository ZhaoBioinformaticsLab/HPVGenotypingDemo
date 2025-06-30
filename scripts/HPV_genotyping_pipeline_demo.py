from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
REFERENCE = "References/HPV_pave_Downloaded_06_28_2025.fa"
BOWTIE2_INDEX_PREFIX = "References/HPV_bowtie2_index"
SRR_LIST = [
    "SRR1186008", "SRR1186009", "SRR1186010",
    "SRR1186011", "SRR1186012", "SRR1186013",
    "SRR1186014", "SRR1186015", "SRR1186016"
]
THREADS = 8
OUTDIR = Path("hpv_typing_results")
OUTDIR.mkdir(exist_ok=True)

def run_cmd(cmd, msg):
    print(f"[RUNNING] {msg}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"[DONE] {msg}")

# Step 1: Index reference for Bowtie2 if not already done
if not Path(BOWTIE2_INDEX_PREFIX + ".1.bt2").exists():
    run_cmd(f"bowtie2-build {REFERENCE} {BOWTIE2_INDEX_PREFIX}", "Building Bowtie2 index")

# Step 2: Trim reads, align with Bowtie2, and process each sample
for srr in SRR_LIST:
    fq_main = Path(srr) / f"{srr}.fastq"
    fq_alt = Path(srr) / f"{srr}_1.fastq"
    fq_raw = fq_main if fq_main.exists() else fq_alt

    if not fq_raw.exists():
        print(f"[SKIP] No FASTQ found for {srr}")
        continue

    fq_clean = OUTDIR / f"{srr}.trimmed.fastq"
    sam = OUTDIR / f"{srr}.sam"
    bam = OUTDIR / f"{srr}.bam"
    sorted_bam = OUTDIR / f"{srr}.sorted.bam"

    run_cmd(f"fastp -i {fq_raw} -o {fq_clean} -q 20 -u 30 -l 50 --thread {THREADS} -h {OUTDIR}/{srr}_fastp.html -j {OUTDIR}/{srr}_fastp.json", f"Trimming {srr}")
    run_cmd(f"bowtie2 -x {BOWTIE2_INDEX_PREFIX} -U {fq_clean} -p {THREADS} -S {sam}", f"Aligning {srr} with Bowtie2")
    run_cmd(f"samtools view -bS {sam} > {bam}", f"Converting to BAM: {srr}")
    run_cmd(f"samtools sort {bam} -o {sorted_bam}", f"Sorting BAM: {srr}")
    run_cmd(f"samtools index {sorted_bam}", f"Indexing BAM: {srr}")

# Step 3: Count mapped reads per reference for each sample
bam_files = list(OUTDIR.glob("*.sorted.bam"))
count_matrix = {}

for bam in bam_files:
    srr = bam.stem.replace(".sorted", "")
    result = subprocess.run(
        f"samtools idxstats {bam}", shell=True,
        stdout=subprocess.PIPE, text=True, check=True
    )
    lines = result.stdout.strip().split("\n")
    count_matrix[srr] = {}
    for line in lines:
        ref, length, mapped, _ = line.split("\t")
        count_matrix[srr][ref] = int(mapped)

# Step 4: Write genotype report
df = pd.DataFrame.from_dict(count_matrix, orient="index").fillna(0).astype(int)
df.to_csv(OUTDIR / "HPV_genotype_report.tsv", sep="\t")

# Generate condensed report (exclude columns with all zeros)
df_condensed = df.loc[:, (df != 0).any(axis=0)]
df_condensed.to_csv(OUTDIR / "HPV_genotype_report_condensed.tsv", sep="\t")

# Step 5: Identify dominant HPV type per sample
dominant_types = df.idxmax(axis=1)
read_counts = df.max(axis=1)
summary = pd.DataFrame({
    "Dominant_HPV_Type": dominant_types,
    "Mapped_Reads": read_counts
})
summary.to_csv(OUTDIR / "HPV_dominant_summary.tsv", sep="\t")
print(summary)
print("\n✅ Summary of dominant HPV types saved to:", OUTDIR / "HPV_dominant_summary.tsv")

# Step 6: Plotting using condensed data
df_plot = pd.read_csv(OUTDIR / "HPV_genotype_report_condensed.tsv", sep="\t", index_col=0)

# Step 7: Write summary text file
summary_file = OUTDIR / "HPV_genotype_summary.txt"
with open(summary_file, "w") as f:
    f.write("HPV Genotyping Summary Report\n")
    f.write("==================================\n")
    f.write(f"Total samples processed: {len(df)}\n\n")

    f.write("Detected HPV Types (non-zero):\n")
    detected_types = df_condensed.columns.tolist()
    for hpv_type in detected_types:
        total_reads = df[hpv_type].sum()
        f.write(f"  - {hpv_type}: {total_reads} mapped reads\n")

    f.write("\nDominant HPV Type per Sample:\n")
    for sample, row in df.iterrows():
        dom_type = row.idxmax()
        dom_reads = row[dom_type]
        f.write(f"  - {sample}: {dom_type} ({dom_reads} reads)\n")

print("📝 Summary file written to:", summary_file)

# Bar plot (stacked)
plt.figure(figsize=(12, 6))
df_sorted = df_plot.loc[df_plot.sum(axis=1).sort_values(ascending=False).index]
df_sorted.plot(kind="bar", stacked=True, colormap="tab20", edgecolor="black")
plt.ylabel("Mapped Reads")
plt.title("HPV Genotype Abundance per Sample (Stacked Bar)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(OUTDIR / "HPV_genotype_barplot.png")
plt.close()

# Heatmap (raw read counts)
plt.figure(figsize=(12, 6))
sns.heatmap(df_plot.T, cmap="viridis", annot=True, fmt="d", cbar_kws={"label": "Read Count"})
plt.title("HPV Genotype Read Count Heatmap")
plt.xlabel("Sample")
plt.ylabel("HPV Genotype")
plt.tight_layout()
plt.savefig(OUTDIR / "HPV_genotype_heatmap.png")
plt.close()

print("✅ Genotyping pipeline completed. Results saved in:", OUTDIR)
