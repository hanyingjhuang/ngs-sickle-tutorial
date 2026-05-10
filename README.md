# NGS pipeline · reads → variants → sickle cell

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hanyingjhuang/ngs-sickle-tutorial/blob/main/notebook.ipynb)

A click-to-run teaching notebook that takes real public sequencing reads end-to-end through a classical short-read variant-calling pipeline and recovers the sickle cell mutation (*HBB* rs334) in a confirmed carrier.

**No install. Real Linux execution.** Click the Colab badge → cells run real `bwa`, `samtools`, `bcftools`, `fastqc`, `fastp` binaries on a free cloud VM.

## Pipeline

```
FASTQ ──▶ FastQC + fastp ──▶ BWA-MEM ──▶ samtools sort/index ──▶ bcftools mpileup/call ──▶ rs334
 raw         QC  trim         align         BAM                     VCF                    HbS
```

| step | tool | output |
| --- | --- | --- |
| 1 | `head` | inspect reads |
| 2 | `fastqc` | quality report |
| 3 | `fastp` | trimmed.fq + stats |
| 4 | `bwa mem` | aligned.sam |
| 5 | `samtools sort/index/flagstat/depth` | aligned.bam + coverage plot |
| 6 | `samtools mpileup` | pileup window |
| 7 | `bcftools mpileup` + `call` | variants.vcf.gz |
| 8 | parse in Python | rs334 → c.20A>T → p.Glu7Val → HbAS |

## Data

- **Reference** · GRCh37 chr11:5,246,001–5,250,000 (the *HBB* gene + flanking), fetched from UCSC. (`data/ref.fa`, 4 KB)
- **Reads** · 438 real Illumina low-coverage reads from **HG02666**, a Gambian (GWD) sample in the [1000 Genomes Project](https://www.internationalgenome.org/), confirmed *HBB* rs334 heterozygote per the 1000G phase 3 VCF. Reads pre-extracted from the public BAM at chr11:5,246,000–5,250,000. (`data/reads.fq`, 97 KB)

To rebuild from scratch:

```bash
samtools view -h "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02666/alignment/HG02666.mapped.ILLUMINA.bwa.GWD.low_coverage.20121211.bam" 11:5246000-5250000 \
  | samtools sort -n - \
  | samtools fastq - > reads.fq
```

## Variant called by the pipeline

```
HBB:2232  T → A   GT=0/1   QUAL=105   DP=7   AD=2,5
```

That's chr11:5,248,232 (GRCh37) — **rs334**, the sickle cell mutation:

| | |
| --- | --- |
| genomic | chr11:5,248,232 T>A |
| HGVS-c  | HBB:c.20A>T (HBB is on the reverse strand) |
| HGVS-p  | p.Glu7Val |
| zygosity | 0/1 — heterozygous (HbAS, sickle cell trait) |
| ClinVar | VCV000015333 — pathogenic |

## Local use

If you don't want Colab, clone and run on any Linux/macOS box that has `bwa`, `samtools`, `bcftools`, `fastp`, `fastqc`, and Jupyter installed:

```bash
git clone https://github.com/hanyingjhuang/ngs-sickle-tutorial
cd ngs-sickle-tutorial
jupyter notebook notebook.ipynb
```

## License

Code: MIT.  
Data: redistributed under the [1000 Genomes data use policy](https://www.internationalgenome.org/data-portal/data-collection) (open use).
