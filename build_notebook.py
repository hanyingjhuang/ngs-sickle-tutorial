"""Build the NGS teaching notebook (Colab-runnable)."""
import json

REPO = "hanyingjhuang/ngs-sickle-tutorial"
RAW = f"https://raw.githubusercontent.com/{REPO}/main"

def md(*lines):
    src = []
    for i, l in enumerate(lines):
        src.append(l + ("\n" if i < len(lines) - 1 else ""))
    return {"cell_type": "markdown", "metadata": {}, "source": src}

def code(*lines):
    src = []
    for i, l in enumerate(lines):
        src.append(l + ("\n" if i < len(lines) - 1 else ""))
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": src}

cells = []

# ===================================================================
#  COVER
# ===================================================================
cells.append(md(
    "# NGS pipeline — from raw reads to a clinical variant",
    "",
    "### Teaching notebook · 生技產業研發 III (Biotech Industry R&D III)",
    "**Taipei Medical University**  ",
    "Author: **Han-Ying Jhuang**, PhD ",
    "Contact: hanyingjhuang@tmu.edu.tw  ",
    "",
    "---",
    "",
    "## What this notebook teaches",
    "",
    "By the end of this notebook you will have run, with your own hands, the entire short-read variant-calling pipeline that hospitals and research labs use every day. You will start with the unprocessed output of an Illumina sequencer (a `FASTQ` file of millions of short DNA reads) and finish by **finding the sickle cell mutation** in a real person's genome.",
    "",
    "We use real public data: 438 sequencing reads taken from sample **HG02666**, a Gambian individual in the [1000 Genomes Project](https://www.internationalgenome.org/), confirmed in the project's official variant calls to be a heterozygous carrier of the sickle cell allele (**rs334**, *HBB* c.20A>T, p.Glu7Val). Our reference is a 4-kilobase slice of human chromosome 11 (GRCh37) covering the *HBB* gene, fetched from the UCSC Genome Browser.",
    "",
    "## Pipeline overview",
    "",
    "```",
    "  raw reads          QC + trim          alignment         sorted BAM",
    "    FASTQ    ─►   FastQC + fastp   ─►   BWA-MEM    ─►    samtools",
    "                                                              │",
    "                                                              ▼",
    "                                  variant call        pileup of bases",
    "                                  bcftools         ◄─    samtools mpileup",
    "                                       │",
    "                                       ▼",
    "                                  rs334 / HbS",
    "                              clinical interpretation",
    "```",
    "",
    "## How to use this notebook",
    "",
    "1. Click **Runtime → Run all** to run every step automatically, *or*",
    "2. Click each cell and press **Shift + Enter** to step through one at a time.",
    "",
    "Every code cell runs *real Linux binaries* (`bwa`, `samtools`, `bcftools`, …) on a free Google cloud VM — exactly the same software you would install on a lab workstation. Nothing here is simulated.",
    "",
    "## Learning objectives",
    "",
    "After completing this notebook you will be able to:",
    "1. Read a FASTQ file and explain its 4-line structure.",
    "2. Run quality control and adapter trimming, and tell when a dataset needs more or less of each.",
    "3. Explain what an aligner does and read a SAM/BAM record.",
    "4. Generate and interpret a pileup at any position in a genome.",
    "5. Call SNVs and indels with `bcftools` and read the resulting VCF.",
    "6. Trace a single nucleotide change to its codon, amino acid, and phenotype.",
    "",
    "---"
))

# ===================================================================
# STEP 0a: install
# ===================================================================
cells.append(md(
    "## 0 · Setting up the environment",
    "",
    "Colab gives us a brand-new Linux machine each time. The first thing we have to do is install the bioinformatics tools we'll need. These are exactly the same command-line programs used in every NGS lab in the world:",
    "",
    "| tool | what it does |",
    "| --- | --- |",
    "| **`bwa`** | Burrows-Wheeler Aligner — maps short reads to a reference genome |",
    "| **`samtools`** | sort, index, view, and pile-up BAM/SAM files |",
    "| **`bcftools`** | call and filter variants from a pile-up |",
    "| **`fastqc`** | classical quality-control report for FASTQ |",
    "| **`fastp`** | fast read trimmer + QC report (a modern alternative to Trimmomatic) |",
    "",
    "The cell below uses `apt-get`, Ubuntu's package manager, to install all five at once."
))

cells.append(code(
    "%%bash",
    "# Install all bioinformatics tools we'll need (≈30 seconds)",
    "apt-get -qq update",
    "apt-get -qq install -y bwa samtools bcftools fastqc fastp 2>&1 | tail -3",
    "",
    "echo",
    "echo \"Installed versions:\"",
    "bwa 2>&1 | grep '^Version'",
    "samtools --version | head -1",
    "bcftools --version | head -1",
    "fastp    --version 2>&1 | head -1",
    "fastqc   --version"
))

cells.append(md(
    "✅ If you see version numbers above, the tools are ready. Notice the version numbers — in real research you record these so others can reproduce your work."
))

# ===================================================================
# STEP 0b: data
# ===================================================================
cells.append(md(
    "## 0b · Downloading the real data",
    "",
    "We need two files to do anything useful:",
    "",
    "1. **A reference genome** (`ref.fa` — a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file). FASTA is the simplest possible sequence format: one `>name` header line followed by the sequence on the lines below. Our reference is a small 4 kb slice of human chromosome 11 covering the **HBB gene** (β-globin), fetched from the UCSC Genome Browser.",
    "",
    "2. **The sequencing reads** (`reads.fq` — a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file). FASTQ is what comes out of an Illumina sequencer: each read is 4 lines (name, sequence, separator, quality). Our reads come from sample **HG02666** in the 1000 Genomes Project — a Gambian individual whose genome was sequenced with low-coverage Illumina chemistry.",
    "",
    "Both files live in the GitHub repo for this lesson. We pull them with `wget`."
))

cells.append(code(
    "%%bash",
    "set -e                                                  # stop on any error",
    "mkdir -p data && cd data",
    f"curl -sLO {RAW}/data/ref.fa                            # 4 kB reference",
    f"curl -sLO {RAW}/data/reads.fq                          # 97 kB of reads",
    "",
    "ls -la                                                  # confirm files arrived",
    "echo",
    "echo \"--- reference ---\"",
    "head -1 ref.fa                                          # FASTA header line",
    "echo \"length: $(grep -v '^>' ref.fa | tr -d '\\n' | wc -c) bp\"",
    "echo",
    "echo \"--- reads ---\"",
    "echo \"read count: $(( $(wc -l < reads.fq) / 4 ))\"      # 4 lines per read"
))

cells.append(md(
    "We now have a 4,000 bp reference and 438 short reads. In a real experiment the FASTQ would be 30–100 GB — millions of reads — but the pipeline is identical. Working with a small subset just makes it run in seconds instead of hours."
))

# ===================================================================
# STEP 1: inspect reads
# ===================================================================
cells.append(md(
    "---",
    "## Step 1 · What does a sequencing read actually look like?",
    "",
    "Before we touch the pipeline, let's just *look* at the data. Each read in a FASTQ file is exactly **4 lines**:",
    "",
    "| line | content | example |",
    "| --- | --- | --- |",
    "| 1 | `@` + read name | `@SRR582169.1324745/1` |",
    "| 2 | the DNA sequence (A, C, G, T, sometimes N) | `TGGCTCTGCCCT…` |",
    "| 3 | `+` (separator) | `+` |",
    "| 4 | the quality string — one character per base | `CEFEHFHHGFGH…` |",
    "",
    "The quality string is the most surprising part. Each character represents the *Phred score* — how confident the sequencer is in that base. Letters early in the alphabet mean low confidence, later letters mean high. Specifically: `Phred = ord(char) - 33`. A character of `'I'` (ASCII 73) means quality 40 (a 1 in 10,000 chance the base is wrong). A character of `'#'` (ASCII 35) means quality 2 — basically a guess.",
    "",
    "Let's look at the first two reads:"
))

cells.append(code(
    "!head -8 data/reads.fq                  # first 2 reads = 8 lines"
))

cells.append(md(
    "Notice the quality strings end with garbage characters like `;` and `<` and `=` — these are low-quality bases at the 3′ end of each read. This is a universal pattern in Illumina sequencing: the chemistry degrades over the length of a read, so the last 10–20 bases are noisier than the first.",
    "",
    "Below is a **visual** version of the same reads. Each base gets a colour (A=green, T=red, C=blue, G=yellow), and a coloured bar height-encodes the quality below it (taller = higher quality)."
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "",
    "with open('data/reads.fq') as f:",
    "    lines = [next(f).rstrip() for _ in range(8)]",
    "",
    "for r in range(2):",
    "    name, seq, _, qual = lines[r*4 : r*4 + 4]",
    "    Q = [ord(c) - 33 for c in qual]",
    "",
    "    fig, ax = plt.subplots(figsize=(14, 1.5))",
    "    cmap = {'A':'#1a7f37', 'T':'#b1272d', 'C':'#0969da', 'G':'#9a6700', 'N':'#777'}",
    "    # Letter row at top",
    "    for i, b in enumerate(seq):",
    "        ax.text(i, 1.4, b, ha='center', va='center',",
    "                color=cmap.get(b,'#000'), family='monospace', fontsize=9)",
    "    # Quality bars below",
    "    bar_colors = ['#1a7f37' if q>=30 else '#9a6700' if q>=20",
    "                  else '#bf6b00' if q>=10 else '#b1272d' for q in Q]",
    "    ax.bar(range(len(Q)), Q, color=bar_colors, width=0.9)",
    "    ax.set_xlim(-0.5, len(seq)-0.5); ax.set_ylim(0, 42)",
    "    ax.set_yticks([0, 20, 30, 40]); ax.set_xticks([])",
    "    ax.set_ylabel('Phred Q')",
    "    ax.set_title(name, fontsize=9, loc='left')",
    "    for s in ['top','right']: ax.spines[s].set_visible(False)",
    "    plt.tight_layout(); plt.show()"
))

cells.append(md(
    "👀 **Look at the right side of each plot.** The bars get shorter and turn yellow / orange / red — those are the low-quality 3′ ends. We will trim them off in Step 3."
))

# ===================================================================
# STEP 2: FastQC
# ===================================================================
cells.append(md(
    "---",
    "## Step 2 · Quality control with FastQC",
    "",
    "Looking at one read tells you nothing about the *whole* dataset. **FastQC** is the standard first-pass QC tool. It reads the entire FASTQ, computes summary statistics (per-base quality, GC content, adapter contamination, sequence duplication, …), and writes an HTML report you can read in a browser.",
    "",
    "FastQC has been the de-facto standard since 2010. Its output is the very first thing every bioinformatician looks at when handed a new dataset."
))

cells.append(code(
    "%%bash",
    "mkdir -p qc",
    "fastqc data/reads.fq -o qc 2>&1 | tail -3        # runs in seconds on small data",
    "ls qc/                                            # what did it produce?"
))

cells.append(md(
    "FastQC produced two files: a `.html` report (the human-readable one) and a `.zip` (raw data, for parsing in a script). Let's display the HTML report inline below:"
))

cells.append(code(
    "from IPython.display import IFrame",
    "IFrame('qc/reads_fastqc.html', width='100%', height=600)"
))

cells.append(md(
    "**How to read this report.** Each section has a tick (✓ pass), warning (⚠), or cross (✗). Don't panic if you see warnings — many are normal for small datasets. Pay attention to:",
    "",
    "- **Per base sequence quality** — should stay green / above Q20. Drops at the end are normal.",
    "- **Per base sequence content** — A/C/G/T should be ~25% each in a random region (it's *not* random for a small targeted region like ours, so warnings are expected).",
    "- **Adapter content** — should be flat at zero. If it climbs at the read ends, you have adapter contamination and need to trim.",
    "",
    "*Discussion question for class:* what biases would you expect from sequencing only the HBB region (4 kb of human DNA), and why does FastQC flag some of them?"
))

# ===================================================================
# STEP 3: trim
# ===================================================================
cells.append(md(
    "---",
    "## Step 3 · Trimming with `fastp`",
    "",
    "FastQC just *reports* problems — it doesn't fix them. To clean up the reads we use a separate tool. Two are widely used:",
    "",
    "- **Trimmomatic** (Java, classical, very configurable)",
    "- **fastp** (C++, fast, modern, also generates its own QC report)",
    "",
    "We use `fastp` here because it's faster and self-documenting. Our settings:",
    "",
    "- `-i` / `-o` — input FASTQ, output trimmed FASTQ",
    "- `-j` / `-h` — write a JSON and HTML report",
    "- `--cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20` — slide a 4-base window from the 3′ end and trim until the average quality in the window stays at Q20 or above. (Sliding-window quality trimming is the standard approach.)"
))

cells.append(code(
    "%%bash",
    "fastp \\",
    "  -i data/reads.fq \\",
    "  -o data/trimmed.fq \\",
    "  -j data/qc.json -h data/qc.html \\",
    "  --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \\",
    "  2>&1 | tail -25                       # show the last lines of fastp's log"
))

cells.append(md(
    "Read the log carefully — it tells you exactly what fastp did:",
    "",
    "- *reads passed filter* / *reads failed*: how many reads survived",
    "- *Q20 / Q30 rate*: percentage of high-quality bases, before vs after",
    "- *adapter trimming*: was any adapter found and removed?",
    "",
    "Now let's plot the per-base quality before and after trimming, using the JSON fastp wrote:"
))

cells.append(code(
    "import json, matplotlib.pyplot as plt",
    "qc = json.load(open('data/qc.json'))",
    "before = qc['read1_before_filtering']['quality_curves']['mean']",
    "after  = qc['read1_after_filtering']['quality_curves']['mean']",
    "",
    "fig, ax = plt.subplots(figsize=(10, 3.5))",
    "# coloured quality bands as the backdrop",
    "ax.axhspan(28, 40, color='#1a7f37', alpha=0.07, label='_high')",
    "ax.axhspan(20, 28, color='#9a6700', alpha=0.10, label='_med')",
    "ax.axhspan(0,  20, color='#b1272d', alpha=0.07, label='_low')",
    "ax.plot(before, '--', color='#b1272d', lw=1.6, label='before trim')",
    "ax.plot(after,  '-',  color='#1a7f37', lw=1.6, label='after trim')",
    "ax.set_xlabel('position along read'); ax.set_ylabel('mean Phred Q')",
    "ax.set_ylim(0, 41); ax.legend(loc='lower left', frameon=False)",
    "for s in ['top','right']: ax.spines[s].set_visible(False)",
    "ax.set_title('Per-base mean quality — before vs. after fastp')",
    "plt.tight_layout(); plt.show()",
    "",
    "# Print a side-by-side stats table",
    "b, a = qc['summary']['before_filtering'], qc['summary']['after_filtering']",
    "print(f\"{'metric':<12} {'before':>14} {'after':>14}\")",
    "print('-' * 42)",
    "print(f\"{'reads':<12} {b['total_reads']:>14,} {a['total_reads']:>14,}\")",
    "print(f\"{'bases':<12} {b['total_bases']:>14,} {a['total_bases']:>14,}\")",
    "print(f\"{'Q30 rate':<12} {b['q30_rate']*100:>13.1f}% {a['q30_rate']*100:>13.1f}%\")",
    "print(f\"{'GC%':<12} {b['gc_content']*100:>13.1f}% {a['gc_content']*100:>13.1f}%\")"
))

cells.append(md(
    "📈 The **green line (after) sits above the red dashed line (before)** at the 3′ end — that's the trimming working. We've thrown away noisy bases. The total number of bases drops, but the *quality* of what remains is higher. Q30 rate climbing means fewer low-quality bases pollute downstream steps."
))

# ===================================================================
# STEP 4: align
# ===================================================================
cells.append(md(
    "---",
    "## Step 4 · Alignment with BWA-MEM",
    "",
    "Each read is a short, anonymous string of A/C/G/T. We don't yet know *where* in the genome it came from. **Alignment** answers that question: for every read, find its best-matching position on the reference.",
    "",
    "**BWA-MEM** is the workhorse aligner for short reads (Burrows-Wheeler Aligner, MEM = Maximal Exact Match algorithm, published by Heng Li in 2009 and updated since). It's used by the 1000 Genomes Project, GATK, and most clinical labs. The MEM algorithm anchors each read on a long stretch that matches the reference exactly, then extends through any mismatches and small indels.",
    "",
    "Two commands:",
    "",
    "1. **`bwa index ref.fa`** — build a Burrows-Wheeler search index of the reference. Run once per reference, takes seconds for our 4 kb (would take ~1 hour for the whole human genome).",
    "2. **`bwa mem ref.fa trimmed.fq > aligned.sam`** — align reads to the indexed reference. Output is **SAM** (Sequence Alignment/Map) format."
))

cells.append(code(
    "%%bash",
    "# Build the index — outputs ref.fa.amb, .ann, .bwt, .pac, .sa",
    "bwa index data/ref.fa 2>&1 | tail -3",
    "",
    "# Align — write SAM to file, route bwa's progress log to a separate file",
    "bwa mem -t 2 data/ref.fa data/trimmed.fq 2> data/bwa.log > data/aligned.sam",
    "",
    "echo",
    "echo '--- bwa log (last lines) ---'",
    "tail -8 data/bwa.log",
    "",
    "echo",
    "echo '--- aligned.sam first 4 lines ---'",
    "head -4 data/aligned.sam | cut -c1-180"
))

cells.append(md(
    "**The SAM format**, briefly. Every line that does *not* start with `@` is one read alignment, with these tab-separated columns:",
    "",
    "| col | name | meaning |",
    "| --- | --- | --- |",
    "| 1 | QNAME | read name |",
    "| 2 | FLAG  | bit flags (mapped? reverse strand? secondary?) |",
    "| 3 | RNAME | reference sequence the read mapped to |",
    "| 4 | POS   | 1-based leftmost position on the reference |",
    "| 5 | MAPQ  | mapping quality (0-60, higher = more confident) |",
    "| 6 | CIGAR | match/mismatch/indel pattern (e.g. `100M` = 100 matches) |",
    "| 9 | TLEN  | template length (insert size for paired reads) |",
    "| 10| SEQ   | the read sequence |",
    "| 11| QUAL  | the read quality string |",
    "",
    "Lines starting with `@` are headers — `@SQ` describes each reference contig, `@PG` records the program that wrote the file. (Provenance!)"
))

# ===================================================================
# STEP 5: sort + index
# ===================================================================
cells.append(md(
    "---",
    "## Step 5 · Sort and index the alignments",
    "",
    "BWA outputs reads in the order it processed them — essentially random along the genome. Every downstream tool (variant callers, browsers, depth calculators) expects reads to be sorted by chromosome and position. We also need an *index* so tools can jump to any genomic region without scanning the whole file.",
    "",
    "Three commands:",
    "",
    "1. **`samtools sort`** — sorts SAM by reference position, outputs compressed binary BAM",
    "2. **`samtools index`** — builds a BAM index (`.bai`) for random access",
    "3. **`samtools flagstat`** — quick alignment statistics"
))

cells.append(code(
    "%%bash",
    "samtools sort  data/aligned.sam -o data/aligned.bam",
    "samtools index data/aligned.bam",
    "",
    "echo '--- flagstat output ---'",
    "samtools flagstat data/aligned.bam"
))

cells.append(md(
    "**Reading flagstat:**",
    "",
    "- *in total* — total reads (including unmapped)",
    "- *mapped* — reads BWA placed somewhere on the reference. Should be ≥95% in a clean run.",
    "- *paired in sequencing* — for paired-end. Our data is single-end so this is 0.",
    "- *duplicates* — PCR duplicates flagged by `MarkDuplicates`. We haven't run that step.",
    "",
    "Now the satisfying visualization: **per-position depth across the reference**. Every position gets one bar — the number of reads covering it. A red dashed line marks the rs334 (sickle cell) site so you can see the data we're going to call variants from."
))

cells.append(code(
    "import subprocess, matplotlib.pyplot as plt",
    "",
    "raw = subprocess.check_output(['samtools', 'depth', '-a', 'data/aligned.bam'], text=True)",
    "rows = [l.split('\\t') for l in raw.strip().split('\\n') if l]",
    "pos = [int(r[1]) for r in rows]",
    "dep = [int(r[2]) for r in rows]",
    "RS334 = 5248232 - 5246000          # in our subset's 1-based coords this is 2232",
    "",
    "fig, ax = plt.subplots(figsize=(11, 3))",
    "ax.fill_between(pos, dep, color='#0a7c7e', alpha=0.4)",
    "ax.plot(pos, dep, color='#0a7c7e', lw=0.9)",
    "ax.axvline(RS334, color='#b1272d', lw=1, ls='--', label='rs334 (HbS)')",
    "ax.set_xlabel('position on HBB region (chr11:5,246,001-5,250,000, GRCh37)')",
    "ax.set_ylabel('read depth')",
    "ax.legend(loc='upper right', frameon=False)",
    "for s in ['top','right']: ax.spines[s].set_visible(False)",
    "ax.set_title('Coverage across the HBB region')",
    "plt.tight_layout(); plt.show()",
    "",
    "print(f'mean depth: {sum(dep)/len(dep):.1f}x')",
    "print(f'positions with ≥1 read: {sum(1 for d in dep if d):,} / {len(dep):,}')",
    "print(f'positions with 0 reads: {sum(1 for d in dep if d==0):,}')"
))

cells.append(md(
    "📊 The coverage is uneven because this is *low-coverage* data (the 1000 Genomes phase 1 strategy: many people, lightly sequenced). High-coverage clinical sequencing would give you a flatter ~30× curve. **Average ~3-5×** is enough to call common variants on a heterozygote — barely. For homozygous calls you ideally want ≥10×."
))

# ===================================================================
# STEP 6: pileup
# ===================================================================
cells.append(md(
    "---",
    "## Step 6 · Pile up reads at every position",
    "",
    "Now we *transpose* the data. Instead of looking at reads, we look at the genome **column by column**: at each reference position, what bases do the reads observe?",
    "",
    "This is what `samtools mpileup` produces. The output looks like:",
    "",
    "```",
    "HBB    2232    T    7    a.aAA.,    DCEHHE,",
    "```",
    "",
    "tab-separated columns: contig, position, **REF base**, depth, **read bases**, **read qualities**.",
    "",
    "The *read bases* string uses a special encoding:",
    "",
    "- `.` = match REF on the **forward** strand",
    "- `,` = match REF on the **reverse** strand",
    "- `A C G T` = mismatch on the forward strand",
    "- `a c g t` = mismatch on the reverse strand",
    "- `^]` = start of a read (followed by mapping quality char)",
    "- `$` = end of a read",
    "- `+N…` / `-N…` = insertion / deletion",
    "",
    "So the example above means: at HBB position 2232, the reference is T, depth is 7. The read bases are `a.aAA.,` — three reads see `A` (in different orientations), four reads see the reference T. **That's a heterozygous variant.**"
))

cells.append(code(
    "%%bash",
    "samtools mpileup -f data/ref.fa data/aligned.bam 2>/dev/null > data/pileup.txt",
    "echo \"total pileup rows: $(wc -l < data/pileup.txt)\"",
    "echo",
    "echo '--- the rows around the sickle position (HBB:2232) ---'",
    "awk '$2 >= 2228 && $2 <= 2236' data/pileup.txt"
))

cells.append(md(
    "**Look at line 2232.** Compare the read-base string to the rules above. How many reads see `A` vs `T`? Is this clearly a variant or noise?",
    "",
    "Below is a **visual pileup** — each column is one reference position, each row in the stack is one read. Bases matching REF are faded, mismatches are highlighted. The red dashed box marks the sickle site."
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "",
    "rows = []",
    "with open('data/pileup.txt') as f:",
    "    for line in f:",
    "        c, p, ref, dp, bases, qual = line.rstrip().split('\\t')",
    "        if 2218 <= int(p) <= 2246:           # window around the sickle site",
    "            rows.append((int(p), ref.upper(), int(dp), bases))",
    "",
    "def parse_bases(refb, s):",
    "    \"\"\"Convert mpileup base string to one base per supporting read.\"\"\"",
    "    out, i = [], 0",
    "    while i < len(s):",
    "        ch = s[i]",
    "        if   ch == '^':              i += 2;  continue   # skip mapq char",
    "        elif ch == '$':              i += 1;  continue",
    "        elif ch in '+-':             # insertion/deletion: skip the bases",
    "            j = i + 1",
    "            while j < len(s) and s[j].isdigit(): j += 1",
    "            n = int(s[i+1:j]); i = j + n; continue",
    "        elif ch in '.,':             out.append(refb)",
    "        elif ch.upper() in 'ACGTN':  out.append(ch.upper())",
    "        elif ch == '*':              out.append('-')",
    "        i += 1",
    "    return out",
    "",
    "cmap = {'A':'#1a7f37','T':'#b1272d','C':'#0969da','G':'#9a6700','N':'#777','-':'#777'}",
    "max_stack = max(r[2] for r in rows)",
    "",
    "fig, ax = plt.subplots(figsize=(13, max_stack*0.22 + 1.5))",
    "for i, (p, ref, dp, bases) in enumerate(rows):",
    "    # REF letter on top",
    "    ax.text(i, max_stack+1, ref, ha='center', va='center',",
    "            color=cmap[ref], family='monospace', fontsize=11, weight='bold')",
    "    # stack of read bases",
    "    for k, b in enumerate(parse_bases(ref, bases)):",
    "        match = (b == ref)",
    "        ax.add_patch(Rectangle((i-0.4, max_stack-1-k), 0.8, 0.85,",
    "                               color=cmap.get(b,'#777'), alpha=0.18 if match else 0.92))",
    "        if not match:",
    "            ax.text(i, max_stack-1-k+0.42, b, ha='center', va='center',",
    "                    color='white', family='monospace', fontsize=8)",
    "",
    "# Highlight the sickle column",
    "if any(p == 2232 for p, *_ in rows):",
    "    RS = next(i for i,(p,*_) in enumerate(rows) if p == 2232)",
    "    ax.add_patch(Rectangle((RS-0.5, -0.5), 1, max_stack+2.2,",
    "                           fill=False, ec='#b1272d', lw=1.5, ls='--'))",
    "",
    "ax.set_xticks(range(len(rows)))",
    "ax.set_xticklabels([str(r[0]) for r in rows], rotation=90, fontsize=8)",
    "ax.set_yticks([]); ax.set_xlim(-1, len(rows)); ax.set_ylim(-1, max_stack+2.5)",
    "ax.set_title('Pileup window — REF row at top, reads stacked below; '",
    "             + 'red box = rs334 (sickle position)', fontsize=10)",
    "for s in ax.spines.values(): s.set_visible(False)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "🔴 **Inside the red box**, at HBB position 2232, you should see a column where roughly half the reads show **`A`** (red) instead of the REF **T**. That's a classic heterozygous SNV — and it's the sickle cell mutation.",
    "",
    "*Eye-balling pileups is exactly how variant callers work, just systematically.*"
))

# ===================================================================
# STEP 7: variant call
# ===================================================================
cells.append(md(
    "---",
    "## Step 7 · Calling variants → VCF",
    "",
    "We've seen *by eye* that there's a variant at HBB:2232. Now let's have a real variant caller find it for us.",
    "",
    "**`bcftools`** does this in two steps:",
    "",
    "1. `bcftools mpileup` — build a BCF (binary VCF) of base likelihoods at every position",
    "2. `bcftools call -mv` — call variants from the likelihoods. `-m` is the multiallelic caller, `-v` means only output positions that are variant.",
    "",
    "The output is a **VCF** (Variant Call Format) file — the universal format for variants. Every line after the header is one variant: chrom, position, REF, ALT, quality, filter, info fields, and per-sample genotype."
))

cells.append(code(
    "%%bash",
    "# Step 7a: pile up base likelihoods → BCF",
    "bcftools mpileup -f data/ref.fa data/aligned.bam -Ou -o data/pile.bcf 2>&1 | tail -2",
    "",
    "# Step 7b: call variants from the likelihoods → compressed VCF",
    "bcftools call -mv -Oz -o data/variants.vcf.gz data/pile.bcf 2>&1 | tail -2",
    "bcftools index -f data/variants.vcf.gz",
    "",
    "echo",
    "echo '--- variants.vcf.gz (skipping ## meta lines) ---'",
    "bcftools view data/variants.vcf.gz 2>/dev/null | grep -v '^##' | head -25"
))

cells.append(md(
    "**The VCF header** (line starting with `#CHROM`) names the columns; **each line below** is one variant call.",
    "",
    "Key fields:",
    "- **POS** — 1-based position",
    "- **REF / ALT** — reference and alternate alleles",
    "- **QUAL** — Phred-scaled call confidence (higher = better)",
    "- **INFO** — `key=value` pairs: `DP` = depth, `MQ` = mean mapping quality, `DP4` = ref-fwd, ref-rev, alt-fwd, alt-rev counts",
    "- **FORMAT** — defines per-sample fields below",
    "- **last column** — per-sample data: `GT` = genotype (`0/0`, `0/1`, `1/1`), `AD` = allele depths",
    "",
    "*Practice question for class:* find a homozygous variant call (genotype `1/1`) in the table above. Why does its INFO field show `DP4=0,0,…` for the reference alleles?"
))

# ===================================================================
# STEP 8: interpret
# ===================================================================
cells.append(md(
    "---",
    "## Step 8 · Interpreting the variant — finding sickle cell",
    "",
    "We started with raw FASTQ and ended with a list of variants. **Now comes the biology.**",
    "",
    "We're specifically looking for **rs334** — the variant that causes sickle cell anaemia. Background:",
    "",
    "- **Position**: chromosome 11, base 5,248,232 (GRCh37 coordinates)",
    "- **In our reference subset** (which starts at 5,246,001), this is HBB position 2,232",
    "- **The change**: the reference T → A on the genomic forward strand",
    "- **The biology**: HBB is on the *reverse* strand of chr11, so on the gene's coding sequence (the cDNA) this is the reverse complement: A → T. In standard HGVS notation we write **`HBB:c.20A>T`**.",
    "- **The protein**: c.20 falls inside codon 7 (`GAG`, glutamic acid). The mutation changes it to `GTG` (valine). HGVS protein notation: **`p.Glu7Val`**.",
    "- **The phenotype**: substituting a polar Glu with a hydrophobic Val on the β-globin surface lets HbS molecules polymerize when oxygen is low. Red blood cells deform into the characteristic sickle shape. Heterozygotes have *sickle cell trait* (HbAS, usually asymptomatic, malaria-protective). Homozygotes have *sickle cell disease* (HbSS, chronic anaemia + vaso-occlusive crises).",
    "",
    "Let's pull rs334 out of our VCF programmatically and translate it:"
))

cells.append(code(
    "import subprocess",
    "",
    "vcf = subprocess.check_output(['bcftools', 'view', 'data/variants.vcf.gz'], text=True)",
    "rows = [l.split('\\t') for l in vcf.splitlines() if l and not l.startswith('#')]",
    "rs334 = next((r for r in rows if int(r[1]) == 2232), None)",
    "",
    "if rs334 is None:",
    "    print('rs334 was not called in this run — try re-running, or '",
    "          + 'check coverage at this position.')",
    "else:",
    "    fmt  = dict(zip(rs334[8].split(':'), rs334[9].split(':')))",
    "    info = dict(kv.split('=', 1) for kv in rs334[7].split(';') if '=' in kv)",
    "    dp   = info.get('DP', '?')",
    "    ad   = fmt.get('AD', '?,?').split(',')",
    "",
    "    print('=' * 64)",
    "    print('  rs334 (sickle cell allele) called from real HG02666 reads')",
    "    print('=' * 64)",
    "    print(f'  position    HBB:{rs334[1]}  =  chr11:5,248,232  (GRCh37)')",
    "    print(f'  REF -> ALT  {rs334[3]} -> {rs334[4]}      (genomic forward strand)')",
    "    print( '  cDNA        c.20 A>T   (HBB on reverse strand; revcomp of T>A)')",
    "    print( '  protein     codon 7  GAG (Glu) -> GTG (Val)  =  p.Glu7Val')",
    "    print(f'  genotype    {fmt[\"GT\"]}     '",
    "          + '(0/1 = heterozygous = sickle cell trait, HbAS)')",
    "    print(f'  depth       {dp}    (REF reads = {ad[0]}, ALT reads = {ad[1]})')",
    "    print(f'  QUAL        {rs334[5]}')",
    "    print()",
    "    print('  → Phenotype: HbAS (sickle cell trait)')",
    "    print('    • Heterozygous carrier; usually asymptomatic.')",
    "    print('    • Confers partial protection against falciparum malaria.')",
    "    print('    • Children of two HbAS carriers have a 1-in-4 risk of HbSS.')",
    "    print()",
    "    print('  → Cross-references:')",
    "    print('    ClinVar VCV000015333 (pathogenic).  OMIM 603903.  dbSNP rs334.')"
))

cells.append(md(
    "🎉 **Congratulations.** You have just walked an Illumina FASTQ from a real human all the way to a clinically meaningful variant call, using exactly the same software a hospital genetics lab would use. The pipeline you just ran is the foundation of clinical genome sequencing, oncology panel sequencing, prenatal carrier screening, and large-scale studies like the 1000 Genomes Project itself."
))

# ===================================================================
# CLOSING
# ===================================================================
cells.append(md(
    "---",
    "## What we just did, in one table",
    "",
    "| step | tool | input → output |",
    "| --- | --- | --- |",
    "| inspect | `head` | reads.fq |",
    "| QC | `fastqc` | qc/reads_fastqc.html |",
    "| trim | `fastp` | reads.fq → trimmed.fq + qc.json |",
    "| align | `bwa index`, `bwa mem` | trimmed.fq + ref.fa → aligned.sam |",
    "| sort/index | `samtools sort/index/flagstat` | aligned.sam → aligned.bam + .bai |",
    "| coverage | `samtools depth` | aligned.bam → depth plot |",
    "| pileup | `samtools mpileup` | aligned.bam → pileup.txt |",
    "| call | `bcftools mpileup`, `bcftools call` | aligned.bam → variants.vcf.gz |",
    "| interpret | parsing in Python | variants.vcf.gz → rs334 → HbS |",
    "",
    "## Try it yourself",
    "",
    "1. **Re-run with a non-carrier sample.** Pick another 1000 Genomes sample (e.g. NA12878, European, almost certainly does not carry rs334) and edit the `wget` URL in the data-download cell. The pipeline will run unchanged but the VCF will not contain a row at HBB:2232.",
    "2. **Try a homozygote.** A few 1000 Genomes samples are HbSS (homozygous for the sickle allele). Their variant call at HBB:2232 will show genotype `1/1` instead of `0/1`, and their `DP4` field will have zero REF reads.",
    "3. **Vary the trimming.** Set `--cut_tail_mean_quality` to 30 (more aggressive) or remove the flag entirely. How does this affect the final variant calls? Could over-trimming *miss* a variant?",
    "4. **Add VEP / SnpEff annotation.** Real clinical pipelines run a final annotator that maps every variant to its consequence (synonymous, missense, splice, …) and looks up its frequency in gnomAD and its clinical record in ClinVar. This is the next step beyond what we did here.",
    "",
    "## Further reading",
    "",
    "- **The 1000 Genomes Project Consortium**, *A global reference for human genetic variation*, Nature 526, 68–74 (2015).",
    "- **Heng Li**, *Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM*, arXiv:1303.3997 (2013).",
    "- **Petr Danecek et al.**, *Twelve years of SAMtools and BCFtools*, GigaScience 10, giab008 (2021).",
    "- **HBB / sickle cell** at OMIM 603903, ClinVar VCV000015333, and dbSNP rs334.",
    "",
    "---",
    "",
    "*Notebook prepared by **Han-Ying Jhuang, PhD** for **生技產業研發 III · Biotech Industry R&D III**, Taipei Medical University. Released under MIT for the code and the [1000 Genomes data use policy](https://www.internationalgenome.org/data-portal/data-collection) for the data. Questions or improvements: [hanyingjhuang@tmu.edu.tw](mailto:hanyingjhuang@tmu.edu.tw).*"
))

# Build & save
nb = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "name": "ngs_sickle_pipeline.ipynb"},
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.10"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}
with open("/tmp/ngs-repo/notebook.ipynb", "w") as f:
    json.dump(nb, f, indent=1)

print(f"wrote {len(cells)} cells")
n_md   = sum(1 for c in cells if c['cell_type']=='markdown')
n_code = sum(1 for c in cells if c['cell_type']=='code')
print(f"  markdown cells: {n_md}")
print(f"  code cells:     {n_code}")
