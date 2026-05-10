"""Build the NGS teaching notebook (Colab-runnable).

Design notes:
- Code cells are source-hidden by default (Colab honors jupyter.source_hidden).
  Students see the result; click to expand the code.
- Each step opens with a short schematic illustration (small matplotlib figure)
  so the concept lands before the code runs.
- FastQC HTML is embedded via <iframe srcdoc> so it renders inside Colab,
  plus a button opens it in a new tab as a fallback.
"""
import json

REPO = "hanyingjhuang/ngs-sickle-tutorial"
RAW = f"https://raw.githubusercontent.com/{REPO}/main"

def md(*lines):
    src = [l + ("\n" if i < len(lines) - 1 else "") for i, l in enumerate(lines)]
    return {"cell_type": "markdown", "metadata": {}, "source": src}

def code(*lines, hide=True):
    """Code cell. Hidden source by default — Colab shows just the title /
    output unless the student clicks to expand."""
    src = [l + ("\n" if i < len(lines) - 1 else "") for i, l in enumerate(lines)]
    cell = {
        "cell_type": "code",
        "metadata": {"jupyter": {"source_hidden": True}} if hide else {},
        "execution_count": None,
        "outputs": [],
        "source": src,
    }
    return cell

cells = []

# ===================================================================
#  COVER
# ===================================================================
cells.append(md(
    "# NGS pipeline — from raw reads to a clinical variant",
    "",
    "**Author** · Han-Ying Jhuang (莊漢英), PhD — Assistant Professor, Taipei Medical University  ",
    "**Contact** · hanyingjhuang@tmu.edu.tw",
    "",
    "---",
    "",
    "## What this notebook teaches",
    "",
    "By the end of this notebook you will have run the entire short-read variant-calling pipeline that hospitals and research labs use every day. You will start with the unprocessed output of an Illumina sequencer (a `FASTQ` file of millions of short DNA reads) and finish by **finding the sickle cell mutation** in a real person's genome.",
    "",
    "We use real public data: 438 sequencing reads from sample **HG02666**, a Gambian individual in the [1000 Genomes Project](https://www.internationalgenome.org/), confirmed in the project's official VCF to be a heterozygous carrier of the sickle cell allele (**rs334**, *HBB* c.20A>T, p.Glu7Val). Our reference is a 4-kilobase slice of human chromosome 11 (GRCh37) covering the *HBB* gene, fetched from the UCSC Genome Browser.",
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
    "## How to use",
    "",
    "1. Click **Runtime → Run all** to run every step automatically, *or*",
    "2. Click each cell and press **Shift + Enter** to step through one at a time.",
    "",
    "## Learning objectives",
    "",
    "1. Read a FASTQ file and explain its 4-line structure.",
    "2. Run quality control and adapter trimming, and tell when a dataset needs more or less.",
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
    "Colab gives us a fresh Linux machine each session. We install the bioinformatics tools with `apt-get`. These are the same command-line programs used in every NGS lab:",
    "",
    "| tool | role |",
    "| --- | --- |",
    "| **`bwa`** | Burrows-Wheeler Aligner — maps short reads to a reference |",
    "| **`samtools`** | sort, index, view, pile-up BAM/SAM |",
    "| **`bcftools`** | call and filter variants |",
    "| **`fastqc`** | classical QC report for FASTQ |",
    "| **`fastp`** | fast trimmer + QC report (modern alternative to Trimmomatic) |"
))

cells.append(code(
    "%%bash",
    "apt-get -qq update",
    "apt-get -qq install -y bwa samtools bcftools fastqc fastp 2>&1 | tail -3",
    "echo",
    "echo \"--- versions installed ---\"",
    "bwa 2>&1 | grep '^Version'",
    "samtools --version | head -1",
    "bcftools --version | head -1",
    "fastp    --version 2>&1 | head -1",
    "fastqc   --version"
))

cells.append(md(
    "✅ If you see version numbers above, the tools are ready."
))

# ===================================================================
# STEP 0b: data
# ===================================================================
cells.append(md(
    "## 0b · Downloading the real data",
    "",
    "We need two files:",
    "",
    "1. **A reference genome** (`ref.fa`) — a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file: one `>name` header followed by the sequence. Ours is a 4 kb slice of human chromosome 11 covering the *HBB* gene.",
    "2. **The sequencing reads** (`reads.fq`) — a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file: 4 lines per read (name · sequence · `+` · quality). Ours: 438 real Illumina reads from sample HG02666."
))

cells.append(code(
    "%%bash",
    "set -e",
    "mkdir -p data && cd data",
    f"curl -sLO {RAW}/data/ref.fa",
    f"curl -sLO {RAW}/data/reads.fq",
    "ls -la",
    "echo",
    "echo \"--- reference ---\"",
    "head -1 ref.fa",
    "echo \"length: $(grep -v '^>' ref.fa | tr -d '\\n' | wc -c) bp\"",
    "echo",
    "echo \"--- reads ---\"",
    "echo \"read count: $(( $(wc -l < reads.fq) / 4 ))\""
))

cells.append(md(
    "We now have a 4,000 bp reference and 438 short reads. In a real experiment the FASTQ would be 30–100 GB; the pipeline is identical, just slower."
))

# ===================================================================
# STEP 1: inspect reads
# ===================================================================
cells.append(md(
    "---",
    "## Step 1 · What does a sequencing read actually look like?"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import FancyBboxPatch",
    "fig, ax = plt.subplots(figsize=(10, 2.6))",
    "ax.set_xlim(0, 10); ax.set_ylim(0, 5); ax.axis('off')",
    "rows = [",
    "    (4, '@SRR582169.1324745/1', 'line 1 — read name (after \"@\")', '#57606a'),",
    "    (3, 'TGGCTCTGCCCTGACTTTTATGCC...',   'line 2 — DNA sequence',           '#1f2328'),",
    "    (2, '+',                              'line 3 — separator',               '#57606a'),",
    "    (1, 'CEFEHFHHGFGHGHFJHHGFGH...',     'line 4 — Phred quality string',    '#57606a'),",
    "]",
    "for y, txt, role, c in rows:",
    "    ax.text(0.3, y, txt, family='monospace', fontsize=11, color=c, va='center')",
    "    ax.annotate(role, xy=(4.5, y), xytext=(7.0, y), fontsize=10, color='#57606a',",
    "                va='center', arrowprops=dict(arrowstyle='->', color='#bbb', lw=0.7))",
    "ax.add_patch(FancyBboxPatch((0.15, 0.5), 4.3, 4.0, boxstyle='round,pad=0.1',",
    "                            ec='#d0d7de', fc='#f6f8fa', lw=0.8))",
    "ax.set_title('FASTQ — one read = 4 lines', loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "Each FASTQ record is **exactly 4 lines** — name, sequence, `+`, and a quality string with one character per base.",
    "",
    "**The quality string** is the most surprising part. Each character is the *Phred score* — how confident the sequencer is in that base. Specifically: `Phred = ord(char) − 33`. `'I'` (ASCII 73) means Q40 (1 in 10,000 chance the base is wrong). `'#'` (ASCII 35) means Q2 — basically a guess.",
    "",
    "Let's look at the first two reads:"
))

cells.append(code(
    "!head -8 data/reads.fq"
))

cells.append(md(
    "Notice the quality strings end with garbage characters like `;<=` — low-quality bases at the 3′ end. This is universal in Illumina sequencing: chemistry degrades over the length of a read.",
    "",
    "Below is a visual version of the first two reads. The base letter under each bar is colored by base; the bar height is the Phred score (taller = higher quality, green = Q≥30, yellow = Q20-29, orange = Q10-19, red = Q<10)."
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
    "    cmap = {'A':'#1a7f37', 'T':'#b1272d', 'C':'#0969da', 'G':'#9a6700', 'N':'#777'}",
    "    bar_colors = ['#1a7f37' if q>=30 else '#9a6700' if q>=20",
    "                  else '#bf6b00' if q>=10 else '#b1272d' for q in Q]",
    "",
    "    fig, ax = plt.subplots(figsize=(14, 2.0))",
    "    ax.bar(range(len(Q)), Q, color=bar_colors, width=0.95)",
    "    ax.set_xlim(-0.5, len(seq)-0.5); ax.set_ylim(0, 42)",
    "    ax.set_yticks([0, 20, 30, 40])",
    "    ax.set_xticks(range(len(seq)))",
    "    ax.set_xticklabels(list(seq), family='monospace', fontsize=8)",
    "    # Color each x-tick label by base",
    "    for tick, b in zip(ax.get_xticklabels(), seq):",
    "        tick.set_color(cmap.get(b, '#000'))",
    "    ax.tick_params(axis='x', length=0, pad=2)   # tighten label-to-bar gap",
    "    ax.set_ylabel('Phred Q')",
    "    ax.set_title(name, fontsize=9, loc='left')",
    "    for s in ['top', 'right']: ax.spines[s].set_visible(False)",
    "    plt.tight_layout(); plt.show()"
))

cells.append(md(
    "👀 Look at the right side of each plot — bars get shorter and turn yellow / orange / red. Those low-quality 3′ ends will be trimmed in Step 3."
))

# ===================================================================
# STEP 2: FastQC
# ===================================================================
cells.append(md(
    "---",
    "## Step 2 · Quality control with FastQC"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "fig, ax = plt.subplots(figsize=(10, 2.5))",
    "ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis('off')",
    "for i in range(5):",
    "    ax.add_patch(Rectangle((0.5, 3.5 - i*0.5), 1.6, 0.35, fc='#cce5ff', ec='#0969da', lw=0.4))",
    "ax.text(1.3, 0.6, 'all reads\\nin reads.fq', ha='center', fontsize=9, color='#57606a')",
    "ax.annotate('', xy=(3.4, 2.2), xytext=(2.3, 2.2),",
    "            arrowprops=dict(arrowstyle='->', color='#0a7c7e', lw=2))",
    "ax.text(2.85, 2.5, 'FastQC', ha='center', fontsize=10, color='#0a7c7e', weight='bold')",
    "tiles = [('per-base Q', 3.4, '✓', '#1a7f37'),",
    "         ('GC content', 2.5, '✓', '#1a7f37'),",
    "         ('adapter',    1.6, '⚠', '#9a6700'),",
    "         ('duplication',0.7, '✓', '#1a7f37')]",
    "for label, y, mark, col in tiles:",
    "    ax.add_patch(Rectangle((4.2, y), 5.0, 0.7, fc='#f6f8fa', ec='#d0d7de', lw=0.5))",
    "    ax.text(4.4, y+0.35, mark, fontsize=14, color=col, va='center')",
    "    ax.text(5.0, y+0.35, label, fontsize=10, color='#1f2328', va='center')",
    "    ax.text(8.8, y+0.35, 'pass' if mark=='✓' else 'warn', fontsize=9, color=col, va='center', ha='right')",
    "ax.set_title('FastQC — read all FASTQs, score each metric', loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "Looking at one read tells you nothing about the *whole* dataset. **FastQC** is the standard first-pass tool. It scans the entire FASTQ and writes an HTML report scoring per-base quality, GC content, adapter contamination, sequence duplication, and more."
))

cells.append(code(
    "%%bash",
    "mkdir -p qc",
    "fastqc data/reads.fq -o qc 2>&1 | tail -3",
    "ls qc/"
))

cells.append(md(
    "FastQC produced a `.html` report. The next cell embeds it inline. If the iframe doesn't render in your browser, click the **Open FastQC report in new tab** button below it."
))

cells.append(code(
    "import base64, html as html_mod",
    "from IPython.display import HTML, display",
    "",
    "with open('qc/reads_fastqc.html') as f:",
    "    fastqc_html = f.read()",
    "",
    "# 1) Inline iframe via srcdoc (self-contained, sandboxed)",
    "escaped = html_mod.escape(fastqc_html, quote=True)",
    "display(HTML(",
    "    f'<iframe srcdoc=\"{escaped}\" '",
    "    f'style=\"width:100%;height:640px;border:1px solid #d0d7de;border-radius:6px;\"></iframe>'",
    "))",
    "",
    "# 2) Open-in-new-tab button (data: URL)",
    "b64 = base64.b64encode(fastqc_html.encode()).decode()",
    "display(HTML(",
    "    f'<a href=\"data:text/html;base64,{b64}\" target=\"_blank\" '",
    "    f'style=\"display:inline-block;padding:8px 16px;background:#0a7c7e;'",
    "    f'color:white;text-decoration:none;border-radius:5px;'",
    "    f'font-family:sans-serif;font-size:13px;margin-top:10px;\">'",
    "    f'📄 Open FastQC report in new tab</a>'",
    "))"
))

cells.append(md(
    "**How to read the report.** Each section shows a tick (✓ pass), warning (⚠), or cross (✗). Don't panic at warnings — many are normal for small datasets. Pay attention to:",
    "",
    "- **Per base sequence quality** — should stay green (Q≥28). Drops at the end are normal.",
    "- **Adapter content** — should be flat at zero. If it climbs at read ends, there's adapter contamination → trim.",
    "- **GC content** — should be roughly bell-shaped around the genome's mean. Bimodal = contamination."
))

# ===================================================================
# STEP 3: trim
# ===================================================================
cells.append(md(
    "---",
    "## Step 3 · Trimming with `fastp`"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "import numpy as np",
    "fig, ax = plt.subplots(figsize=(10, 2.6))",
    "ax.set_xlim(-3, 70); ax.set_ylim(0, 6); ax.axis('off')",
    "n = 50",
    "Q = 38 - (np.linspace(0, 1, n)**1.6)*32",
    "colors = ['#1a7f37' if q>=30 else '#9a6700' if q>=20 else '#bf6b00' if q>=10 else '#b1272d' for q in Q]",
    "for i, c in enumerate(colors):",
    "    ax.add_patch(Rectangle((i, 4), 0.95, 0.7, fc=c, ec='none', alpha=0.85))",
    "ax.text(-1.5, 4.35, 'read', ha='right', fontsize=10, color='#57606a', va='center')",
    "cut = n",
    "for i in range(n-4, -1, -1):",
    "    if np.mean(Q[i:i+4]) >= 20:",
    "        cut = i+4; break",
    "ax.plot([cut, cut], [3.7, 5.1], color='#b1272d', lw=2, ls='--')",
    "ax.text(cut, 5.4, '✂  cut here', ha='center', fontsize=10, color='#b1272d')",
    "for i in range(cut):",
    "    ax.add_patch(Rectangle((i, 1.6), 0.95, 0.7, fc=colors[i], ec='none', alpha=0.85))",
    "ax.text(-1.5, 1.95, 'kept', ha='right', fontsize=10, color='#1a7f37', va='center')",
    "ax.text(n+1, 4.35, 'fastp slides a 4-base window\\nfrom 3′; trims when mean Q < 20',",
    "        fontsize=9, color='#57606a', va='center')",
    "ax.set_title('fastp — sliding-window 3′ quality trim', loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "FastQC just *reports* problems — it doesn't fix them. **fastp** trims adapters, slides a 4-base window from the 3′ end, and cuts when the window's mean quality drops below Q20."
))

cells.append(code(
    "%%bash",
    "fastp \\",
    "  -i data/reads.fq \\",
    "  -o data/trimmed.fq \\",
    "  -j data/qc.json -h data/qc.html \\",
    "  --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \\",
    "  2>&1 | tail -20"
))

cells.append(md(
    "Now plot the per-base mean quality before and after, from the JSON report fastp wrote:"
))

cells.append(code(
    "import json, matplotlib.pyplot as plt",
    "qc = json.load(open('data/qc.json'))",
    "before = qc['read1_before_filtering']['quality_curves']['mean']",
    "after  = qc['read1_after_filtering']['quality_curves']['mean']",
    "",
    "fig, ax = plt.subplots(figsize=(10, 3.5))",
    "ax.axhspan(28, 40, color='#1a7f37', alpha=0.07)",
    "ax.axhspan(20, 28, color='#9a6700', alpha=0.10)",
    "ax.axhspan(0,  20, color='#b1272d', alpha=0.07)",
    "ax.plot(before, '--', color='#b1272d', lw=1.6, label='before trim')",
    "ax.plot(after,  '-',  color='#1a7f37', lw=1.6, label='after trim')",
    "ax.set_xlabel('position along read'); ax.set_ylabel('mean Phred Q')",
    "ax.set_ylim(0, 41); ax.legend(loc='lower left', frameon=False)",
    "for s in ['top','right']: ax.spines[s].set_visible(False)",
    "ax.set_title('Per-base mean quality')",
    "plt.tight_layout(); plt.show()",
    "",
    "b, a = qc['summary']['before_filtering'], qc['summary']['after_filtering']",
    "print(f\"{'metric':<12} {'before':>14} {'after':>14}\")",
    "print('-' * 42)",
    "print(f\"{'reads':<12} {b['total_reads']:>14,} {a['total_reads']:>14,}\")",
    "print(f\"{'bases':<12} {b['total_bases']:>14,} {a['total_bases']:>14,}\")",
    "print(f\"{'Q30 rate':<12} {b['q30_rate']*100:>13.1f}% {a['q30_rate']*100:>13.1f}%\")",
    "print(f\"{'GC%':<12} {b['gc_content']*100:>13.1f}% {a['gc_content']*100:>13.1f}%\")"
))

cells.append(md(
    "📈 The green line (after) sits above the red dashed line (before) at the 3′ end — that's the trim. Total bases drop, but the *quality* of what remains is higher."
))

# ===================================================================
# STEP 4: align
# ===================================================================
cells.append(md(
    "---",
    "## Step 4 · Alignment with BWA-MEM"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "import numpy as np",
    "fig, ax = plt.subplots(figsize=(10, 3.5))",
    "ax.set_xlim(0, 70); ax.set_ylim(0, 8); ax.axis('off')",
    "ax.add_patch(Rectangle((2, 5.7), 56, 0.8, fc='#e2e6ea', ec='#57606a', lw=0.5))",
    "ax.text(1.5, 6.1, 'REF', ha='right', fontsize=10, color='#57606a')",
    "ax.text(30, 7.0, 'reference (chr11 — HBB region)', ha='center', fontsize=9, color='#57606a')",
    "np.random.seed(2)",
    "n_reads = 8",
    "starts = sorted(np.random.choice(range(2, 50), n_reads, replace=False))",
    "lanes = [0,1,0,1,0,1,0,1]",
    "for s, lane in zip(starts, lanes):",
    "    ax.add_patch(Rectangle((s, 4.6 - lane*0.7), 8, 0.55, fc='#cce5ff', ec='#0969da', lw=0.4))",
    "ax.annotate('', xy=(30, 5.5), xytext=(30, 4.8),",
    "            arrowprops=dict(arrowstyle='->', color='#0a7c7e', lw=2))",
    "ax.text(34, 5.0, 'BWA-MEM', fontsize=10, color='#0a7c7e', weight='bold')",
    "for s, lane in zip(starts, lanes):",
    "    ax.add_patch(Rectangle((s, 2.3 - lane*0.7), 8, 0.55, fc='#cce5ff', ec='#0969da', lw=0.4))",
    "ax.text(1.5, 2.5, 'reads', ha='right', fontsize=10, color='#57606a')",
    "ax.add_patch(Rectangle((starts[2]+5, 2.3 - lanes[2]*0.7), 0.7, 0.55, fc='#b1272d', ec='none'))",
    "ax.text(60, 1.7, '◀ a mismatch (potential variant)', fontsize=9, color='#b1272d', va='center')",
    "ax.set_title('BWA-MEM — for every read, find best position on REF',",
    "             loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "Each read is a short anonymous string. **Alignment** finds where on the reference it came from.",
    "",
    "**BWA-MEM** is the workhorse aligner for short reads (Burrows-Wheeler Aligner, MEM = Maximal Exact Match). Used by 1000 Genomes, GATK, and most clinical labs.",
    "",
    "Two commands: `bwa index` builds a search index of the reference (run once); `bwa mem` aligns reads."
))

cells.append(code(
    "%%bash",
    "bwa index data/ref.fa 2>&1 | tail -3",
    "bwa mem -t 2 data/ref.fa data/trimmed.fq 2> data/bwa.log > data/aligned.sam",
    "echo",
    "echo '--- bwa log (last lines) ---'",
    "tail -8 data/bwa.log",
    "echo",
    "echo '--- aligned.sam first 4 lines ---'",
    "head -4 data/aligned.sam | cut -c1-180"
))

cells.append(md(
    "**SAM format**: lines starting with `@` are headers; the rest are alignments. Every alignment row has these tab-separated columns:",
    "",
    "| col | name | meaning |",
    "| --- | --- | --- |",
    "| 1 | QNAME | read name |",
    "| 2 | FLAG  | bit flags (mapped? reverse strand?) |",
    "| 3 | RNAME | reference contig |",
    "| 4 | POS   | leftmost mapping position (1-based) |",
    "| 5 | MAPQ  | mapping quality (0–60, higher = more confident) |",
    "| 6 | CIGAR | match/mismatch/indel pattern, e.g. `100M` |",
    "| 10| SEQ   | the read sequence |",
    "| 11| QUAL  | the read quality string |"
))

# ===================================================================
# STEP 5: sort + index
# ===================================================================
cells.append(md(
    "---",
    "## Step 5 · Sort and index the alignments"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "import numpy as np",
    "fig, ax = plt.subplots(figsize=(10, 2.7))",
    "ax.set_xlim(-2, 60); ax.set_ylim(0, 5); ax.axis('off')",
    "np.random.seed(7)",
    "starts_random = np.random.choice(range(2, 48), 8, replace=False)",
    "starts_sorted = sorted(starts_random)",
    "ax.text(-1.5, 4.0, 'unsorted SAM', fontsize=10, color='#57606a')",
    "for i, s in enumerate(starts_random):",
    "    ax.add_patch(Rectangle((s, 3.5 - (i%4)*0.25), 5, 0.18, fc='#cce5ff', ec='#0969da', lw=0.4))",
    "ax.annotate('', xy=(30, 2.4), xytext=(30, 3.1),",
    "            arrowprops=dict(arrowstyle='->', color='#0a7c7e', lw=2))",
    "ax.text(33, 2.7, 'samtools sort + index', fontsize=10, color='#0a7c7e', weight='bold')",
    "ax.text(-1.5, 1.8, 'sorted BAM (+ .bai)', fontsize=10, color='#57606a')",
    "for i, s in enumerate(starts_sorted):",
    "    ax.add_patch(Rectangle((s, 1.3 - (i%4)*0.25), 5, 0.18, fc='#cce5ff', ec='#0969da', lw=0.4))",
    "ax.set_title('Sort reads by position so callers can stream through quickly',",
    "             loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "BWA outputs reads in random order. Variant callers, browsers, and depth tools all expect them sorted by position. We also need an index for random access into any region of the genome."
))

cells.append(code(
    "%%bash",
    "samtools sort  data/aligned.sam -o data/aligned.bam",
    "samtools index data/aligned.bam",
    "echo '--- flagstat ---'",
    "samtools flagstat data/aligned.bam"
))

cells.append(md(
    "Coverage across the reference, with rs334 (sickle site) marked:"
))

cells.append(code(
    "import subprocess, matplotlib.pyplot as plt",
    "raw = subprocess.check_output(['samtools', 'depth', '-a', 'data/aligned.bam'], text=True)",
    "rows = [l.split('\\t') for l in raw.strip().split('\\n') if l]",
    "pos = [int(r[1]) for r in rows]",
    "dep = [int(r[2]) for r in rows]",
    "RS334 = 5248232 - 5246000",
    "fig, ax = plt.subplots(figsize=(11, 3))",
    "ax.fill_between(pos, dep, color='#0a7c7e', alpha=0.4)",
    "ax.plot(pos, dep, color='#0a7c7e', lw=0.9)",
    "ax.axvline(RS334, color='#b1272d', lw=1, ls='--', label='rs334 (HbS)')",
    "ax.set_xlabel('position on HBB region (chr11:5,246,001-5,250,000, GRCh37)')",
    "ax.set_ylabel('read depth')",
    "ax.legend(loc='upper right', frameon=False)",
    "for s in ['top','right']: ax.spines[s].set_visible(False)",
    "ax.set_title('Coverage across HBB')",
    "plt.tight_layout(); plt.show()",
    "print(f'mean depth: {sum(dep)/len(dep):.1f}x')",
    "print(f'positions with ≥1 read: {sum(1 for d in dep if d):,} / {len(dep):,}')"
))

cells.append(md(
    "📊 Coverage is uneven because this is *low-coverage* sequencing (the 1000 Genomes phase-1 design: many people, lightly sequenced). Average ~3–5× is barely enough for heterozygous calls."
))

# ===================================================================
# STEP 6: pileup
# ===================================================================
cells.append(md(
    "---",
    "## Step 6 · Pile up reads at every position"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "fig, ax = plt.subplots(figsize=(10, 3.6))",
    "ax.set_xlim(-1, 13); ax.set_ylim(0, 7); ax.axis('off')",
    "REF = 'TCTCCTTAGAGT'",
    "cmap = {'A':'#1a7f37', 'T':'#b1272d', 'C':'#0969da', 'G':'#9a6700', '·':'#aaa'}",
    "for i, b in enumerate(REF):",
    "    ax.add_patch(Rectangle((i, 5.7), 0.92, 0.7, fc='#e2e6ea', ec='#d0d7de', lw=0.4))",
    "    ax.text(i+0.46, 6.05, b, ha='center', va='center', family='monospace', fontsize=11, color=cmap[b])",
    "ax.text(-0.5, 6.05, 'REF', ha='right', va='center', fontsize=10, color='#57606a')",
    "stack = ['·····A···   ', '·····T···   ', '·····A···   ', '·····A···   ', '·····T···   ']",
    "for k, row in enumerate(stack):",
    "    for i, b in enumerate(row):",
    "        if b == ' ': continue",
    "        is_match = b == '·'",
    "        ax.add_patch(Rectangle((i, 4.6-k*0.6), 0.92, 0.55,",
    "                               fc=cmap.get(b,'#aaa') if not is_match else '#f6f8fa',",
    "                               ec='#d0d7de', lw=0.3))",
    "        if not is_match:",
    "            ax.text(i+0.46, 4.6-k*0.6+0.27, b, ha='center', va='center', family='monospace',",
    "                    fontsize=11, color='white', weight='bold')",
    "ax.text(-0.5, 4.0, 'reads', ha='right', va='center', fontsize=10, color='#57606a')",
    "ax.add_patch(Rectangle((4.96, 1.4), 1, 5.2, fill=False, ec='#b1272d', lw=1.4, ls='--'))",
    "ax.text(5.46, 0.9, 'mixed bases here = variant', ha='center', fontsize=9, color='#b1272d')",
    "ax.set_title('mpileup — column-wise view of every reference position',",
    "             loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "We *transpose* the data. Instead of looking at reads, we look at the genome **column by column**: at each reference position, what bases do the reads observe?",
    "",
    "Output format (tab-separated): **contig · position · REF · depth · read bases · read qualities**.",
    "",
    "The read-bases string uses a special encoding: `.` = match REF on forward strand · `,` = match REF on reverse strand · `A C G T` = mismatch (forward) · `a c g t` = mismatch (reverse) · `^]` = read start · `$` = read end · `+N…` / `-N…` = insertion / deletion."
))

cells.append(code(
    "%%bash",
    "samtools mpileup -f data/ref.fa data/aligned.bam 2>/dev/null > data/pileup.txt",
    "echo \"total pileup rows: $(wc -l < data/pileup.txt)\"",
    "echo",
    "echo '--- rows around the sickle position (HBB:2232) ---'",
    "awk '$2 >= 2228 && $2 <= 2236' data/pileup.txt"
))

cells.append(md(
    "**Look at line 2232.** How many reads see `A` vs `T`? That's the variant we're looking for.",
    "",
    "Below — visual pileup. Each column is one reference position; each row of the stack is one read. Bases matching REF are faded; mismatches are highlighted. Red dashed box = sickle site."
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle",
    "rows = []",
    "with open('data/pileup.txt') as f:",
    "    for line in f:",
    "        c, p, ref, dp, bases, qual = line.rstrip().split('\\t')",
    "        if 2218 <= int(p) <= 2246:",
    "            rows.append((int(p), ref.upper(), int(dp), bases))",
    "",
    "def parse_bases(refb, s):",
    "    out, i = [], 0",
    "    while i < len(s):",
    "        ch = s[i]",
    "        if   ch == '^':              i += 2;  continue",
    "        elif ch == '$':              i += 1;  continue",
    "        elif ch in '+-':",
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
    "fig, ax = plt.subplots(figsize=(13, max_stack*0.22 + 1.5))",
    "for i, (p, ref, dp, bases) in enumerate(rows):",
    "    ax.text(i, max_stack+1, ref, ha='center', va='center',",
    "            color=cmap[ref], family='monospace', fontsize=11, weight='bold')",
    "    for k, b in enumerate(parse_bases(ref, bases)):",
    "        match = (b == ref)",
    "        ax.add_patch(Rectangle((i-0.4, max_stack-1-k), 0.8, 0.85,",
    "                               color=cmap.get(b,'#777'), alpha=0.18 if match else 0.92))",
    "        if not match:",
    "            ax.text(i, max_stack-1-k+0.42, b, ha='center', va='center',",
    "                    color='white', family='monospace', fontsize=8)",
    "if any(p == 2232 for p, *_ in rows):",
    "    RS = next(i for i,(p,*_) in enumerate(rows) if p == 2232)",
    "    ax.add_patch(Rectangle((RS-0.5, -0.5), 1, max_stack+2.2,",
    "                           fill=False, ec='#b1272d', lw=1.5, ls='--'))",
    "ax.set_xticks(range(len(rows)))",
    "ax.set_xticklabels([str(r[0]) for r in rows], rotation=90, fontsize=8)",
    "ax.set_yticks([]); ax.set_xlim(-1, len(rows)); ax.set_ylim(-1, max_stack+2.5)",
    "ax.set_title('Pileup window — REF row at top, reads stacked below; '",
    "             + 'red box = rs334', fontsize=10)",
    "for s in ax.spines.values(): s.set_visible(False)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "🔴 At HBB:2232 (red box), roughly half the reads show **A** instead of REF **T** — a classic heterozygous SNV."
))

# ===================================================================
# STEP 7: variant call
# ===================================================================
cells.append(md(
    "---",
    "## Step 7 · Calling variants → VCF"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle, FancyBboxPatch",
    "fig, ax = plt.subplots(figsize=(10, 3.0))",
    "ax.set_xlim(0, 14); ax.set_ylim(0, 6); ax.axis('off')",
    "ax.add_patch(Rectangle((1, 4.6), 0.7, 0.6, fc='#e2e6ea', ec='#d0d7de'))",
    "ax.text(1.35, 4.9, 'T', ha='center', va='center', family='monospace', fontsize=12, color='#b1272d', weight='bold')",
    "ax.text(0.6, 4.9, 'REF', ha='right', va='center', fontsize=9, color='#57606a')",
    "stack = ['T','A','T','A','A']",
    "cmap = {'A':'#1a7f37','T':'#b1272d'}",
    "for k, b in enumerate(stack):",
    "    ax.add_patch(Rectangle((1, 3.8 - k*0.55), 0.7, 0.5,",
    "                           fc=cmap[b] if b!='T' else '#f6f8fa', ec='#d0d7de', lw=0.3))",
    "    if b != 'T':",
    "        ax.text(1.35, 3.8-k*0.55+0.25, b, ha='center', va='center', family='monospace', color='white', weight='bold')",
    "    else:",
    "        ax.text(1.35, 3.8-k*0.55+0.25, '·', ha='center', va='center', family='monospace', color='#aaa')",
    "ax.text(1.35, 0.6, 'pileup column\\nat HBB:2232', ha='center', fontsize=9, color='#57606a')",
    "ax.annotate('', xy=(5.5, 3.5), xytext=(2.6, 3.5),",
    "            arrowprops=dict(arrowstyle='->', color='#0a7c7e', lw=2))",
    "ax.text(4.05, 3.9, 'bcftools', ha='center', fontsize=10, color='#0a7c7e', weight='bold')",
    "ax.add_patch(FancyBboxPatch((6.0, 2.6), 7.5, 1.7, boxstyle='round,pad=0.1',",
    "                            fc='#f6f8fa', ec='#d0d7de'))",
    "ax.text(6.3, 3.85, 'CHROM  POS    REF  ALT  QUAL  GT', family='monospace', fontsize=9, color='#57606a')",
    "ax.text(6.3, 3.30, 'HBB    2232   T    A    105   0/1', family='monospace', fontsize=11, color='#1f2328')",
    "ax.text(9.75, 0.6, 'one row of variants.vcf', ha='center', fontsize=9, color='#57606a')",
    "ax.set_title('bcftools — pileup likelihoods → variant calls (VCF)',",
    "             loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "We've seen the variant by eye. Now let `bcftools` call it formally. Two steps:",
    "",
    "1. `bcftools mpileup` — base likelihoods at every position (BCF format)",
    "2. `bcftools call -mv` — call variant sites only",
    "",
    "The output is a **VCF** (Variant Call Format) — the universal variant file."
))

cells.append(code(
    "%%bash",
    "bcftools mpileup -f data/ref.fa data/aligned.bam -Ou -o data/pile.bcf 2>&1 | tail -2",
    "bcftools call -mv -Oz -o data/variants.vcf.gz data/pile.bcf 2>&1 | tail -2",
    "bcftools index -f data/variants.vcf.gz",
    "echo",
    "echo '--- variants.vcf.gz (## meta lines hidden) ---'",
    "bcftools view data/variants.vcf.gz 2>/dev/null | grep -v '^##' | head -25"
))

cells.append(md(
    "Each non-header line is one variant. Key fields: **POS** (1-based), **REF/ALT**, **QUAL**, **INFO** (`DP` = depth, `MQ` = mean MAPQ, `DP4` = ref-fwd, ref-rev, alt-fwd, alt-rev), **FORMAT** + per-sample (`GT` = `0/0`/`0/1`/`1/1`)."
))

# ===================================================================
# STEP 8: interpret
# ===================================================================
cells.append(md(
    "---",
    "## Step 8 · Interpreting the variant — finding sickle cell"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle, FancyBboxPatch",
    "fig, ax = plt.subplots(figsize=(11, 3.6))",
    "ax.set_xlim(0, 14); ax.set_ylim(0, 7); ax.axis('off')",
    "cmap = {'A':'#1a7f37','T':'#b1272d','C':'#0969da','G':'#9a6700'}",
    "ax.text(0.3, 5.7, 'REF cDNA codon 7:', fontsize=10, color='#57606a')",
    "for i, b in enumerate('GAG'):",
    "    ax.add_patch(Rectangle((4 + i*0.7, 5.4), 0.65, 0.7, fc=cmap[b], alpha=0.2, ec=cmap[b]))",
    "    ax.text(4 + i*0.7 + 0.32, 5.75, b, ha='center', va='center', family='monospace', fontsize=12, color=cmap[b], weight='bold')",
    "ax.annotate('', xy=(7.5, 5.75), xytext=(6.4, 5.75), arrowprops=dict(arrowstyle='->', color='#57606a'))",
    "ax.add_patch(FancyBboxPatch((7.7, 5.4), 2.4, 0.7, boxstyle='round,pad=0.05', fc='#f6f8fa', ec='#d0d7de'))",
    "ax.text(8.9, 5.75, 'Glu (E)', ha='center', va='center', fontsize=11, color='#1f2328')",
    "ax.text(11.0, 5.75, 'normal β-globin', va='center', fontsize=10, color='#57606a')",
    "ax.annotate('', xy=(8.0, 4.0), xytext=(8.0, 5.2), arrowprops=dict(arrowstyle='->', color='#b1272d', lw=2))",
    "ax.text(8.3, 4.6, 'rs334 — single base A→T (cDNA)', fontsize=9, color='#b1272d')",
    "ax.text(0.3, 3.3, 'ALT cDNA codon 7:', fontsize=10, color='#57606a')",
    "for i, b in enumerate('GTG'):",
    "    fc = '#b1272d' if i == 1 else cmap[b]",
    "    alpha = 0.5 if i == 1 else 0.2",
    "    ax.add_patch(Rectangle((4 + i*0.7, 3.0), 0.65, 0.7, fc=fc, alpha=alpha, ec=cmap[b]))",
    "    ax.text(4 + i*0.7 + 0.32, 3.35, b, ha='center', va='center', family='monospace', fontsize=12, color=cmap[b], weight='bold')",
    "ax.annotate('', xy=(7.5, 3.35), xytext=(6.4, 3.35), arrowprops=dict(arrowstyle='->', color='#57606a'))",
    "ax.add_patch(FancyBboxPatch((7.7, 3.0), 2.4, 0.7, boxstyle='round,pad=0.05', fc='#f6f8fa', ec='#b1272d', lw=1.5))",
    "ax.text(8.9, 3.35, 'Val (V)', ha='center', va='center', fontsize=11, color='#b1272d', weight='bold')",
    "ax.text(11.0, 3.35, 'HbS β-globin', va='center', fontsize=10, color='#b1272d')",
    "ax.add_patch(FancyBboxPatch((0.3, 0.4), 13.4, 1.7, boxstyle='round,pad=0.05',",
    "                            fc='#fff5f5', ec='#b1272d', lw=1.0))",
    "ax.text(0.7, 1.55, '🩸 Phenotype:', fontsize=11, color='#b1272d', weight='bold')",
    "ax.text(0.7, 1.05,",
    "        'het (HbAS) → sickle cell trait (carrier, malaria-protective)   ·   '",
    "        'hom (HbSS) → sickle cell anaemia',",
    "        fontsize=10, color='#1f2328')",
    "ax.set_title('rs334 — the molecular basis of sickle cell disease', loc='left', color='#57606a', fontsize=11)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "We started with raw FASTQ and ended with a VCF. **Now the biology.**",
    "",
    "We're looking for **rs334**:",
    "- chr11:5,248,232 (GRCh37) → in our subset, HBB:2,232",
    "- T → A on the genomic forward strand. *HBB* is on the reverse strand, so on the cDNA: A → T (HGVS `c.20A>T`)",
    "- That falls in codon 7 (`GAG`, glutamic acid). The mutation makes it `GTG` (valine). HGVS protein: `p.Glu7Val`",
    "- Polar Glu → hydrophobic Val on the β-globin surface lets HbS polymerize when oxygen is low; red blood cells deform into the characteristic sickle shape",
    "- Heterozygotes have **sickle cell trait** (HbAS, usually asymptomatic, malaria-protective). Homozygotes have **sickle cell anaemia** (HbSS)"
))

cells.append(code(
    "import subprocess",
    "vcf = subprocess.check_output(['bcftools', 'view', 'data/variants.vcf.gz'], text=True)",
    "rows = [l.split('\\t') for l in vcf.splitlines() if l and not l.startswith('#')]",
    "rs334 = next((r for r in rows if int(r[1]) == 2232), None)",
    "",
    "if rs334 is None:",
    "    print('rs334 was not called in this run — check coverage at HBB:2232.')",
    "else:",
    "    fmt  = dict(zip(rs334[8].split(':'), rs334[9].split(':')))",
    "    info = dict(kv.split('=', 1) for kv in rs334[7].split(';') if '=' in kv)",
    "    dp   = info.get('DP', '?')",
    "    ad   = fmt.get('AD', '?,?').split(',')",
    "    print('=' * 64)",
    "    print('  rs334 (sickle cell allele) called from real HG02666 reads')",
    "    print('=' * 64)",
    "    print(f'  position    HBB:{rs334[1]}  =  chr11:5,248,232  (GRCh37)')",
    "    print(f'  REF -> ALT  {rs334[3]} -> {rs334[4]}      (genomic forward strand)')",
    "    print( '  cDNA        c.20 A>T   (HBB on reverse strand; revcomp of T>A)')",
    "    print( '  protein     codon 7  GAG (Glu) -> GTG (Val)  =  p.Glu7Val')",
    "    print(f'  genotype    {fmt[\"GT\"]}     (0/1 = heterozygous = sickle cell trait, HbAS)')",
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
    "🎉 You have just walked an Illumina FASTQ from a real human all the way to a clinically meaningful variant call, using the same software a hospital genetics lab would use."
))

# ===================================================================
# CLOSING
# ===================================================================
cells.append(md(
    "---",
    "## Pipeline summary",
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
    "1. **Run with a non-carrier sample.** Pick another 1000 Genomes sample (e.g. NA12878, European — almost certainly no rs334) and edit the `curl` URL in the data-download cell. The pipeline runs unchanged but the VCF will not contain a row at HBB:2232.",
    "2. **Try a homozygote.** A few 1000 Genomes samples are HbSS. Their variant call at HBB:2232 will show `1/1` instead of `0/1`, and `DP4` will have zero REF reads.",
    "3. **Vary the trim aggressiveness.** Set `--cut_tail_mean_quality` to 30 (more aggressive). How do the final variant calls change? Could over-trimming *miss* a variant?",
    "4. **Add VEP / SnpEff.** Real clinical pipelines run an annotator that maps every variant to its consequence (synonymous / missense / splice / …) and looks up gnomAD frequencies and ClinVar significance.",
    "",
    "## Further reading",
    "",
    "- 1000 Genomes Project Consortium · *A global reference for human genetic variation*, Nature 526, 68–74 (2015)",
    "- Heng Li · *Aligning sequence reads with BWA-MEM*, arXiv:1303.3997 (2013)",
    "- Petr Danecek et al. · *Twelve years of SAMtools and BCFtools*, GigaScience 10, giab008 (2021)",
    "- HBB / sickle cell · OMIM 603903, ClinVar VCV000015333, dbSNP rs334",
    "",
    "---",
    "",
    "*Notebook by **Han-Ying Jhuang (莊漢英), PhD** — Assistant Professor, Taipei Medical University. Code released under MIT; data redistributed under the [1000 Genomes data use policy](https://www.internationalgenome.org/data-portal/data-collection). Questions or improvements: [hanyingjhuang@tmu.edu.tw](mailto:hanyingjhuang@tmu.edu.tw).*"
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

n_md   = sum(1 for c in cells if c['cell_type']=='markdown')
n_code = sum(1 for c in cells if c['cell_type']=='code')
print(f"wrote {len(cells)} cells  ({n_md} markdown, {n_code} code)")
