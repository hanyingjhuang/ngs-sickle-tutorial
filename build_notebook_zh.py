"""簡化版中文教學筆記本（繁體中文・台灣用語）。

設計原則：
- 程式碼保留英文
- Markdown 解說與 print() 輸出使用繁體中文（台灣用語）
- matplotlib 圖表內的標籤保留英文，避免 Colab CJK 字型對齊問題
- 比英文版更精簡：每步驟僅留 1–2 句說明
"""
import json

REPO = "hanyingjhuang/ngs-sickle-tutorial"
RAW = f"https://raw.githubusercontent.com/{REPO}/main"

def md(*lines):
    src = [l + ("\n" if i < len(lines) - 1 else "") for i, l in enumerate(lines)]
    return {"cell_type": "markdown", "metadata": {}, "source": src}

def code(*lines, hide=True):
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
#  封面
# ===================================================================
cells.append(md(
    "# NGS 變異分析流程 — 從原始定序資料到臨床變異",
    "",
    "**作者** · 莊漢英 助理教授・臺北醫學大學  ",
    "**聯絡** · hanyingjhuang@tmu.edu.tw",
    "",
    "---",
    "",
    "## 這份筆記能學到什麼",
    "",
    "本筆記從 Illumina 定序儀的原始輸出（FASTQ 檔）開始，一路操作到最後一步——在真實受試者的基因組中**找出鐮刀型細胞貧血症的致病變異**。",
    "",
    "資料完全公開：使用 [千人基因組計畫](https://www.internationalgenome.org/) 中的 **HG02666**（甘比亞籍個體）讀段，已知為鐮刀型基因 (rs334) 異型合子帶因者；參考序列為人類 GRCh37 第 11 號染色體上 *HBB* 基因附近的 4 kb 片段。",
    "",
    "## 流程",
    "",
    "```",
    "  原始讀段          品管 + 修剪          比對             已排序 BAM",
    "    FASTQ    ─►   FastQC + fastp   ─►   BWA-MEM    ─►    samtools",
    "                                                              │",
    "                                                              ▼",
    "                                  變異呼叫              鹼基堆疊圖",
    "                                  bcftools         ◄─    samtools mpileup",
    "                                       │",
    "                                       ▼",
    "                                  rs334 / HbS",
    "                                  臨床判讀",
    "```",
    "",
    "## 使用方式",
    "",
    "點擊 **執行階段 → 全部執行 (Runtime → Run all)**，或逐格按 **Shift + Enter**。",
    "",
    "**程式碼預設摺疊**——學生看到的是結果，點左側細條才會展開原始指令。每個步驟跑的都是 *真正的 Linux 二進位執行檔*（`bwa`、`samtools`、`bcftools`…），不是模擬。",
    "",
    "**註**：圖內標籤保留英文（避免 Colab 預設字型缺中文字符），其餘解說與輸出皆為中文。"
))

# ===================================================================
# Step 0
# ===================================================================
cells.append(md(
    "## 0 · 環境設定",
    "",
    "Colab 每次給的是全新 Linux 機器，先用 `apt-get` 把工具裝上。"
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
    "✅ 看到版本號就代表工具已準備好。"
))

# ===================================================================
# Step 0b
# ===================================================================
cells.append(md(
    "## 0b · 下載真實資料",
    "",
    "需要兩個檔案：**參考基因組** (`ref.fa`) 與**定序讀段** (`reads.fq`)，都從 GitHub 直接抓。"
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

# ===================================================================
# Step 1
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 1 · 一條定序讀段長什麼樣？"
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
    "ax.set_title('FASTQ — one read = 4 lines', loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "FASTQ 一條讀段固定 **4 行**：名稱、序列、`+`、品質字串。",
    "",
    "品質字串的每個字元都代表那個鹼基的 Phred 分數（`Phred = ord(字元) − 33`）。`'I'` (Q40) 表示 1 萬分之 1 的錯誤率；`'#'` (Q2) 幾乎是亂猜。"
))

cells.append(code(
    "!head -8 data/reads.fq"
))

cells.append(md(
    "讀段尾端常出現 `;<=` 這類雜訊字元——Illumina 化學反應到讀段末端會劣化，3′ 端品質變差。"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "with open('data/reads.fq') as f:",
    "    lines = [next(f).rstrip() for _ in range(8)]",
    "for r in range(2):",
    "    name, seq, _, qual = lines[r*4 : r*4 + 4]",
    "    Q = [ord(c) - 33 for c in qual]",
    "    cmap = {'A':'#1a7f37', 'T':'#b1272d', 'C':'#0969da', 'G':'#9a6700', 'N':'#777'}",
    "    bar_colors = ['#1a7f37' if q>=30 else '#9a6700' if q>=20",
    "                  else '#bf6b00' if q>=10 else '#b1272d' for q in Q]",
    "    fig, ax = plt.subplots(figsize=(14, 2.0))",
    "    ax.bar(range(len(Q)), Q, color=bar_colors, width=0.95)",
    "    ax.set_xlim(-0.5, len(seq)-0.5); ax.set_ylim(0, 42)",
    "    ax.set_yticks([0, 20, 30, 40])",
    "    ax.set_xticks(range(len(seq)))",
    "    ax.set_xticklabels(list(seq), family='monospace', fontsize=8)",
    "    for tick, b in zip(ax.get_xticklabels(), seq):",
    "        tick.set_color(cmap.get(b, '#000'))",
    "    ax.tick_params(axis='x', length=0, pad=2)",
    "    ax.set_ylabel('Phred Q')",
    "    ax.set_title(name, fontsize=9, loc='left')",
    "    for s in ['top', 'right']: ax.spines[s].set_visible(False)",
    "    plt.tight_layout(); plt.show()"
))

cells.append(md(
    "👀 注意右側——長條變短、由黃轉橘紅，這就是要在步驟 3 修掉的低品質尾端。"
))

# ===================================================================
# Step 2
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 2 · 用 FastQC 做品管"
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
    "         ('adapter',    1.6, '!', '#9a6700'),",
    "         ('duplication',0.7, '✓', '#1a7f37')]",
    "for label, y, mark, col in tiles:",
    "    ax.add_patch(Rectangle((4.2, y), 5.0, 0.7, fc='#f6f8fa', ec='#d0d7de', lw=0.5))",
    "    ax.text(4.4, y+0.35, mark, fontsize=14, color=col, va='center')",
    "    ax.text(5.0, y+0.35, label, fontsize=10, color='#1f2328', va='center')",
    "    ax.text(8.8, y+0.35, 'pass' if mark=='✓' else 'warn', fontsize=9, color=col, va='center', ha='right')",
    "ax.set_title('FastQC — read all FASTQs, score each metric', loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "FastQC 是業界標準的品管工具：掃過整份 FASTQ，產出 HTML 報告，逐項給通過 / 警告 / 失敗。"
))

cells.append(code(
    "%%bash",
    "mkdir -p qc",
    "fastqc data/reads.fq -o qc 2>&1 | tail -3",
    "ls qc/"
))

cells.append(md(
    "下方第一個格子內嵌 FastQC 的 HTML 報告。如果無法在內嵌畫面中顯示，請按下方按鈕另開分頁查看。"
))

cells.append(code(
    "import base64, html as html_mod",
    "from IPython.display import HTML, display",
    "with open('qc/reads_fastqc.html') as f:",
    "    fastqc_html = f.read()",
    "escaped = html_mod.escape(fastqc_html, quote=True)",
    "display(HTML(",
    "    f'<iframe srcdoc=\"{escaped}\" '",
    "    f'style=\"width:100%;height:640px;border:1px solid #d0d7de;border-radius:6px;\"></iframe>'",
    "))",
    "b64 = base64.b64encode(fastqc_html.encode()).decode()",
    "display(HTML(",
    "    f'<a href=\"data:text/html;base64,{b64}\" target=\"_blank\" '",
    "    f'style=\"display:inline-block;padding:8px 16px;background:#0a7c7e;'",
    "    f'color:white;text-decoration:none;border-radius:5px;'",
    "    f'font-family:sans-serif;font-size:13px;margin-top:10px;\">'",
    "    f'📄 另開分頁顯示 FastQC 報告</a>'",
    "))"
))

cells.append(md(
    "報告中要特別注意 **每位置品質**（應在綠色區）、**接頭污染**（應為零）、**GC 含量**（應呈鐘形分布）。"
))

# ===================================================================
# Step 3
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 3 · 用 fastp 修剪低品質尾端"
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
    "ax.text(cut, 5.4, 'cut here', ha='center', fontsize=10, color='#b1272d')",
    "for i in range(cut):",
    "    ax.add_patch(Rectangle((i, 1.6), 0.95, 0.7, fc=colors[i], ec='none', alpha=0.85))",
    "ax.text(-1.5, 1.95, 'kept', ha='right', fontsize=10, color='#1a7f37', va='center')",
    "ax.text(n+1, 4.35, \"slide a 4-bp window from 3';\\ncut where mean Q < 20\",",
    "        fontsize=9, color='#57606a', va='center')",
    "ax.set_title(\"fastp — sliding-window 3' quality trim\", loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "FastQC 只報告問題；**fastp** 才會動手修剪。它從 3′ 端滑動 4-鹼基視窗，平均品質低於 Q20 就剪掉。"
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
    "從 fastp 寫的 JSON 把修剪前後的每位置平均品質畫出來："
))

cells.append(code(
    "import json, matplotlib.pyplot as plt",
    "qc = json.load(open('data/qc.json'))",
    "before = qc['read1_before_filtering']['quality_curves']['mean']",
    "after  = qc['read1_after_filtering']['quality_curves']['mean']",
    "fig, ax = plt.subplots(figsize=(10, 3.5))",
    "ax.axhspan(28, 40, color='#1a7f37', alpha=0.07)",
    "ax.axhspan(20, 28, color='#9a6700', alpha=0.10)",
    "ax.axhspan(0,  20, color='#b1272d', alpha=0.07)",
    "ax.plot(before, '--', color='#b1272d', lw=1.6, label='before trim')",
    "ax.plot(after,  '-',  color='#1a7f37', lw=1.6, label='after trim')",
    "ax.set_xlabel('position along read'); ax.set_ylabel('mean Phred Q')",
    "ax.set_ylim(0, 41); ax.legend(loc='lower left', frameon=False)",
    "for s in ['top','right']: ax.spines[s].set_visible(False)",
    "ax.set_title('per-base mean quality')",
    "plt.tight_layout(); plt.show()",
    "b, a = qc['summary']['before_filtering'], qc['summary']['after_filtering']",
    "print(f\"{'指標':<10} {'修剪前':>12} {'修剪後':>12}\")",
    "print('-' * 40)",
    "print(f\"讀段數    {b['total_reads']:>12,} {a['total_reads']:>12,}\")",
    "print(f\"鹼基數    {b['total_bases']:>12,} {a['total_bases']:>12,}\")",
    "print(f\"Q30 比率  {b['q30_rate']*100:>11.1f}% {a['q30_rate']*100:>11.1f}%\")",
    "print(f\"GC%      {b['gc_content']*100:>11.1f}% {a['gc_content']*100:>11.1f}%\")"
))

cells.append(md(
    "📈 綠線（修剪後）在 3′ 端高於紅虛線（修剪前）——雜訊鹼基被剪掉了。鹼基總量變少，但留下來的品質更高。"
))

# ===================================================================
# Step 4
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 4 · 用 BWA-MEM 比對到參考序列"
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
    "starts = sorted(np.random.choice(range(2, 50), 8, replace=False))",
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
    "ax.text(60, 1.7, '* mismatch (potential variant)', fontsize=9, color='#b1272d', va='center')",
    "ax.set_title('BWA-MEM — find best position on REF for each read',",
    "             loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "每條讀段都是來歷不明的短字串。**比對 (alignment)** 就是替每條讀段找出它在參考序列上最對得起來的位置。**BWA-MEM** 是業界主流的短讀段比對工具。",
    "",
    "兩個指令：`bwa index` 替參考建索引（只做一次），`bwa mem` 真正執行比對。"
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
    "**SAM 格式**：`@` 開頭的是檔頭，其他每行是一筆比對結果。重要欄位：`POS`（比對位置）、`MAPQ`（比對品質 0–60）、`CIGAR`（一致 / 不一致 / 插入缺失模式）、`SEQ` / `QUAL`（序列與品質）。"
))

# ===================================================================
# Step 5
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 5 · 排序與索引 BAM"
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
    "ax.set_title('sort by position so callers can stream through quickly',",
    "             loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "BWA 輸出順序混亂；變異呼叫工具與基因瀏覽器都要求 BAM 依位置排序，並產生索引 (`.bai`) 才能隨機存取。"
))

cells.append(code(
    "%%bash",
    "samtools sort  data/aligned.sam -o data/aligned.bam",
    "samtools index data/aligned.bam",
    "echo '--- flagstat ---'",
    "samtools flagstat data/aligned.bam"
))

cells.append(md(
    "全區段的覆蓋深度（紅虛線標出 rs334 鐮刀位點）："
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
    "ax.set_title('coverage across HBB')",
    "plt.tight_layout(); plt.show()",
    "print(f'平均深度：{sum(dep)/len(dep):.1f}x')",
    "print(f'有讀段覆蓋的位置：{sum(1 for d in dep if d):,} / {len(dep):,}')"
))

cells.append(md(
    "📊 深度不均是因為這是**低覆蓋**定序（千人計畫第一期策略：人多但每人讀得淺）。平均 3–5× 勉強夠呼叫雜合變異。"
))

# ===================================================================
# Step 6
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 6 · 鹼基堆疊（pileup）"
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
    "stack = ['*****A***   ', '*****T***   ', '*****A***   ', '*****A***   ', '*****T***   ']",
    "for k, row in enumerate(stack):",
    "    for i, b in enumerate(row):",
    "        if b == ' ': continue",
    "        is_match = b == '*'",
    "        ax.add_patch(Rectangle((i, 4.6-k*0.6), 0.92, 0.55,",
    "                               fc=cmap.get(b,'#aaa') if not is_match else '#f6f8fa',",
    "                               ec='#d0d7de', lw=0.3))",
    "        if not is_match:",
    "            ax.text(i+0.46, 4.6-k*0.6+0.27, b, ha='center', va='center', family='monospace',",
    "                    fontsize=11, color='white', weight='bold')",
    "ax.text(-0.5, 4.0, 'reads', ha='right', va='center', fontsize=10, color='#57606a')",
    "ax.add_patch(Rectangle((4.96, 1.4), 1, 5.2, fill=False, ec='#b1272d', lw=1.4, ls='--'))",
    "ax.text(5.46, 0.9, 'mismatched bases here = variant', ha='center', fontsize=9, color='#b1272d')",
    "ax.set_title('mpileup — column-wise view of every reference position',",
    "             loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "把資料**轉置**：不再看每條讀段，而是逐個參考位置看「這個位置上所有讀段觀察到什麼鹼基？」",
    "",
    "輸出格式（tab 分隔）：**染色體・位置・REF・深度・觀察到的鹼基・品質字串**。鹼基欄位的編碼：`.` 正股一致、`,` 反股一致、`A C G T` 正股不一致、`a c g t` 反股不一致。"
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
    "**看 2232 那一行。** 多少條讀段是 `A`、多少是 `T`？這就是要找的變異。",
    "",
    "下方視覺化堆疊圖：每欄是一個參考位置，每列是一條讀段。一致鹼基淡化處理，不一致鹼基用顏色標出；紅虛線框就是鐮刀位點。"
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
    "ax.set_title('pileup window — REF row at top, reads stacked below; '",
    "             + 'red box = rs334', fontsize=11)",
    "for s in ax.spines.values(): s.set_visible(False)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "🔴 在 HBB:2232（紅框）這欄，大約一半讀段顯示 **A**、一半 **T** ——典型的**雜合 SNV**。"
))

# ===================================================================
# Step 7
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 7 · 呼叫變異 → VCF"
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
    "        ax.text(1.35, 3.8-k*0.55+0.25, '*', ha='center', va='center', family='monospace', color='#aaa')",
    "ax.text(1.35, 0.6, 'pileup column\\nat HBB:2232', ha='center', fontsize=9, color='#57606a')",
    "ax.annotate('', xy=(5.5, 3.5), xytext=(2.6, 3.5),",
    "            arrowprops=dict(arrowstyle='->', color='#0a7c7e', lw=2))",
    "ax.text(4.05, 3.9, 'bcftools', ha='center', fontsize=10, color='#0a7c7e', weight='bold')",
    "ax.add_patch(FancyBboxPatch((6.0, 2.6), 7.5, 1.7, boxstyle='round,pad=0.1',",
    "                            fc='#f6f8fa', ec='#d0d7de'))",
    "ax.text(6.3, 3.85, 'CHROM  POS    REF  ALT  QUAL  GT', family='monospace', fontsize=9, color='#57606a')",
    "ax.text(6.3, 3.30, 'HBB    2232   T    A    105   0/1', family='monospace', fontsize=11, color='#1f2328')",
    "ax.text(9.75, 0.6, 'one row of variants.vcf', ha='center', fontsize=9, color='#57606a')",
    "ax.set_title('bcftools — pileup likelihoods to variant calls (VCF)',",
    "             loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "用肉眼看到變異後，交給 `bcftools` 正式呼叫。先 `mpileup` 算機率，再 `call -mv` 挑出真正變異的位置。輸出是 **VCF** 檔——通用的變異格式。"
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
    "每行非檔頭資料 = 一個變異。重點欄位：**POS**（位置）、**REF/ALT**（參考/替代等位）、**QUAL**（品質）、**INFO** 內的 `DP`（深度）、**GT**（基因型 `0/0` / `0/1` / `1/1`）。"
))

# ===================================================================
# Step 8
# ===================================================================
cells.append(md(
    "---",
    "## 步驟 8 · 判讀變異——找到鐮刀型變異"
))

cells.append(code(
    "import matplotlib.pyplot as plt",
    "from matplotlib.patches import Rectangle, FancyBboxPatch",
    "fig, ax = plt.subplots(figsize=(11, 3.6))",
    "ax.set_xlim(0, 14); ax.set_ylim(0, 7); ax.axis('off')",
    "cmap = {'A':'#1a7f37','T':'#b1272d','C':'#0969da','G':'#9a6700'}",
    "ax.text(0.3, 5.7, 'REF cDNA codon 7:', fontsize=10, color='#57606a')",
    "for i, b in enumerate('GAG'):",
    "    ax.add_patch(Rectangle((4.5 + i*0.7, 5.4), 0.65, 0.7, fc=cmap[b], alpha=0.2, ec=cmap[b]))",
    "    ax.text(4.5 + i*0.7 + 0.32, 5.75, b, ha='center', va='center', family='monospace', fontsize=12, color=cmap[b], weight='bold')",
    "ax.annotate('', xy=(8.0, 5.75), xytext=(6.9, 5.75), arrowprops=dict(arrowstyle='->', color='#57606a'))",
    "ax.add_patch(FancyBboxPatch((8.2, 5.4), 2.4, 0.7, boxstyle='round,pad=0.05', fc='#f6f8fa', ec='#d0d7de'))",
    "ax.text(9.4, 5.75, 'Glu (E)', ha='center', va='center', fontsize=11, color='#1f2328')",
    "ax.text(11.0, 5.75, 'normal beta-globin', va='center', fontsize=10, color='#57606a')",
    "ax.annotate('', xy=(8.0, 4.0), xytext=(8.0, 5.2), arrowprops=dict(arrowstyle='->', color='#b1272d', lw=2))",
    "ax.text(8.3, 4.6, 'rs334 -- single base A>T (cDNA)', fontsize=9, color='#b1272d')",
    "ax.text(0.3, 3.3, 'ALT cDNA codon 7:', fontsize=10, color='#57606a')",
    "for i, b in enumerate('GTG'):",
    "    fc = '#b1272d' if i == 1 else cmap[b]",
    "    alpha = 0.5 if i == 1 else 0.2",
    "    ax.add_patch(Rectangle((4.5 + i*0.7, 3.0), 0.65, 0.7, fc=fc, alpha=alpha, ec=cmap[b]))",
    "    ax.text(4.5 + i*0.7 + 0.32, 3.35, b, ha='center', va='center', family='monospace', fontsize=12, color=cmap[b], weight='bold')",
    "ax.annotate('', xy=(8.0, 3.35), xytext=(6.9, 3.35), arrowprops=dict(arrowstyle='->', color='#57606a'))",
    "ax.add_patch(FancyBboxPatch((8.2, 3.0), 2.4, 0.7, boxstyle='round,pad=0.05', fc='#f6f8fa', ec='#b1272d', lw=1.5))",
    "ax.text(9.4, 3.35, 'Val (V)', ha='center', va='center', fontsize=11, color='#b1272d', weight='bold')",
    "ax.text(11.0, 3.35, 'HbS beta-globin', va='center', fontsize=10, color='#b1272d')",
    "ax.add_patch(FancyBboxPatch((0.3, 0.4), 13.4, 1.7, boxstyle='round,pad=0.05',",
    "                            fc='#fff5f5', ec='#b1272d', lw=1.0))",
    "ax.text(0.7, 1.55, 'Phenotype:', fontsize=11, color='#b1272d', weight='bold')",
    "ax.text(0.7, 1.05,",
    "        'het (HbAS) -> sickle cell trait (carrier, malaria-protective)   |   '",
    "        'hom (HbSS) -> sickle cell anaemia',",
    "        fontsize=10, color='#1f2328')",
    "ax.set_title('rs334 -- molecular basis of sickle cell disease',",
    "             loc='left', color='#57606a', fontsize=12)",
    "plt.tight_layout(); plt.show()"
))

cells.append(md(
    "現在進入**生物學**。我們要找的就是 **rs334**：",
    "- chr11:5,248,232 (GRCh37) → 在我們的子集中即 HBB:2,232",
    "- 基因組正股上 T → A。*HBB* 基因位於反股，因此在 cDNA 上是 A → T（HGVS `c.20A>T`）",
    "- 落在第 7 密碼子（`GAG`，麩胺酸 Glu），突變後變成 `GTG`（纈胺酸 Val）。HGVS 蛋白質：`p.Glu7Val`",
    "- 極性 Glu → 疏水 Val 出現在 β-球蛋白表面，缺氧時 HbS 會聚合，紅血球變形成鐮刀狀",
    "- 雜合 = **鐮刀型基因特徵**（HbAS，多無症狀、有抗瘧保護）；同合 = **鐮刀型細胞貧血症**（HbSS）"
))

cells.append(code(
    "import subprocess",
    "vcf = subprocess.check_output(['bcftools', 'view', 'data/variants.vcf.gz'], text=True)",
    "rows = [l.split('\\t') for l in vcf.splitlines() if l and not l.startswith('#')]",
    "rs334 = next((r for r in rows if int(r[1]) == 2232), None)",
    "",
    "if rs334 is None:",
    "    print('rs334 was not called in this run.')",
    "else:",
    "    fmt  = dict(zip(rs334[8].split(':'), rs334[9].split(':')))",
    "    info = dict(kv.split('=', 1) for kv in rs334[7].split(';') if '=' in kv)",
    "    dp   = info.get('DP', '?')",
    "    ad   = fmt.get('AD', '?,?').split(',')",
    "    print('=' * 64)",
    "    print('  從 HG02666 真實讀段呼叫到的 rs334 (鐮刀型變異)')",
    "    print('=' * 64)",
    "    print(f'  位置        HBB:{rs334[1]}  =  chr11:5,248,232  (GRCh37)')",
    "    print(f'  REF -> ALT  {rs334[3]} -> {rs334[4]}      (基因組正股)')",
    "    print( '  cDNA        c.20 A>T   (HBB 在反股，反向互補後)')",
    "    print( '  蛋白質      第 7 密碼子 GAG (Glu) -> GTG (Val) = p.Glu7Val')",
    "    print(f'  基因型      {fmt[\"GT\"]}     (0/1 = 雜合 = 鐮刀型基因特徵 HbAS)')",
    "    print(f'  覆蓋深度    {dp}    (REF 讀段 = {ad[0]}, ALT 讀段 = {ad[1]})')",
    "    print(f'  品質        {rs334[5]}')",
    "    print()",
    "    print('  → 表現型: HbAS (鐮刀型基因特徵)')",
    "    print('    • 雜合帶因者，多數無症狀。')",
    "    print('    • 對惡性瘧原蟲有部分保護作用。')",
    "    print('    • 兩名 HbAS 結婚，子女 1/4 機率為 HbSS (鐮刀型細胞貧血症)。')",
    "    print()",
    "    print('  → 資料庫:')",
    "    print('    ClinVar VCV000015333 (致病).  OMIM 603903.  dbSNP rs334.')"
))

cells.append(md(
    "🎉 你剛剛從一份真實 Illumina FASTQ，一路操作到具有臨床意義的變異呼叫——使用的就是醫院遺傳實驗室在跑的同一套軟體。"
))

# ===================================================================
# 結尾
# ===================================================================
cells.append(md(
    "---",
    "## 流程總表",
    "",
    "| 步驟 | 工具 | 輸入 → 輸出 |",
    "| --- | --- | --- |",
    "| 檢視 | `head` | reads.fq |",
    "| QC | `fastqc` | qc/reads_fastqc.html |",
    "| 修剪 | `fastp` | reads.fq → trimmed.fq |",
    "| 比對 | `bwa index`, `bwa mem` | trimmed.fq + ref.fa → aligned.sam |",
    "| 排序索引 | `samtools sort/index` | aligned.sam → aligned.bam + .bai |",
    "| 堆疊 | `samtools mpileup` | aligned.bam → pileup.txt |",
    "| 呼叫 | `bcftools mpileup`/`call` | aligned.bam → variants.vcf.gz |",
    "| 判讀 | Python | variants.vcf.gz → rs334 → HbS |",
    "",
    "## 延伸練習",
    "",
    "1. **換成非帶因者**：把資料下載 cell 中的 `curl` URL 改成另一個千人計畫樣本（例如 NA12878，歐洲人，幾乎不可能帶 rs334），整條流程不變，但 VCF 不會出現 HBB:2232 的紀錄。",
    "2. **試同合子**：少數樣本是 HbSS，他們的 HBB:2232 變異會顯示 `1/1`，`DP4` 中 REF 讀段會是 0。",
    "",
    "---",
    "",
    "*筆記：**莊漢英 助理教授**・臺北醫學大學。程式碼以 MIT 授權；資料依[千人基因組計畫資料使用政策](https://www.internationalgenome.org/data-portal/data-collection)再散布。意見回饋歡迎來信 [hanyingjhuang@tmu.edu.tw](mailto:hanyingjhuang@tmu.edu.tw).*"
))

# Build & save
nb = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "name": "notebook_zh.ipynb"},
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.10"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}
with open("/tmp/ngs-repo/notebook_zh.ipynb", "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

n_md   = sum(1 for c in cells if c['cell_type']=='markdown')
n_code = sum(1 for c in cells if c['cell_type']=='code')
print(f"wrote {len(cells)} cells  ({n_md} markdown, {n_code} code)")
