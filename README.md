## LoopSeq analysis repo

**Expected install time:** < 1 hr for code, 1-2 hr for dataset.  
**Expected computing time:** < 2 hr when samples are run in parallel.

This repository contains scripts + reference FASTAs used to analyze LoopSeq long-read FASTA files, summarize results across samples, and generate the plots committed in this repo.

### What’s in this repo

- **`scripts/analyze_aa_mutations.py`**: main per-sample analysis (filters reads, anchors by gene-start, matches to reference, calls AA mutations, writes `sample_results.xlsx` + plots).
- **`scripts/enumerate_sample_keys.py`**: reads `sample_keys.xlsx` and prints one TSV line per sample describing which FASTA to analyze and which reference/gene-start to use.
- **`scripts/build_all_summary.py`**: aggregates all `aa_analysis_results_*/*/sample_results.xlsx` into `all_summary.xlsx`.
- **`reference_fastas/`**: reference sequences (`sp23_*`, `sp24_*`) and gene-start FASTAs.
- **`sample_keys.xlsx`**: workbook that maps samples to reference type and drives batch enumeration.
- **`aa_analysis_results_*/`**: per-sample output directories (one per analyzed sample).

### Data organization (expected layout)

This repo is designed to live at a project root that also contains your FASTA data (often large and typically ignored by git).

Minimal structure:

```text
all_files/
  reference_fastas/
    sp23_dna.fa
    sp23_gene_start.fa
    sp24_dna.fa
    sp24_gene_start.fa
  sample_keys.xlsx
  scripts/
    analyze_aa_mutations.py
    enumerate_sample_keys.py
    build_all_summary.py
  aa_analysis_results_<sample_id>_cluster_1based/
    sample_results.xlsx
    mutation_table.csv
    reference_used.txt
    mutation_set_distribution.pdf
    dna_mutation_distance_distribution.pdf
    observed_vs_expected_frequencies.pdf
    ...
  samples_1-24/            # optional: large input tree (ignored in git)
  samples_25-48/           # optional: large input tree (ignored in git)
  initial_set/             # optional: large input tree (ignored in git)
```

### How sample IDs map to output folders

Per-sample results are written to directories named like:

```text
aa_analysis_results_<sample_id>_cluster_1based/
```

Examples:
- `aa_analysis_results_1-48_C4_n27_cluster_1based/`
- `aa_analysis_results_initial_set_r2_A1_cluster_1based/`

The `<sample_id>` string is also used in `all_summary.xlsx`.

### How to prepare FASTA inputs

`scripts/analyze_aa_mutations.py` expects:
- **Input reads** in FASTA format (`--input_fasta`)
- A **reference FASTA** (`--reference_fasta`) matching the sample’s intended template (sp23 or sp24)
- A **gene-start sequence** (`--gene_start`) used to anchor the subread
- A **gene length** (`--gene_length`) for the extracted region

You can organize FASTA inputs anywhere on disk, but the default workflow assumes:
- For samples 1–24: input FASTAs live under `samples_1-24/**/...trimmed.fa`
- For samples 25–48: input FASTAs live under `samples_25-48/**/...trimmed.fa`
- For “initial_set”: input FASTAs live under `initial_set/**/..._trimmed.fa`

Those conventions are implemented in `scripts/enumerate_sample_keys.py`. If your FASTA naming differs, either:
- update that script to match your filesystem, or
- bypass it and run `analyze_aa_mutations.py` directly (example below).

### Step 1 — Enumerate samples from `sample_keys.xlsx`

This prints one TSV row per sample with:
`sample_name  input_fasta  reference_fasta  gene_length  gene_start`

```bash
python3 scripts/enumerate_sample_keys.py /path/to/all_files > samples.tsv
```

### Step 2 — Run the per-sample analysis

Run a single sample directly:

```bash
python3 scripts/analyze_aa_mutations.py \
  --input_fasta "/abs/path/to/sample_trimmed.fa" \
  --reference_fasta "reference_fastas/sp23_dna.fa" \
  --gene_start "$(python3 - <<'PY'
from pathlib import Path
p=Path('reference_fastas/sp23_gene_start.fa')
seq=[]
for line in p.read_text().splitlines():
    if line.startswith('>'): continue
    seq.append(line.strip())
print(''.join(seq).upper())
PY
)" \
  --gene_length 1233 \
  --output_dir "aa_analysis_results_example_cluster_1based"
```

Batch-running multiple samples:
- Use `samples.tsv` as the source of truth.
- For each row, pass the corresponding fields to `analyze_aa_mutations.py`.
- (This repo also includes `scripts/submit_batch_analysis.sh` which can be adapted for your scheduler.)

### Per-sample outputs you can expect

Inside each `aa_analysis_results_*` directory, the analysis writes:
- **`sample_results.xlsx`**
  - `mutation_set_categories`
  - `plotted_mutation_sets`
  - `dna_distance_distribution`
- **Plots**
  - `mutation_set_distribution.(pdf|png)`
  - `dna_mutation_distance_distribution.(pdf|png)`
  - `observed_vs_expected_frequencies.(pdf|png)`
- **Tables**
  - `mutation_table.csv`
- **Metadata**
  - `reference_used.txt`

### Step 3 — Build `all_summary.xlsx`

Aggregate all per-sample summaries:

```bash
python3 scripts/build_all_summary.py /path/to/all_files
```

This writes:
- `all_summary.xlsx` → sheet `all_samples`
- `all_summary.xlsx` → sheet `plotted_reads_gt_1000` (subset with `total_plotted_reads > 1000`)

The summary workbook also includes metadata columns (e.g. `timepoint`, `replicate`, `Starting point`) used by the plotting scripts.

### Plotting scripts committed in this repo

- **3D mutation-set distributions** (grouped by replicate + starting point):

```bash
python3 scripts/plot_grouped_3d_mutation_set_distributions.py --base_dir /path/to/all_files
```

Outputs in `plots_3d_mutation_set_distributions/`.

- **3D mutation-distance vs time**:

```bash
python3 scripts/plot_3d_mutation_distance_vs_time.py --base_dir /path/to/all_files
```

Outputs in `plots_3d_mutation_distance/`.

- **2D jitter plot of mutation distance vs time**:

```bash
python3 scripts/plot_violin_mutation_distance_vs_time.py --base_dir /path/to/all_files
```

Outputs in `plots_mutation_distance_violin/`.

### Demo (downsampled FASTA + expected output)

This repo includes a small, downsampled demo FASTA and the corresponding expected analysis outputs so you can quickly verify your setup without processing the full (large) input FASTA files.

Downsampling logic (for reproducibility):
- We start from the full `replicate = 1`, `Starting point = WT/E` input FASTA for cluster `1-48_C4_n27` (`L1ba3f0d1_sample_BROAD_27_contig_list_trimmed.fa`).
- From the full-sample workbook `aa_analysis_results_1-48_C4_n27_cluster_1based/sample_results.xlsx`, sheet `plotted_mutation_sets`, we take genotypes with `read_count > 100`, sort by `read_count` descending, and keep the top **10** labels (fewer if not enough genotypes pass the filter).
- We keep **every** read that maps to the same mutation-set label as one of those chosen genotypes (same alignment and mutation-calling logic as `scripts/analyze_aa_mutations.py`).
- For this sample the plotted sheet lists only **five** genotypes, all with `read_count > 100`, so all five are included and the demo FASTA contains every mapped read for those five (no further subsampling).

The chosen genotype labels are listed in `demo/rep1_startWT_E_demo/demo_selected_genotypes.txt`.

#### Demo assets included

The demo lives under:

```text
demo/rep1_startWT_E_demo/
  input_demo_plotted_genotypes_all_reads.fa
  demo_selected_genotypes.txt
  output/
    sample_results.xlsx
    mutation_set_distribution.pdf
    dna_mutation_distance_distribution.pdf
    observed_vs_expected_frequencies.pdf
    mutation_table.csv
    summary_statistics.txt
    reference_used.txt
```

The expected output files above were produced with:
- `replicate = 1`
- `Starting point = WT/E`
- `reference_fasta = reference_fastas/sp24_dna.fa`

#### Run the demo (recompute outputs)

1. Extract `gene_start` from the matching reference gene-start FASTA.
2. Re-run the analysis on the included downsampled demo FASTA.

```bash
GENE_START="$(python3 - <<'PY'
from pathlib import Path
p = Path('reference_fastas/sp24_gene_start.fa')
seq = []
in_seq = False
for line in p.read_text().splitlines():
    line = line.strip()
    if not line:
        continue
    if line.startswith('>'):
        if seq:
            break
        in_seq = True
        continue
    if in_seq:
        seq.append(line)
print(''.join(seq).upper())
PY
)"

python3 scripts/analyze_aa_mutations.py \
  --input_fasta demo/rep1_startWT_E_demo/input_demo_plotted_genotypes_all_reads.fa \
  --reference_fasta reference_fastas/sp24_dna.fa \
  --gene_start "$GENE_START" \
  --gene_length 1233 \
  --output_dir demo/rep1_startWT_E_demo/output_rerun \
  --sample_keys sample_keys.xlsx
```

#### Compare your outputs to the expected output

Fast check: compare SHA-256 hashes of the expected files.

```bash
sha256sum demo/rep1_startWT_E_demo/output_rerun/*.pdf demo/rep1_startWT_E_demo/output_rerun/*.xlsx \
  demo/rep1_startWT_E_demo/output_rerun/*.csv demo/rep1_startWT_E_demo/output_rerun/*.txt \
  2>/dev/null || true
```

Expected hashes for the shipped outputs:

```text
3c868b1b3330e779bb7eef26deb15ec30209256a49fb0d1c5ab24e2fc313d8ef  demo/rep1_startWT_E_demo/output/sample_results.xlsx
13100b253c4fc9578f85ddd4023e00109a30bfa577b51f3beef54ba8f2103fca  demo/rep1_startWT_E_demo/output/mutation_set_distribution.pdf
96da2c7608b6182b2832c5221ae219aa4e4555d169b5c3a74985bcbd2977a580  demo/rep1_startWT_E_demo/output/dna_mutation_distance_distribution.pdf
b19797f0d2e01264ae7c2cc196970af988cd0299ab3d9dea9da44dfbbcf9b298  demo/rep1_startWT_E_demo/output/observed_vs_expected_frequencies.pdf
fec515483786529b4721f20216b322234d80a27f7f9fc771f791b84dc3e2fc1e  demo/rep1_startWT_E_demo/output/mutation_table.csv
8e5ff1ccbd05caab95b04279686a7c2e39f476ef018dfdc6ca3807f478d585c5  demo/rep1_startWT_E_demo/output/mutation_pairs_raw_data.csv
6dd08eb3bb5bc716316df260431660e9858dc74144d39fe3e5aa27c9f98827df  demo/rep1_startWT_E_demo/output/summary_statistics.txt
37379bd8331fb75cfa290b20837e300f79430837f60cde2be875cc04c0b4eb1a  demo/rep1_startWT_E_demo/output/reference_used.txt
```

If hashes differ (e.g. due to matplotlib/OS/font differences), compare `sample_results.xlsx` numerically:
- sheet `mutation_set_categories` (key summary columns)
- sheet `plotted_mutation_sets` (mutation set → read count)

### Notes on git / large files

This project commonly includes very large input directories (`samples_1-24/`, `samples_25-48/`, `initial_set/`) and many derived outputs. `.gitignore` is set up to keep large raw inputs and obvious output artifacts out of version control, but in this repo some plot outputs are committed intentionally.

