# Reproduction code — *Pathway-Level Profiling Reduces Model Optimism in Chemoresistance Modeling for Colorectal Cancer Treated with Fluoropyrimidine-Based Adjuvant Chemotherapy: A Nested Cross-Validation Benchmark*

This repository holds the analysis code and the derived result tables for the
manuscript. It is deliberately **not** a data repository: all expression input is
public and is re-downloaded by the scripts named below.

## Status of deposit

This code is being released alongside a manuscript under review. The repository
is public, so the caveats below are part of the release, not internal notes.

## KNOWN ISSUES

### 1. Two of five cohorts were analysed on a linear expression scale — resolved 2026-09-24

`scripts/03_deg_analysis.R` contains the assertion

```r
# Microarray data is already log2-transformed (range ~4-16)
```

and applies no check. That is true for GSE39582, GSE104645 and GSE72970 and
**false for GSE28702 and GSE69657**, whose stored matrices have maxima of
405,930 and 26,153 respectively (measured 2026-09-22 on the deposited objects;
see `verify/scale_check.txt`). limma therefore fitted a linear-scale design for
those two cohorts, and the stored tables carry `AveExpr` of 574 and 1794 versus
6.5 and 8.0 for the log2 cohorts, with |log2FC| reaching 3,106 and 3,905.

Consequences, in order of severity:

* The `logFC` values for GSE28702 and GSE69657 are not interpretable as fold changes.
* The **p-values and directions** for those two cohorts come from a test on an
  untransformed, heteroscedastic scale dominated by the most abundant
  transcripts, so they are not comparable to the other three cohorts.
* The five-cohort Stouffer meta-analysis therefore did **not** have a defensible
  basis as deposited. **Resolved 2026-09-24:** re-running GSE28702 and GSE69657 on
  log2 scale (same limma core and same `annotate_*` group rules as
  `scripts/03_deg_analysis.R`) and applying BH FDR < 0.05 **without** the IQR probe
  filter yields **nine genes** (LDLRAD4, PRDM2, PCGF5, BAG5, APOL6, STON2, WDFY3,
  SLC25A28, KIAA1257). This nine-gene set replaces the 253-gene signature computed
  before the scale correction; the manuscript reports the nine and cites the 253
  only as the superseded pre-correction artifact.
* Within 36.7% of genes (3,895/10,618) in the meta-analysis table, at least one
  cohort contributes an |logFC| above 5.

### 2. The meta-analysis statistic reproduces, but its generating script is absent

`results_tables/meta_analysis/deg_meta_results.csv` (Supplementary Table S1) is
the five-cohort Stouffer output, but the only meta-analysis code in this
repository is Fisher's method in `scripts/03_deg_analysis.R`, writing a
different file. The Stouffer statistic itself is recoverable: re-deriving it from
the stored per-cohort p-values, with weights `sqrt(n_cohort)` and the sign taken
from `logFC`, reproduces the recorded values (LDLRAD4 z = -5.45 against a
recorded -5.446; WIPI2 -4.74 against -4.741) and yields 286 genes at
BH FDR < 0.05, of which 243 (85%) appear in the published table.

What is **not** recoverable is the input path. That reconstruction consumes the
`DEG_<cohort>_mapped.csv` tables; run instead through the pipeline in
`scripts/03_deg_analysis.R` as committed, the same five-cohort meta-analysis
yields single-digit gene counts. The difference is the probe filter described in
item 5, which the published tables evidently did not use. So treat
`meta_analysis/deg_meta_*` as **a result whose arithmetic is reproducible and
whose preparation is not**.

### 3. The per-probe DEG tables lost their identifier column

`results_tables/DEG_<cohort>_limma.csv` begin at `logFC`; the probe ID was
dropped when the tables were written (`write.csv(..., row.names = FALSE)`), so
individual rows cannot be traced back to a probe. A `Gene` column is present but,
because `scripts/03_deg_analysis.R` sets `res$Gene <- rownames(res)` on a
probe-indexed fit, it holds probe identifiers rather than symbols; the symbols
live in a separate `gene_symbol` column, which only the `_mapped.csv` files have
(23,520 distinct symbols for each GPL570 cohort, 19,565 for GPL6480).

### 4. The pathway panel is fixed and literature-defined

`results_tables/pathway_activity/` scores 44 gene sets — 26 KEGG legacy and 18
Hallmark — hard-coded in `scripts/14_ssgsea_pathway_scoring.R`. They were
enumerated before any cohort was scored and no set was added, dropped or
re-ranked using expression data here. The list, per-set gene counts and the
check that all 44 fall inside the `minSize = 10 / maxSize = 500` window are in
Supplementary Table S16.

### 5. The Methods' probe filter postdates the tables it is supposed to describe

The Methods state that an interquartile-range filter reduced the GPL570 probe set
from 54,675 to roughly 19,000 before differential testing. The filter **is**
implemented, in `analyze_deg_limma()` in `scripts/03_deg_analysis.R`:

```r
# Filter low-expression genes (v6.0 fix: IQR-based filter for microarray)
probe_iqr  <- apply(expr_mat, 1, IQR, na.rm = TRUE)
probe_med  <- apply(expr_mat, 1, median, na.rm = TRUE)
overall_med <- median(expr_mat, na.rm = TRUE)
keep <- probe_med >= overall_med & probe_iqr >= 0.5
```

but the deposited tables were not produced by it. The file dates on the analysis
project settle it: `DEG_<cohort>_limma.csv` and the 1 GB workspace image are
timestamped 2026-05-14, `meta_analysis/deg_meta_*.csv` 2026-05-17, and
`03_deg_analysis.R` was last modified 2026-06-02 — nineteen and sixteen days
later. The re-run wrapper added alongside it (`scripts/run_03_deg.R`) writes to a
scratch directory, and no filtered output was ever promoted into `results/tables`.
So this is a stale-artifact problem, not an invented method: the code describes
what the analysis should do, the shipped tables show what it did at the time.
**Resolution (2026-09-24):** the finalized manuscript discloses this IQR filter as
*considered during audit but not applied*, because it discards more than half of
all GPL570 probes, and the scale-corrected nine-gene meta-analysis result in
issue 1 is computed without it. Applying the filter would reduce the nine genes
to zero (smallest q = 0.052), which is precisely why it is not applied.

The tables are unfiltered, and demonstrably so. All four GPL570 cohorts carry
exactly **54,675** rows (GSE104645, GPL6480, carries all **41,093**), and their
mapping outputs are **identically 23,520 mapped and 8,893 unmapped rows** despite
sample counts of 164, 83, 124 and 30 — because `scripts/06_probe_to_gene_mapping.R`
reads each `DEG_<cohort>_limma.csv` and maps it onto the platform annotation, so
identical counts mean the input probe sets were identical and unfiltered.
Re-running the same limma design with the filter active gives cohort-specific
probe counts of 18,846, 23,251, 19,944 and 15,923, and collapses the five-cohort
meta-analysis from 286 genes to single digits. Regenerating the DEG and
meta-analysis tables with the committed script — which would also carry the scale
repair of items 1 and 6 — is the coherent way to close this, and it is not yet
done. This is the largest open discrepancy between the manuscript text and the
deposited artifacts.

### 6. Quantified: what the scale repair alone is worth

Keeping the code, samples, filter, mapping and Stouffer definition fixed and
switching only the log2 transform for GSE28702 and GSE69657, the five-cohort
signature moves from **4 genes to 9** (BH FDR < 0.05): 4 genes are significant
either way (BAG5, LDLRAD4, PCGF5, PRDM2), 5 are added, none are lost, Jaccard
0.444. The input-scale defect is therefore real but mild in effect, and it is not
the origin of the reported signature. See `verify/` and the note in item 5.

## Data sources

Public, re-downloadable; nothing proprietary or patient-identifiable is
included.

| Accession | Platform | Role |
|---|---|---|
| GSE39582 | GPL570 | discovery (adjuvant, fluoropyrimidine-based) |
| GSE104645 | GPL6480 | validation (restricted to oxaliplatin-containing regimens) |
| GSE28702 | GPL570 | validation (mFOLFOX6) |
| GSE72970 | GPL570 | validation (regimen unrestricted; 71% irinotecan-containing) |
| GSE69657 | GPL570 | validation (neoadjuvant FOLFOX4) |
| GSE156915 | GPL29069 | CMS / DDIR / MSI / mutation annotation |
| GSE132465 | — | single-cell validation |
| TCGA-COADREAD | RNA-seq | external application of the frozen nomogram |
| MSigDB v2024.1 | h.all, c2.cp.kegg_legacy | gene sets |

## Environment

R 4.6.0 and Python 3.14. Package versions actually used are recorded in
`MANIFEST_gdc_download.txt` (GDC inputs) and in `verify/session_info.txt`.

`limma`, `edgeR`, `Biobase` and `GSVA` are **not vendored here** and were not
present in the R library at the time of this deposit, which is why
reproducing the DEG layer requires installing them first (see `install_deps.R`).

## Layout

```
scripts/            numbered analysis scripts, run in order; 00_run_all.R orchestrates
  01-07b            data intake, DEG, WGCNA, candidate pool, enrichment
  09-12             gene-level ensemble fingerprint and its validation
  13                ComBat batch correction (QC + fingerprint transfer only)
  14-19             ssGSEA panel, pathway modelling, meta-analysis enrichment
  21-39             PRS, nomogram, nested cross-validation, external validation
results_tables/     derived tables only (no clinical contact data, no matrices)
verify/             checks backing the KNOWN ISSUES above
```

## What is deliberately excluded

* raw and processed expression data, and the GDC tree — re-download instead;
* `results_tables/*_clinical_data.csv` — these were parsed from GEO series
  matrices and carry the original submitters' names and e-mail addresses;
* `scripts/node_modules/` and `__pycache__/`;
* the 170 MB combined expression matrix.

TCGA participant barcodes are retained: they are pseudonymous, already public,
and required to reproduce the TCGA analyses.

Machine-specific absolute paths were replaced by `/path/to/xelox_project` and
`/home/user` before release, so set `PROJECT_ROOT` accordingly before running.

## v2.1 (2026-09-24): the corrected run is now part of the deposit

The repository previously carried only the May-2026 pipeline, so a reviewer could
not reproduce the gene-level numbers the manuscript now reports. This revision
publishes the corrected analysis alongside the original:

| path | what it is |
|---|---|
| `scripts/reanalysis/rerun_deg3.R` | per-cohort limma re-run with GSE28702 and GSE69657 log2(x+1)-transformed, `annotate_*` group definitions copied verbatim from `02_extract_xelox_groups.R` |
| `scripts/reanalysis/meta_final.R`, `meta_export.R`, `meta_calibrate.R`, `sig9_table.R` | Stouffer recombination weighted by sqrt(n analysed per cohort), BH correction, SI Table S1 export, nine-gene table |
| `results_tables/reanalysis_2026-09/deg_meta_corrected_noIQRfilter.csv` | **the table behind the manuscript**: 9 genes at `stouffer_fdr < 0.05` (LDLRAD4, PRDM2, PCGF5, BAG5, APOL6, STON2, WDFY3, SLC25A28, KIAA1257) |
| `results_tables/reanalysis_2026-09/deg_meta_corrected_IQRfiltered.csv` | same run with the IQR probe filter additionally applied: 0 genes at FDR < 0.05, smallest q = 0.052 |
| `results_tables/reanalysis_2026-09/s1_*.csv`, `signature_9_percohort.csv` | Supplementary Table S1 direction matrix, per-cohort effect sizes and the nine-gene per-cohort table |
| `results_tables/reanalysis_2026-09/agg_lasso_nested.csv`, `complexity_matched_nested.csv`, `perfold_*.csv` | nested cross-validation and complexity-matched comparison behind 0.0704 vs 0.3671, the 81 % figure and Figure 5 |
| `results_tables/reanalysis_2026-09/external_validation_TCGA.csv` | frozen-nomogram TCGA validation (C-index 0.451, calibration slope -0.450) |
| `results_tables/reanalysis_2026-09/*_log.txt` | run logs of the above, with input sizes and timestamps |
| `results_tables/_superseded/deg_meta_results_pre_log2fix.csv` | the pre-correction meta-analysis (253 genes); referenced in the manuscript only as the superseded count |

`scripts/03_deg_analysis.R` keeps its original behaviour and carries a STATUS
note saying which two defects make it inconsistent with the manuscript, so the
historical pipeline and the reported numbers stay distinguishable.

Machine-specific absolute paths in this revision were again normalised to
`/path/to/xelox_project`, `/path/to/xelox_reanalysis`, `/path/to/Rlibs`, `/tmp`
and `/home/user`; result tables were not rewritten.
