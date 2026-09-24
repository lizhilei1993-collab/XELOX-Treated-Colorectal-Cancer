deg_meta_results_pre_log2fix.csv is the pre-correction meta-analysis: GSE28702 and
GSE69657 were entered on their raw linear scale, which is why it reports 253 genes at
BH FDR < 0.05 (stouffer_fdr). Those 253 genes are used nowhere in the manuscript.
The corrected tables are in ../reanalysis_2026-09/: deg_meta_corrected_noIQRfilter.csv
is the one the paper reports (9 genes at FDR < 0.05); deg_meta_corrected_IQRfiltered.csv
shows the effect of additionally applying the IQR probe filter (0 genes, smallest
q = 0.052), which the Methods state was considered and not applied.
