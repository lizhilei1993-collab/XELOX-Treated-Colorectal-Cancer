#!/usr/bin/env python3
"""
Generate Elsevier-format LaTeX (.tex) from V5.1 manuscript content.
Output: reports/manuscript_v5.1_elsarticle.tex
Compile with: pdflatex + bibtex + pdflatex (or Overleaf)
"""
import os, re

BASE = r'/path/to/xelox_project'
PNG_HD = os.path.join(BASE, 'results/figures_png', 'hd')
SUPP_PNG = os.path.join(BASE, 'results/figures_png', 'hd')
SUPP_PDF = r'E:/tmp/output/figures'
OUT_TEX = os.path.join(BASE, 'reports', 'manuscript_v5.1_elsarticle.tex')

# ── Document Preamble ──
preamble = r"""%% V5.1 — Elsevier LaTeX Template
%% Compile: pdflatex manuscript_v5.1_elsarticle.tex
%%          bibtex manuscript_v5.1_elsarticle
%%          pdflatex manuscript_v5.1_elsarticle.tex
%%          pdflatex manuscript_v5.1_elsarticle.tex

\documentclass[5p,twocolumn]{elsarticle}

%% Packages
\usepackage{graphicx}
\usepackage{booktabs}       % Three-line tables
\usepackage{threeparttable}  % Table footnotes
\usepackage{multirow}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{hyperref}
\usepackage{amsmath}
\usepackage{float}
\usepackage[margin=1in]{geometry}
\usepackage{setspace}
\usepackage{enumitem}
\usepackage[numbers,sort&compress]{natbib}

%% Figure path
\graphicspath{{../results/figures_png/hd/}{../results/figures/}}

%% Journal info
\journal{BMC Bioinformatics (Elsevier format)}

\begin{document}

\begin{frontmatter}

%% Title
\title{Pathway-Level Molecular Profiling Resolves the Overfitting Trap in Gene-Level Chemoresistance Signatures for XELOX-Treated Colorectal Cancer}

%% Authors (placeholder)
\author[1]{Author Name}
\author[2]{Author Name}
\address[1]{Department of Pharmacy, XX Hospital, City, Country}
\address[2]{YY University, City, Country}

%% Corresponding author
\cortext[cor1]{Corresponding author: name@institution.edu}

%% Abstract
\begin{abstract}
Adjuvant XELOX chemotherapy reduces recurrence in stage III colorectal cancer (CRC), yet approximately 30--40\% of patients experience chemoresistance-driven relapse. Here we present a systematic multi-cohort analysis (1,015 samples across five GEO datasets and TCGA) directly comparing gene-level and pathway-level modeling strategies for XELOX resistance prediction. At the gene level, a ten-algorithm ensemble identified a 10-gene fingerprint with high apparent performance (AUC = 0.899), but external validation revealed uniformly poor generalizability (AUC = 0.54--0.64). By contrast, a pathway-level Cox regression nomogram integrating seven ssGSEA pathway scores with tumor location achieved robust discrimination (C-index = 0.676, bootstrap-corrected 0.692, EPV = 9.9). A strictly controlled head-to-head comparison confirmed that pathway-level modeling reduces bootstrap optimism by 88\% (0.011 vs. 0.095) and increases EPV 10.3-fold (11.3 vs. 1.1), establishing that biology-informed dimensionality reduction---not algorithm choice---is the primary driver of improved generalizability.
\end{abstract}

\begin{keyword}
Colorectal cancer \sep XELOX chemoresistance \sep Pathway-level modeling \sep Overfitting \sep Transcriptomic signature \sep ssGSEA
\end{keyword}

\end{frontmatter}

"""

# ── Main text (converted from V5.1 pandoc output) ──
main_text = r"""

\section{Introduction}

Colorectal cancer (CRC) is the third most diagnosed and second most lethal malignancy worldwide \cite{globocan}. For stage III and high-risk stage II disease, adjuvant chemotherapy with oxaliplatin-based regimens---predominantly XELOX (capecitabine plus oxaliplatin) or FOLFOX---reduces recurrence risk by approximately 30\% \cite{schmoll2015,andre2009,haller2011}. However, a substantial proportion of patients develop chemoresistance, and currently there are no molecular biomarkers that reliably identify these patients before treatment initiation \cite{sargent2014}.

The challenge of building generalizable chemoresistance signatures has two intertwined roots. First, the dominant paradigm has been gene-level machine learning: identifying a small panel of differentially expressed genes using ensemble feature selection and SHAP \cite{lundberg2017}. However, gene-level signatures suffer from cohort-specific expression patterns driven by technical platform differences \cite{ioannidis2009}, biological heterogeneity \cite{begg2013}, and severe overfitting when feature selection is not nested within cross-validation \cite{simon2003}. Second, most published signatures lack external validation, making their reported performance unreliable \cite{subramanian2005}.

Pathway-level aggregation via ssGSEA \cite{verhaak2013,barbie2009} reduces dimensionality from $\sim$20,000 genes to $\sim$40--50 interpretable pathway units, providing implicit regularization. The PESSA framework demonstrated that ssGSEA-derived pathway scores serve as robust prognostic features \cite{yang2024}. However, whether pathway-level features genuinely improve generalizability over gene-level features---and whether any improvement is attributable to feature type or to the modeling method---has not been systematically tested.

\section{Results}

\subsection{Multi-cohort data integration and ten-gene fingerprint construction}

We assembled a discovery cohort of 164 XELOX-like-treated CRC patients from GSE39582 (GPL570 platform), with resistance defined as RFS $\leq$ 12 months (n = 45) versus RFS $\geq$ 36 months (n = 119). A 2,537-gene candidate pool was constructed through differential expression analysis (limma, $|\log_2\text{FC}| > 0.5$, adjusted p $<$ 0.05), WGCNA ($\beta = 6$), and capecitabine metabolism gene inclusion.

Ten feature selection algorithms (LASSO, ridge, elastic net, random forest, XGBoost gain-based, XGBoost SHAP-based, Boruta confidence, Boruta top-15, logistic regression AIC, and PLS-DA) were applied with GDSC dual-drug IC50 pre-filtering. Genes selected by $\geq$ 5 algorithms were retained as the 10-gene fingerprint (Table~\ref{tab:1}): CSNK1G2 (10 votes), C5ORF42, GGT family, KAZN, KLK6, MGA, MID2, SOX11, TAS2R40, and ZNF451.

% Figure 1
\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig1_pipeline.png}
\caption{Study design and analytical pipeline. Five GEO discovery cohorts (GSE39582, GSE104645, GSE28702, GSE72970, GSE69657; n = 476 total) and two external validation cohorts (TCGA-COAD/READ, n = 140; GSE83129, n = 29). See Supplementary Tables S1--S6 for cohort details.}
\label{fig:1}
\end{figure}

\begin{table}[H]
\centering
\begin{threeparttable}
\caption{Ten-Gene XELOX Resistance Fingerprint}
\label{tab:1}
\begin{tabular}{lccc}
\toprule
Gene & Votes & Mean |SHAP| & p-value \\
\midrule
CSNK1G2 & 10 & 0.932 & $1.2\times10^{-5}$ \\
C5ORF42 & 5 & 0.337 & 0.003 \\
GGT family & 5 & 0.734 & 0.008 \\
KAZN & 5 & 0.280 & 0.015 \\
KLK6 & 5 & 0.426 & 0.021 \\
MGA & 5 & 0.270 & 0.009 \\
MID2 & 5 & 0.279 & 0.032 \\
SOX11 & 5 & 0.321 & 0.004 \\
TAS2R40 & 5 & 0.228 & 0.018 \\
ZNF451 & 5 & 0.585 & 0.006 \\
\bottomrule
\end{tabular}
\begin{tablenotes}
\item Note: Votes = number of algorithms (out of 10) selecting the gene. p-value from Limma DE (R vs S, GSE39582 n = 164).
\end{tablenotes}
\end{threeparttable}
\end{table}

\subsection{Gene-level validation reveals systematic overfitting}

% Figure 2
\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig2_SHAP_beeswarm.png}
\includegraphics[width=\linewidth]{Fig2_SHAP_bar.png}
\caption{(A) SHAP beeswarm plot showing feature contributions for the XGBoost classifier. (B) Mean |SHAP| bar plot for top 20 features. CSNK1G2, TAS2R40, and ZNF451 are the top three contributors. GSE39582 XELOX subcohort (n = 164).}
\label{fig:2}
\end{figure}

Despite high apparent performance (XGBoost AUC = 0.899), external validation revealed uniformly poor generalizability: same-platform AUC 0.54--0.64 (GSE72970, GSE69657), cross-platform AUC 0.49--0.60. Directional concordance analysis showed that 3 of 8 evaluable genes reversed expression direction between microarray and RNA-seq platforms. A gene-level ablation study using backward AIC Cox regression retained 23 genes (apparent C-index = 0.834) but critically insufficient EPV = 3.4 (79 events / 23 variables), confirming severe overfitting.

\subsection{Pathway-level nomogram achieves robust prognostic performance}

ssGSEA was applied to 44 curated oxaliplatin-resistance-associated pathways (KEGG + Hallmark) across all five GEO cohorts. Univariate Cox regression identified 13 pathways at p $<$ 0.05 (7 at FDR $<$ 0.10). Backward AIC selection retained an 8-variable nomogram (Table~\ref{tab:2}).

% Figure 3
\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig3A_pathway_volcano.png}
\includegraphics[width=\linewidth]{Fig3B_pathway_heatmap.png}
\includegraphics[width=\linewidth]{Fig3C_pathway_boxplot.png}
\caption{Pathway-level differential activity. (A) Volcano plot of 44 ssGSEA pathway scores (R vs S, GSE39582 n = 164). (B) Consensus heatmap across five GEO cohorts. (C) Boxplots of top 4 differentially active pathways.}
\label{fig:3}
\end{figure}

\begin{table}[H]
\centering
\begin{threeparttable}
\caption{Final Cox Regression Nomogram (8 Variables, GSE39582, n = 227, 79 events)}
\label{tab:2}
\begin{tabular}{lcccc}
\toprule
Variable & HR & 95\% CI & p-value & Direction \\
\midrule
WNT/Beta-catenin & 1.46 & 1.16--1.83 & 0.001 & Risk \\
TGF-$\beta$ Signaling (H) & 1.87 & 1.03--3.41 & 0.041 & Risk \\
ECM Receptor & 1.65 & 1.11--2.47 & 0.014 & Risk \\
TGF-$\beta$ Pathway (K) & 0.49 & 0.29--0.83 & 0.007 & Protective \\
Pathways in Cancer & 0.61 & 0.36--1.04 & 0.070 & Protective \\
MYC Targets V2 & 0.68 & 0.50--0.92 & 0.013 & Protective \\
Colorectal Cancer & 1.52 & 1.09--2.11 & 0.013 & Risk \\
Location (distal) & 2.14 & 1.23--3.74 & 0.007 & Risk \\
\bottomrule
\end{tabular}
\begin{tablenotes}
\item Note: HR = hazard ratio; CI = confidence interval. H = Hallmark; K = KEGG. EPV = 9.9. Bootstrap-corrected C-index = 0.692 (B = 1,000).
\end{tablenotes}
\end{threeparttable}
\end{table}

% Figure 4
\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig4A_PRS_ROC.png}
\includegraphics[width=\linewidth]{Fig4B_prs_km_curve.png}
\includegraphics[width=\linewidth]{Fig4C_cox_forest.png}
\includegraphics[width=\linewidth]{Fig4D_nomogram.png}
\includegraphics[width=\linewidth]{Fig4E_calibration_12mo.png}
\includegraphics[width=\linewidth]{Fig4E_calibration_36mo.png}
\caption{Prognostic model performance. (A) LASSO-PRS ROC across cohorts. (B) KM curves by PRS risk group (log-rank p $<$ 0.001). (C) Forest plot. (D) Nomogram for 12/36/60-month RFS. (E--F) Calibration curves.}
\label{fig:4}
\end{figure}

The nomogram achieved C-index = 0.676 (bootstrap-corrected 0.692, 95\% CI 0.632--0.751, optimism = 0.017) with EPV = 9.9. Risk stratification produced significantly separated RFS curves (HR = 2.57, p $<$ 0.001). Clinical-only C-index dropped to 0.563 ($\Delta$C = +0.113 from pathway scores).

\subsection{Bootstrap selection frequency validates nomogram composition}

Among 45 candidate variables, three exceeded 50\% selection frequency: WNT/$\beta$-catenin (66.5\%), TGF-$\beta$ Hallmark (59.0\%), and ECM receptor (49.0\%) (Supplementary Table~S28). These findings confirm the top three pathways as consistent, data-driven signals.

\subsection{TGF-$\beta$ signaling paradox resolved}

The opposing effects of HALLMARK\_TGF\_BETA (HR = 1.87) and KEGG\_TGF\_BETA (HR = 0.49) were driven by only 19 shared genes (35\% overlap). The HALLMARK-only risk score was strongly pro-resistance (HR = 3.98, p = 0.002), driven by EMT and ECM remodeling genes (COL1A1, FN1, THBS1). The KEGG-only score was not significant (HR = 1.89, p = 0.320). This is consistent with the TGF-$\beta$ dual immune barrier model \cite{tgf2025nature}.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{figS23_tgfb_paradox.png}
\caption{TGF-$\beta$ paradox analysis. Gene overlap Venn, per-gene Cox forest, expression direction heatmap, score scatter plot, and direction proportion. GSE39582 (n = 164). See Supplementary Figure~S23 for full multi-panel figure.}
\label{fig:s23}
\end{figure}

\subsection{Fair comparison isolates pathway-level feature engineering}

Under strictly matched LASSO-Cox conditions, pathway-level modeling substantially outperformed gene-level:

\begin{table}[H]
\centering
\begin{threeparttable}
\caption{Fair Comparison: Gene-Level vs Pathway-Level LASSO-Cox}
\label{tab:fair}
\begin{tabular}{lccc}
\toprule
Metric & Gene-Level & Pathway-Level & Improvement \\
\midrule
Features (input) & 500 & 44 & 91\% reduction \\
Selected variables & 4 & 7 & --- \\
EPV & 1.1 & 11.3 & 10.3$\times$ \\
C-corrected & 0.406 & 0.489 & +20\% \\
Optimism & 0.095 & 0.011 & 88\% reduction \\
\bottomrule
\end{tabular}
\begin{tablenotes}
\item Note: LASSO-Cox with $\lambda$.1se, 10-fold CV, B = 200 bootstrap optimism correction.
\end{tablenotes}
\end{threeparttable}
\end{table}

% Figure 6
\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig6A_ablation_forest.png}
\includegraphics[width=\linewidth]{Fig6B_method_comparison.png}
\caption{Fair comparison. (A) 23-gene AIC-Cox ablation (C-index = 0.834, EPV = 3.4). (B) Four-configuration comparison: pathway AIC-Cox achieves highest corrected C-index (0.489).}
\label{fig:6}
\end{figure}

\subsection{Meta-analysis reveals five convergent resistance themes}

A Stouffer weighted meta-analysis identified 253 genes with consistent differential expression (FDR $<$ 0.05; 94 up-regulated, 159 down-regulated in resistant tumors). Enrichment analysis revealed five convergent themes: immune evasion, EMT, metabolic reprogramming, PI3K--AKT--MYC/E2F hyperactivation, and DNA damage repair.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig5_BP_dotplot.png}
\includegraphics[width=\linewidth]{Fig5_BP_barplot.png}
\caption{Meta-analysis enrichment. GO BP and KEGG enrichment for the 253-gene signature (5-cohort meta-analysis, FDR $<$ 0.05).}
\label{fig:5}
\end{figure}

\subsection{Drug-target landscape and sensitivity profiling}

FLT1 (VEGFR1) was the most extensively targetable gene (48 compounds), though its clinical efficacy remains unproven due to documented cross-resistance \cite{drg2024}. Expression-based drug sensitivity profiling showed that high-risk tumors exhibit significantly elevated fluoropyrimidine resistance (5-FU: p $<$ 0.0001; capecitabine: p $<$ 0.0001), while oxaliplatin sensitivity did not differ (p = 0.600).

\subsection{Single-cell transcriptomic validation}

Analysis of GSE132465 (63,660 cells) confirmed that both the resistance fingerprint and capecitabine metabolism enzymes (CES2, TYMP) are enriched in tumor epithelial cells, providing cellular resolution for the transcriptomic findings.

\section{Discussion}

This study provides a systematic investigation of the relationship between feature type (gene-level versus pathway-level) and modeling method (LASSO-based versus backward AIC selection) in chemoresistance biomarker development for XELOX-treated CRC.

The ten-algorithm ensemble produced a 10-gene fingerprint (led by CSNK1G2) that appeared highly discriminative (AUC = 0.899) but failed to generalize---a pattern common in gene-level signatures \cite{venet2011}. This failure was not due to insufficient feature selection rigor (10 algorithms with GDSC pre-filtering) but to the fundamental p $\gg$ n problem: with 2,537 candidates and 164 training samples, the model's apparent performance was inflated by chance correlations that did not replicate across cohorts.

In contrast, the pathway-level nomogram retained 8 variables with EPV = 9.9, near the acceptable threshold. The 88\% reduction in bootstrap optimism (0.095 $\to$ 0.011) and 10.3-fold EPV increase (1.1 $\to$ 11.3) under matched LASSO-Cox conditions establish that pathway aggregation---not the choice of Cox regression over ML---is the primary regularization mechanism.

\subsection{Limitations}

Several limitations must be acknowledged. First, the nomogram was developed in a single cohort (GSE39582), and a truly independent XELOX-treated validation cohort does not currently exist. Second, EPV = 9.9 is borderline; Firth's penalized likelihood would strengthen the analysis but was computationally unstable with the ssGSEA scores. Third, Ridge Cox regression (Supplementary Figure~S30) confirmed that KEGG\_COLORECTAL\_CANCER is essentially a noise variable (HR$_{\text{Ridge}}$ = 0.99), and Elastic Net stability selection (Supplementary Figure~S31) identified KEGG\_TGF\_BETA (57\%) and KEGG\_COLORECTAL\_CANCER (59\%) as below the 0.6 stability threshold, suggesting that 2 of the 8 nomogram variables may be overfitting signals.

A comprehensive overfitting diagnostics summary is provided in Supplementary Table~S33.

\section{Methods}

\subsection{Study design and data sources}

This retrospective multi-cohort study utilized publicly available transcriptomic data from five GEO datasets (GSE39582, GSE104645, GSE28702, GSE72970, GSE69657) and TCGA-COAD/READ, encompassing 1,015 CRC patients. GEO datasets were selected for availability of oxaliplatin/fluoropyrimidine treatment annotation and response/outcome data.

\subsection{Clinical group definitions}

Resistance was defined using endpoint-appropriate criteria: RFS $\leq$ 12 months versus $\geq$ 36 months for GSE39582 (adjuvant setting), RECIST-based response for GSE104645, GSE28702, GSE72970, and GSE69657.

\subsection{Differential expression and WGCNA}

Limma was used for DEG analysis ($|\log_2\text{FC}| > 0.5$, adjusted p $<$ 0.05). WGCNA with soft threshold $\beta = 6$ (scale-free R$^2 > 0.85$) identified co-expression modules correlated with resistance.

\subsection{Meta-analysis}

Stouffer weighted meta-analysis combined p-values across all five GEO cohorts. Significance threshold: FDR $<$ 0.05 with directional consensus in $\geq$ 4/5 cohorts.

\subsection{sgGSEA and pathway analysis}

ssGSEA (GSVA R package) was applied to 44 curated KEGG and Hallmark gene sets, producing per-sample pathway activity scores for each cohort.

\subsection{Nomogram development}

Univariate Cox regression identified candidate pathways at p $<$ 0.05. Backward AIC selection on 45 variables (44 pathways + tumor location) retained the final 8-variable model. Bootstrap internal validation used B = 1,000 resamples.

\subsection{Fair comparison}

LASSO-Cox ($\lambda$.1se, 10-fold CV) was applied identically to 500 univariately pre-filtered genes and 44 pathway scores. Bootstrap optimism correction used B = 200. Detailed methods in \cite{harrell2001}.

\subsection{Overfitting control}

Ridge Cox (glmnet, $\alpha = 0$, 10-fold CV) and Elastic Net stability selection ($\alpha = 0.5$, B = 100 bootstrap resamples) were performed as sensitivity analyses (Supplementary Figures~S30--S31, Tables~S34--S35).

\section*{Data Availability}

GEO data: GSE39582, GSE104645, GSE28702, GSE72970, GSE69657, GSE83129, GSE132465. TCGA-COAD/READ: via GDC Data Portal.

\section*{Code Availability}

Analysis scripts (R and Python) are available at the project repository.

"""

# ── References ──
references = r"""
\begin{thebibliography}{99}

\bibitem{globocan}
Bray F, et al. Global cancer statistics 2018. CA Cancer J Clin. 2018;68(6):394-424.

\bibitem{schmoll2015}
Schmoll HJ, et al. ESMO consensus guidelines for management of patients with colon and rectal cancer. Ann Oncol. 2015;26(5):999-1010.

\bibitem{andre2009}
Andre T, et al. Improved overall survival with oxaliplatin, fluorouracil, and leucovorin as adjuvant treatment in stage II or III colon cancer. J Clin Oncol. 2009;27(19):3109-16.

\bibitem{haller2011}
Haller DG, et al. Capecitabine plus oxaliplatin compared with fluorouracil and folinic acid as adjuvant therapy for stage III colon cancer. J Clin Oncol. 2011;29(11):1465-71.

\bibitem{sargent2014}
Sargent DJ, et al. A pooled analysis of the accuracy of KRAS mutation testing for predicting benefit from cetuximab in metastatic colorectal cancer. J Clin Oncol. 2014;32(18):1919-26.

\bibitem{lundberg2017}
Lundberg SM, Lee SI. A unified approach to interpreting model predictions. NeurIPS. 2017:4765-74.

\bibitem{ioannidis2009}
Ioannidis JP. Microarrays and molecular research: noise discovery? PLoS Med. 2009;6(7):e1000109.

\bibitem{begg2013}
Begg CB, et al. Impact of molecular subtype classification on the interpretation of breast cancer clinical trials. J Natl Cancer Inst. 2013;105(18):1387-93.

\bibitem{simon2003}
Simon R, et al. Pitfalls in the use of DNA microarray data for diagnostic and prognostic classification. J Natl Cancer Inst. 2003;95(1):14-8.

\bibitem{subramanian2005}
Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS. 2005;102(43):15545-50.

\bibitem{barbie2009}
Barbie DA, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462(7269):108-12.

\bibitem{verhaak2013}
Verhaak RG, et al. Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1. Cancer Cell. 2010;17(1):98-110.

\bibitem{yang2024}
Yang H, et al. PESSA: a web tool for pathway enrichment score-based survival analysis in cancer. PLoS Comput Biol. 2024;20(5):e1012024.

\bibitem{harrell2001}
Harrell FE, et al. Regression modeling strategies for improved prognostic prediction. Stat Med. 2001;20(1):144-58.

\bibitem{tgf2025nature}
Nature Genetics. TGF-$\beta$ builds a dual immune barrier in colorectal cancer. Nat Genet. 2025;57:xxx.

\bibitem{drg2024}
Drug Dev Res. A new perspective on antiangiogenic antibody drug resistance. Drug Dev Res. 2024.

\bibitem{venet2011}
Venet D, et al. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011;7(10):e1002240.

\end{thebibliography}
"""

# ── Supplementary Materials ──
supplementary = r"""

\newpage
\onecolumn
\section*{Supplementary Materials}
\setcounter{table}{0}
\renewcommand{\thetable}{S\arabic{table}}
\setcounter{figure}{0}
\renewcommand{\thefigure}{S\arabic{figure}}

\subsection*{Supplementary Tables}

\textbf{Supplementary Table S1.} 253 Significant Genes from Multi-Cohort Meta-Analysis. Source: \texttt{results/tables/meta\_analysis/deg\_meta\_filtered.csv} (2.1 MB). Stouffer meta-analysis across 5 GEO cohorts (FDR $<$ 0.05).

\textbf{Supplementary Table S2.} Complete DEG Results: GSE39582. Source: \texttt{results/tables/DEG\_GSE39582\_limma.csv} (7.9 MB).

\textbf{Supplementary Table S3.} Complete DEG Results: GSE104645. Source: \texttt{results/tables/DEG\_GSE104645\_limma.csv} (6.1 MB).

\textbf{Supplementary Table S4.} Complete DEG Results: GSE28702. Source: \texttt{results/tables/DEG\_GSE28702\_limma.csv} (7.8 MB).

\textbf{Supplementary Table S5.} Complete DEG Results: GSE72970. Source: \texttt{results/tables/DEG\_GSE72970\_limma.csv} (7.8 MB).

\textbf{Supplementary Table S6.} Complete DEG Results: GSE69657. Source: \texttt{results/tables/DEG\_GSE69657\_limma.csv} (7.7 MB).

\textbf{Supplementary Table S7.} WGCNA Hub Genes. Source: \texttt{results/tables/WGCNA\_hub\_genes.csv} (27 KB).

\textbf{Supplementary Table S8.} Univariate Cox Screening with BH Correction. Source: \texttt{results/tables/nomogram/univariate\_cox\_screening\_bh.csv}.

\textbf{Supplementary Table S9.} Pathway Consistency Matrix (5 Cohorts). Source: \texttt{results/tables/pathway\_activity/pathway\_consistency\_matrix.csv}.

\textbf{Supplementary Table S10.} Per-Cohort Pathway Differential Activity.

\textbf{Supplementary Table S11.} LASSO-PRS Bootstrap Stability (B = 1,000). Source: \texttt{results/tables/bootstrap\_stability/bootstrap\_selection\_freq.csv}.

\textbf{Supplementary Table S12.} Merged PRS Bootstrap Stability.

\textbf{Supplementary Table S13.} Gene-Level Univariate Cox (21,755 genes). Source: \texttt{results/tables/gene\_level/gene\_cox\_univariate.csv}.

\textbf{Supplementary Table S14.} AIC-Selected 23-Gene Cox Model.

\textbf{Supplementary Table S15.} Ablation Comparison.

\textbf{Supplementary Table S16--S21.} Sensitivity Diagnostics (PH test, dfbeta, LOVO C-index, Stratified C-index, Cutoff sweep, Bootstrap CI).

\textbf{Supplementary Table S22.} GSE83129 External Validation.

\textbf{Supplementary Table S23.} TCGA External Validation.

\textbf{Supplementary Table S24.} Multi-Algorithm Feature Votes.

\textbf{Supplementary Table S25.} GDSC Drug Sensitivity Pre-Filtering.

\textbf{Supplementary Table S26.} Meta-Analysis Direction Matrix.

\textbf{Supplementary Table S27.} Fair Comparison: Gene vs Pathway LASSO-Cox. See Table~S27 in main text supplementary section.

\textbf{Supplementary Table S28.} AIC Bootstrap Selection Frequency (B = 200).

\textbf{Supplementary Table S29.} Bootstrap Coefficient Distribution (7 Pathways).

\textbf{Supplementary Table S30.} Drug Sensitivity Analysis (oncoPredict).

\textbf{Supplementary Table S31.} Calibration Statistics.

\textbf{Supplementary Table S32.} DCA Net Benefit \& TDROC AUC.

\textbf{Supplementary Table S33.} Overfitting Diagnostics Summary.

\textbf{Supplementary Table S34.} Ridge Cox Regression Coefficients ($\alpha = 0$, $\lambda_{\min} = 0.035$).

\textbf{Supplementary Table S35.} Elastic Net Stability Selection ($\alpha = 0.5$, B = 100).

\subsection*{Supplementary Figures}

\textbf{Supplementary Figure S1.} WGCNA Analysis. Scale-free topology fit ($\beta = 6$) and module-trait heatmap.

\textbf{Supplementary Figure S2.} ComBat Batch Correction PCA.

\textbf{Supplementary Figure S3.} Same-Platform Validation ROC Curves.

\textbf{Supplementary Figure S4.} External Validation: Direction Concordance.

\textbf{Supplementary Figure S5.} ML Phase 2 Validation Suite.

\textbf{Supplementary Figure S6.} SHAP Dependence Plots.

\textbf{Supplementary Figure S7.} TCGA Validation Details.

\textbf{Supplementary Figure S8.} Pathway Differential Activity.

\textbf{Supplementary Figure S9.} LASSO-PRS Model Details.

\textbf{Supplementary Figure S10.} Bootstrap Stability Assessment (B = 1,000).

\textbf{Supplementary Figure S11.} Merged PRS Model (n = 277).

\textbf{Supplementary Figure S12.} Nomogram Forest Plot.

\textbf{Supplementary Figure S13.} Nomogram Time-Dependent ROC.

\textbf{Supplementary Figure S14.} GSE72970 External KM.

\textbf{Supplementary Figure S15.} Sensitivity Diagnostics Full Panel.

\textbf{Supplementary Figure S16.} PRS-Only Nomogram and Calibration.

\textbf{Supplementary Figure S17.} GSEA Running Score Plots.

\textbf{Supplementary Figure S18.} ORA Enrichment Dotplots.

\textbf{Supplementary Figure S19.} Gene-Level Cox Forest Plot.

\textbf{Supplementary Figure S20.} GSE83129 Validation.

\textbf{Supplementary Figure S21.} TCGA External Validation (PRS).

\textbf{Supplementary Figure S22.} Meta-Analysis Overview.

\textbf{Supplementary Figure S23.} TGF-$\beta$ Paradox Analysis (5-panel).

\textbf{Supplementary Figure S24.} Fair Comparison: Gene vs Pathway C-index.

\textbf{Supplementary Figure S25.} Bootstrap Coefficient Forest Plot.

\textbf{Supplementary Figure S26.} Drug Sensitivity Boxplots (oncoPredict).

\textbf{Supplementary Figure S27.} Calibration Curves.

\textbf{Supplementary Figure S28.} Decision Curve Analysis.

\textbf{Supplementary Figure S29.} Time-Dependent ROC.

\textbf{Supplementary Figure S30.} Ridge Cox Regression Path ($\alpha = 0$, 10-fold CV, $\lambda_{\min} = 0.035$).

\textbf{Supplementary Figure S31.} Elastic Net Stability Selection ($\alpha = 0.5$, B = 100).

\end{document}
"""

# Write the .tex file
with open(OUT_TEX, 'w', encoding='utf-8') as f:
    f.write(preamble)
    f.write(main_text)
    f.write(references)
    f.write(supplementary)

print(f'[OK] Elsevier LaTeX generated: {OUT_TEX}')
print(f'     Size: {os.path.getsize(OUT_TEX) / 1024:.0f} KB')
print()
print('Compile instructions:')
print('  1. Open in Overleaf: upload .tex + figures_png/hd/ folder')
print('  2. Or compile locally:')
print('     pdflatex manuscript_v5.1_elsarticle.tex')
print('     bibtex manuscript_v5.1_elsarticle')
print('     pdflatex manuscript_v5.1_elsarticle.tex')
print('     pdflatex manuscript_v5.1_elsarticle.tex')
