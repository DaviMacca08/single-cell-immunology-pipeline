Trajectory Analysis of Intestinal T-Cell Compartments in Crohn’s Disease
================
Bioinformatics Analysis Service — Davide Maccarrone
2026-09-02

# Executive Summary

This report describes the trajectory inference analysis performed on the
intraepithelial lymphocyte (IEL) and lamina propria lymphocyte (LPL)
compartments of the Jaeger et al. 2021 dataset, using Slingshot as the
sole trajectory inference engine. Unlike the pseudobulk DESeq2,
CellChat, and LIANA analyses conducted elsewhere in this project, the
trajectory analysis does **not** attempt a formal statistical comparison
of Crohn’s disease (CD) versus Control along pseudotime. This choice was
made after diagnostic inspection of both compartments revealed
condition-associated cell-type composition imbalance so severe that most
cell states are effectively condition-exclusive, making formal
cross-condition trajectory comparisons poorly informative for many
lineages and highly sensitive to donor composition, given the strong
condition-associated cell-state imbalance and the n = 4 donor design.

Instead, Slingshot was fit once per compartment on the pooled
(condition-agnostic) dataset, producing a single descriptive map of cell
states. Disease condition is then overlaid onto this map purely
descriptively — via condiments’ `imbalance_score`, pseudotime density by
condition, and compositional profiles along pseudotime bins — without
any formal hypothesis test. This is a deliberate and explicitly
justified methodological departure from the condiments/tradeSeq
framework introduced in earlier drafts of this analysis, and is
documented as such throughout.

Two substantive findings emerged during quality control that materially
affect interpretation, and one of them directly resolves an open
question from the bulk RNA-seq cross-validation report:

1.  **IEL** — a shared trajectory trunk connects CD4-labeled and
    CD8-labeled cell states without a clean, monotonic co-receptor
    reprogramming gradient on canonical markers (CD4, CD8A, CD8B, RUNX3,
    ZBTB7B); this trunk should be read as a region of transcriptional
    similarity rather than a literal CD4→CD8αα developmental conversion.
2.  **LPL** — several CD-associated populations, including the FOXP3+
    activated Treg population whose increased relative abundance in CD
    had been flagged as an unresolved discrepancy with Jaeger et al.
    2021) in the bulk RNA-seq report, are contributed almost exclusively
          by a single CD donor, not fully resolved by Harmony
          integration. This finding provides a donor-level explanation
          for the observed discrepancy in this dataset: the apparent
          increase in LPL Treg abundance is driven almost entirely by a
          single CD donor and therefore cannot be interpreted as a
          generalizable disease-associated shift.

# Analytical Workflow

The analytical pipeline was applied independently to each compartment
(IEL and LPL), following an identical structure:

- SingleCellExperiment (SCE) construction from Seurat v5 objects, with
  raw counts pulled from the joined, whole-object layer to avoid
  multi-layer ambiguity, and explicit assignment of `counts`/`logcounts`
  assays.
- Condition-imbalance diagnostic (`condiments::imbalance_score`)
  computed on the full compartment UMAP, prior to any trajectory
  fitting.
- Removal of non-lymphoid / technically negligible populations
  (contaminant epithelial and myeloid clusters) from the cells used for
  trajectory fitting.
- Slingshot lineage inference (`getLineages`) on the PCA embedding,
  using cell-type labels as cluster input, an explicitly verified
  biological root cluster, and `omega = TRUE` to avoid forcing
  geometrically distant, biologically unrelated clusters onto a single
  connected tree.
- Principal curve fitting (`getCurves`), pseudotime and curve-weight
  extraction (`slingPseudotime`, `slingCurveWeights`).
- Descriptive, condition-agnostic characterization of pseudotime:
  per-lineage visualization, condition density overlays, and
  compositional profiles along pseudotime bins, plus descriptive
  marker-trend plots (binned means, no model fit).
- No condiments statistical tests (`topologyTest`, `progressionTest`,
  `differentiationTest`, `fateSelectionTest`) and no tradeSeq model
  fitting (`fitGAM`, `associationTest`, `startVsEndTest`,
  `conditionTest`) were run as part of the final reported workflow — see
  the rationale below.

## Methodological Rationale: Why No Formal CD vs Control Trajectory Testing

The condiments framework (Roux de Bézieux et al., *Nat Commun* 2024) is
designed to formalize the comparison of a shared trajectory across two
conditions: does topology differ, does pseudotime progression differ
along a shared lineage, does fate selection differ at a branch point.
This framework presupposes that most or all trajectory states are
meaningfully populated by both conditions. Diagnostic inspection of both
IEL and LPL (see per-compartment sections below) showed that this
precondition is not met: the large majority of annotated cell states in
both compartments are condition-exclusive or near-exclusive (≥ 90–99% of
cells from a single condition), a compositional pattern already
established and characterized elsewhere in this project via the
CellChat/LIANA validity-tier framework (fully comparable / imbalanced /
CD-emergent / Control-dominant).

Applying condiments’ or tradeSeq’s condition-aware tests to trajectory
branches that are condition-exclusive would not test whether the same
biological process proceeds at different kinetics between conditions; it
would merely re-derive, with an inferential p-value, the compositional
imbalance already documented via simple contingency tables. Combined
with the project’s n = 4 donor design (2 CD / 2 Control per compartment,
never modeled with a donor-blocking term in any analysis layer of this
project), the effective replication available for a genuine
cross-condition kinetic comparison is negligible on most branches. For
this reason, condition is treated throughout this analysis as a
**descriptive covariate** overlaid on a single, condition-agnostic
trajectory map, rather than as an experimental factor in a formal
statistical model — the same logic already applied to CD-emergent and
Control-dominant populations in the CellChat/LIANA communication
analysis of this project.

A single formal check — `condiments::topologyTest_multipleSamples`,
which accounts for donor-of-origin via a `Samples` argument rather than
treating every cell as an independent replicate — was considered as an
optional, one-time documented consistency check rather than a decision
gate, but was ultimately not required to justify the descriptive
approach given how unambiguous the compositional evidence already was in
both compartments.

# Results — IEL Compartment

## Composition Imbalance

Cell-type composition in IEL is severely imbalanced by condition (global
ratio: 44% CD / 56% Control). The large majority of annotated
populations are condition-exclusive or near-exclusive: CD-dominant
examples include FOXP3+ IL2RA+ Treg (79 CD / 0 Control), Terminal
effector CD8 T cells (529 / 0), TH17-like CD4 T cells (1619 / 2), and
Cytotoxic γδ T cells (2113 / 11); Control-dominant examples include
IL7R+ memory-like CD8 T cells (0 / 1481), γδ T cells / NK-like cytotoxic
T cells (1 / 2287), and “Naive CD4 T cells” (4 / 1935 — not to be
confused with the marker-verified “Naïve CD4” root cluster below, a
separate, correctly-spelled label in the annotation).

<div style="display:flex; gap:10px; align-items:flex-start;">

<img src="../../results/plots/trajectory_IEL/ImbalanceScore_IEL.png"
style="width:50.0%" /> <img
src="../../results/plots/trajectory_IEL/Imbalance_Score_fitted_IEL.png"
style="width:50.0%" />

</div>

## Root Cluster Selection

Two candidate labels (“Naïve CD4” and “Naïve / central memory CD4”) were
both confirmed as genuinely quiescent/naive by mean expression of CCR7,
TCF7, SELL, and LEF1 across all IEL clusters, prior to selecting a root.
“Naïve CD4” scored highest on all four markers (CCR7 1.40, TCF7 1.19,
SELL 1.56, LEF1 0.66) and was selected as the Slingshot root
(`start.clus`). A third, superficially similar label (“Naive CD4 T
cells”) was explicitly ruled out despite its name, showing near-zero
CCR7/SELL (0.009 / 0.001) — a reminder that cluster labels were not
taken at face value for this decision. No IEL population showed a
comparable naive-like profile for the CD8 lineage, a limitation carried
into the topology discussion below.

## Trajectory Topology

Slingshot (`getLineages`, PCA embedding, `omega = TRUE`,
`start.clus = "Naïve CD4"`) recovered 7 lineages, all sharing “Naïve
CD4” as their root — consistent with, though not independently
re-confirmed via an explicit `start.given` metadata check the way it was
for LPL (see below). `omega = TRUE` was explicitly confirmed active
(`omega_scale = 1.5`) via `metadata(lineages_iel)$slingParams`. Lineages
1–2 share a long common trunk that alternates between CD4-labeled and
CD8-labeled states (Naïve CD4 → Naïve/CM CD4 → Activated memory CD4 →
IL7R+ memory-like CD8 → Effector mixed → CD39+ tissue-resident CD8 →
TH17-like innate CD4 → Cycling → TH17-like CD4 → terminal CD8/Treg).
This trunk persisted unchanged after `omega = TRUE` was applied,
indicating the geometric distance between these CD4- and CD8-labeled
clusters in PCA space does not exceed the omega threshold (1.5 × median
MST edge length) — i.e. this adjacency is a property of the data, not an
artifact the omega correction could have removed.

![](../../results/plots/trajectory_IEL/IEL_trunk.png)

This adjacency was investigated against the documented phenomenon of
local CD4→CD8αα lineage re-commitment in gut intraepithelial lymphocytes
(Mucida et al., *Nat Immunol* 2013, 14:281–289; Reis et al., *Nat
Immunol* 2013, 14:271–280; reviewed in Park, Moon & Lee, *BMB Reports*
2016, 49:11–17), which would predict a monotonic decline in CD4/ZBTB7B
and a monotonic rise in CD8A/CD8B/RUNX3 along the trunk. The observed
pattern does **not** support this: CD4 is uniformly low throughout the
trunk (a known technical limitation of CD4 mRNA detection in scRNA-seq,
not by itself informative), ZBTB7B remains near-zero throughout without
ever being elevated at the presumptive root, and CD8A/RUNX3 rise through
the middle of the trunk but decline toward the terminal state rather
than continuing to increase. The trunk is therefore reported as a region
of shared transcriptional similarity — consistent with the general
convergence of gut IEL states on partial CD8αα/RUNX3 programs regardless
of nominal CD4/CD8 identity — rather than as evidence of a literal,
ordered lineage-conversion process.

## Pseudotime

<div style="display:flex; gap:10px; align-items:flex-start;">

<img src="../../results/plots/trajectory_IEL/Pseuotime_overall_IEL.png"
style="width:32.0%" /> <img
src="../../results/plots/trajectory_IEL/pseudotime_IEL_lineage_1.png"
style="width:32.0%" /> <img
src="../../results/plots/trajectory_IEL/pseudotime_IEL_lineage_7.png"
style="width:32.0%" />

</div>

Pseudotime gradients are spatially coherent across all 7 lineages, with
no evidence of embedding breakdown. Lineages 5–7 are short (3–4 nodes)
and terminate in small, spatially isolated populations; their pseudotime
ranges, while wide, are compressed into narrow UMAP regions and should
be interpreted with more caution than the longer lineages.

## Descriptive Condition Overlay

Condition is overlaid on the IEL trajectory purely descriptively, via
pseudotime density by condition and a compositional profile (proportion
CD) across ten pseudotime bins per lineage (figures produced by the
pipeline as `Pseudotime_distribution_condition.png` and
`CompositionProfile.png` for the IEL compartment; see project figure
archive — these particular files were superseded on-disk by same-named
LPL exports and are referenced here by content rather than re-embedded).

Lineages 1, 2, and 3 show a consistent, non-monotonic (“U-shaped”)
compositional pattern: CD dominates the earliest pseudotime bin (an
inherited effect of the CD-skewed “Naïve CD4” root, 678 CD / 8 Control),
Control dominates the intermediate bins, and CD dominates again toward
the terminal bins. Lineage 7 (terminating in CXCL13+ TFH-like T cells,
94.9% Control) is the most extreme case: CD is confined almost entirely
to the first pseudotime bin, consistent with its near-total Control
dominance in the underlying cell-type table.

![](../../results/plots/trajectory_IEL/MarkerTrends_along_pseudotime.png)

## IEL — Compartment-Specific Limitations

- No naive-like CD8 population was identified in IEL by canonical
  markers, meaning the Slingshot root is necessarily CD4-lineage; any
  CD8 branch is therefore rooted through the shared trunk rather than
  through an independently verified CD8 origin.
- The shared CD4/CD8 trunk (Lineages 1–2) should not be read as evidence
  of a literal reprogramming trajectory; see marker verification above.
- Lineages 5–7 are short and terminate in small, spatially compressed
  populations; pseudotime within them is more susceptible to noise than
  in the longer lineages.
- As throughout this project, all compositional and descriptive
  condition comparisons are limited by the n = 4 donor design (2 CD / 2
  Control), with no donor-blocking term used in any analytical layer.

# Results — LPL Compartment

## Composition Imbalance

As in IEL, the large majority of LPL populations are condition-exclusive
or near-exclusive. Notable examples: FOXP3+ activated Treg (2205 CD / 0
Control), NK-like cytotoxic T cells (296 / 0), TRM-like CD8 cytotoxic T
cells (204 / 0), Activated CD4 memory T cells (5578 / 1), Activated CD4
T cells stress response (15 / 2628), and T follicular helper (24 /
2000).

<div style="display:flex; gap:10px; align-items:flex-start;">

<img src="../../results/plots/trajectory_LPL/Condition_umap_LPL.png"
style="width:50.0%" />
<img src="../../results/plots/trajectory_LPL/ImbalanceScore_LPL.png"
style="width:50.0%" />

</div>

## Root Cluster Selection

“Naive CD4 T cells” was used as the Slingshot root (`start.clus`),
confirmed via `start.given = TRUE`. Unlike IEL, this selection was
**not** independently cross-checked against canonical naive markers
(CCR7/TCF7/SELL/LEF1) prior to fitting; it is used here because it is
the only unambiguously naive-labeled population in the LPL annotation
set and its label did not require the same disambiguation that IEL’s
near-duplicate naive labels did. This is noted as a lighter-touch
verification than was applied in IEL.

## Trajectory Topology

Slingshot (`omega = TRUE`) recovered 11 lineages. Ten of these (Lineages
1–10) share “Naive CD4 T cells” as a common root and remain connected in
a single tree; `omega` correctly isolated only Lineage 11 (root and end
both “Exhaustion-like CD4 memory T cells”, n = 36 cells total, 3 CD / 33
Control) as an independent sub-tree, which should be treated as
non-conclusive given its small size and spatial isolation.

<div style="display:flex; gap:10px; align-items:flex-start;">

<img src="../../results/plots/trajectory_LPL/lineages_LPL.png"
style="width:50.0%" /> <img
src="../../results/plots/trajectory_LPL/pseudotime_LPL_lineage_11.png"
style="width:50.0%" />

</div>

## Donor Effect on Lineages 1–10 (Critical Finding)

Visual inspection of the LPL UMAP revealed two spatially disconnected
point clouds (“islands”), with the smaller, right-hand island almost
entirely composed of CD cells. Because `omega = TRUE` did not separate
this island into an independent sub-tree, Lineages 1–10 extend across
both islands via PCA-space adjacency that is not visually apparent in
the UMAP projection.

Donor-level inspection resolved this decisively: of the four donors,
donor **1818 (CD)** contributes 15,247 of its 15,264 total cells (99.9%)
to the right-hand UMAP island, while the other three donors are almost
entirely confined to the left-hand island. Cross-referencing with
cell-type composition shows this single donor accounts for the
overwhelming majority of several of the compartment’s most CD-skewed
populations:

| Cell type | Cells from donor 1818 | Total cells (compartment) |
|:---|:--:|:--:|
| FOXP3+ activated Treg | 2204 | 2205 |
| Activated CD4 memory T cells | 5578 | 5579 |
| NK-like cytotoxic T cells | 296 | 296 |
| TRM-like CD8 cytotoxic T cells | 204 | 204 |

These populations are therefore contributed almost exclusively by a
single donor, not resolved by Harmony integration, rather than being
shared, reproducible CD-associated states. With n = 2 CD donors in this
compartment, a population contributed at this level by one donor is
statistically indistinguishable from a donor-specific effect.

## Connection to the Bulk RNA-seq Report: The Treg Discrepancy, Resolved

The bulk RNA-seq cross-validation report (GSE193677 analysis) flagged an
unresolved discrepancy: this project’s pseudobulk/scRNA-seq
compositional analysis found an increased relative Treg proportion in CD
(LPL: 0.12 CD vs 0.05 Control), directly contradicting Jaeger et al.
(2021), who reported reduced Treg proportions in inflamed CD tissue.
That report explicitly declined to resolve the discrepancy, attributing
it tentatively to possible differences in Treg subtype definition or
sampling between studies.

The donor-level finding above provides a direct, mechanistic resolution:
**2204 of the 2205 FOXP3+ activated Treg cells in this project’s LPL
compartment originate from a single CD donor (1818)**. The apparent
CD-associated increase in Treg proportion is therefore not
distinguishable from a donor-specific effect and should not be
interpreted as a generalizable, disease-associated compositional shift.
This finding should be incorporated into the manuscript’s discussion of
the Treg discrepancy in place of the previously open-ended language.

## Pseudotime

<div style="display:flex; gap:10px; align-items:flex-start;">

<img src="../../results/plots/trajectory_LPL/Pseudotime_overall_LPL.png"
style="width:32.0%" />

<img
src="../../results/plots/trajectory_LPL/pseudotime_LPL_lineage_1.png"
style="width:32.0%" />

<img
src="../../results/plots/trajectory_LPL/pseudotime_LPL_lineage_8.png"
style="width:32.0%" />

</div>

## Descriptive Condition Overlay

<div style="display:flex; gap:10px; align-items:flex-start;">

<img
src="../../results/plots/trajectory_LPL/Imbalance_Score_fitted_LPL.png"
style="width:32.0%" /> <img
src="../../results/plots/trajectory_LPL/Pseudotime_distribution_condition.png"
style="width:32.0%" />
<img src="../../results/plots/trajectory_LPL/CompositionProfile.png"
style="width:32.0%" />

</div>

## LPL — Compartment-Specific Limitations

- A subset of Lineages 1–10 (those reaching the donor-1818-dominated
  UMAP region, including at minimum the lineages terminating in FOXP3+
  activated Treg, NK-like cytotoxic T cells, TRM-like CD8 cytotoxic T
  cells, and passing through Activated CD4 memory T cells) cannot be
  distinguished from a single-donor effect and must not be reported as
  generalizable CD-associated trajectory findings.
- Root selection (“Naive CD4 T cells”) was not independently
  marker-verified as it was in IEL; this is a lower-rigor step than the
  equivalent IEL decision.
- Lineage 11 is based on 36 cells total and is non-conclusive.
- As in IEL, all findings are subject to the project-wide n = 4 donor
  constraint (2 CD / 2 Control), compounded here by the donor-specific
  effect identified above.

# Project-Wide Limitations and Considerations

## Donor sample size

Both IEL and LPL compartments are limited to n = 4 total donors (2 CD /
2 Control). No donor-blocking or paired design has been used in any
analytical layer of this project, including this trajectory analysis.
Condiments’ donor-aware statistical tests
(`topologyTest_multipleSamples`, `progressionTest_multipleSamples`,
`fateSelectionTest_multipleSamples`) were considered but not adopted as
the primary analytical framework, given that the compositional imbalance
documented above makes most branches ineligible for meaningful
cross-condition testing regardless of the test’s donor-awareness.

Because trajectory inference and descriptive pseudotime summaries were
performed at the cell level, donor-to-donor variability cannot be
estimated reliably with only two donors per condition. Consequently, the
inferred topology and pseudotime distributions should be interpreted as
dataset-level descriptive structures rather than donor-generalized
disease effects.

## No formal trajectory inference method comparison

Slingshot was used as the sole trajectory inference method. An
independent topology cross-validation using a second trajectory
inference framework was not performed and therefore represents a
methodological limitation of the analysis.

## tradeSeq not applied

tradeSeq (`fitGAM` and downstream tests) was attempted on the IEL
compartment but was not computationally tractable on standard hardware
even after `nknots` selection via `evaluateK` (k range truncated to 3–7
for feasibility; the AIC-relative elbow fell clearly at k = 4, within
the tested range) and gene-set filtering strategies; a single `fitGAM`
run was projected at approximately 1 day 19 hours and was abandoned. No
gene-level, model-based differential expression along pseudotime is
reported for either compartment. Descriptive, non-model-based
marker-trend plots (binned mean expression) are provided as a
lower-rigor substitute for the biological questions that tradeSeq would
otherwise have addressed.

![](../../results/plots/trajectory_IEL/evaluateK_plot.png)

## Interpretive tier framework

Consistent with the validity-tier framework already applied to
CellChat/LIANA and pseudobulk DESeq2 elsewhere in this project (fully
comparable / imbalanced / CD-emergent / Control-dominant), no trajectory
finding in this report should be read as a differential comparison
between CD and Control. All condition-related observations are
descriptive characterizations of where each condition’s cells sit on a
condition-agnostic map, not statistically tested differences in
trajectory kinetics.

# Reproducibility

All analyses were performed in R using Seurat v5 (`LayerData`-based
assay access), SingleCellExperiment, Slingshot, and condiments
(`imbalance_score` only). A fixed random seed (`set.seed(1234)`) was
used throughout. Complete session information for both compartments is
stored as `sessionInfo_Trajectory_IEL.txt` and
`sessionInfo_Trajectory_LPL.txt` in the project’s logs directory;
**exact package version numbers should be inserted here from those logs
prior to manuscript submission**, rather than restated from memory in
this report.

| Component            | Version            |
|:---------------------|:-------------------|
| R                    | 4.6.0              |
| Seurat               | v5 (LayerData API) |
| SingleCellExperiment | 1.34.0             |
| slingshot            | 2.20.0             |
| condiments           | 1.20.0             |

Software versions used in the trajectory analysis (placeholders to be
completed from saved session logs).

# Conclusion

This trajectory analysis provides a condition-agnostic, descriptive map
of T-cell states in both the IEL and LPL compartments of the Jaeger et
al. dataset, deliberately avoiding formal CD-versus-Control statistical
testing given the severity of condition-associated compositional
imbalance and the project’s n = 4 donor constraint. Two findings arising
from careful quality control materially strengthen the project’s overall
interpretive rigor rather than merely qualifying it: the IEL shared
CD4/CD8 trunk is shown, through direct marker verification, not to
support a literal reprogramming narrative; and the LPL donor-effect
finding provides a concrete, donor-level explanation to a discrepancy
with the source literature that had previously been left as an open
question in the bulk RNA-seq report. Both findings should be carried
forward explicitly into the manuscript rather than treated as
analysis-stage-only caveats.
