Cell–Cell Communication Analysis (CellChat v2) of Intestinal IEL and LPL
Compartments in Crohn’s Disease
================
Davide Maccarrone
2026-08-19

# Executive Summary

This report presents the cell–cell communication analysis of
intraepithelial (IEL) and lamina propria (LPL) lymphocyte compartments
in Crohn’s disease (CD) versus control intestinal tissue,using CellChat
v2 (Jin et al., 2024) as an orthogonal, network-level complement to the
pseudobulk differential expression and GSEA analyses reported
separately. Where GSEA characterizes *which programs* are active within
a cell type, CellChat characterizes *who is signaling to whom*, and
through which ligand–receptor pairs — a distinct and complementary axis
of biological information.

Two independent CellChat objects were built per compartment (Control and
CD), never merged prior to inference, and subsequently aligned
(`liftCellChat()`) for differential comparison. Across both
compartments, the analysis converges on a small number of
well-supported, literature-concordant findings:

- **A denser, structurally more complex communication network in CD** in
  both compartments, consistent with active mucosal inflammation, though
  the mechanism differs: in IEL the CD network is *diffuse* (more
  pathways, lower average edge strength — partly a function of more cell
  populations clearing the detection threshold in CD), whereas in LPL
  the CD network is *proportionally more intense* (more interactions and
  more total strength, ~1.6–1.7× Control).
- **MHC-II signaling shows a disease-associated reorganization involving
  regulatory T-cell populations in both compartments.**, but through
  population-specific rather than generic mechanisms: in IEL,
  CD-emergent Treg/macrophage/TH17-like populations dominate the CD-side
  receiver profile; in LPL, the receiving Treg subset is *replaced*
  between conditions (Activated effector Treg in Control → FOXP3⁺
  activated Treg in CD), not merely intensified.
- **LIGHT (TNFSF14) signaling toward FOXP3⁺ IL2RA⁺ Treg via LTβR, with a
  broader TNFRSF14 (HVEM) axis toward multiple populations (IEL)**, and
  **SIRPG–CD47 T–T costimulatory signaling (LPL)** are identified as
  compartment-specific, biologically coherent, literature-supported
  axes.
- **CD103 (ITGAE/ITGB7)–E-cadherin** emerges as the most specific
  ligand–receptor pair recovered in the IEL compartment, mechanistically
  consistent with the established role of this axis in intraepithelial
  lymphocyte retention.
- The single most important methodological caveat governing the entire
  analysis is **severe, condition-dependent cell type imbalance**,
  formalized here as an explicit validity classification per
  compartment. All findings should be interpreted through this
  framework.

This report consolidates and reframes findings already produced and
internally reviewed in the IEL and LPL CellChat pipelines (scripts
`13–16` and `17–20` respectively).

# Analytical Rationale

## Why cell–cell communication inference as a complementary method

Pseudobulk DESeq2 and GSEA describe transcriptional programs *within* a
cell type. They cannot, by construction, describe the *directionality*
of intercellular signaling — which population is the likely source of a
ligand, which population expresses the cognate receptor, and how that
sender–receiver structure reorganizes with disease. CellChat v2 infers
this structure from ligand–receptor co-expression using a curated
database of literature-supported interactions and provides both
quantitative summaries (aggregated network strength/count, information
flow per pathway) and cell-type-resolved output.

## Pipeline architecture

- Two fully independent CellChat objects per compartment
  (`createCellChat` → `subsetData` →
  `identifyOverExpressedGenes/Interactions` → `computeCommunProb`
  (trimean method, package default) →
  `filterCommunication(min.cells = 10)` → `computeCommunProbPathway` →
  `aggregateNet`), never merged prior to inference.
- `liftCellChat()` was used to align cell type levels across conditions
  for comparison; `netAnalysis_computeCentrality()` was re-run after
  lifting because lifted matrices differ in dimension from the
  originals.
- Reproducibility: `netEmbedding(..., umap.method = "uwot")` with
  `set_seed(1234)` immediately before each stochastic block (functional
  similarity embedding, `selectK`, `identifyCommunicationPatterns`).
- IFN-II, IL17, and TGFb signaling are absent from CellChat results
  across both compartments and both conditions. TNF signaling is absent
  in IEL (both conditions) but is detected in LPL (both Control and CD).
  This partial absence should not be interpreted as biological absence
  of these cytokine programs — it most plausibly reflects the detection
  threshold of the trimean-based default method, combined with
  compartment-specific differences in which cell populations express the
  relevant ligand-receptor pairs at sufficient frequency to pass
  filtering.

## Scope

Unlike the GSEA analysis, CellChat operates directly on single-cell
expression and therefore includes every annotated cell type in each
compartment, subject to the `min.cells = 10` per-condition filter
applied independently to each condition’s network.

> ⚠️ **Interpretation rule**
>
> Cell types that fail the `min.cells = 10` threshold in one condition
> cannot support a standard quantitative CD-vs-Control comparison for
> that population. The resulting outputs must be labeled as
> within-condition characterization, imbalanced comparison, or fully
> comparable comparison.

# Methodological Framework: Cell Type Comparison Validity

`filterCommunication(min.cells = 10)` excludes any cell type with fewer
than 10 cells in a given condition from that condition’s network;
`liftCellChat()` subsequently re-inserts it as an all-zero row/column
when aligning objects for comparison. Both compartments show
substantial, condition-dependent imbalance in cell type composition.
This imbalance is the single most important lens through which every
result in this report must be read.

Tier boundaries are set on the absolute cell count on the weaker side of
each population (not on the CD:Control ratio): fewer than 10 cells on
the weaker side aligns with the `filterCommunication(min.cells = 10)`
threshold and places a population in the emergent/dominant tier
(within-condition characterization only); 10–~30 cells on the weaker
side is treated as imbalanced (comparison possible but indicative);
above ~30 cells on both sides is treated as fully comparable.

## IEL compartment

| Tier | Cell types (n CD / n Control) | Interpretation |
|----|----|----|
| **Fully comparable** | GZMK⁺ effector memory CD8 T cells (183/686); Effector T cell mixed state (213/356); CXCL13⁺ TFH-like T cells (75/1394); Cycling T cells (207/35) | Fully interpretable CD-vs-Control differential comparison |
| **Imbalanced** | CD39⁺ tissue-resident CD8 T cells (2551/20); Cytotoxic γδ T cells (2113/11); NK-like cytotoxic T cells (249/10) | Comparison possible but Control side is statistically weak; treat as indicative |
| **CD-emergent** | FOXP3⁺ IL2RA⁺ Treg (79/0); CST3⁺LYZ⁺ macrophages (35/1); TH17-like CD4 T cells (1619/2); TH17-like innate-associated CD4 T cells (832/6); Terminal effector CD8 T cells (529/0); Naïve CD4 (678/8); Epithelial cells/enterocytes (20/6) | **Within-CD characterization only** |
| **Control-dominant** | γδ/NK-like cytotoxic T cells (1/2287); Naive CD4 T cells (4/1935); TH17 effector (4/1508); IL7R⁺ memory-like CD8 T cells (0/1481); Naïve/central memory CD4 (5/1238); Activated memory CD4 (4/983) | **Within-Control characterization only** |

## LPL compartment

| Tier | Cell types (n CD / n Control) | Interpretation |
|----|----|----|
| **Fully comparable** | Activated cytotoxic CD8 T cells (5650/2002); Naive CD4 T cells (2284/3756); Activated effector Treg (873/617); CD4 TRM-like activated memory T cells (753/249) | Fully interpretable CD-vs-Control differential comparison |
| **Imbalanced (moderate)** | γδ / innate cytotoxic T cells (486/57); Cycling T cells (230/23) | Comparison possible; one side weaker — indicative estimates |
| **CD-heavy** | Th17-like activated CD4 T cells (4118/53); Innate-like/TRM CD8 (2194/16); IFN-activated CD8 cytotoxic T cells (1474/47); Activated CD4 memory T cells (5578/1) | Characterization predominantly in CD |
| **Control-heavy** | Activated CD4 T cells stress response (15/2628); T follicular helper (24/2000) | Characterization predominantly in Control |
| **CD-exclusive** | FOXP3⁺ activated Treg (2205/0); NK-like cytotoxic T cells (296/0); TRM-like CD8 cytotoxic T cells (204/0) | **Within-CD characterization only** |
| **Excluded** | Epithelial/stromal contamination (36/0) | Declared contamination; excluded from pathway deep-dives and hierarchy plots |

> Structural note

> Unlike IEL, the LPL compartment contains **no myeloid or epithelial
> populations** — it is entirely T-cell-composed. This is a structural
> property of the compartment and explains why LPL lacks direct
> analogues to the IEL macrophage- or epithelial-centered pathways MIF
> and CDH1.

# Global Network Architecture

| Metric | IEL (CD vs Control) | LPL (CD vs Control) |
|----|----|----|
| Interaction count | ~3× more pathways detected in CD | 2208 vs 1301 (~1.7×) |
| Interaction strength | Lower average strength per edge in CD (diffuse network) | 161.0 vs 101.4 (~1.6×) |
| Communication-pattern complexity | Shared K across conditions | Outgoing K=4 Control vs K=8 CD; Incoming K=7 vs K=7 |
| Overall pattern | More diffuse, less concentrated network in CD | Proportionally more intense network in CD |

`compareInteractions()`, `netVisual_diffInteraction()`, and
`netVisual_heatmap()` operate on the full aggregated network and are
therefore valid at the network level without depending on the
statistical power of a single cell type. Functional similarity embedding
and `rankSimilarity` likewise describe pathway-level reorganization.

**Figure 1.** Global interaction count and strength, Control versus CD,
for both compartments.

![](../../results/plots/CellChat_IEL/merge/CellChat_IEL_compareInteractions_count.png)
![](../../results/plots/CellChat_IEL/merge/CellChat_IEL_compareInteractions_weight.png)
![](../../results/plots/CellChat_LPL/merge/CellChat_LPL_compareInteractions_count.png)
![](../../results/plots/CellChat_LPL/merge/CellChat_LPL_compareInteractions_weight.png)

A compartment-level asymmetry is worth noting explicitly: the larger
apparent pathway count in CD IEL reflects a combination of genuine
transcriptional/functional heterogeneity and the near-mechanical
consequence of more distinct cell populations clearing the
`min.cells = 10` threshold in CD.

# IEL Compartment: Pathway-Level Findings

Six pathways were selected for deep-dive interpretation: **MHC-I,
MHC-II, LIGHT, MIF, GALECTIN, CDH1**.

The `celltypes_of_interest` set comprised seven populations: FOXP3⁺
IL2RA⁺ Treg, CST3⁺LYZ⁺ macrophages, TH17-like CD4 T cells, CD39⁺
tissue-resident CD8 T cells, Epithelial cells / enterocytes, GZMK⁺
effector memory CD8 T cells, and Effector T cell mixed state.

By cell-count validity tier (Control vs. CD; see Section 3.1): FOXP3⁺
IL2RA⁺ Treg (0/79), TH17-like CD4 T cells (2/1619), and CST3⁺LYZ⁺
macrophages (1/35) are condition-emergent and are interpreted
within-condition only. CD39⁺ tissue-resident CD8 T cells (20/2551) is
severely imbalanced toward CD and is not treated as a genuine shared
population despite being nominally present in both conditions.
Epithelial cells / enterocytes (6/20) have low absolute counts in both
conditions, limiting statistical power for either within- or
cross-condition comparison. Only GZMK⁺ effector memory CD8 T cells
(686/183) and Effector T cell mixed state (356/213) have cell numbers in
both conditions sufficient to support a genuine quantitative
cross-condition comparison.

## MHC-I

Ubiquitous signaling pattern in both conditions, as expected for a
broadly expressed pathway. Beyond simple homogenization, the
receiver/sender structure itself shifts between conditions: CST3⁺LYZ⁺
macrophages and Epithelial cells/enterocytes were not detected as active
participants in the inferred MHC-I network in Control (zero interactions
as either source or target) and emerge as active senders only in CD.
Conversely, CD39⁺ tissue-resident CD8 T cells displays a strong increase
in inferred MHC-I incoming signaling in CD, consistent with its
expansion and increased representation within the CD network. This
pattern is compatible with inflammatory IFN-driven remodeling described
in intestinal inflammation, and should be read as a structural
reorganization rather than a simple intensification of a shared axis.
The hierarchy plots below use CD39⁺ tissue-resident CD8 T cells and
Cytotoxic γδ T cells as receivers: Cytotoxic γδ T cells is the one
target consistently populated in both conditions, while CD39⁺ TRM CD8
becomes substantial only in CD (near-zero in Control).

**Figure 6.** MHC-I circle plots (Control and CD) and hierarchy plots
(Control and CD).

#### Circle plots

**Control**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CTRL_MHC-I_circle.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CD_MHC-I_circle.png)

#### Hierarchy plots

**Control**

![](../../results/plots/CellChat_IEL/merge/pathways/hierarchy/CellChat_IEL_hierarchy_MHC_I_Control.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_IEL/merge/pathways/hierarchy/CellChat_IEL_hierarchy_MHC_I_CD.png)

## MHC-II — *receiver population remodeling rather than simple amplification*

The inferred MHC-II communication network shows distinct receiver
profiles between Control and CD.

  
**Control:** MHC-II signaling is predominantly associated with TH17
effector cells, a strongly Control-enriched population (1508 Control / 4
CD).

  
**CD:** the inferred MHC-II network involves FOXP3⁺ IL2RA⁺ Treg,
TH17-like innate-associated CD4 T cells, and CST3⁺LYZ⁺ macrophages,
reflecting the emergence of CD-associated receiver populations. Because
several of these populations are strongly condition-skewed or absent in
one condition, this result should not be interpreted as quantitative
amplification of a conserved Treg-centered circuit. Instead, the data
support a remodeling of the cellular context in which MHC-II
communication is detected during CD.

**Literature context.** MHC-II presentation toward regulatory and
TH17-like CD4 T cell populations is consistent with CD40-axis biology
described in intestinal inflammation, where CD103⁺ dendritic cells and
macrophage-lineage antigen-presenting cells shape the local Treg/TH17
balance (Barthels et al. 2017; Wang et al. 2026). The population
substitution observed here — a different receiving Treg-associated
population in CD versus Control, rather than intensification of one
shared population — should be read as a shift in which cell type
occupies this signaling role under inflammation, not as direct evidence
of the CD40 mechanism itself, which CellChat’s ligand-receptor database
does not model explicitly for this pathway.

**Figure 7.** MHC-II circle, bubble, and hierarchy plots.

#### Circle plots

**Control**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CTRL_MHC-II_circle.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CD_MHC-II_circle.png)

### Bubble plot

![](../../results/plots/CellChat_IEL/merge/pathways/bubble/CellChat_IEL_bubble_MHC_II.png)

### Hierarchy plots

**Control**

![](../../results/plots/CellChat_IEL/merge/pathways/hierarchy/CellChat_IEL_hierarchy_MHC_II_Control.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_IEL/merge/pathways/hierarchy/CellChat_IEL_hierarchy_MHC_II_CD.png)

## LIGHT (TNFSF14) — *within-CD characterization*

Within the CD LIGHT network, TH17-like CD4 T cells show the highest
inferred outgoing contribution, followed by CST3⁺LYZ⁺ macrophages. The
LIGHT–LTBR axis includes FOXP3⁺ IL2RA⁺ Treg among the inferred receiver
populations, while the broader LIGHT–TNFRSF14 (HVEM) axis reaches a wide
set of receiver populations including CXCL13⁺ TFH-like T cells,
TH17-like innate-associated CD4 T cells, Naïve CD4, Cycling T cells,
epithelial cells, and CST3⁺LYZ⁺ macrophages themselves (self-signaling).
This should be presented as a within-CD characterization rather than as
a confirmatory differential pathway, given the CD-emergent status of the
Treg receiver population.

**Figure 8.** LIGHT circle plot and expression plots.

#### Circle plot

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CD_LIGHT_circle.png)

#### Gene expression

![](../../results/plots/CellChat_IEL/merge/pathways/violin/CellChat_IEL_geneExpr_LIGHT.png)

**Literature context.** LIGHT/TNFSF14 is associated with IBD
susceptibility and has been reported as elevated in Crohn’s disease. Its
receptors include HVEM and LTβR, with potentially divergent biological
consequences. The detected LIGHT–LTβR signaling involving
macrophage-associated populations in IEL should therefore not be
described simply as pro-inflammatory.

## MIF — *mixed, with CD-associated network remodeling*

The inferred MIF communication pattern differs substantially between
Control and CD. In Control, MIF involves mainly Control-enriched
populations, limiting direct interpretation of a quantitative
cross-condition comparison. In CD, the inferred MIF network becomes
broader, with increased participation of CST3⁺LYZ⁺ macrophages and
additional contribution from GZMK⁺ effector memory CD8 T cells.
CST3⁺LYZ⁺ macrophages represent a prominent node within the
CD-associated MIF network; however, this population is strongly
CD-enriched (35 CD / 1 Control) and should therefore be interpreted
primarily as a within-CD characterization. In contrast, GZMK⁺ effector
memory CD8 T cells are present in both conditions (183 CD / 686 Control)
and represent the most suitable population for direct CD-versus-Control
interpretation of MIF-associated communication changes. Overall, MIF
should be interpreted as a condition-associated network remodeling event
rather than a simple pathway activation comparison.

**Figure 9.** MIF circle plots.

**Control**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CTRL_MIF_circle.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CD_MIF_circle.png)

## GALECTIN — *CD-emergent, epithelial-centered communication pattern*

GALECTIN signaling was detected in the CD IEL network and showed an
epithelial-centered communication structure. In the circle plot,
Epithelial cells / enterocytes represent the dominant visual hub,
including a prominent self-signaling component, indicating autocrine
epithelial GALECTIN communication within the inferred network.

FOXP3⁺ IL2RA⁺ Treg cells participate in the GALECTIN network but do not
represent the dominant signaling hub based on the network visualization.
Given that both epithelial and Treg populations are strongly
condition-dependent in this dataset, this result should be interpreted
as a within-CD communication pattern rather than a differential pathway
activation claim.

**Figure 10.** GALECTIN signaling in IEL under Crohn’s disease (CD).

#### Circle plot

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CD_GALECTIN_circle.png)

#### Bubble plot

![](../../results/plots/CellChat_IEL/merge/pathways/bubble/CellChat_IEL_bubble_GALECTIN.png)

## CDH1 (E-cadherin) — *epithelial-associated CDH1 communication pattern*

CDH1 signaling in the IEL CD network shows a prominent
epithelial-centered structure, with Epithelial cells / enterocytes
displaying the strongest visual hub and autocrine component in the
circle plot.

Although FOXP3⁺ IL2RA⁺ Treg cells contribute to the inferred CDH1
communication network, the current visualization does not support
describing Treg as the exclusive or dominant sender population. The
canonical epithelial CDH1–CD103 (ITGAE/ITGB7) interaction remains
biologically plausible, but sender assignment should follow the complete
CellChat communication table rather than visual interpretation alone.

**Figure 11.** CDH1 (E-cadherin) signaling in IEL under Crohn’s disease
(CD).

#### Circle plot

![](../../results/plots/CellChat_IEL/pathways/CellChat_IEL_CD_CDH1_circle.png)

#### Chord plot

![](../../results/plots/CellChat_IEL/merge/pathways/chord/CellChat_IEL_chord_CDH1_CD.png)

#### Hierarchy plot

![](../../results/plots/CellChat_IEL/merge/pathways/hierarchy/CellChat_IEL_hierarchy_CDH1_CD.png)

> ## Biological interpretation

> CD103–E-cadherin is an established retention mechanism for
> intraepithelial lymphocytes and CD103⁺ Tregs. In this dataset it
> should be presented as a strong, mechanistically anchored
> **CD-population characterization**, not as a differential
> CD-vs-Control claim.

## Signaling-role scatter plots and communication patterns

The signaling-role scatter plots summarize CellChat-inferred changes in
incoming and outgoing communication strength for selected IEL
populations. These outputs provide a network-level overview of pathway
redistribution between conditions but should not be interpreted as
direct measurements of biological activation.

Interpretation is constrained by the strong cell-type imbalance
described in Section 3.1. CD-emergent populations (including FOXP3⁺
IL2RA⁺ Treg, CST3⁺LYZ⁺ macrophages, and TH17-like CD4 T cells) primarily
represent CD-network characterization rather than quantitative
CD-versus-Control comparisons. More balanced populations, such as GZMK⁺
effector memory CD8 T cells and effector T-cell mixed states, provide
more interpretable comparative information.

**Figure 12.** CellChat signaling-role scatter plots for selected IEL
populations.

#### FOXP3⁺ IL2RA⁺ Treg

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_FOXP3_IL2RA_Treg.png)

#### CST3⁺LYZ⁺ macrophages

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_CST3_LYZ_macrophages.png)

#### TH17-like CD4 T cells

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_TH17_like_CD4_T_cells.png)

#### CD39⁺ tissue-resident CD8 T cells

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_CD39_tissue_resident_CD8_T_cells.png)

#### Epithelial cells / Enterocytes

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_Epithelial_cells_enterocytes.png)

#### GZMK⁺ effector memory CD8 T cells

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_GZMK_effector_memory_CD8_T_cells.png)

#### Effector T cell (mixed state)

![](../../results/plots/CellChat_IEL/merge/ct_interest/CellChat_IEL_signalingChanges_Effector_T_cell_mixed_state.png)

**Figure 13.** NMF-derived communication patterns in the IEL network.

NMF decomposition was performed independently for outgoing and incoming
communication programs. Four outgoing patterns (K = 4) and six incoming
patterns (K = 6) were identified. Pattern identities are
condition-specific and should not be interpreted as directly
corresponding between Control and CD. Differences in pattern composition
therefore reflect global communication-network reorganization rather
than a direct comparison of individual latent factors.

#### Outgoing communication patterns

##### Control

**River plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_Control_communicationPatterns_outgoing_river.png)

**Dot plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_Control_communicationPatterns_outgoing_dot.png)

##### Crohn’s disease (CD)

**River plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_CD_communicationPatterns_outgoing_river.png)

**Dot plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_CD_communicationPatterns_outgoing_dot.png)

#### Incoming communication patterns

##### Control

**River plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_Control_communicationPatterns_incoming_river.png)

**Dot plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_Control_communicationPatterns_incoming_dot.png)
\### Crohn’s disease (CD)

**River plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_CD_communicationPatterns_incoming_river.png)

**Dot plot**

![](../../results/plots/CellChat_IEL/merge/pathways/pattern/CellChat_IEL_CD_communicationPatterns_incoming_dot.png)

In the CD incoming dot plot, FOXP3⁺ Treg is the row with the broadest
pathway association. This is a robust description of the Treg-CD
profile, but not a differential claim because Treg is not measurable in
Control.

# LPL Compartment: Pathway-Level Findings

Six pathways were selected: **MHC-I, MHC-II, SIRP, GALECTIN, IL2,
LIGHT**. ADGRE and IL10 were discarded after direct verification via
`subsetCommunication()` showed too few interactions to justify
inclusion.

| Pathway | Real sender | Interactions (Control / CD) | Robustness |
|----|----|---:|----|
| MHC-I | Broad, near-ubiquitous | \>600 rows/condition | Solid |
| MHC-II | Broad in CD, restricted in Control | 83 / 146 | Solid for T follicular helper; descriptive for Treg axis |
| SIRP | Reorganized sender/receiver set | 17 / 24 | Solid; narrate as reorganization |
| GALECTIN | IFN-activated CD8 cytotoxic T cells | 0 / 38 | Moderate, CD-exclusive |
| IL2 | Innate-like/TRM CD8 | 0 / 11 | Weak, receptor-nuance dependent |
| LIGHT | Activated CD4 T cells stress response | 0 / 8 | Fragile; n=15 CD sender cells |

## MHC-I

MHC-I represents a broadly detected signaling pathway in both LPL
conditions. The dominant inferred receiver population is **Activated
cytotoxic CD8 T cells**, a fully comparable population between CD and
Control (5650 and 2002 cells, respectively).

Because this population is present at high abundance in both conditions,
the MHC-I communication pattern represents one of the more robust
pathway-level comparisons in the LPL compartment. The signaling-role
scatter highlights a marked redistribution of incoming MHC-I-associated
signaling toward **Activated effector Treg**, another population that is
quantitatively represented in both conditions (873 CD / 617 Control).

These results indicate a reorganization of MHC-I-associated
communication within the LPL immune network, rather than direct evidence
of altered MHC-I biological activity. CellChat infers ligand–receptor
communication potential and does not measure antigen presentation
efficiency or functional T-cell activation.

**Figure 14.** MHC-I signaling in LPL under Control and Crohn’s disease
(CD), together with the signaling-role change observed for Activated
effector Treg.

#### Circle plots

**Control**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CTRL_MHC-I_circle.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CD_MHC-I_circle.png)

#### Signaling-role change

![](../../results/plots/CellChat_LPL/merge/ct_interest/CellChat_LPL_signalingChanges_Activated_effector_Treg.png)

## MHC-II — *population substitution rather than simple intensification*

MHC-II signaling shows a condition-dependent redistribution of inferred
receiver populations rather than a simple increase of a shared
Treg-associated circuit.

In Control, MHC-II communication is primarily associated with **T
follicular helper cells**, with additional contribution from **Activated
effector Treg** and **Cycling T cells**.

In CD, the receiving profile is reorganized, with **FOXP3⁺ activated
Treg** emerging as a major inferred receiver population, together with
contribution from **IFN-activated CD8 cytotoxic T cells**. Because
FOXP3⁺ activated Treg cells are absent in Control (2205 CD / 0 Control),
this represents a CD-specific network characterization rather than a
quantitative increase of a conserved Treg pathway.

Activated effector Treg, although present in both conditions (873 CD /
617 Control), is no longer among the dominant inferred MHC-II receiver
populations in CD.

Therefore, the main interpretation is a **change in the cellular
identity of the MHC-II receiving compartment**, rather than
amplification of the same Treg signaling circuit across conditions.

**Figure 15.** MHC-II signaling in LPL under Control and Crohn’s disease
(CD), showing the condition-specific receiving profiles and the inferred
signaling hierarchy.

#### Circle plots

**Control**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CTRL_MHC-II_circle.png)
**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CD_MHC-II_circle.png)

#### Signaling hierarchy

**Control**

![](../../results/plots/CellChat_LPL/merge/pathways/hierarchy/CellChat_LPL_hierarchy_MHC_II_Control.png)

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/merge/pathways/hierarchy/CellChat_LPL_hierarchy_MHC_II_CD.png)

## SIRP (SIRPG–CD47) — *network reorganization*

SIRP signaling is detected in both LPL conditions (17 interactions in
Control and 24 in CD) but shows one of the strongest changes in
pathway-level organization between conditions.

Functional similarity analysis indicates that shared pathways occupy a
relatively coherent region of the inferred communication landscape,
whereas pathways with condition-specific representation are positioned
separately within the embedding space. This should be interpreted as a
redistribution of communication-network structure rather than evidence
of a newly activated biological pathway.

Among evaluated pathways, SIRP displays the largest rank-based
displacement between Control and CD, indicating a marked change in its
relative communication pattern across conditions.

Given the T-cell-dominated composition of the LPL compartment, the
observed SIRP pattern is most consistent with altered
SIRPG–CD47-associated immune cell communication rather than classical
myeloid checkpoint signaling.

**Figure 16.** SIRP signaling and pathway-level functional similarity in
LPL under Control and Crohn’s disease (CD).

#### Circle plots

**Control**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CTRL_SIRP_circle.png)
**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CD_SIRP_circle.png)

#### Functional similarity embedding

![](../../results/plots/CellChat_LPL/merge/CellChat_LPL_embeddingPairwise_functional.png)

#### Rank similarity

![](../../results/plots/CellChat_LPL/merge/CellChat_LPL_rankSimilarity_functional.png)

## GALECTIN — *CD-associated, within-CD characterized*

GALECTIN signaling was detected in the CD LPL CellChat network through
38 inferred ligand–receptor interactions, whereas no GALECTIN
interactions were detected in Control.

Direct communication inference should be distinguished from network
topology: the circle plot highlights T follicular helper cells as a
major network hub, including prominent self-signaling, whereas sender
assignment should be derived from the underlying ligand–receptor
communication table rather than visual centrality alone.

Given the strong condition-dependent cellular composition of LPL and the
absence of detected GALECTIN interactions in Control, this result should
be interpreted as a CD-associated communication pattern rather than
evidence of disease-specific GALECTIN induction.

**Figure 17.** GALECTIN signaling in LPL under Crohn’s disease (CD),
showing the inferred CD-exclusive communication pattern and
ligand–receptor interactions.

#### Circle plot

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CD_GALECTIN_circle.png)

## IL2 — *CD-associated, receptor-context dependent*

IL2 signaling was detected only in the CD CellChat network, with 11
inferred ligand–receptor interactions and no detected IL2 interactions
in Control. The only detected sender population was Innate-like/TRM CD8
T cells.

Because this population is strongly imbalanced between conditions (2194
CD cells versus 16 Control cells), the result should be interpreted as a
CD-associated communication pattern rather than a quantitative increase
of IL2 signaling. CellChat supports the presence of an IL2 communication
program within the CD cellular landscape but cannot establish
disease-specific induction.

The inferred receptor usage was not restricted to a canonical
Treg-specific configuration. Predicted IL2 receptor complexes included
IL2RA-containing and IL2RB/IL2RG-containing configurations across
different receiver populations, suggesting that the inferred IL2
signaling landscape involves multiple cellular contexts rather than a
single high-affinity IL2–Treg axis.

**Figure 18.** IL2 signaling in LPL under Crohn’s disease (CD), showing
the inferred CD-exclusive communication network.

#### Circle plot

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CD_IL2_circle.png)

## LIGHT — *fragile, population-limited*

LIGHT signaling was detected only in the CD CellChat network, with the
only inferred sender population being Activated CD4 T cells stress
response. This population is extremely imbalanced between conditions (15
CD cells versus 2628 Control cells), limiting the interpretation of this
interaction.

The only inferred ligand–receptor pair was TNFSF14–TNFRSF14 (HVEM),
representing a distinct receptor context from the IEL LIGHT signaling
pattern, where TNFSF14–LTβR interactions were detected. Given the very
small number of CD sender cells and the strong population imbalance,
this result should be considered hypothesis-generating rather than
evidence of disease-associated LIGHT induction.

**Figure 19.** LIGHT (TNFSF14) signaling in LPL under Crohn’s disease
(CD), showing the inferred LIGHT–HVEM communication axis.

#### Circle plot

**Crohn’s disease (CD)**

![](../../results/plots/CellChat_LPL/pathways/CellChat_LPL_CD_LIGHT_circle.png)

## Communication patterns and cell-type signaling changes

Non-negative matrix factorization (NMF) was used to identify
higher-order communication programs in the LPL compartment. Control and
CD required different outgoing pattern numbers (Control K=4; CD K=8),
whereas incoming communication patterns showed the same decomposition
size (K=7 in both conditions).

Because the outgoing decomposition was performed with different K
values, individual outgoing pattern indices cannot be interpreted as
one-to-one equivalents across conditions. The outgoing results should
therefore be interpreted as differences in the complexity and
organization of inferred signaling programs rather than direct
pattern-to-pattern comparisons.

Incoming patterns, where the same K was selected in both conditions,
allow a more direct qualitative comparison of the overall organization
of receiver-associated signaling programs, while still remaining
dependent on the underlying cell-type composition and CellChat inference
framework.

**Figure 20.** NMF-derived communication patterns in LPL under Control
and Crohn’s disease (CD). Outgoing communication programs are shown as
river plots, while incoming communication programs are summarized using
dot plots. Control outgoing and incoming decompositions use K = 4 and K
= 7, respectively; CD uses K = 8 and K = 7.

#### Outgoing communication patterns

##### Control

**River plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_Control_communicationPatterns_outgoing_river.png)

**Dot plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_Control_communicationPatterns_outgoing_dot.png)

##### Crohn’s disease (CD)

**River plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_CD_communicationPatterns_outgoing_river.png)

**Dot plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_CD_communicationPatterns_outgoing_dot.png)

#### Incoming communication patterns

##### Control

**River plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_Control_communicationPatterns_incoming_river.png)

**Dot plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_Control_communicationPatterns_incoming_dot.png)

##### Crohn’s disease (CD)

**River plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_CD_communicationPatterns_incoming_river.png)

**Dot plot**

![](../../results/plots/CellChat_LPL/merge/pathways/pattern/CellChat_LPL_CD_communicationPatterns_incoming_dot.png)

# Cross-Validation with GSEA

The GSEA finding of **IL2-STAT5 signaling up in CD within LPL-Treg**
finds a directional counterpart in CellChat because IL2 is a pathway
gained in CD. The convergence is **partial, not clean**: receptor usage
is not specifically high-affinity/Treg-selective and the sender is
Innate-like/TRM CD8 rather than Treg.

This should therefore be presented as **directional convergence between
two independent methods**, not as precise mechanistic validation of the
same molecular circuit.

# Known Technical Issues (CellChat v2) — Methods Note

The following behaviors were diagnosed empirically during the LPL
analysis:

1.  **`netVisual_bubble()` errors** when a pathway is absent from the
    second dataset in `comparison`. The comparison index must therefore
    reflect the dataset in which the pathway is actually present.
2.  **`netVisual_bubble()` can silently return an empty plot** when the
    true sender is excluded by `sources.use` / `targets.use`. This
    affected LIGHT, IL2, and GALECTIN in LPL; unrestricted plots were
    therefore regenerated.
3.  **Hierarchy plots** require at least two receiver elements. SIRP,
    IL2, and GALECTIN were excluded from hierarchy visualization where a
    second receiver could not be assigned without biological
    arbitrariness. MHC-I and MHC-II retained hierarchy plots.

# Comparison with External Literature

| Finding | Independent literature | Concordance |
|----|----|----|
| MHC-II / CD40-axis signaling toward Treg and TH17-like populations | CD40 signaling and intestinal Treg/TH17 biology; Barthels et al. 2017; Wang et al. 2026 | Mechanistically consistent |
| LIGHT–LTβR toward Treg/macrophages, IEL | LIGHT association with IBD and LTβR biology | Consistent, with dual-role nuance |
| CDH1(E-cadherin)–CD103, IEL | Established retention mechanism for IEL and CD103⁺ Tregs | Direct mechanistic support |
| SIRPG–CD47 reorganization, LPL | SIRPG/CD47 T-cell costimulation literature | Mechanistically plausible |
| IEL γδ T-cell compositional reduction in CD | Jaeger et al. 2021; Song et al. 2024 | Source-paper concordance + independent replication |
| LP CD8⁺/cytotoxic expansion in CD | Jaeger et al. 2021 | Direct source-paper concordance |

> Correction note

> An earlier internal draft cited a separate “Chen et al. (2021, Nat.
> Commun., PMID 33771991)” as an independent replication cohort. Direct
> verification showed that PMID 33771991 corresponds to Jaeger et
> al. (2021) itself. Therefore, those findings are corroborated by the
> source publication, not by an independent cohort. The genuinely
> independent replication identified for the γδ-IEL reduction is Song et
> al. (2024).

# Cross-Compartment Synthesis

1.  **Network complexity is elevated in CD in both compartments, but
    through different mechanisms.** IEL shows a diffuse expansion (more
    pathways, lower per-edge strength, partly threshold-driven); LPL
    shows a proportional expansion (more interactions and more strength
    together).
2.  **MHC-II is the most consistent disease-associated hub across both
    compartments**, but neither compartment shows a simple
    intensification of a stable circuit. The dominant receiving
    population changes.
3.  **Regulatory T-cell biology is a recurring node of interest** across
    MHC-II, GALECTIN, CDH1 in IEL and MHC-I, MHC-II, IL2 in LPL, with
    centrality, signaling-role and NMF outputs converging on increased
    involvement of regulatory T-cell-associated populations. picture.
4.  **LPL lacks the myeloid/epithelial signaling axes present in IEL**,
    including the MIF-macrophage and CDH1-epithelial analogues.
5.  **TNF/IFNG/IL17/TGFb are largely undetected in CellChat in both
    compartments** (TNF is detected in LPL only), despite prominence in
    GSEA/DESeq2. This is most plausibly a detection-threshold issue and
    must not be interpreted as biological absence.

# Limitations

1.  **Cell-type imbalance** is the dominant methodological constraint
    and governs the validity tier of every claim.
2.  **IEL-Cytotoxic-TRM-like / LPL macropopulation heterogeneity:** the
    IEL cytotoxic TRM-like population aggregates CD8 αβ and γδ lineages
    that could not be statistically separated.
3.  **Detection-threshold artifacts:** absence of IFN-II/IL17/TGFb (and
    of TNF in IEL) should not be interpreted as biological absence.
4.  **Asymmetry in `celltypes_of_interest`:** the IEL deep-dive set is
    weighted toward CD-emergent populations, whereas LPL explicitly
    balances comparable and CD-relevant populations.
5.  **LIGHT (LPL)** rests on an extremely small sender population (n=15)
    and should be treated as hypothesis-generating.
6.  **No independent re-derivation in this report:** quantitative values
    are taken from previously reviewed pipeline outputs.
7.  **LPL session continuity:** scripts 17/18 were not re-verified
    line-by-line in the session that produced the LPL findings log;
    their description follows the established pipeline architecture.
    Script 20 was fully verified.

# Reproducibility

## Environment

| Package  | Version               |
|----------|-----------------------|
| R        | 4.6.0                 |
| Seurat   | 5.5.0                 |
| CellChat | v2 (Jin et al., 2024) |

## Reference scripts

- `13_CellChat_IEL_CTRL.r`
- `14_CellChat_IEL_CD.r`
- `15_CellChat_IEL_Merge.r`
- `16_CellChat_IEL_Comparison.r`
- `17_CellChat_LPL_CTRL.r`
- `18_CellChat_LPL_CD.r`
- `19_CellChat_LPL_Merge.r`
- `20_CellChat_LPL_Comparison.r`

## Reproducibility settings

Random seed: `1234`, called immediately before each stochastic block
(functional similarity embedding, `selectK`,
`identifyCommunicationPatterns`).

`netEmbedding(..., umap.method = "uwot")` was used as the R-native UMAP
backend, avoiding the Python/reticulate dependency.

# Deliverables

- Per-condition, per-compartment CellChat objects and lifted/merged
  comparison objects.
- Global network summary plots.
- RankNet information-flow plots.
- Pathway-specific bubble, circle, gene-expression, chord and hierarchy
  plots.
- Signaling-role scatter plots.
- Communication-pattern recognition outputs.
- Functional similarity embeddings and `rankSimilarity`.
- This consolidated Quarto report.

> ## Recommended manuscript figures

> The strongest candidates identified in the original report are:

> Differential network, both compartments.  
> RankNet, both compartments.  
> MHC-II population substitution, IEL and LPL.  
> CDH1–CD103 axis, IEL.  
> MHC-I → Activated effector Treg, LPL.  
> Functional similarity embedding + `rankSimilarity`, LPL.  

# Conclusion

Cell–cell communication analysis of the IEL and LPL compartments
provides a network-level complement to pseudobulk/GSEA findings,
converging on a small set of well-supported, literature-concordant
observations: a denser and structurally more complex communication
network in CD in both compartments; a recurring, mechanistically nuanced
role for MHC-II and Treg-centered signaling; compartment-specific axes
(LIGHT–LTβR toward Treg in IEL; SIRPG–CD47 reorganization in LPL); and a
transcriptomic support for a CDH1–CD103 interaction consistent with
known IEL retention biology.

The explicit, prospectively applied cell-type validity framework is an
integral part of the analytical design and should be retained in any
future CellChat extension or cross-validation with LIANA.

# Appendix A: Glossary of Abbreviations

| Term           | Definition                                     |
|----------------|------------------------------------------------|
| CD             | Crohn’s disease                                |
| IEL            | Intraepithelial lymphocytes                    |
| LPL            | Lamina propria lymphocytes                     |
| TRM            | Tissue-resident memory (T cell)                |
| APC            | Antigen-presenting cell                        |
| NMF            | Non-negative matrix factorization              |
| L–R            | Ligand–receptor                                |
| MHC-I / MHC-II | Major histocompatibility complex, class I / II |
| GSEA           | Gene Set Enrichment Analysis                   |

# Appendix B: References

1.  Jaeger N., et al. (2021). *Single-cell analyses of Crohn’s disease
    tissues reveal intestinal intraepithelial T cells heterogeneity and
    altered subset distributions*. Nature Communications, 12, 1921.
2.  Jin S., Plikus M.V., Nie Q. (2024). CellChat for systematic analysis
    of cell-cell communication from single-cell transcriptomics. *Nature
    Protocols*, 20, 180–219.
    <https://doi.org/10.1038/s41596-024-01045-4> communication from
    single-cell transcriptomics\*. bioRxiv 2023.11.05.565674.
3.  Barthels C., et al. (2017). *CD40-signalling abrogates induction of
    RORγt⁺ Treg cells by intestinal CD103⁺ DCs and causes fatal
    colitis*. Nature Communications, 8, 14715.
4.  Wang M., et al. (2026). *Fusobacterium nucleatum* drives
    CD40-mediated dendritic cell activation and Th17/Treg imbalance to
    exacerbate intestinal inflammation in Crohn’s disease. Frontiers in
    Immunology.
5.  Fan et al. (2021). Elevated levels of the cytokine LIGHT in Crohn’s
    disease. medRxiv.
6.  Herro R., et al. (2018) and related works on TNFSF14/LTβR in
    intestinal inflammation resolution.
7.  Piccio L., et al. (2005). Adhesion of human T cells to
    antigen-presenting cells through SIRPβ2–CD47 interaction
    costimulates T-cell proliferation. Blood, 105(6), 2421–2427.
8.  Yang et al. (2025). Integrin CD103 expression in naive CD8⁺ T cells
    promotes cytokine-driven acquisition of memory phenotype and
    effector function. Immunity.
9.  Roberts et al. (2021). Integrin αE(CD103)β7 in epithelial cancer.
10. Song et al. (2024). Hyperactivation and enhanced cytotoxicity of
    reduced CD8⁺ γδ T cells in the intestine of patients with Crohn’s
    disease correlates with disease activity. BMC Immunology.
11. Kong et al. (2024). Dysregulation of γδ intraepithelial lymphocytes
    precedes Crohn’s disease-like ileitis. Science Immunology.

# Appendix C: Figure Index

|  \# | Compartment | Main output                            | Section |
|----:|-------------|----------------------------------------|--------:|
|   1 | Both        | Global interaction count/strength      |       4 |
|   2 | Both        | Differential network / heatmap         |       4 |
|   3 | IEL         | RankNet                                |       4 |
|   4 | LPL         | Functional similarity + rankSimilarity |  4, 6.3 |
|   5 | LPL         | RankNet                                |       4 |
|   6 | IEL         | MHC-I circle + hierarchy               |     5.1 |
|   7 | IEL         | MHC-II population substitution         |     5.2 |
|   8 | IEL         | LIGHT–LTβR axis                        |     5.3 |
|   9 | IEL         | MIF                                    |     5.4 |
|  10 | IEL         | GALECTIN                               |     5.5 |
|  11 | IEL         | CDH1–CD103                             |     5.6 |
|  12 | IEL         | Signaling-role scatter                 |     5.7 |
|  13 | IEL         | NMF communication patterns             |     5.7 |
|  14 | LPL         | MHC-I / Activated effector Treg        |     6.1 |
|  15 | LPL         | MHC-II population substitution         |     6.2 |
|  16 | LPL         | SIRP reorganization                    |     6.3 |
|  17 | LPL         | GALECTIN                               |     6.4 |
|  18 | LPL         | IL2 receptor usage                     |     6.5 |
|  19 | LPL         | LIGHT, n=15 sender cells               |     6.6 |
|  20 | LPL         | NMF communication patterns             |     6.7 |
