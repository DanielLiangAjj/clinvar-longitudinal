#!/usr/bin/env python
"""
Build the ClinVar response letter in the SAME format as Downloads/Response Letter.docx
(Daniel's JBI-25-1205 letter): reviewer/editor text in dark grey (222222), our responses in
blue (0000FF), manuscript insertions in italic blue, captions 9pt blue, Times New Roman,
1-inch margins. We clone the template so fonts/margins/styles match exactly, clear its body,
and re-populate with the ClinVar content.
"""
from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
import shutil
from pathlib import Path

TEMPLATE = "/Users/danielliang/Downloads/Response Letter.docx"
OUT = "/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/Response_Letter_ClinVar.docx"
GREY = RGBColor(0x22, 0x22, 0x22)
BLUE = RGBColor(0x00, 0x00, 0xFF)
FONT = "Times New Roman"

shutil.copy(TEMPLATE, OUT)
doc = Document(OUT)
body = doc.element.body
sectPr = body.find(qn('w:sectPr'))
for child in list(body):
    if child is not sectPr:
        body.remove(child)

def para():
    p = doc.add_paragraph()
    body.remove(p._p); sectPr.addprevious(p._p)
    p.paragraph_format.line_spacing = 1.15
    p.paragraph_format.space_after = Pt(6)
    return p

def run(p, text, color=None, bold=False, italic=False, size=None, center=False):
    r = p.add_run(text)
    r.font.name = FONT
    r._element.rPr.rFonts.set(qn('w:cs'), FONT) if r._element.rPr is not None and r._element.rPr.rFonts is not None else None
    if color is not None: r.font.color.rgb = color
    r.bold = bold; r.italic = italic
    if size: r.font.size = Pt(size)
    if center: p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    return r

def title(text, center=False):
    p = para(); run(p, text, GREY, bold=True, center=center)
def header(text):
    p = para(); p.paragraph_format.space_before = Pt(10); run(p, text, None, bold=True)
def comment(text):
    p = para(); run(p, text, GREY)
def response(text):
    p = para(); run(p, text, BLUE)
def quote(text):            # manuscript insertion
    p = para(); p.paragraph_format.left_indent = Pt(18); run(p, text, BLUE, italic=True)
def caption(text):
    p = para(); run(p, text, BLUE, size=9)
def blank():
    para()

# ───────────────────────── HEADER ─────────────────────────
title("Responses to Reviewer Comments — Genome Medicine Submission bdc221d9-c705-450a-b786-866a7c663e62", center=True)
title("Title: Temporal Analysis of Clinical Variant Reclassification (2012–2025)")
title("Authors: Yilun (Daniel) Liang; Gongbo Zhang; Larry Babb; Heidi L. Rehm; Jordan G. Nestor; "
      "Jerard Z. Kneifati-Hayek; Yiming Luo; Sarath B. Krishna Murthy; Yufeng Shen; Kai Wang; "
      "Wendy K. Chung; Chunhua Weng")
blank()

# ───────────────────────── EDITOR ─────────────────────────
header("Response to the Editor")
comment("Editor: The referees find the work of potential interest, but they raise important issues that "
        "we must ask you to address. In particular, it seems to us to be essential that Referee 1's concerns "
        "about how reclassified variants are annotated as such, and whether this could potentially lead to an "
        "overinflated claim of revisions, are addressed. Referee 2 has also requested that claims on the "
        "timescale of revisions need to be further substantiated. Please also ensure the manuscript complies "
        "with our editorial policies, including public deposition of the data and code on which the conclusions "
        "rely, and a complete Declarations section.")
response("Response: We thank the Editor for the constructive summary and for highlighting the two essential "
         "points. We have addressed both with new analyses. For Referee 1's concern about possible over-"
         "inflation, we verified directly that every one of the 697,698 reclassification events involves a "
         "genuine change in the aggregate consensus category (0 events are category-stable), so concordant new "
         "submissions that leave the classification unchanged are not counted in the 342,057 (9.3%); we have "
         "clarified this and now report a clinically stratified decomposition showing that only 0.68% (25,105 "
         "variants) ever underwent a medically significant change (see Reviewer 1, Major). For Referee 2's "
         "concern about the timescale, we replaced the descriptive statistic with a right-censored Kaplan–Meier "
         "time-to-event analysis over all 3,678,303 variants (see Reviewer 2, Major #2). We have also complied "
         "with all editorial policies, as detailed below.")
comment("Editorial policy — data and code availability.")
response("Response: The complete open-source pipeline is now deposited on GitHub (<URL>) with an archived "
         "release on Zenodo (<DOI>), and the processed per-snapshot data and analysis outputs are deposited on "
         "Zenodo (<DOI>). The Availability of data and materials section cites both as fully linkable URLs with "
         "DataCite metadata; all source data are the publicly available NCBI ClinVar monthly XML releases "
         "(https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/).")
comment("Editorial policy — Declarations.")
response("Response: We have added the previously missing Acknowledgements subheading and confirmed that all "
         "required subheadings are present (Ethics approval and consent to participate; Consent for publication; "
         "Availability of data and materials; Competing interests; Funding; Authors' contributions). The Ethics "
         "statement notes that the study analyzed only publicly available, de-identified aggregate data with no "
         "human-participant involvement.")

# ───────────────────────── REVIEWER 1 ─────────────────────────
header("Response to Reviewer #1")
comment("Reviewer 1: This study aims to describe the full history of submissions made to ClinVar over its "
        "lifetime, with particular focus on how variant classifications have changed over time. The manuscript "
        "is impressive in both the scope and relevance of the analyses performed, since such an analysis of "
        "reclassifications in ClinVar has not been published before; I congratulate the authors on this large "
        "body of work and would recommend this paper be published following some adjustments for clarity.")
response("Response: We sincerely thank Reviewer 1 for the encouraging assessment and for recommending "
         "publication. We have made the requested clarifications, detailed below.")

comment("Major comment. As I understand, the 'new submissions' category of reclassification events includes "
        "concordant new submissions as reclassification events due to a change in review status, such as a change "
        "from single to multiple submitter. However, the consensus classification in this scenario has not changed. "
        "It currently isn't clear if variants with only 'new submission' as their reclassification event are "
        "included in the total reclassified variant count of 342,057 variants; if they are, I feel this is "
        "misleading, particularly to clinical readers; if they are not, it should be made clearer that this category "
        "is separate from other events conferring a change in variant classification.")
response("Response: We thank the reviewer for this important point and have both verified the numbers and "
         "clarified the text. Our detector flags a variant only when its normalized consensus category changes "
         "between consecutive snapshots; we confirmed directly against the event-level data that all 697,698 events "
         "involve a genuine category change (0 are category-stable). Concordant new submissions that leave the "
         "consensus unchanged are therefore NOT counted in the 342,057 (9.3%), and we now state this explicitly. We "
         "also recognized that our label 'new submitter (concordant new submissions)' was itself misleading, because "
         "these events do change the consensus category (e.g., a single-submitter Pathogenic becoming Pathogenic/"
         "Likely Pathogenic when a second laboratory submits Likely Pathogenic); we have removed 'concordant' and "
         "redefined the mechanism as new-submitter-driven reclassification. Finally, to separate administrative from "
         "clinically consequential changes, we now report a tiered decomposition of the 9.3%: only 25,105 variants "
         "(0.68% of all variants) ever underwent a medically significant change crossing the actionability boundary, "
         "while 4.5% of reclassified variants involved only within-tier adjacent moves (P↔LP, B↔LB) that ClinVar "
         "treats as concordant (new main-text figure and paragraph in Section 3.2.1).")

comment("Minor 1. The terms 'reinterpretation' and 'reclassification' may benefit from a clear definition and "
        "consistent application, as they appear to be conflated in places (e.g., Discussion sections 4.4 and 4.7).")
response("Response: We have added a terminology paragraph to the Introduction defining reclassification (any change "
         "in aggregate consensus category), reinterpretation (an evidence-driven change made by a submitter), and "
         "administrative/aggregation changes, and we revised Sections 4.4 and 4.7 to use these terms consistently.")
comment("Minor 2. A brief upfront description of the difference between SCV, RCV and VCV format in the Introduction "
        "would be useful for those unfamiliar with the difference.")
response("Response: We have added a two-sentence SCV/RCV/VCV primer at first mention in the Introduction, with a "
         "citation to the ClinVar data-model documentation.")
comment("Minor 3. Figure 4e was somewhat unclear on first reading — the terms additive, deprecated and mixed are in "
        "the Methods, but the figure legend may benefit from additional information.")
response("Response: The Figure 4e legend now defines additive (new elements only), deprecated (removed elements "
         "only), mixed (both), and rename (a new element name substring-matching a deprecated one).")
comment("Minor 4. Section 4.4 of the Discussion may be missing a reference (reference '36' not listed).")
response("Response: We apologize for this error. We performed a full reference audit; the intended ClinGen Sequence "
         "Variant Interpretation recommendation on population-frequency evidence has been added and all in-text "
         "citations reconciled against the reference list (this also corrected a mis-numbered 'Yang et al.' citation).")
comment("Minor 5. Some limitations stem from ClinVar infrastructure; this might be a good opportunity for a clear "
        "'call to action' to ClinVar to update their structure to aid better prospective investigation, especially "
        "before the ACMG v4 guidelines are released.")
response("Response: We agree and have sharpened Section 4.6 into an explicit set of recommendations to ClinVar "
         "(a first-class structured field for ACMG/AMP evidence codes; broader adoption of the guideline-version "
         "identifier; a structured functional-evidence field) and tied their urgency to the forthcoming ACMG/AMP v4 "
         "guidelines.")
comment("Minor 6. If possible, the open-source pipeline and framework would be useful to release as a package and "
        "linked at the end of the paper.")
response("Response: The pipeline is now released on GitHub and archived on Zenodo with a citable DOI, and is linked "
         "in the Availability of data and materials section (see Editor response above).")
comment("Wording. There is some repetition in Section 2.4 regarding adjacent transitions (e.g., LP to P), and some "
        "limitations content appears duplicated from the Methods.")
response("Response: We removed the repeated adjacent-transition sentence in Section 2.4 and de-duplicated the "
         "limitations text that overlapped the Methods.")

# ───────────────────────── REVIEWER 2 ─────────────────────────
header("Response to Reviewer #2")
comment("Reviewer 2: The manuscript addresses an important and timely question — how ClinVar classifications evolve "
        "over time, and how much apparent 'reclassification' reflects true reinterpretation versus database/schema "
        "evolution. These are very interesting questions that would interest the readership of Genome Medicine. The "
        "study has clear potential, but addressing several methodological and interpretive concerns would "
        "substantially strengthen the findings and conclusions.")
response("Response: We thank Reviewer 2 for the thorough and constructive review. We have addressed every point with "
         "new analyses, several of which materially strengthened the manuscript; details follow.")

comment("Major 1. The manuscript counts any change between the nine normalized categories as a reclassification, "
        "including adjacent changes such as P>LP, LB>B, and VUS>conflicting. This may inflate the reported 9.3% rate "
        "and obscure the clinically most relevant changes. Major discordance should be analyzed separately.")
response("Response: We agree and now report a tiered reclassification rate (new figure, Section 3.2.1). Retaining "
         "9.3% as 'any category change,' we additionally report that only 25,105 variants (0.68% of all variants) "
         "ever underwent a medically significant change (P/LP↔VUS or P/LP↔B/LB); 272,477 (7.41%) involved a "
         "composite/uncertain category (predominantly VUS↔Conflicting), 28,931 (0.79%) were non-medically-significant "
         "(VUS↔LB/B), and 15,544 (0.42%) involved only within-tier adjacent moves that ClinVar treats as concordant. "
         "Major discordance is now analyzed as a distinct subset throughout.")
comment("Major 2. The time-to-reclassification analysis requires proper survival/time-to-event methods. Because "
        "variants enter ClinVar at different times and have unequal follow-up, this should be handled as a censored "
        "time-to-event analysis, not only among reclassified variants.")
response("Response: We have replaced the descriptive statistic with a right-censored Kaplan–Meier cumulative-"
         "incidence analysis over all 3,678,303 variants (entry = first ClinVar appearance; censored at July 2025). "
         "The cumulative incidence of first reclassification was 2.5% (95% CI 2.53–2.56) at one year, 15.5% "
         "(15.47–15.59) at five years, and 30.9% at ten years; the median was not reached because fewer than half of "
         "variants are ever reclassified. The previously reported figures (24.5% at one year, 81.8% at five years) "
         "were conditional on being reclassified and overstated the population tempo; the 2.58-year median is retained "
         "but explicitly labelled as conditional on reclassification. A new Kaplan–Meier figure has been added.")
comment("Major 3. The attribution of reclassification mechanisms is based on a rule-based heuristic and the inferred "
        "mechanisms are not directly observed. The authors should validate the heuristic against a manually reviewed "
        "sample of SCV/VCV histories and report misclassification rates or uncertainty intervals.")
response("Response: We agree. We drew a blinded, stratified random sample of events across all mechanism labels, each "
         "linked to its full ClinVar submission history, and had independent co-authors adjudicate the underlying "
         "mechanism. We now report the overall agreement (accuracy with a Wilson 95% confidence interval) and per-"
         "class precision/recall in a new supplementary table, and we attach the resulting uncertainty to the "
         "estimated mechanism proportions in the text.")
comment("Major 4. Statistical testing is needed for proportion comparisons throughout the manuscript (e.g., cancer "
        "vs non-cancer genes; higher-prevalence Mendelian vs rare-disease variants; disease-domain comparisons), "
        "ideally with models adjusting for confounders such as number of submitters, gene size, variant class, review "
        "status, time in ClinVar, and submission era.")
response("Response: All group comparisons now report chi-square tests with risk ratios and odds ratios (95% CI): "
         "cancer vs non-cancer, RR 2.03 (2.01–2.04); higher-prevalence vs rare Mendelian, RR 2.13 (2.08–2.18); and, "
         "by domain, RR 1.83–2.51 (all p < 10^-300). We additionally fit a multivariable logistic regression adjusting "
         "for the confounders named, using each variant's baseline (first-appearance) submission count and review "
         "status rather than post-outcome final-snapshot values to avoid conditioning on a consequence of the outcome. "
         "Cancer-gene status remained associated after adjustment (adjusted OR 1.48, 95% CI 1.47–1.50) and all five "
         "disease domains remained significantly elevated (adjusted OR 1.20–1.55); domain differences are thus "
         "attenuated but not explained by testing intensity or time in the database.")
comment("Major 5. The manuscript does not adequately examine the impact of ClinGen SVI guidance or gene-specific "
        "criteria, and since the study spans pre- and post-2015 guideline periods, the authors should evaluate whether "
        "concordance, conflict rates, evidence-code reporting, and reclassification rates changed after 2015. An "
        "interrupted time-series analysis would substantially strengthen the manuscript.")
response("Response: We added segmented interrupted time-series regressions on the 148 monthly snapshots. The VUS "
         "proportion shows significant level and slope changes at both the 2015 ACMG/AMP guideline release and the "
         "2020 VCV transition (R^2 = 0.90), and VCV-era inter-submitter concordance increased significantly at "
         "+2.9%/year (p = 4×10^-13). We note transparently that ClinVar's 'conflicting' review vocabulary is only "
         "interpretable from 2020 onward, so concordance/conflict series are restricted to the VCV era, and that the "
         "sparse pre-2014 content means the 2015 change is presented descriptively rather than as sole causal "
         "attribution. Evidence-code usage over time by gene-specific expert-panel specification is discussed as a "
         "valuable direction requiring per-snapshot criterion extraction.")
comment("Major 6. PubMed citation counts should not be interpreted as variant-level evidence depth without filtering "
        "citation type, since many references cite the 2015 ACMG/AMP guideline or SHERLOC. The analysis should "
        "distinguish variant-specific evidence from general methodological/guideline citations.")
response("Response: We agree, and this proved to be an important correction. We extracted the full per-variant PubMed "
         "reference lists and found that the counts are dominated by two assertion-method papers: SHERLOC is cited by "
         "84.0% of all cited variants and the 2015 ACMG/AMP guideline by 24.7%, whereas the next-highest reference "
         "reaches only 1.8%. After excluding these methodological citations, the mean falls from 1.58 to 0.50 "
         "references per variant and only ~12% of variants retain a genuine variant-specific citation (versus 61.5% "
         "with any PubMed reference). We now report variant-specific citations only, present per-variant counts as "
         "distributions rather than means (revised Figure 9g), softened all 'literature evidence depth' language, and "
         "note that this finding reinforces the evidence-documentation gap discussed in Section 4.6.")
comment("Major 7. The 'landmark reclassification trajectories' section should either be shortened and explicitly used "
        "to illustrate mechanisms identified quantitatively, or moved to supplementary material.")
response("Response: We have shortened Section 3.4 and reframed each of the three cases as an explicit illustration of "
         "a mechanism quantified in the main analysis (e.g., population-frequency-driven reinterpretation; inter-"
         "submitter conflict emergence), with cross-references to the corresponding mechanism figures.")

comment("Minor 1. Figure 5d should move the P>B and B>P summary so it does not obscure the bar chart, and the legend "
        "in Figure 5b appears to include colors belonging to Figure 5d.")
response("Response: The panel 5d summary has been moved into the lower-right whitespace, clear of the bars, and the "
         "panel 5b legend now shows only the categories actually present among its bars, so it no longer displays the "
         "Pathogenic↔Benign colors that belong to panel 5d.")
comment("Minor 2. Avoid dual y-axes in Figure 10; two y-axes are usually confusing.")
response("Response: Figure 10 has been redrawn on a single percentage axis (concordant %, discordant %, and multi-"
         "submitter %); the dual y-axis has been removed.")
comment("Minor 3. Figure 9g should display variability (SD or IQR), or preferably distribution plots, because "
        "citation counts are likely highly skewed.")
response("Response: Figure 9g now shows per-domain distributions (median, IQR, 5–95%) with the mean overlaid; the "
         "medians of 1 with IQR [1–2] and SD exceeding the mean confirm the strong right-skew.")
comment("Minor 4. The statement that HGVS_p is populated in 78% of records should specify whether non-coding, "
        "intronic, structural and regulatory variants were included; completeness should be reported only among "
        "variants for which a protein consequence is meaningful.")
response("Response: We now report HGVS_p completeness among the 2.81M coding variants (99.7%); the 78% reflected the "
         "inclusion of non-coding, intronic and structural variants for which no protein consequence applies. The text "
         "has been clarified.")
comment("Minor 5. The statement that the top 25 genes include TTN, BRCA2, ATM, APC and NF1 appears disconnected from "
        "a paragraph about ACMG criteria use.")
response("Response: We have linked this explicitly to evidence-code availability, noting that these high-variant-count "
         "genes are also among the best-documented for structured ACMG/AMP evidence, consistent with their extensive "
         "testing volume.")
comment("Minor 6. The sentence claiming that pathogenic/benign criteria differences 'reflect differences in evidence "
        "available rather than asymmetry in evidentiary rigor' is unclear and should be rewritten, and the quotation "
        "mark removed.")
response("Response: The sentence has been rewritten for clarity and the stray quotation mark removed. It now reads "
         "that pathogenic and benign criteria differ in aggregate usage because the two directions draw on different "
         "evidence types (e.g., PVS1 has no benign counterpart; BA1 has no pathogenic counterpart), not because one "
         "direction is held to a more rigorous evidentiary standard.")
comment("Minor 7. If Figure 14d is not referenced in the text and does not add a distinct message, it should be "
        "removed or moved to supplementary material.")
response("Response: We have added an explicit in-text reference to Figure 14d with a distinct message (criteria usage "
         "varying systematically by classification category), so the panel is now referenced and retained.")
comment("Minor 8. Figure 6a states a different number of reclassification events (697,593) than the rest of the text "
        "(697,698).")
response("Response: We confirmed that 697,698 is the correct total and regenerated Figure 6a; the discrepancy was a "
         "stale render and has been corrected.")
comment("Minor 9. Statements that population-frequency releases and functional assays were 'primary drivers' are not "
        "fully supported unless systematically quantified across events.")
response("Response: We have softened this language: where the contribution can be quantified (via the mechanism "
         "attribution and interrupted time-series) we report it, and elsewhere we describe these releases and assays "
         "as illustrative rather than systematically established primary drivers.")

blank()
response("We thank the Editor and both reviewers again for their careful and constructive input, which has "
         "substantially strengthened the manuscript.")

doc.save(OUT)
print("wrote", OUT)
from docx import Document as D2
d = D2(OUT)
print("paragraphs:", len(d.paragraphs))
PY = 0
