# -*- coding: utf-8 -*-
"""
Build the final manual-review workbook for the mechanism-heuristic validation
(Reviewer 2, Major point #3).

Inputs
------
  outputs/validation_sample_for_manual_review.csv   the 48 sampled events + heuristic labels
  validation_scv_histories.csv                      per-submitter (SCV) records at the two
                                                    snapshots bracketing each event (server run)

Output
------
  Mechanism_Validation_Review_FINAL.xlsx   six tabs:
     1. Start here            one-screen protocol
     2. Review form           48 events, blinded, three fill-in columns
     3. SCV evidence          per-event prev-vs-curr submission diff (the observed record)
     4. Label definitions     the eight allowed mechanism labels
     5. Heuristic labels      the algorithm's guess — deliberately on its own tab (open last)
     6. Raw SCV data          unmodified server output, for auditing

Design notes
------------
* BLINDING: the heuristic's label is not on the review form. Reviewers record their own judgment
  first; the "agree" column is a formula that resolves itself once a mechanism is entered. This
  lets us state in the response letter that adjudication was blinded to the heuristic's output.
* NO SUGGESTED ANSWERS: the "What changed" column reports only observable facts (which SCVs were
  added/removed/changed). It never names or hints at a mechanism label.
* Every event carries an explicit evidence-quality flag so reviewers know when the workbook alone
  is sufficient and when they must open the ClinVar page.

Re-run this script after a refreshed validation_scv_histories.csv arrives; it regenerates cleanly.
"""
import csv, collections
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Border, Side, Alignment
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation
from openpyxl.formatting.rule import FormulaRule
from openpyxl.workbook.defined_name import DefinedName

REV    = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses")
SAMPLE = REV / "outputs" / "validation_sample_for_manual_review.csv"
OUT    = REV / "Mechanism_Validation_Review_FINAL.xlsx"

# Prefer the gap-filled histories from 24_complete_scv_histories_local.py when they exist;
# fall back to the raw server output otherwise.
_COMPLETE = REV / "validation_scv_histories_COMPLETE.csv"
SCV = _COMPLETE if _COMPLETE.exists() else REV / "validation_scv_histories.csv"

SENTINEL = {"SNAPSHOT_FILE_MISSING", "VARIANT_NOT_IN_SNAPSHOT"}
NOT_CAPTURED = "(not captured)"

# ── palette ───────────────────────────────────────────────────────────────────
NAVY   = PatternFill("solid", fgColor="1F4E79")
SLATE  = PatternFill("solid", fgColor="44546A")
GREY   = PatternFill("solid", fgColor="F2F2F2")
INPUT  = PatternFill("solid", fgColor="FFF2CC")
ADDED  = PatternFill("solid", fgColor="E2EFDA")
REMOVD = PatternFill("solid", fgColor="FCE4E4")
CHANGD = PatternFill("solid", fgColor="FFF0CC")
BLUE   = PatternFill("solid", fgColor="DCE6F1")
WARN   = PatternFill("solid", fgColor="FFE5E5")

WHITE  = Font(color="FFFFFF", bold=True)
BOLD   = Font(bold=True)
DIM    = Font(color="595959")
TITLE  = Font(bold=True, size=14, color="1F4E79")
SMALL  = Font(size=9, color="808080")
THIN   = Border(*[Side(style="thin", color="D9D9D9")] * 4)
CENTER = Alignment(horizontal="center", vertical="center")
WRAP   = Alignment(vertical="top", wrap_text=True)
WRAPC  = Alignment(vertical="center", wrap_text=True)

LABELS = ["reinterpretation", "new_submitter", "conflict_emerged", "conflict_resolved",
          "batch_event", "likely_schema_driven", "likely_knowledge_driven",
          "transition_ambiguous"]

CANON = {"pathogenic": "Pathogenic", "likely pathogenic": "Likely pathogenic",
         "benign": "Benign", "likely benign": "Likely benign",
         "uncertain significance": "Uncertain significance",
         "risk factor": "risk factor", "not provided": "not provided"}


def norm_cls(s):
    s = (s or "").strip()
    return CANON.get(s.lower(), s)


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


# ── load ──────────────────────────────────────────────────────────────────────
sample = read_csv(SAMPLE)
scv_rows = [r for r in read_csv(SCV) if r["scv_id"] not in SENTINEL]

hist = collections.defaultdict(lambda: {"prev": {}, "curr": {}})
for r in scv_rows:
    hist[r["review_id"]][r["role"]][r["scv_id"]] = r


def event_diff(rid, exp_prev=None, exp_curr=None):
    """Return the observed prev->curr submission diff for one event. Facts only.

    exp_prev/exp_curr are the submission counts recorded when the sample was drawn (from the
    release XML). If the evidence source yields fewer rows than that, a submission is missing
    here that existed in the release — say so loudly, because a reviewer who cannot see it may
    conclude nothing happened when in fact a laboratory joined or withdrew.
    """
    p, c = hist[rid]["prev"], hist[rid]["curr"]
    added, removed, shared = sorted(set(c) - set(p)), sorted(set(p) - set(c)), sorted(set(p) & set(c))
    changed = [k for k in shared
               if norm_cls(p[k]["classification"]).lower() != norm_cls(c[k]["classification"]).lower()
               and p[k]["classification"].strip() and c[k]["classification"].strip()]
    total = len(p) + len(c)
    withcls = sum(1 for d in (p, c) for v in d.values() if v["classification"].strip())
    quality = "Full" if total and withcls == total else ("None" if withcls == 0 else "Partial")

    short = []
    for role, have, exp in (("earlier", len(p), exp_prev), ("later", len(c), exp_curr)):
        if exp is not None and have < exp:
            short.append(f"{exp - have} of the {exp} submissions on record at the {role} snapshot "
                         f"could not be recovered")
    if short:
        quality = "Discrepancy"
    return {"prev": p, "curr": c, "added": added, "removed": removed, "shared": shared,
            "changed": changed, "quality": quality, "short": short}


def label_of(scv_id, side):
    r = side.get(scv_id)
    if not r:
        return ""
    sub = (r["submitter"] or "").strip()
    return sub if sub else f"SCV {scv_id}"


def plur(n, one, many=None):
    return f"{n} {one if n == 1 else (many or one + 's')}"


def change_summary(d):
    """Plain factual sentence describing what the record shows. Never names a mechanism."""
    bits = []
    if d["added"]:
        bits.append(plur(len(d["added"]), "submission") + " added")
    if d["removed"]:
        bits.append(plur(len(d["removed"]), "submission") + " removed")
    if not d["added"] and not d["removed"]:
        bits.append("no submissions added or removed")
    n = len(d["shared"])
    if n:
        if d["changed"]:
            bits.append(f"{len(d['changed'])} of the {n} carried-over "
                        f"{'submission' if n == 1 else 'submissions'} changed classification")
        elif d["quality"] == "None":
            bits.append(f"{plur(n, 'submission')} carried over, "
                        f"{'its' if n == 1 else 'their'} "
                        f"{'call' if n == 1 else 'calls'} not captured — check ClinVar")
        else:
            bits.append(f"{plur(n, 'submission')} carried over with unchanged classification")
    out = "; ".join(bits) + "."
    if d.get("short"):
        out += "  NOTE: " + "; ".join(d["short"]) + " — see the ClinVar page."
    return out


QUALITY_NOTE = {
    "Full":    "Full — workbook shows every submitter's call at both snapshots.",
    "Partial": "Partial — some per-submitter calls missing; confirm on ClinVar.",
    "None":    "Membership only — which submissions were added/removed is reliable, "
               "but no per-submitter classifications were captured. Open the ClinVar page.",
    "Discrepancy": "INCOMPLETE — a submission on record could not be recovered, so the list "
                   "below is missing at least one row. Judge this event from the ClinVar page, "
                   "not from this workbook.",
}

wb = Workbook()

# ══ Tab 1 · Start here ════════════════════════════════════════════════════════
ws = wb.active
ws.title = "Start here"
ws.sheet_view.showGridLines = False
ws.column_dimensions["A"].width = 3
ws.column_dimensions["B"].width = 112

LINES = [
    ("title", "Manual validation of the ClinVar reclassification-mechanism heuristic"),
    ("sub",   "Genome Medicine revision · Reviewer 2, Major point #3 · 48 events"),
    ("gap",   ""),
    ("h",     "What we are asking you to do"),
    ("p",     "Our manuscript assigns every reclassification event a mechanism (a submitter reinterpreted the "
              "variant, a new laboratory joined, a conflict emerged, and so on) using a rule-based heuristic that "
              "reads only metadata: dates, submission counts, and review-status wording. Reviewer 2 asked us to "
              "check that heuristic against human judgment. You are the human judgment."),
    ("p",     "For each of the 48 events, decide which single mechanism best explains why the variant's consensus "
              "classification changed, and record it. We then compare your answers to the algorithm's and report "
              "its accuracy with a confidence interval."),
    ("cal",   "We are testing the algorithm, not you. Disagreements are the most valuable rows in this file — "
              "please do not try to guess what the algorithm said."),
    ("gap",   ""),
    ("h",     "How to work through it"),
    ("n1",    "Go to the 'Review form' tab. Each row is one event: a variant whose consensus classification "
              "changed between two monthly ClinVar snapshots (Prev date to Curr date)."),
    ("n2",    "Read the 'What the record shows' column. This is the observed submission-level history — which "
              "submissions were added, removed, or changed their call between the two snapshots. The 'SCV evidence' "
              "tab has the same information broken out submitter by submitter."),
    ("n3",    "Open the ClinVar link when you need more (submitter comments, evidence, cited literature). The "
              "'Evidence' column tells you when the workbook alone is not enough."),
    ("n4",    "Enter your answer in the three yellow columns: the mechanism (dropdown), your confidence, and an "
              "optional note. Definitions of the eight mechanisms are on the 'Label definitions' tab."),
    ("n5",    "The 'Agrees with algorithm' column fills in by itself once you choose a mechanism. That is your "
              "only feedback — the algorithm's actual label sits on the last tab, and we would rather you did not "
              "open it until you have finished."),
    ("gap",   ""),
    ("h",     "Two things worth knowing"),
    ("b",     "Judge from the record, not from the rules. The point is to catch cases where metadata misleads — "
              "a withdrawal and a new submission that cancel out in the counts, a bulk reprocessing that looks "
              "like an individual update, a real evidence change during the 2020 format switch. If nothing fits, "
              "pick the closest label and say why in the notes."),
    ("b",     "Confidence is not a penalty. Mark Low freely when the record is genuinely unclear; we report "
              "accuracy both overall and restricted to confident calls, and 'the record does not say' is a "
              "legitimate finding we can report."),
    ("gap",   ""),
    ("h",     "Effort"),
    ("p",     "Roughly 2 to 4 minutes per event, so about 2 to 3 hours in total. The 48 rows split cleanly "
              "between reviewers — just agree who takes which review_id range and keep the notes column populated."),
    ("gap",   ""),
    ("h",     "When you are done"),
    ("p",     "Save the workbook (keep the .xlsx format) and send it back to Yilun. The scoring script reads the "
              "'Review form' tab directly and produces the per-class precision and recall table and the accuracy "
              "confidence interval that go into the revision."),
    ("gap",   ""),
    ("note",  "Internal working document. Submitter names appear here because you need them to judge; if any part "
              "of this file is ever published as supplementary material, they must first be replaced with "
              "'Laboratory A/B/C…' in line with the manuscript's anonymisation policy."),
]

r = 2
for kind, text in LINES:
    c = ws.cell(r, 2, text)
    if kind == "title":
        c.font = TITLE
    elif kind == "sub":
        c.font = Font(size=10, color="808080")
    elif kind == "h":
        c.font = Font(bold=True, size=11, color="1F4E79")
    elif kind in ("p", "b", "n1", "n2", "n3", "n4", "n5"):
        prefix = {"n1": "1. ", "n2": "2. ", "n3": "3. ", "n4": "4. ", "n5": "5. ",
                  "b": "•  "}.get(kind, "")
        c.value = prefix + text
        c.font = Font(size=10)
    elif kind == "cal":
        c.font = Font(size=10, bold=True, color="1F4E79")
        c.fill = BLUE
    elif kind == "note":
        c.font = Font(size=9, italic=True, color="808080")
    c.alignment = Alignment(vertical="top", wrap_text=True)
    est = max(1, (len(str(c.value)) // 105) + 1)
    ws.row_dimensions[r].height = 15 * est if kind != "title" else 24
    r += 1

# ══ Tab 2 · Review form ═══════════════════════════════════════════════════════
wf = wb.create_sheet("Review form")
wf.sheet_view.showGridLines = False

COLS = [("review_id", 9), ("Gene", 11), ("Variation ID", 12), ("ClinVar link", 15),
        ("Prev date", 10), ("Curr date", 10), ("Classification before", 15),
        ("Classification after", 15), ("Subs before", 8), ("Subs after", 8),
        ("What the record shows (observed)", 46), ("Evidence", 34),
        ("YOUR ANSWER: mechanism", 24), ("Confidence", 12), ("Notes", 34),
        ("Agrees with algorithm", 13)]

wf.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(COLS))
b = wf.cell(1, 1, "Fill in the three yellow columns. Definitions are on the 'Label definitions' tab. "
                  "The last column fills in by itself.")
b.font = Font(bold=True, color="FFFFFF")
b.fill = SLATE
b.alignment = Alignment(horizontal="left", vertical="center")
wf.row_dimensions[1].height = 20

for j, (h, w) in enumerate(COLS, 1):
    c = wf.cell(2, j, h)
    c.fill = INPUT if 13 <= j <= 15 else NAVY
    c.font = Font(bold=True, color="7F6000") if 13 <= j <= 15 else WHITE
    c.border, c.alignment = THIN, Alignment(horizontal="center", vertical="center", wrap_text=True)
    wf.column_dimensions[get_column_letter(j)].width = w
wf.row_dimensions[2].height = 30

for i, s in enumerate(sample):
    rw = i + 3
    rid = s["review_id"]
    d = event_diff(rid, int(s["Prev_Num_Subs"]), int(s["Curr_Num_Subs"]))
    vals = [rid, s["Gene"], s["Variation_ID"], "ClinVar page",
            s["Prev_Date"], s["Curr_Date"], s["Prev_Category"], s["Curr_Category"],
            s["Prev_Num_Subs"], s["Curr_Num_Subs"],
            change_summary(d), QUALITY_NOTE[d["quality"]], "", "", "", ""]
    for j, v in enumerate(vals, 1):
        c = wf.cell(rw, j, v)
        c.border = THIN
        c.alignment = WRAP if j in (11, 12, 15) else WRAPC
        if 13 <= j <= 15:
            c.fill = INPUT
        elif j == 16:
            c.fill = GREY
            c.alignment = CENTER
        else:
            c.fill = GREY
            c.font = DIM
    link = wf.cell(rw, 4)
    link.hyperlink = s["clinvar_link"]
    link.font = Font(color="0563C1", underline="single")
    if d["quality"] != "Full":
        wf.cell(rw, 12).fill = WARN
        wf.cell(rw, 12).font = Font(color="9C0006", size=9)
    else:
        wf.cell(rw, 12).font = Font(color="375623", size=9)
    wf.cell(rw, 11).font = Font(size=9, color="404040")
    # auto-agreement: resolves only after the reviewer commits an answer.
    # Case-insensitive on purpose — a pasted value bypasses the dropdown, and casing
    # should never be the reason a row is scored as a disagreement.
    wf.cell(rw, 16).value = (
        f'=IF(TRIM(M{rw})="","",'
        f'IF(LOWER(TRIM(M{rw}))=LOWER(\'Heuristic labels\'!B{rw}),"Y","N"))')
    wf.row_dimensions[rw].height = 46

last = len(sample) + 2
wf.freeze_panes = "C3"
wf.auto_filter.ref = f"A2:P{last}"

# ══ Tab 3 · SCV evidence ══════════════════════════════════════════════════════
we = wb.create_sheet("SCV evidence")
we.sheet_view.showGridLines = False
EW = [14, 13, 40, 24, 24, 26, 26, 14]
for j, w in enumerate(EW, 1):
    we.column_dimensions[get_column_letter(j)].width = w

we.merge_cells(start_row=1, start_column=1, end_row=1, end_column=8)
c = we.cell(1, 1, "Observed submission-level record for each event, taken from the ClinVar release XML at the two "
                  "snapshots. 'Added' and 'Removed' are always reliable; a blank classification means that field "
                  "was not captured for that snapshot and must be checked on ClinVar.")
c.font = Font(bold=True, color="FFFFFF")
c.fill = SLATE
c.alignment = Alignment(vertical="center", wrap_text=True)
we.row_dimensions[1].height = 30

r = 3
for s in sample:
    rid = s["review_id"]
    d = event_diff(rid, int(s["Prev_Num_Subs"]), int(s["Curr_Num_Subs"]))
    pd_, cd_ = s["Prev_Date"], s["Curr_Date"]

    we.merge_cells(start_row=r, start_column=1, end_row=r, end_column=8)
    h = we.cell(r, 1, f"{rid}  ·  {s['Gene']}  ·  Variation {s['Variation_ID']}      "
                      f"{pd_} → {cd_}      {s['Prev_Category']} → {s['Curr_Category']}      "
                      f"submissions {s['Prev_Num_Subs']} → {s['Curr_Num_Subs']}")
    h.fill, h.font = NAVY, WHITE
    h.alignment = Alignment(vertical="center")
    we.row_dimensions[r].height = 20
    r += 1

    if d["quality"] != "Full":
        we.merge_cells(start_row=r, start_column=1, end_row=r, end_column=8)
        n = we.cell(r, 1, "   " + QUALITY_NOTE[d["quality"]])
        n.fill, n.font = WARN, Font(size=9, color="9C0006", italic=True)
        r += 1

    hdr = ["Status", "SCV ID", "Submitter", f"Classification {pd_}", f"Classification {cd_}",
           f"Review status {pd_}", f"Review status {cd_}", "Last evaluated"]
    for j, t in enumerate(hdr, 1):
        c = we.cell(r, j, t)
        c.fill, c.font, c.border = GREY, Font(bold=True, size=9), THIN
        c.alignment = CENTER
    r += 1

    order = ([("Added", k) for k in d["added"]] + [("Removed", k) for k in d["removed"]] +
             [("Changed call", k) for k in d["changed"]] +
             [("Carried over", k) for k in d["shared"] if k not in d["changed"]])
    for status, k in order:
        pr, cu = d["prev"].get(k), d["curr"].get(k)
        src = cu or pr
        row = [status, k, (src["submitter"] or "").strip() or "(not captured)",
               norm_cls(pr["classification"]) if pr else "—",
               norm_cls(cu["classification"]) if cu else "—",
               (pr["review_status"] or "").strip() if pr else "—",
               (cu["review_status"] or "").strip() if cu else "—",
               (src["last_evaluated"] or "").strip()]
        row = [v if str(v).strip() else NOT_CAPTURED for v in row]
        fill = {"Added": ADDED, "Removed": REMOVD, "Changed call": CHANGD}.get(status, GREY)
        for j, v in enumerate(row, 1):
            c = we.cell(r, j, v)
            c.fill, c.border = fill, THIN
            c.alignment = WRAPC
            c.font = BOLD if (j == 1 and status != "Carried over") else Font(size=9)
            if v == NOT_CAPTURED:
                c.font = SMALL
        r += 1
    r += 1

# ══ Tab 4 · Label definitions ═════════════════════════════════════════════════
wd = wb.create_sheet("Label definitions")
wd.sheet_view.showGridLines = False
for j, w in enumerate([26, 96], 1):
    wd.column_dimensions[get_column_letter(j)].width = w

DEFS = [
    ("section", "Use these five for events outside the 2020 format-transition window", ""),
    ("reinterpretation",
     "An existing submitter changed its own classification on new evidence or updated criteria, and the "
     "consensus moved with it. The set of submitters is usually unchanged."),
    ("new_submitter",
     "A new laboratory added a submission and that shifted the consensus, without creating a formal conflict."),
    ("conflict_emerged",
     "The variant became 'conflicting interpretations' — a submitter's call now disagrees across tiers "
     "(P/LP versus VUS/B/LB) with the others."),
    ("conflict_resolved",
     "A previously conflicting variant reached consensus; the review status moved away from 'conflicting'."),
    ("batch_event",
     "The change is part of a mass update rather than an individual clinical decision — a bulk reprocessing "
     "or a large institutional submission in which many variants moved together. The affected snapshots are "
     "2013-12, 2014-01, 2014-02, 2018-06 and 2020-01."),
    ("section", "Use these three only when the Curr date falls in 2020-01 to 2020-06", ""),
    ("note",
     "In this window ClinVar switched its public export from the RCV format to the VCV format, which changed "
     "how submissions are aggregated. A consensus category can therefore appear to change with no new clinical "
     "evidence behind it. Decide whether the record shows a real evidence change or an artefact of the switch.", ""),
    ("likely_schema_driven",
     "The change looks like an artefact of the RCV-to-VCV aggregation switch — submission counts or review-status "
     "wording shift at the format boundary with no genuine new interpretation."),
    ("likely_knowledge_driven",
     "A genuine evidence-based update that happened to coincide with the format switch; the underlying "
     "submissions or evidence really did change."),
    ("transition_ambiguous",
     "The record does not let you tell schema from evidence. Use it rather than guessing."),
]

wd.merge_cells(start_row=1, start_column=1, end_row=1, end_column=2)
c = wd.cell(1, 1, "The eight mechanisms — enter exactly one per event")
c.font, c.fill = Font(bold=True, size=12, color="FFFFFF"), NAVY
c.alignment = Alignment(vertical="center")
wd.row_dimensions[1].height = 22

r = 3
for item in DEFS:
    kind, *rest = item
    if kind == "section":
        wd.merge_cells(start_row=r, start_column=1, end_row=r, end_column=2)
        c = wd.cell(r, 1, rest[0])
        c.font, c.fill = Font(bold=True, color="1F4E79"), BLUE
        c.alignment = Alignment(vertical="center")
        wd.row_dimensions[r].height = 20
    elif kind == "note":
        wd.merge_cells(start_row=r, start_column=1, end_row=r, end_column=2)
        c = wd.cell(r, 1, rest[0])
        c.font = Font(size=9, italic=True, color="595959")
        c.alignment = WRAP
        wd.row_dimensions[r].height = 42
    else:
        a = wd.cell(r, 1, kind)
        a.font, a.fill, a.border, a.alignment = Font(bold=True, name="Consolas", size=10), GREY, THIN, WRAPC
        bcell = wd.cell(r, 2, rest[0])
        bcell.border, bcell.alignment, bcell.font = THIN, WRAP, Font(size=10)
        wd.row_dimensions[r].height = 15 * max(1, len(rest[0]) // 92 + 1) + 8
    r += 1

r += 1
wd.merge_cells(start_row=r, start_column=1, end_row=r, end_column=2)
c = wd.cell(r, 1, "If none of the eight fits well, choose the closest, mark your confidence Low, and explain in "
                  "the Notes column. A documented poor fit is more useful to us than a forced label.")
c.font, c.fill = Font(size=10, bold=True, color="1F4E79"), BLUE
c.alignment = WRAP
wd.row_dimensions[r].height = 32

# Contiguous vocabulary block -> named range 'Mechanisms'. Driving the dropdown from real cells
# (rather than an inline string) keeps the list editable, dodges Excel's 255-character limit on
# inline validation lists, and gives the paste-guard something to COUNTIF against.
r += 2
c = wd.cell(r, 1, "Allowed values")
c.font, c.fill, c.alignment = Font(bold=True, color="1F4E79"), BLUE, Alignment(vertical="center")
c2 = wd.cell(r, 2, "This list drives the dropdown on the Review form. Do not edit it.")
c2.font, c2.fill = Font(size=9, italic=True, color="595959"), BLUE
vocab_first = r + 1
for k, lab in enumerate(LABELS):
    c = wd.cell(vocab_first + k, 1, lab)
    c.font, c.border, c.alignment = Font(name="Consolas", size=10), THIN, WRAPC
vocab_last = vocab_first + len(LABELS) - 1
wb.defined_names.add(DefinedName(
    "Mechanisms", attr_text=f"'Label definitions'!$A${vocab_first}:$A${vocab_last}"))

# ── validation + guards on the Review form (needs 'Mechanisms' to exist first) ────────────────
dv_m = DataValidation(type="list", formula1="Mechanisms", allow_blank=True,
                      showErrorMessage=True, errorTitle="Use a listed mechanism",
                      error="Pick one of the eight labels from the dropdown "
                            "(they are defined on the 'Label definitions' tab).")
dv_c = DataValidation(type="list", formula1='"High,Medium,Low"', allow_blank=True,
                      showErrorMessage=True, errorTitle="Confidence",
                      error="Choose High, Medium, or Low.")
wf.add_data_validation(dv_m)
wf.add_data_validation(dv_c)
dv_m.add(f"M3:M{last}")
dv_c.add(f"N3:N{last}")

# Paste-guard: Excel silently accepts pasted values that violate a dropdown, so flag anything
# outside the vocabulary in red rather than letting it reach the scorer unnoticed.
wf.conditional_formatting.add(
    f"M3:M{last}",
    FormulaRule(formula=[f'AND(M3<>"",COUNTIF(Mechanisms,M3)=0)'],
                fill=PatternFill("solid", fgColor="FFC7CE"),
                font=Font(color="9C0006", bold=True), stopIfTrue=False))
# Quiet visual feedback on the auto-agreement column.
wf.conditional_formatting.add(
    f"P3:P{last}",
    FormulaRule(formula=['P3="Y"'], fill=PatternFill("solid", fgColor="E2EFDA"),
                font=Font(color="375623", bold=True)))
wf.conditional_formatting.add(
    f"P3:P{last}",
    FormulaRule(formula=['P3="N"'], fill=PatternFill("solid", fgColor="FFF0CC"),
                font=Font(color="9C6500", bold=True)))

# ══ Tab 5 · Heuristic labels ══════════════════════════════════════════════════
wh = wb.create_sheet("Heuristic labels")
wh.sheet_view.showGridLines = False
for j, w in enumerate([10, 26, 62], 1):
    wh.column_dimensions[get_column_letter(j)].width = w
wh.merge_cells(start_row=1, start_column=1, end_row=1, end_column=3)
c = wh.cell(1, 1, "The algorithm's answers — please finish the Review form before reading this tab.")
c.font, c.fill = Font(bold=True, color="FFFFFF"), PatternFill("solid", fgColor="9C0006")
c.alignment = Alignment(vertical="center")
wh.row_dimensions[1].height = 22
for j, t in enumerate(["review_id", "heuristic_label", "signals the algorithm used"], 1):
    c = wh.cell(2, j, t)
    c.fill, c.font, c.border, c.alignment = NAVY, WHITE, THIN, CENTER
# row alignment with the review form is load-bearing: the agree formula reads B<same row>
for i, s in enumerate(sample):
    rw = i + 3
    for j, v in enumerate([s["review_id"], s["heuristic_label"], s["heuristic_signals"]], 1):
        c = wh.cell(rw, j, v)
        c.border, c.font, c.alignment = THIN, Font(size=9), WRAPC
wh.freeze_panes = "A3"

# ══ Tab 6 · Raw SCV data ══════════════════════════════════════════════════════
wr = wb.create_sheet("Raw SCV data")
wr.sheet_view.showGridLines = False
raw = read_csv(SCV)
hdr = list(raw[0].keys()) if raw else []
for j, h in enumerate(hdr, 1):
    c = wr.cell(1, j, h)
    c.fill, c.font, c.border, c.alignment = NAVY, WHITE, THIN, CENTER
    wr.column_dimensions[get_column_letter(j)].width = {"submitter": 42, "condition": 34,
                                                        "review_status": 30}.get(h, 14)
for i, row in enumerate(raw, 2):
    for j, h in enumerate(hdr, 1):
        c = wr.cell(i, j, row[h])
        c.border, c.font, c.alignment = THIN, Font(size=9), WRAPC
wr.freeze_panes = "A2"
wr.auto_filter.ref = f"A1:{get_column_letter(len(hdr))}{len(raw) + 1}"

wb.save(OUT)

# ── report ────────────────────────────────────────────────────────────────────
q = collections.Counter(event_diff(s["review_id"], int(s["Prev_Num_Subs"]),
                                   int(s["Curr_Num_Subs"]))["quality"] for s in sample)
print(f"wrote {OUT.name}")
print(f"  Review form   : {len(sample)} events, inputs M/N/O, agreement auto-computed in P")
print(f"  SCV evidence  : {len(scv_rows)} submission rows across {len(sample)} events")
print(f"  evidence quality — Full {q['Full']}, Partial {q['Partial']}, Membership only {q['None']}")
