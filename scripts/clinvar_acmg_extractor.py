#!/usr/bin/env python3
"""
ClinVar ACMG/AMP Assertion Criteria Extractor
==============================================
Standalone script that reads v5 checkpoint files and performs comprehensive
two-tier extraction of ACMG/AMP evidence criteria from submission data.

Tier 1 — Explicit codes: Enhanced regex capturing strength modifiers
         (PVS1_Strong, PM2_Supporting, PM2:moderate, etc.)
Tier 2 — NLP-inferred codes: ~50 keyword/phrase patterns mapped to ACMG
         criteria with negation detection.

Requirements:
    Python 3.7+  (no external packages needed)

Usage:
    python3 clinvar_acmg_extractor.py \
        --checkpoint-dir ./clinvar_output/checkpoints \
        --output acmg_criteria_detailed.csv
"""

import os
import sys
import re
import csv
import gzip
import pickle
import time
import argparse
from pathlib import Path
from collections import Counter, defaultdict


# ═══════════════════════════════════════════════════════════════════════
# ACMG/AMP CRITERIA DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

# Canonical ACMG/AMP criteria codes
ALL_ACMG_CODES = [
    "PVS1",
    "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1",
    "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
]

# Strength classes for grouping
STRENGTH_CLASS = {}
for code in ALL_ACMG_CODES:
    if code.startswith("PVS"):
        STRENGTH_CLASS[code] = "Very Strong Pathogenic"
    elif code.startswith("PS"):
        STRENGTH_CLASS[code] = "Strong Pathogenic"
    elif code.startswith("PM"):
        STRENGTH_CLASS[code] = "Moderate Pathogenic"
    elif code.startswith("PP"):
        STRENGTH_CLASS[code] = "Supporting Pathogenic"
    elif code.startswith("BA"):
        STRENGTH_CLASS[code] = "Stand-alone Benign"
    elif code.startswith("BS"):
        STRENGTH_CLASS[code] = "Strong Benign"
    elif code.startswith("BP"):
        STRENGTH_CLASS[code] = "Supporting Benign"


# ═══════════════════════════════════════════════════════════════════════
# TIER 1: ENHANCED EXPLICIT CODE EXTRACTION
# ═══════════════════════════════════════════════════════════════════════

# Base ACMG code pattern — matches PVS1, PS1, PM2, PP3, BA1, BS1, BP4, etc.
_BASE_CODE = r'(?:PVS[1-9]|PS[1-9]\d?|PM[1-9]\d?|PP[1-9]\d?|BA[1-9]|BS[1-9]\d?|BP[1-9]\d?)'

# Enhanced pattern capturing optional strength modifiers:
#   PVS1_Strong, PM2_Supporting, PS3:moderate, BP4 (supporting), PM2_P, etc.
_STRENGTH_SUFFIX = (
    r'(?:'
    r'[_:\s]*'
    r'(?:very[_ ]?strong|strong|moderate|supporting|stand[_ ]?alone'
    r'|VeryStrong|Strong|Moderate|Supporting|Standalone'
    r'|_VS|_S|_M|_P|_SA'
    r')'
    r')?'
)

# Full pattern: word-boundary code + optional strength modifier
EXPLICIT_PATTERN = re.compile(
    r'\b(' + _BASE_CODE + r')' + _STRENGTH_SUFFIX,
    re.IGNORECASE,
)

# Pattern for parenthetical strength: "PM2 (moderate)"
PAREN_STRENGTH_PATTERN = re.compile(
    r'\b(' + _BASE_CODE + r')\s*\('
    r'(very\s*strong|strong|moderate|supporting|stand[- ]?alone)'
    r'\)',
    re.IGNORECASE,
)


def extract_explicit_codes(text):
    """Tier 1: Extract explicit ACMG codes from text.

    Returns a set of normalized base codes (e.g., 'PM2', 'PVS1').
    """
    if not text:
        return set()
    codes = set()
    # Standard pattern
    for m in EXPLICIT_PATTERN.finditer(text):
        code = m.group(1).upper()
        # Normalize to canonical form
        code = _normalize_code(code)
        if code:
            codes.add(code)
    # Parenthetical pattern
    for m in PAREN_STRENGTH_PATTERN.finditer(text):
        code = m.group(1).upper()
        code = _normalize_code(code)
        if code:
            codes.add(code)
    return codes


def _normalize_code(raw_code):
    """Normalize an extracted code to canonical form."""
    raw_code = raw_code.upper().strip()
    # Match against known codes
    for canonical in ALL_ACMG_CODES:
        if raw_code == canonical:
            return canonical
    # Check if it's a valid pattern but not in our predefined list
    # (future ACMG codes like PS5, PM7 etc.)
    m = re.match(r'^(PVS|PS|PM|PP|BA|BS|BP)(\d+)$', raw_code)
    if m:
        return raw_code
    return None


# ═══════════════════════════════════════════════════════════════════════
# TIER 2: NLP-INFERRED CODE EXTRACTION
# ═══════════════════════════════════════════════════════════════════════

# Negation detection window (characters before match to check for negation)
NEGATION_WINDOW = 40
NEGATION_PATTERN = re.compile(
    r'\b(no|not|without|absence of|lack of|lacking|negative for|ruled out|'
    r'unlikely|excluded|does not|did not|cannot|insufficient|failed to)\b',
    re.IGNORECASE,
)

# NLP patterns: list of (compiled_regex, inferred_code, description)
NLP_PATTERNS = []


def _add_nlp(code, patterns, desc=""):
    """Register NLP patterns for a given ACMG code."""
    for pat_str in patterns:
        compiled = re.compile(pat_str, re.IGNORECASE)
        NLP_PATTERNS.append((compiled, code, desc or pat_str))


# PVS1 — Null variant in a gene where LOF is a known mechanism
_add_nlp("PVS1", [
    r'\bnull\s+variant\b',
    r'\bloss[\s-]*of[\s-]*function\b',
    r'\btruncating\s+(?:variant|mutation)\b',
    r'\bnonsense[\s-]*mediated\s+decay\b',
    r'\bNMD\b',
    r'\bframeshift\s+(?:variant|mutation|truncat)',
    r'\bpremature\s+(?:stop|termination|truncation)\b',
    r'\binitiation\s+codon\b',
    r'\bcanonical\s+splice\s+site\b',
], "Null variant / LOF")

# PS1 — Same amino acid change as established pathogenic
_add_nlp("PS1", [
    r'\bsame\s+amino\s+acid\s+change\b',
    r'\bpreviously\s+established\s+pathogenic\b',
    r'\bknown\s+pathogenic\s+(?:variant|mutation)\s+at\s+(?:the\s+)?same\b',
    r'\bidentical\s+(?:amino\s+acid|missense)\s+(?:change|substitution)\b',
], "Same AA change as pathogenic")

# PS2 — De novo (confirmed)
_add_nlp("PS2", [
    r'\bconfirmed?\s+de\s*novo\b',
    r'\bde\s*novo\s*\(?\s*(?:maternity|paternity)\s+(?:and\s+(?:maternity|paternity)\s+)?confirmed\b',
    r'\bde\s*novo\b(?=.*\b(?:paternity|maternity|parentage)\s+confirmed)',
], "Confirmed de novo")

# PS3 — Functional studies supportive of damaging effect
_add_nlp("PS3", [
    r'\bfunctional\s+stud(?:y|ies)\b',
    r'\bfunctional\s+assay\b',
    r'\bin\s+vitro\b',
    r'\bin\s+vivo\b',
    r'\bwell[\s-]*established\s+functional\b',
    r'\bfunctional\s+(?:analysis|characterization|evidence|data|impact|effect)\b',
    r'\bfunctional\s+protein\s+(?:analysis|assay)\b',
    r'\benzymatic\s+(?:activity|assay)\b',
    r'\btranscriptional?\s+(?:activity|assay)\b',
    r'\bknockout\s+(?:mouse|mice|model)\b',
    r'\bsplice\s+assay\b',
    r'\bminigene\b',
    r'\bRNA\s+(?:analysis|studies|assay)\b',
], "Functional studies")

# PS4 — Prevalence in affected significantly increased
_add_nlp("PS4", [
    r'\bprevalence\s+(?:in\s+affected\s+)?significantly\s+increased\b',
    r'\bodds\s+ratio\b',
    r'\bcase[\s-]*control\b',
    r'\bstatistically\s+significant\s+(?:increase|association|enrichment)\b',
    r'\bsignificantly\s+(?:more\s+frequent|enriched)\s+in\s+(?:affected|cases|patients)\b',
], "Prevalence / case-control")

# PM1 — Mutational hot spot / critical functional domain
_add_nlp("PM1", [
    r'\bmutational?\s+hot\s*spot\b',
    r'\bcritical\s+(?:domain|region|residue)\b',
    r'\bfunctional\s+domain\b',
    r'\bwell[\s-]*established\s+(?:functional\s+)?domain\b',
    r'\bactive\s+site\b',
    r'\bbinding\s+(?:site|domain|pocket)\b',
    r'\bcatalytic\s+(?:site|domain|residue)\b',
], "Hot spot / critical domain")

# PM2 — Absent from controls / population databases
_add_nlp("PM2", [
    r'\babsent\s+from\s+(?:controls|gnomAD|ExAC|ESP|databases?|population)\b',
    r'\bnot\s+(?:found|present|observed|seen)\s+in\s+(?:gnomAD|ExAC|ESP|controls|population)\b',
    r'\bnot\s+in\s+(?:gnomAD|ExAC)\b',
    r'\blow\s+(?:allele\s+)?frequency\b(?=.*(?:gnomAD|control|population|ExAC))',
    r'\b(?:allele\s+)?frequency\s+(?:of\s+)?0(?:\.0+)?\s+in\b',
    r'\bzero\s+(?:allele\s+)?frequency\b',
    r'\brare\s+(?:in|across)\s+(?:all\s+)?population\b',
    r'\babsent\s+(?:in|from)\s+(?:all\s+)?population\s+databases?\b',
], "Absent from controls")

# PM3 — In trans with pathogenic variant (recessive)
_add_nlp("PM3", [
    r'\bin\s+trans\s+with\s+(?:a\s+)?pathogenic\b',
    r'\bcompound\s+heterozygous\b',
    r'\bcompound\s+het\b',
    r'\btrans\s+(?:configuration|phase)\b',
    r'\bbiallelic\b',
], "In trans / compound het")

# PM4 — Protein length change (in-frame)
_add_nlp("PM4", [
    r'\bin[\s-]*frame\s+(?:deletion|insertion|indel)\b',
    r'\bprotein\s+length\s+chang(?:e|ing)\b',
    r'\bstop[\s-]*loss\b',
], "In-frame / protein length change")

# PM5 — Novel missense at known pathogenic position
_add_nlp("PM5", [
    r'\bnovel\s+missense\s+(?:at|in)\s+(?:a\s+)?(?:known\s+)?pathogenic\s+(?:position|residue)\b',
    r'\bdifferent\s+(?:amino\s+acid|missense)\s+(?:change|substitution)\s+at\s+(?:the\s+)?same\s+(?:position|residue|codon)\b',
], "Novel missense at pathogenic position")

# PM6 — Assumed de novo
_add_nlp("PM6", [
    r'\bassumed\s+de\s*novo\b',
    r'\bde\s*novo\b(?=.*\b(?:not\s+confirmed|without\s+(?:parental|parentage)|unconfirmed))',
], "Assumed de novo")

# PP1 — Co-segregation with disease
_add_nlp("PP1", [
    r'\bco[\s-]*segregat(?:ion|es|ing)\b',
    r'\bsegregat(?:es|ion|ing)\s+with\s+(?:the\s+)?disease\b',
    r'\bsegregat(?:es|ion|ing)\s+(?:with|in)\s+(?:affected|family|families)\b',
    r'\bfamilial\s+segregation\b',
], "Co-segregation")

# PP2 — Missense in gene with low rate of benign missense
_add_nlp("PP2", [
    r'\bmissense\s+(?:variant|mutation)\s+in\s+(?:a\s+)?gene\s+(?:with|where)\s+(?:a\s+)?low\s+rate\b',
    r'\bgenetic?\s+constraint\b',
], "Missense in constrained gene")

# PP3 — Computational / in silico evidence supports deleterious
_add_nlp("PP3", [
    r'\bcomputational\s+evidence\b',
    r'\bin\s+silico\b',
    r'\bREVEL\b',
    r'\bCADD\b',
    r'\bpredicted\s+(?:deleterious|damaging|pathogenic)\b',
    r'\bPolyphen[\s-]*2?\b',
    r'\bSIFT\b(?=.*(?:damaging|deleterious))',
    r'\bMutationTaster\b',
    r'\bmultiple\s+(?:lines\s+of\s+)?computational\s+evidence\b',
    r'\bbioinformatic\b',
    r'\balgorithm(?:s|ic)?\s+(?:predict|support|agree)\b',
], "Computational / in silico")

# PP4 — Patient phenotype highly specific for gene
_add_nlp("PP4", [
    r'\bpatient\s+(?:\'s\s+)?phenotype\s+(?:is\s+)?highly\s+specific\b',
    r'\bphenotype\s+(?:is\s+)?(?:highly\s+)?specific\s+for\b',
    r'\bspecific\s+phenotype\b(?=.*(?:single\s+gene|genetic\s+etiology))',
], "Phenotype specificity")

# PP5 — Reputable source / expert panel
_add_nlp("PP5", [
    r'\breputable\s+source\b',
    r'\bClinGen\b',
    r'\bexpert\s+panel\b(?=.*(?:classif|report|assert|review))',
], "Reputable source")

# BA1 — Allele frequency > 5%
_add_nlp("BA1", [
    r'\ballele\s+frequency\s*(?:>|greater\s+than)\s*(?:0\.0?5|5\s*%)',
    r'\btoo\s+common\b',
    r'\bMAF\s*(?:>|greater\s+than)\s*(?:0\.0?5|5\s*%)',
    r'\bpopulation\s+frequency\s*(?:>|greater\s+than)\s*(?:0\.0?5|5\s*%)',
    r'\bcommon\s+(?:polymorphism|variant)\b(?=.*(?:frequency|MAF|AF|gnomAD))',
], "Allele frequency > 5%")

# BS1 — Allele frequency greater than expected for disorder
_add_nlp("BS1", [
    r'\ballele\s+frequency\s+(?:is\s+)?greater\s+than\s+expected\b',
    r'\bfrequency\s+(?:too\s+high|higher\s+than\s+expected)\b',
    r'\b(?:allele\s+)?frequency\s+(?:is\s+)?inconsistent\s+with\b',
    r'\bMAF\s+(?:is\s+)?(?:too\s+high|greater\s+than\s+expected)\b',
], "AF greater than expected")

# BS2 — Observed in healthy adult (dominant) or healthy homozygous (recessive)
_add_nlp("BS2", [
    r'\bobserved\s+in\s+(?:a\s+)?healthy\s+(?:adult|individual|control)\b',
    r'\bhomozygous\s+(?:in\s+)?(?:a\s+)?healthy\b',
    r'\bunaffected\s+(?:individual|adult|carrier)s?\s+(?:with|carrying)\b',
], "Observed in healthy individual")

# BS3 — Functional studies show no damaging effect
_add_nlp("BS3", [
    r'\bfunctional\s+stud(?:y|ies)\s+show\s+no\s+damaging\b',
    r'\bnormal\s+protein\s+function\b',
    r'\bno\s+(?:functional\s+)?(?:impact|effect)\s+(?:on|in)\s+(?:protein\s+)?function\b',
    r'\bfunctional(?:ly)?\s+(?:normal|neutral|benign)\b',
    r'\bno\s+(?:significant\s+)?(?:effect|impact)\s+on\s+(?:splicing|expression|activity)\b',
], "Functional studies benign")

# BS4 — Lack of segregation in affected family members
_add_nlp("BS4", [
    r'\black\s+of\s+segregation\b',
    r'\bdoes\s+not\s+segregate\b',
    r'\bnon[\s-]*segregation\b',
    r'\bfails?\s+to\s+segregate\b',
    r'\bsegregation\s+(?:analysis\s+)?(?:does\s+)?not\s+support\b',
], "Lack of segregation")

# BP1 — Missense in gene where only truncating cause disease
_add_nlp("BP1", [
    r'\bmissense\s+(?:variant|mutation|change)\s+in\s+(?:a\s+)?gene\s+(?:where|for\s+which)\s+(?:only\s+)?truncating\b',
], "Missense where only truncating")

# BP3 — In-frame in repetitive region without known function
_add_nlp("BP3", [
    r'\bin[\s-]*frame\s+(?:deletion|insertion)\s+in\s+(?:a\s+)?repetitive\s+region\b',
], "In-frame in repetitive region")

# BP4 — Computational evidence suggests benign
_add_nlp("BP4", [
    r'\bcomputational\s+evidence\s+suggest(?:s|ing)\s+(?:no\s+impact|benign|tolerated|neutral)\b',
    r'\bpredicted\s+(?:benign|tolerated|neutral)\b',
    r'\bin\s+silico\b(?=.*(?:benign|tolerated|neutral|no[\s-]impact))',
    r'\btolerated\s+by\s+(?:SIFT|PolypPhen|prediction)\b',
], "Computational evidence benign")

# BP5 — Found in case with alternate molecular basis
_add_nlp("BP5", [
    r'\balternate\s+molecular\s+basis\b',
    r'\balternative\s+(?:genetic\s+)?(?:cause|explanation|diagnosis)\b',
], "Alternate molecular basis")

# BP6 — Reputable source reports benign
_add_nlp("BP6", [
    r'\breputable\s+source\s+(?:recently\s+)?report(?:s|ed)\s+(?:as\s+)?benign\b',
], "Reputable source benign")

# BP7 — Synonymous with no predicted splice impact
_add_nlp("BP7", [
    r'\bsynonymous\s+(?:variant\s+)?(?:with|and)\s+no\s+(?:predicted\s+)?splic(?:e|ing)\s+impact\b',
    r'\bsilent\s+(?:change|variant|mutation)\b',
    r'\bsynonymous\b(?=.*(?:no\s+(?:effect|impact)\s+on\s+splic|does\s+not\s+affect\s+splic))',
    r'\bsynonymous\b(?=.*\bnot\s+(?:near|at)\s+(?:a\s+)?(?:splice|exon[\s-]intron)\b)',
], "Synonymous / no splice impact")


def _has_negation(text, match_start):
    """Check if a negation word appears within NEGATION_WINDOW chars before match."""
    window_start = max(0, match_start - NEGATION_WINDOW)
    window = text[window_start:match_start]
    return bool(NEGATION_PATTERN.search(window))


def extract_inferred_codes(text):
    """Tier 2: Infer ACMG codes from natural language patterns.

    Returns a set of inferred codes (e.g., 'PM2', 'PS3').
    Applies negation detection to avoid false positives.
    """
    if not text:
        return set()
    codes = set()
    for pattern, code, _desc in NLP_PATTERNS:
        for m in pattern.finditer(text):
            if not _has_negation(text, m.start()):
                codes.add(code)
                break  # One match per pattern is enough
    return codes


# ═══════════════════════════════════════════════════════════════════════
# CHECKPOINT LOADING
# ═══════════════════════════════════════════════════════════════════════

def find_latest_checkpoint(checkpoint_dir):
    """Find the most recent VCV checkpoint file."""
    ckpt_dir = Path(checkpoint_dir)
    if not ckpt_dir.exists():
        return None

    # Look for checkpoint files
    candidates = []
    for f in ckpt_dir.iterdir():
        if f.name.endswith('.checkpoint.pkl.gz'):
            candidates.append(f)
        elif f.name.endswith('.checkpoint.pkl'):
            candidates.append(f)

    if not candidates:
        return None

    # Prefer VCV checkpoints, sort by name (date-based naming)
    vcv = [c for c in candidates if 'VCV' in c.name or 'vcv' in c.name]
    pool = vcv if vcv else candidates
    pool.sort(key=lambda p: p.name)
    return pool[-1]  # Latest by name


def load_checkpoint(path):
    """Load a checkpoint file (gzipped or plain pickle)."""
    path = Path(path)
    print(f"Loading checkpoint: {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
    t0 = time.time()

    if path.name.endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            data = pickle.load(f)
    else:
        with open(path, 'rb') as f:
            data = pickle.load(f)

    elapsed = time.time() - t0
    print(f"  Loaded in {elapsed:.1f}s")

    # Handle different checkpoint formats
    if isinstance(data, dict):
        if 'data' in data:
            snapshot_date = data.get('snapshot_date', 'unknown')
            variants = data['data']
            print(f"  Snapshot date: {snapshot_date}")
            print(f"  Variants: {len(variants):,}")
            return snapshot_date, variants
        else:
            # The dict itself is the variant data
            print(f"  Variants: {len(data):,}")
            return 'unknown', data
    else:
        print(f"  WARNING: Unexpected checkpoint format: {type(data)}")
        return 'unknown', {}


# ═══════════════════════════════════════════════════════════════════════
# CLASSIFICATION CATEGORIZATION
# ═══════════════════════════════════════════════════════════════════════

def categorize_classification(cls_str):
    """Categorize a classification string."""
    if not cls_str:
        return "Other"
    c = cls_str.lower().strip()
    if "conflicting" in c:
        return "Conflicting"
    if "uncertain" in c or "vus" in c:
        return "VUS"
    if "pathogenic/likely" in c or "pathogenic, likely" in c:
        return "P/LP"
    if "benign/likely" in c or "benign, likely" in c:
        return "B/LB"
    if "likely pathogenic" in c:
        return "LP"
    if "pathogenic" in c:
        return "P"
    if "likely benign" in c:
        return "LB"
    if "benign" in c:
        return "B"
    return "Other"


# ═══════════════════════════════════════════════════════════════════════
# MAIN EXTRACTION PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def process_variants(variants, output_path):
    """Process all variants and extract ACMG criteria."""
    print(f"\nProcessing {len(variants):,} variants ...")
    t0 = time.time()

    # Statistics counters
    total_variants = 0
    total_submissions = 0
    submissions_with_comments = 0
    comments_with_explicit = 0
    comments_with_inferred = 0
    comments_with_any = 0
    submissions_with_any_codes = 0

    explicit_code_freq = Counter()
    inferred_code_freq = Counter()
    all_code_freq = Counter()
    code_by_category = defaultdict(Counter)  # category -> code -> count

    rows = []

    for variant_id, variant in variants.items():
        total_variants += 1
        if total_variants % 100_000 == 0:
            elapsed = time.time() - t0
            rate = total_variants / elapsed
            print(f"  {total_variants:,} variants processed "
                  f"({rate:.0f}/s, {elapsed:.0f}s elapsed)")

        # Extract variant-level info
        vcv_id = variant.get('vcv_id', '') or ''
        gene = variant.get('gene_symbol', '') or ''
        classification = variant.get('classification', '') or ''
        category = categorize_classification(classification)

        submissions = variant.get('submissions', []) or []
        num_subs = len(submissions)
        num_comments = 0

        variant_explicit = set()
        variant_inferred = set()

        for sub in submissions:
            total_submissions += 1

            # Collect all text fields to search
            comment = sub.get('comment', '') or ''
            sub_classification = sub.get('classification', '') or ''
            existing_acmg = sub.get('acmg_codes', []) or []

            # Build combined text for NLP search
            text_fields = []
            if comment:
                text_fields.append(comment)
                num_comments += 1
                submissions_with_comments += 1
            if sub_classification:
                text_fields.append(sub_classification)

            combined_text = ' '.join(text_fields)

            # Tier 1: Explicit codes from all text
            sub_explicit = extract_explicit_codes(combined_text)

            # Also re-normalize existing acmg_codes from v5
            for code_str in existing_acmg:
                normalized = extract_explicit_codes(code_str)
                sub_explicit.update(normalized)

            # Tier 2: NLP-inferred codes from combined text
            sub_inferred = extract_inferred_codes(combined_text)
            # Remove codes already found explicitly
            sub_inferred -= sub_explicit

            # Track per-submission stats
            has_explicit = bool(sub_explicit)
            has_inferred = bool(sub_inferred)
            has_any = has_explicit or has_inferred

            if comment:
                if has_explicit:
                    comments_with_explicit += 1
                if has_inferred:
                    comments_with_inferred += 1
                if has_any:
                    comments_with_any += 1

            if has_any:
                submissions_with_any_codes += 1

            variant_explicit.update(sub_explicit)
            variant_inferred.update(sub_inferred)

            # Update frequency counters
            for code in sub_explicit:
                explicit_code_freq[code] += 1
                all_code_freq[code] += 1
            for code in sub_inferred:
                inferred_code_freq[code] += 1
                all_code_freq[code] += 1

        # Remove inferred that overlap with explicit at variant level
        variant_inferred -= variant_explicit
        all_codes = variant_explicit | variant_inferred

        # Track per-category code usage
        for code in all_codes:
            code_by_category[category][code] += 1

        # Build row
        rows.append({
            'Variation_ID': variant_id,
            'VCV_ID': vcv_id,
            'Gene': gene,
            'Classification': classification,
            'Category': category,
            'Num_Submissions': num_subs,
            'Num_Comments': num_comments,
            'Explicit_Codes': ';'.join(sorted(variant_explicit)) if variant_explicit else '',
            'Inferred_Codes': ';'.join(sorted(variant_inferred)) if variant_inferred else '',
            'All_Codes': ';'.join(sorted(all_codes)) if all_codes else '',
            'Num_Explicit': len(variant_explicit),
            'Num_Inferred': len(variant_inferred),
            'Num_All_Codes': len(all_codes),
        })

    elapsed = time.time() - t0
    print(f"  Done: {total_variants:,} variants in {elapsed:.0f}s")

    # ── Write CSV ────────────────────────────────────────────────────
    print(f"\nWriting {output_path} ...")
    fieldnames = [
        'Variation_ID', 'VCV_ID', 'Gene', 'Classification', 'Category',
        'Num_Submissions', 'Num_Comments',
        'Explicit_Codes', 'Inferred_Codes', 'All_Codes',
        'Num_Explicit', 'Num_Inferred', 'Num_All_Codes',
    ]
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written {len(rows):,} rows")

    # ── Build statistics dict for visualization ──────────────────────
    stats = {
        'total_variants': total_variants,
        'total_submissions': total_submissions,
        'submissions_with_comments': submissions_with_comments,
        'comments_with_explicit': comments_with_explicit,
        'comments_with_inferred': comments_with_inferred,
        'comments_with_any': comments_with_any,
        'submissions_with_any_codes': submissions_with_any_codes,
        'explicit_code_freq': dict(explicit_code_freq),
        'inferred_code_freq': dict(inferred_code_freq),
        'all_code_freq': dict(all_code_freq),
        'code_by_category': {cat: dict(counts) for cat, counts in code_by_category.items()},
    }

    # Write stats JSON alongside CSV
    import json
    stats_path = Path(output_path).with_suffix('.json')
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"  Statistics saved to {stats_path.name}")

    # ── Print console statistics ─────────────────────────────────────
    _print_statistics(stats)

    return stats


def _print_statistics(stats):
    """Print comprehensive extraction statistics to console."""
    print()
    print("=" * 70)
    print("ACMG/AMP CRITERIA EXTRACTION STATISTICS")
    print("=" * 70)

    tv = stats['total_variants']
    ts = stats['total_submissions']
    sc = stats['submissions_with_comments']
    ce = stats['comments_with_explicit']
    ci = stats['comments_with_inferred']
    ca = stats['comments_with_any']
    sa = stats['submissions_with_any_codes']

    print(f"\n{'Overview':}")
    print(f"  Total variants:              {tv:>12,}")
    print(f"  Total submissions:           {ts:>12,}")
    print(f"  Submissions with comments:   {sc:>12,}  ({sc/max(ts,1)*100:5.1f}%)")
    print(f"  Submissions with any codes:  {sa:>12,}  ({sa/max(ts,1)*100:5.1f}%)")

    print(f"\n{'Comment-level extraction rates':}")
    print(f"  Comments with explicit codes: {ce:>11,}  ({ce/max(sc,1)*100:5.1f}% of comments)")
    print(f"  Comments with inferred codes: {ci:>11,}  ({ci/max(sc,1)*100:5.1f}% of comments)")
    print(f"  Comments with ANY codes:      {ca:>11,}  ({ca/max(sc,1)*100:5.1f}% of comments)")
    print(f"  Uplift from NLP inference:    {ci-0:>11,}  "
          f"(+{(ca-ce)/max(ce,1)*100:.1f}% over explicit-only)")

    # Per-criteria frequency table
    print(f"\n{'Per-criteria frequency (top 30)':}")
    print(f"  {'Code':<8} {'Explicit':>10} {'Inferred':>10} {'Total':>10}  {'Source'}")
    print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10}  {'-'*20}")

    ef = stats['explicit_code_freq']
    inf = stats['inferred_code_freq']
    af = stats['all_code_freq']

    top_codes = sorted(af.keys(), key=lambda c: af[c], reverse=True)[:30]
    for code in top_codes:
        e = ef.get(code, 0)
        i = inf.get(code, 0)
        t = af.get(code, 0)
        src = "explicit" if e > 0 and i == 0 else "inferred" if e == 0 and i > 0 else "both"
        print(f"  {code:<8} {e:>10,} {i:>10,} {t:>10,}  {src}")

    # Strength class summary
    print(f"\n{'Strength class summary':}")
    class_totals = defaultdict(int)
    for code, count in af.items():
        cls = STRENGTH_CLASS.get(code, "Other")
        class_totals[cls] += count

    for cls in ["Very Strong Pathogenic", "Strong Pathogenic", "Moderate Pathogenic",
                "Supporting Pathogenic", "Stand-alone Benign", "Strong Benign",
                "Supporting Benign"]:
        ct = class_totals.get(cls, 0)
        print(f"  {cls:<28} {ct:>10,}")

    print()
    print("=" * 70)


# ═══════════════════════════════════════════════════════════════════════
# CLI ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="ClinVar ACMG/AMP Criteria Extractor — "
                    "Extract criteria from checkpoint files"
    )
    parser.add_argument("--checkpoint-dir", "-d",
                        default="./clinvar_output/checkpoints",
                        help="Directory containing checkpoint files")
    parser.add_argument("--checkpoint-file", "-f",
                        help="Specific checkpoint file to use (overrides --checkpoint-dir)")
    parser.add_argument("--output", "-o",
                        default="acmg_criteria_detailed.csv",
                        help="Output CSV path")
    args = parser.parse_args()

    print("=" * 70)
    print("ClinVar ACMG/AMP Assertion Criteria Extractor")
    print("=" * 70)
    print(f"  Tier 1: Enhanced explicit regex extraction")
    print(f"  Tier 2: NLP-inferred extraction (~{len(NLP_PATTERNS)} patterns)")
    print()

    # Find and load checkpoint
    if args.checkpoint_file:
        ckpt_path = Path(args.checkpoint_file)
        if not ckpt_path.exists():
            print(f"ERROR: Checkpoint file not found: {ckpt_path}")
            sys.exit(1)
    else:
        ckpt_path = find_latest_checkpoint(args.checkpoint_dir)
        if not ckpt_path:
            print(f"ERROR: No checkpoint files found in {args.checkpoint_dir}")
            print(f"  Expected *.checkpoint.pkl.gz files")
            sys.exit(1)

    snapshot_date, variants = load_checkpoint(ckpt_path)

    if not variants:
        print("ERROR: No variant data found in checkpoint")
        sys.exit(1)

    # Process
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    process_variants(variants, output_path)

    print(f"\nOutput files:")
    print(f"  CSV:   {output_path}")
    print(f"  Stats: {output_path.with_suffix('.json')}")
    print("Done.")


if __name__ == "__main__":
    main()
