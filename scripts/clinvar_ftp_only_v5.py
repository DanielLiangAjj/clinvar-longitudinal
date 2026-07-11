#!/usr/bin/env python3
"""
ClinVar FTP-Only Analysis Tool v5.1
====================================

FIXES from v5.0:
- Fixed XML iterparse memory leak (root element children never cleared)
- Fixed Phase 3: each checkpoint now loaded only ONCE instead of once-per-batch
- Streaming CSV output instead of accumulating all VariantFTPRecord objects in RAM
- Graceful handling of missing checkpoints (warning instead of error + continue)
- Fixed oscillation detection bug (was checking deduplicated list, never triggered)
- String interning to reduce memory usage of repeated classification/review strings
- Explicit gc.collect() after freeing large checkpoint data

Author: Generated for ClinVar Research Project
Date: 2026-02
"""

import os
import sys
import json
import gzip
import gc
import pickle
import logging
import argparse
import time
import re
import csv
import traceback
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Set, Iterator
from dataclasses import dataclass, field, asdict
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import xml.etree.ElementTree as ET


# ============================================================================
# LOGGING - File + Console with proper formatting
# ============================================================================

def setup_logging(output_dir: Path, name: str = "clinvar") -> logging.Logger:
    """Setup logging to both file and console."""
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{name}_{timestamp}.log"

    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Clear existing handlers
    logger.handlers = []

    # File handler - detailed logging
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(fh)

    # Console handler - info and above
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(ch)

    logger.info(f"Logging to: {log_file}")
    return logger


# ============================================================================
# DATA STRUCTURES (kept for backward compatibility with existing pickle files)
# ============================================================================

@dataclass
class SubmissionRecord:
    """Individual submission (SCV) from a snapshot."""
    scv_id: str = ""
    submitter: str = ""
    submitter_id: str = ""
    classification: str = ""
    review_status: str = ""
    last_evaluated: str = ""
    date_submitted: str = ""
    condition: str = ""
    condition_id: str = ""
    origin: str = ""
    affected_status: str = ""
    method: str = ""
    pmids: List[str] = field(default_factory=list)
    acmg_codes: List[str] = field(default_factory=list)
    comment: str = ""


@dataclass
class SnapshotClassification:
    """Classification data from a single snapshot for a variant."""
    snapshot_date: str = ""
    classification: str = ""
    review_status: str = ""
    num_submissions: int = 0
    submitters: List[str] = field(default_factory=list)


@dataclass
class VariantFTPRecord:
    """Complete variant record parsed from FTP files."""
    variation_id: str = ""
    vcv_id: str = ""
    allele_id: str = ""
    dbsnp_id: str = ""
    variant_type: str = ""
    gene_symbol: str = ""
    gene_id: str = ""
    hgvs_c: str = ""
    hgvs_p: str = ""
    hgvs_g: str = ""
    molecular_consequence: str = ""
    chromosome: str = ""
    start: int = 0
    stop: int = 0
    ref_allele: str = ""
    alt_allele: str = ""
    assembly: str = ""
    current_classification: str = ""
    current_review_status: str = ""
    classification_history: List[SnapshotClassification] = field(default_factory=list)
    first_classification: str = ""
    first_classification_date: str = ""
    classification_changed: bool = False
    change_type: str = ""
    num_changes: int = 0
    had_oscillation: bool = False
    submissions: List[SubmissionRecord] = field(default_factory=list)
    all_pmids: List[str] = field(default_factory=list)
    all_acmg_codes: List[str] = field(default_factory=list)
    all_submitters: List[str] = field(default_factory=list)
    total_submissions: int = 0
    date_created: str = ""
    date_last_updated: str = ""


# ============================================================================
# SNAPSHOT PARSING FUNCTIONS
# ============================================================================

def parse_vcv_snapshot_full(snapshot_path: Path, logger=None) -> Dict[str, Dict]:
    """Parse VCV format snapshot with FULL extraction.

    FIX v5.1: Uses root.clear() to prevent iterparse memory leak.
    Raises on mid-parse errors so partial results are NOT saved as checkpoints.
    """
    results = {}
    count = 0
    log = logger.info if logger else print

    log(f"Parsing VCV snapshot: {snapshot_path.name}")
    start_time = time.time()

    with gzip.open(snapshot_path, 'rb') as f:
        context = ET.iterparse(f, events=('start', 'end'))
        root = None

        for event, elem in context:
            if event == 'start':
                if root is None:
                    root = elem
                continue

            # event == 'end'
            if elem.tag == 'VariationArchive':
                record = _extract_vcv_record(elem)
                if record and record.get('variation_id'):
                    results[record['variation_id']] = record
                    count += 1

                    if count % 100000 == 0:
                        elapsed = time.time() - start_time
                        rate = count / elapsed if elapsed > 0 else 0
                        log(f"  [{snapshot_path.name}] Parsed {count:,} variants ({rate:.0f}/s)")

                # Clear root's children to free memory (standard iterparse idiom)
                if root is not None:
                    root.clear()

    elapsed = time.time() - start_time
    log(f"  [{snapshot_path.name}] Completed: {count:,} variants in {elapsed:.1f}s ({count/max(elapsed,1):.0f}/s)")

    return results


def _safe_int(val, default=0):
    """Safely convert a value to int without raising."""
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


def _extract_vcv_record(elem) -> Optional[Dict]:
    """Extract all data from a VariationArchive element.

    IMPORTANT: Only returns None if variation_id is genuinely missing.
    Individual field extraction failures are caught per-field so that
    a single bad field doesn't drop the entire variant.
    """
    vcv_id = elem.get('Accession', '')
    variation_id = elem.get('VariationID') or ''

    if not variation_id and vcv_id:
        match = re.search(r'VCV0*(\d+)', vcv_id)
        if match:
            variation_id = match.group(1)

    if not variation_id:
        return None

    record = {
        'variation_id': variation_id,
        'vcv_id': vcv_id,
        'date_created': elem.get('DateCreated', ''),
        'date_updated': elem.get('DateLastUpdated', ''),
        'version': elem.get('Version', '1'),
    }

    # Classification — critical field, protected separately
    try:
        classification = None
        review_status = None

        for cls_path in ['.//Classifications/GermlineClassification',
                         './/Classifications/SomaticClinicalImpact',
                         './/Classifications/OncogenicityClassification',
                         './/GermlineClassification',
                         './/ClinicalSignificance',
                         './/Interpretations/Interpretation']:
            cls_elem = elem.find(cls_path)
            if cls_elem is not None:
                desc = cls_elem.find('Description')
                if desc is not None and desc.text:
                    classification = desc.text.strip()
                review = cls_elem.find('ReviewStatus')
                if review is not None and review.text:
                    review_status = review.text.strip()
                if classification:
                    break

        # Fallback: search for ReviewStatus at broader paths (fixes VCV 2020-2024 format
        # where ReviewStatus is at InterpretedRecord level, not inside Interpretation)
        if not review_status:
            for rs_path in ['.//InterpretedRecord/ReviewStatus',
                            './/ClassifiedRecord/ReviewStatus',
                            './/ReviewStatus']:
                rs_elem = elem.find(rs_path)
                if rs_elem is not None and rs_elem.text:
                    review_status = rs_elem.text.strip()
                    break

        record['classification'] = classification or ''
        record['review_status'] = review_status or ''
    except Exception:
        record['classification'] = ''
        record['review_status'] = ''

    # Variant info — non-critical metadata, failures should not drop the variant
    try:
        simple_allele = elem.find('.//SimpleAllele')
        if simple_allele is not None:
            record['allele_id'] = simple_allele.get('AlleleID', '')
            record['variant_type'] = simple_allele.get('Type', '')

            gene_elem = simple_allele.find('.//Gene')
            if gene_elem is not None:
                record['gene_symbol'] = gene_elem.get('Symbol', '')
                record['gene_id'] = gene_elem.get('GeneID', '')

            for hgvs in simple_allele.findall('.//HGVSlist/HGVS'):
                hgvs_type = hgvs.get('Type', '')
                # Try both old schema (NucleotideChange) and new schema (NucleotideExpression/Expression)
                nuc_text = None
                nuc_change = hgvs.find('NucleotideChange')
                if nuc_change is not None and nuc_change.text:
                    nuc_text = nuc_change.text
                else:
                    nuc_expr = hgvs.find('NucleotideExpression')
                    if nuc_expr is not None:
                        expr = nuc_expr.find('Expression')
                        if expr is not None and expr.text:
                            nuc_text = expr.text
                        elif nuc_expr.text:
                            nuc_text = nuc_expr.text

                if nuc_text:
                    if ('coding' in hgvs_type.lower() or ':c.' in nuc_text) and not record.get('hgvs_c'):
                        record['hgvs_c'] = nuc_text
                    elif ('genomic' in hgvs_type.lower() or ':g.' in nuc_text) and not record.get('hgvs_g'):
                        record['hgvs_g'] = nuc_text

                prot_text = None
                prot_change = hgvs.find('ProteinChange')
                if prot_change is not None and prot_change.text:
                    prot_text = prot_change.text
                else:
                    prot_expr = hgvs.find('ProteinExpression')
                    if prot_expr is not None:
                        expr = prot_expr.find('Expression')
                        if expr is not None and expr.text:
                            prot_text = expr.text
                        elif prot_expr.text:
                            prot_text = prot_expr.text

                if prot_text and not record.get('hgvs_p'):
                    record['hgvs_p'] = prot_text

            mol_cons = simple_allele.find('.//MolecularConsequence')
            if mol_cons is not None:
                record['molecular_consequence'] = mol_cons.get('Type', '')

            for loc in simple_allele.findall('.//Location/SequenceLocation'):
                if loc.get('Assembly') == 'GRCh38':
                    record['chromosome'] = loc.get('Chr', '')
                    record['start'] = _safe_int(loc.get('start', 0) or 0)
                    record['stop'] = _safe_int(loc.get('stop', 0) or 0)
                    record['ref_allele'] = loc.get('referenceAllele', '')
                    record['alt_allele'] = loc.get('alternateAllele', '')
                    record['assembly'] = 'GRCh38'
                    break

            for xref in simple_allele.findall('.//XRefList/XRef'):
                if xref.get('DB') == 'dbSNP':
                    record['dbsnp_id'] = 'rs' + xref.get('ID', '')
                    break
    except Exception:
        pass  # Variant still included with whatever metadata was extracted

    # Submissions
    try:
        submissions = []
        for assertion in elem.findall('.//ClinicalAssertion'):
            sub = _extract_submission(assertion)
            if sub:
                submissions.append(sub)

        record['submissions'] = submissions
        record['num_submissions'] = len(submissions)

        all_pmids = set()
        all_acmg = set()
        all_submitters = set()

        for sub in submissions:
            all_pmids.update(sub.get('pmids', []))
            all_acmg.update(sub.get('acmg_codes', []))
            if sub.get('submitter'):
                all_submitters.add(sub['submitter'])

        record['all_pmids'] = list(all_pmids)
        record['all_acmg_codes'] = list(all_acmg)
        record['all_submitters'] = list(all_submitters)
    except Exception:
        record.setdefault('submissions', [])
        record.setdefault('num_submissions', 0)

    return record


def _extract_submission(assertion_elem, logger=None) -> Optional[Dict]:
    """Extract submission data from a ClinicalAssertion element."""
    try:
        sub = {
            'scv_id': assertion_elem.get('ID', ''),
            'date_submitted': assertion_elem.get('DateCreated', ''),
            'date_updated': assertion_elem.get('DateLastUpdated', ''),
        }

        submitter = assertion_elem.find('.//ClinVarSubmissionID')
        if submitter is not None:
            sub['submitter'] = submitter.get('submitter', '')
            sub['submitter_id'] = submitter.get('submitterOrgID', '')

        # Search for per-submission classification under multiple XML paths:
        #   .//Classification     — VCV 2024+ format
        #   .//ClinicalSignificance — RCV (pre-2020) format
        #   .//Interpretation     — VCV 2020-2023 format (Description child)
        for cls_path in ['.//Classification', './/ClinicalSignificance', './/Interpretation']:
            cls_elem = assertion_elem.find(cls_path)
            if cls_elem is not None:
                desc = cls_elem.find('Description')
                if desc is not None and desc.text:
                    sub['classification'] = desc.text.strip()
                review = cls_elem.find('ReviewStatus')
                if review is not None and review.text:
                    sub['review_status'] = review.text.strip()
                date_eval = cls_elem.get('DateLastEvaluated')
                if date_eval:
                    sub['last_evaluated'] = date_eval
                if sub.get('classification'):
                    break

        trait = assertion_elem.find('.//TraitSet/Trait')
        if trait is not None:
            name = trait.find('Name/ElementValue')
            if name is not None and name.text:
                sub['condition'] = name.text.strip()
            for xref in trait.findall('.//XRef'):
                if xref.get('DB') in ['MedGen', 'OMIM', 'Orphanet']:
                    sub['condition_id'] = f"{xref.get('DB')}:{xref.get('ID')}"
                    break

        observed = assertion_elem.find('.//ObservedIn')
        if observed is not None:
            sample = observed.find('Sample')
            if sample is not None:
                origin = sample.find('Origin')
                if origin is not None and origin.text:
                    sub['origin'] = origin.text.strip()
                affected = sample.find('AffectedStatus')
                if affected is not None and affected.text:
                    sub['affected_status'] = affected.text.strip()

        method = assertion_elem.find('.//Method/MethodType')
        if method is not None and method.text:
            sub['method'] = method.text.strip()

        pmids = []
        for citation in assertion_elem.findall('.//Citation'):
            for id_elem in citation.findall('ID'):
                if id_elem.get('Source') == 'PubMed' and id_elem.text:
                    pmids.append(id_elem.text.strip())
        sub['pmids'] = pmids

        acmg_codes = []
        # Check both assertion method and comment fields for ACMG codes
        for attr in assertion_elem.findall('.//AttributeSet/Attribute'):
            if attr.get('Type') in ('AssertionMethod', 'AssertionMethod', 'AssertionMethodDocument'):
                text = attr.text or ''
                codes = re.findall(
                    r'\b(PVS[1-9]|PS[1-9]\d?|PM[1-9]\d?|PP[1-9]\d?|BA[1-9]|BS[1-9]\d?|BP[1-9]\d?)\b',
                    text.upper()
                )
                acmg_codes.extend(codes)
        for comment_elem in assertion_elem.findall('.//Comment'):
            if comment_elem.text:
                codes = re.findall(
                    r'\b(PVS[1-9]|PS[1-9]\d?|PM[1-9]\d?|PP[1-9]\d?|BA[1-9]|BS[1-9]\d?|BP[1-9]\d?)\b',
                    comment_elem.text.upper()
                )
                acmg_codes.extend(codes)
        sub['acmg_codes'] = list(set(acmg_codes))

        comment = assertion_elem.find('.//Comment')
        if comment is not None and comment.text:
            sub['comment'] = comment.text.strip()[:500]

        return sub

    except Exception as e:
        scv = assertion_elem.get('ID', 'unknown')
        if logger:
            logger.warning(f"Dropped RCV submission {scv}: {e}")
        return None


def parse_rcv_snapshot_full(snapshot_path: Path, logger=None) -> Dict[str, Dict]:
    """Parse RCV format snapshot with FULL extraction.

    FIX v5.1: Uses root.clear() to prevent iterparse memory leak.
    Raises on mid-parse errors so partial results are NOT saved as checkpoints.
    """
    variant_records: Dict[str, List[Dict]] = defaultdict(list)
    count = 0
    log = logger.info if logger else print

    log(f"Parsing RCV snapshot: {snapshot_path.name}")
    start_time = time.time()

    with gzip.open(snapshot_path, 'rb') as f:
        context = ET.iterparse(f, events=('start', 'end'))
        root = None

        for event, elem in context:
            if event == 'start':
                if root is None:
                    root = elem
                continue

            # event == 'end'
            if elem.tag == 'ClinVarSet':
                record = _extract_rcv_record(elem)
                if record and record.get('variation_id'):
                    variant_records[record['variation_id']].append(record)
                    count += 1

                    if count % 100000 == 0:
                        elapsed = time.time() - start_time
                        rate = count / elapsed if elapsed > 0 else 0
                        log(f"  [{snapshot_path.name}] Parsed {count:,} RCV records ({rate:.0f}/s)")

                # Clear root's children to free memory
                if root is not None:
                    root.clear()

    # Aggregate by variation ID
    results = {}
    for var_id, records in variant_records.items():
        records.sort(key=lambda x: (_review_weight(x.get('review_status', '')),
                                    x.get('date_updated', '')), reverse=True)
        best = records[0]

        all_submissions = []
        all_pmids = set()
        all_acmg = set()
        all_submitters = set()

        for r in records:
            all_submissions.extend(r.get('submissions', []))
            all_pmids.update(r.get('all_pmids', []))
            all_acmg.update(r.get('all_acmg_codes', []))
            all_submitters.update(r.get('all_submitters', []))

        best['submissions'] = all_submissions
        best['num_submissions'] = len(all_submissions)
        best['all_pmids'] = list(all_pmids)
        best['all_acmg_codes'] = list(all_acmg)
        best['all_submitters'] = list(all_submitters)

        results[var_id] = best

    elapsed = time.time() - start_time
    log(f"  [{snapshot_path.name}] Completed: {count:,} RCV records -> {len(results):,} variants in {elapsed:.1f}s")

    return results


def _extract_rcv_record(elem) -> Optional[Dict]:
    """Extract data from a ClinVarSet element (RCV format).

    Handles both MeasureSet (simple variants) and GenotypeSet (compound het /
    complex genotype variants) so no variant type is silently dropped.
    """
    ref_assertion = elem.find('.//ReferenceClinVarAssertion')
    if ref_assertion is None:
        return None

    # Try direct-child MeasureSet first, then GenotypeSet for compound het.
    # IMPORTANT: must NOT use './/MeasureSet' (descendant search) because
    # GenotypeSet contains nested MeasureSets — that would match the nested
    # one instead of falling through to the GenotypeSet path.
    measure_set = ref_assertion.find('MeasureSet')
    genotype_set = None
    if measure_set is None:
        genotype_set = ref_assertion.find('GenotypeSet')
        if genotype_set is None:
            return None
        variation_id = genotype_set.get('ID')
    else:
        variation_id = measure_set.get('ID')

    if not variation_id:
        return None

    record = {
        'variation_id': variation_id,
        'rcv_id': ref_assertion.get('ID', ''),
        'date_created': ref_assertion.get('DateCreated', ''),
        'date_updated': ref_assertion.get('DateLastUpdated', ''),
    }

    # Classification
    try:
        cls_elem = ref_assertion.find('ClinicalSignificance')
        if cls_elem is not None:
            desc = cls_elem.find('Description')
            if desc is not None and desc.text:
                record['classification'] = desc.text.strip()
            review = cls_elem.find('ReviewStatus')
            if review is not None and review.text:
                record['review_status'] = review.text.strip()
    except Exception:
        pass

    # Collect all Measure elements (from MeasureSet or nested in GenotypeSet)
    try:
        if measure_set is not None:
            all_measures = measure_set.findall('.//Measure')
        else:
            all_measures = genotype_set.findall('.//MeasureSet//Measure')
            if not record.get('variant_type'):
                record['variant_type'] = genotype_set.get('Type', '')

        for measure in all_measures:
            if not record.get('variant_type'):
                record['variant_type'] = measure.get('Type', '')

            if not record.get('gene_symbol'):
                for gene in measure.findall('.//MeasureRelationship/Gene'):
                    symbol = gene.get('Symbol', '')
                    if symbol:
                        record['gene_symbol'] = symbol
                        record['gene_id'] = gene.get('GeneID', '')
                        break

            for attr in measure.findall('.//AttributeSet/Attribute'):
                attr_type = attr.get('Type', '')
                if 'HGVS' in attr_type and attr.text:
                    if ('coding' in attr_type.lower() or ':c.' in attr.text) and not record.get('hgvs_c'):
                        record['hgvs_c'] = attr.text
                    elif ('protein' in attr_type.lower() or ':p.' in attr.text) and not record.get('hgvs_p'):
                        record['hgvs_p'] = attr.text
                    elif ('genomic' in attr_type.lower() or ':g.' in attr.text) and not record.get('hgvs_g'):
                        record['hgvs_g'] = attr.text

            if not record.get('chromosome'):
                for loc in measure.findall('.//SequenceLocation'):
                    if loc.get('Assembly') == 'GRCh38':
                        record['chromosome'] = loc.get('Chr', '')
                        record['start'] = _safe_int(loc.get('start', 0) or 0)
                        record['stop'] = _safe_int(loc.get('stop', 0) or 0)
                        record['assembly'] = 'GRCh38'
                        break

            if not record.get('dbsnp_id'):
                for xref in measure.findall('.//XRef'):
                    if xref.get('DB') == 'dbSNP':
                        record['dbsnp_id'] = 'rs' + xref.get('ID', '')
                        break
    except Exception:
        pass  # Keep variant with whatever metadata was extracted

    # Submissions
    try:
        submissions = []
        for assertion in elem.findall('.//ClinVarAssertion'):
            sub = _extract_rcv_submission(assertion)
            if sub:
                submissions.append(sub)

        record['submissions'] = submissions
        record['num_submissions'] = len(submissions)

        all_pmids = set()
        all_acmg = set()
        all_submitters = set()

        for sub in submissions:
            all_pmids.update(sub.get('pmids', []))
            all_acmg.update(sub.get('acmg_codes', []))
            if sub.get('submitter'):
                all_submitters.add(sub['submitter'])

        record['all_pmids'] = list(all_pmids)
        record['all_acmg_codes'] = list(all_acmg)
        record['all_submitters'] = list(all_submitters)
    except Exception:
        record.setdefault('submissions', [])
        record.setdefault('num_submissions', 0)

    return record


def _extract_rcv_submission(assertion_elem, logger=None) -> Optional[Dict]:
    """Extract submission from ClinVarAssertion element."""
    try:
        sub = {
            'scv_id': assertion_elem.get('ID', ''),
            'date_submitted': assertion_elem.get('DateCreated', ''),
            'date_updated': assertion_elem.get('DateLastUpdated', ''),
        }

        submitter = assertion_elem.find('.//ClinVarSubmissionID')
        if submitter is not None:
            sub['submitter'] = submitter.get('submitter', '')

        cls_elem = assertion_elem.find('.//ClinicalSignificance')
        if cls_elem is not None:
            desc = cls_elem.find('Description')
            if desc is not None and desc.text:
                sub['classification'] = desc.text.strip()
            review = cls_elem.find('ReviewStatus')
            if review is not None and review.text:
                sub['review_status'] = review.text.strip()

        pmids = []
        for citation in assertion_elem.findall('.//Citation'):
            for id_elem in citation.findall('ID'):
                if id_elem.get('Source') == 'PubMed' and id_elem.text:
                    pmids.append(id_elem.text.strip())
        sub['pmids'] = pmids

        acmg_codes = []
        for comment in assertion_elem.findall('.//Comment'):
            if comment.text:
                codes = re.findall(
                    r'\b(PVS[1-9]|PS[1-9]\d?|PM[1-9]\d?|PP[1-9]\d?|BA[1-9]|BS[1-9]\d?|BP[1-9]\d?)\b',
                    comment.text.upper()
                )
                acmg_codes.extend(codes)
        sub['acmg_codes'] = list(set(acmg_codes))

        return sub

    except Exception as e:
        scv = assertion_elem.get('ID', 'unknown')
        if logger:
            logger.warning(f"Dropped RCV submission {scv}: {e}")
        return None


def _review_weight(review_status: str) -> int:
    """Convert review status to weight for sorting."""
    status = (review_status or '').lower()
    if 'practice guideline' in status:
        return 4
    elif 'expert panel' in status:
        return 3
    elif 'multiple submitters' in status:
        return 2
    elif 'criteria provided' in status:
        return 1
    return 0


def categorize_classification(classification: str) -> str:
    """Categorize classification into P/LP/P_LP/VUS/LB/B/B_LB/Conflicting/Other.

    FIX v5.1: "Pathogenic/Likely pathogenic" and "Benign/Likely benign" are
    now separate categories (P/LP, B/LB) instead of being lumped into LP/LB.
    """
    if not classification:
        return "Other"

    cls_lower = classification.lower().strip()

    if 'conflicting' in cls_lower:
        return "Conflicting"
    if 'uncertain' in cls_lower or 'vus' in cls_lower:
        return "VUS"
    # Check combined forms BEFORE individual ones
    if 'pathogenic/likely' in cls_lower or 'pathogenic, likely' in cls_lower:
        return "P/LP"
    if 'benign/likely' in cls_lower or 'benign, likely' in cls_lower:
        return "B/LB"
    if 'likely pathogenic' in cls_lower:
        return "LP"
    if 'pathogenic' in cls_lower:
        return "P"
    if 'likely benign' in cls_lower:
        return "LB"
    if 'benign' in cls_lower:
        return "B"

    return "Other"


# ============================================================================
# CHECKPOINT MANAGER
# ============================================================================

class CheckpointManager:
    """Manages checkpoint files for parsed snapshots."""

    def __init__(self, checkpoint_dir: Path, logger: logging.Logger):
        self.checkpoint_dir = checkpoint_dir
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger

    def get_checkpoint_path(self, snapshot_name: str) -> Path:
        """Get checkpoint file path for a snapshot."""
        # Remove .xml.gz extension, add .checkpoint.pkl.gz
        base_name = snapshot_name.replace('.xml.gz', '')
        return self.checkpoint_dir / f"{base_name}.checkpoint.pkl.gz"

    def has_checkpoint(self, snapshot_name: str) -> bool:
        """Check if checkpoint exists for a snapshot."""
        checkpoint_path = self.get_checkpoint_path(snapshot_name)
        return checkpoint_path.exists() and checkpoint_path.stat().st_size > 100

    def save_checkpoint(self, snapshot_name: str, data: Dict[str, Dict], snapshot_date: str):
        """Save parsed snapshot data to checkpoint."""
        checkpoint_path = self.get_checkpoint_path(snapshot_name)
        try:
            checkpoint_data = {
                'snapshot_date': snapshot_date,
                'snapshot_name': snapshot_name,
                'num_variants': len(data),
                'data': data
            }
            with gzip.open(checkpoint_path, 'wb') as f:
                pickle.dump(checkpoint_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            self.logger.info(f"  Saved checkpoint: {checkpoint_path.name} ({len(data):,} variants)")
        except Exception as e:
            self.logger.error(f"Failed to save checkpoint {checkpoint_path}: {e}")

    def load_checkpoint(self, snapshot_name: str) -> Optional[Tuple[str, Dict[str, Dict]]]:
        """Load checkpoint data. Returns (snapshot_date, data) or None."""
        checkpoint_path = self.get_checkpoint_path(snapshot_name)
        if not checkpoint_path.exists():
            return None
        try:
            with gzip.open(checkpoint_path, 'rb') as f:
                checkpoint_data = pickle.load(f)
            return checkpoint_data['snapshot_date'], checkpoint_data['data']
        except Exception as e:
            self.logger.error(f"Failed to load checkpoint {checkpoint_path}: {e}")
            return None

    def list_checkpoints(self) -> List[str]:
        """List all available checkpoints (filtering out corrupt/tiny ones)."""
        return [f.name.replace('.checkpoint.pkl.gz', '.xml.gz')
                for f in self.checkpoint_dir.glob('*.checkpoint.pkl.gz')
                if f.stat().st_size > 100]


# ============================================================================
# WORKER FUNCTION FOR MULTIPROCESSING
# ============================================================================

def parse_snapshot_worker(args: Tuple[Path, str, Path]) -> Tuple[str, str, int, Optional[str]]:
    """
    Worker function for parallel snapshot parsing with checkpointing.

    Returns: (snapshot_date, snapshot_name, num_variants, error_message or None)
    """
    snapshot_path, format_type, checkpoint_dir = args
    snapshot_name = snapshot_path.name

    # Extract date from filename
    match = re.search(r'(\d{4}-\d{2})', snapshot_name)
    snapshot_date = match.group(1) if match else snapshot_path.stem

    try:
        # Parse based on format
        if format_type == "rcv":
            results = parse_rcv_snapshot_full(snapshot_path)
        else:
            results = parse_vcv_snapshot_full(snapshot_path)

        # Save checkpoint
        checkpoint_path = checkpoint_dir / f"{snapshot_name.replace('.xml.gz', '')}.checkpoint.pkl.gz"
        checkpoint_data = {
            'snapshot_date': snapshot_date,
            'snapshot_name': snapshot_name,
            'num_variants': len(results),
            'data': results
        }
        with gzip.open(checkpoint_path, 'wb') as f:
            pickle.dump(checkpoint_data, f, protocol=pickle.HIGHEST_PROTOCOL)

        return snapshot_date, snapshot_name, len(results), None

    except Exception as e:
        error_msg = f"{e}\n{traceback.format_exc()}"
        return snapshot_date, snapshot_name, 0, error_msg


# ============================================================================
# MAIN ANALYZER
# ============================================================================

class ClinVarFTPAnalyzerV5:
    """Memory-efficient ClinVar analyzer with checkpointing."""

    def __init__(self, cache_dir: str, output_dir: str):
        self.cache_dir = Path(cache_dir)
        self.snapshot_dir = self.cache_dir / "ftp_snapshots"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Setup logging
        self.logger = setup_logging(self.output_dir)

        # Setup checkpoint manager
        self.checkpoint_dir = self.output_dir / "checkpoints"
        self.checkpoint_mgr = CheckpointManager(self.checkpoint_dir, self.logger)

        self.RCV_YEARS = [2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
        self.VCV_YEARS = [2020, 2021, 2022, 2023, 2024, 2025, 2026]

    def find_cached_snapshots(self, start_year: int = 2012,
                               end_year: int = 2027) -> List[Tuple[Path, str]]:
        """Find all cached snapshot files."""
        snapshots = []

        if not self.snapshot_dir.exists():
            self.logger.error(f"Snapshot directory not found: {self.snapshot_dir}")
            return snapshots

        for f in sorted(self.snapshot_dir.iterdir()):
            if not f.name.endswith('.xml.gz'):
                continue
            if f.stat().st_size < 1000:
                continue

            if 'ClinVarFullRelease' in f.name:
                format_type = 'rcv'
            elif 'ClinVarVariationRelease' in f.name or 'ClinVarVCVRelease' in f.name:
                format_type = 'vcv'
            else:
                continue

            match = re.search(r'(\d{4})', f.name)
            if match:
                year = int(match.group(1))
                if start_year <= year < end_year:
                    snapshots.append((f, format_type))

        def extract_date(item):
            match = re.search(r'(\d{4}-\d{2})', item[0].name)
            return match.group(1) if match else "9999-99"

        snapshots.sort(key=extract_date)
        return snapshots

    def run(self, start_year: int = 2012, end_year: int = 2027,
            max_workers: int = 4, batch_size: int = 500000) -> Dict:
        """Run analysis with checkpointing and memory-efficient processing."""

        start_time = time.time()

        self.logger.info("=" * 70)
        self.logger.info("ClinVar FTP-Only Analysis v5.1 (Memory-Efficient)")
        self.logger.info("=" * 70)
        self.logger.info(f"Year range: {start_year}-{end_year}")
        self.logger.info(f"Snapshot directory: {self.snapshot_dir}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Checkpoint directory: {self.checkpoint_dir}")
        self.logger.info(f"Workers: {max_workers}")

        try:
            # =================================================================
            # PHASE 1: Find snapshots
            # =================================================================
            self.logger.info("\n" + "=" * 70)
            self.logger.info("PHASE 1: Finding Cached Snapshots")
            self.logger.info("=" * 70)

            snapshots = self.find_cached_snapshots(start_year, end_year)

            if not snapshots:
                self.logger.error("No snapshots found!")
                return None

            self.logger.info(f"Found {len(snapshots)} snapshots")

            # Check for existing checkpoints
            existing_checkpoints = set(self.checkpoint_mgr.list_checkpoints())
            to_parse = [(p, f) for p, f in snapshots if p.name not in existing_checkpoints]
            already_parsed = [(p, f) for p, f in snapshots if p.name in existing_checkpoints]

            self.logger.info(f"Already parsed (checkpoints exist): {len(already_parsed)}")
            self.logger.info(f"Need to parse: {len(to_parse)}")

            # =================================================================
            # PHASE 2: Parse snapshots (with checkpointing)
            # =================================================================
            self.logger.info("\n" + "=" * 70)
            self.logger.info("PHASE 2: Parsing Snapshots")
            self.logger.info("=" * 70)

            if to_parse:
                self.logger.info(f"Parsing {len(to_parse)} snapshots with {max_workers} workers...")

                # Prepare args for workers
                parse_args = [(p, f, self.checkpoint_dir) for p, f in to_parse]

                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    futures = {executor.submit(parse_snapshot_worker, args): args for args in parse_args}

                    completed = 0
                    for future in as_completed(futures):
                        completed += 1
                        try:
                            snapshot_date, snapshot_name, num_variants, error = future.result()
                            if error:
                                self.logger.error(f"Error parsing {snapshot_name}: {error}")
                            else:
                                self.logger.info(f"[{completed}/{len(to_parse)}] Parsed {snapshot_name}: {num_variants:,} variants")
                        except Exception as e:
                            self.logger.error(f"Worker exception: {e}\n{traceback.format_exc()}")
            else:
                self.logger.info("All snapshots already parsed (using checkpoints)")

            # =================================================================
            # PHASE 3: Build history & stream results (memory-efficient)
            # =================================================================
            self.logger.info("\n" + "=" * 70)
            self.logger.info("PHASE 3: Building History & Writing Results")
            self.logger.info("=" * 70)

            stats = self._build_and_stream_results(snapshots)

            elapsed = time.time() - start_time
            self.logger.info("\n" + "=" * 70)
            self.logger.info(f"COMPLETE in {timedelta(seconds=int(elapsed))}")
            if stats:
                self.logger.info(f"Total variants: {stats['summary']['total_variants']:,}")
                self.logger.info(f"Changed: {stats['summary']['variants_with_changes']:,} "
                                 f"| Stable: {stats['summary']['variants_stable']:,}")
            self.logger.info("=" * 70)

            return stats

        except Exception as e:
            self.logger.error(f"FATAL ERROR: {e}\n{traceback.format_exc()}")
            raise

    def _build_and_stream_results(self, snapshots: List[Tuple[Path, str]]) -> Dict:
        """
        Memory-efficient Phase 3+4+5: load each checkpoint ONCE, build
        lightweight history in memory, then stream all output to CSV files.

        FIX v5.1: Replaces the old approach that reloaded every checkpoint
        for every batch of variants, causing massive RAM usage.
        """

        # String interning cache - classification and review_status strings
        # repeat across millions of variants; interning saves ~50% memory.
        _str_cache: Dict[str, str] = {}

        def _intern(s: str) -> str:
            if not s:
                return ''
            if s not in _str_cache:
                _str_cache[s] = s
            return _str_cache[s]

        # Lightweight data structures:
        #   variant_histories: var_id -> list of (date, classification, review_status, num_submissions)
        #   variant_latest:    var_id -> dict with only the fields needed for CSV output
        variant_histories: Dict[str, list] = {}
        variant_latest: Dict[str, dict] = {}

        # ---- Single pass: load each checkpoint exactly ONCE ----
        self.logger.info("Loading checkpoints (single pass)...")
        loaded_count = 0
        skipped_count = 0

        for snapshot_path, _ in snapshots:
            result = self.checkpoint_mgr.load_checkpoint(snapshot_path.name)
            if result is None:
                self.logger.warning(f"Skipping {snapshot_path.name} - no checkpoint available")
                skipped_count += 1
                continue

            snapshot_date, data = result
            loaded_count += 1
            self.logger.info(
                f"  [{loaded_count}] Loading {snapshot_path.name}: "
                f"{len(data):,} variants"
            )

            for var_id, record in data.items():
                # Append lightweight history entry
                cls_str = _intern(record.get('classification', ''))
                rev_str = _intern(record.get('review_status', ''))

                if var_id not in variant_histories:
                    variant_histories[var_id] = []

                variant_histories[var_id].append((
                    snapshot_date,
                    cls_str,
                    rev_str,
                    record.get('num_submissions', 0),
                ))

                # Update latest record fields.
                # - "State" fields (classification, review, submissions): always
                #   take the latest snapshot value (last wins).
                # - "Metadata" fields (gene, hgvs, coords, etc.): only overwrite
                #   if the new value is non-empty, so we don't lose good data
                #   when a later snapshot fails to extract a field.
                new_fields = {
                    'gene_symbol': record.get('gene_symbol', ''),
                    'gene_id': record.get('gene_id', ''),
                    'variant_type': record.get('variant_type', ''),
                    'molecular_consequence': record.get('molecular_consequence', ''),
                    'hgvs_c': record.get('hgvs_c', ''),
                    'hgvs_p': record.get('hgvs_p', ''),
                    'dbsnp_id': record.get('dbsnp_id', ''),
                    'allele_id': record.get('allele_id', ''),
                    'chromosome': record.get('chromosome', ''),
                    'start': record.get('start', 0),
                    'stop': record.get('stop', 0),
                    'assembly': record.get('assembly', ''),
                    'date_created': record.get('date_created', ''),
                    'date_updated': record.get('date_updated', ''),
                }

                if var_id not in variant_latest:
                    # First time seeing this variant — store everything
                    new_fields['classification'] = cls_str
                    new_fields['review_status'] = rev_str
                    new_fields['num_submissions'] = record.get('num_submissions', 0)
                    new_fields['num_pmids'] = len(record.get('all_pmids', []))
                    new_fields['num_acmg'] = len(record.get('all_acmg_codes', []))
                    variant_latest[var_id] = new_fields
                else:
                    existing = variant_latest[var_id]
                    # Always update state fields (latest snapshot = current state)
                    existing['classification'] = cls_str
                    existing['review_status'] = rev_str
                    existing['num_submissions'] = record.get('num_submissions', 0)
                    existing['num_pmids'] = len(record.get('all_pmids', []))
                    existing['num_acmg'] = len(record.get('all_acmg_codes', []))
                    # Only overwrite metadata if new value is non-empty
                    for key, val in new_fields.items():
                        if val:
                            existing[key] = val

            # Free checkpoint data immediately
            del data
            gc.collect()

        total_variants = len(variant_histories)
        self.logger.info(
            f"Loaded {loaded_count} checkpoints "
            f"({skipped_count} skipped), "
            f"{total_variants:,} unique variants"
        )

        if total_variants == 0:
            self.logger.error("No variant data found in any checkpoint!")
            return None

        # ---- Stream output to CSV files ----
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")

        summary_path = self.output_dir / f"summary_{ts}.csv"
        changed_path = self.output_dir / f"changed_variants_{ts}.csv"
        history_path = self.output_dir / f"classification_history_{ts}.csv"
        stats_path = self.output_dir / f"statistics_{ts}.json"

        CSV_HEADER = [
            'Variation_ID', 'VCV_ID', 'Gene', 'Variant_Type', 'Molecular_Consequence',
            'Current_Classification', 'Current_Category', 'Review_Status',
            'First_Classification', 'First_Date', 'Classification_Changed',
            'Change_Type', 'Num_Changes', 'Had_Oscillation',
            'Num_Submissions', 'Num_PMIDs', 'Num_ACMG_Codes',
            'HGVS_c', 'HGVS_p', 'dbSNP', 'Chromosome', 'Start', 'Stop', 'Assembly'
        ]
        HISTORY_HEADER = [
            'Variation_ID', 'VCV_ID', 'Gene', 'Snapshot_Date',
            'Classification', 'Category', 'Review_Status', 'Num_Submissions'
        ]

        # Stats accumulators
        changed_count = 0
        stable_count = 0
        oscillation_count = 0
        classification_dist = Counter()
        change_type_dist = Counter()
        gene_dist = Counter()
        submission_total = 0
        submission_max = 0
        submission_list_for_median = []
        processed = 0

        self.logger.info("Writing output files (streaming)...")

        with open(summary_path, 'w', newline='', encoding='utf-8') as sf, \
             open(changed_path, 'w', newline='', encoding='utf-8') as cf, \
             open(history_path, 'w', newline='', encoding='utf-8') as hf:

            summary_writer = csv.writer(sf)
            changed_writer = csv.writer(cf)
            history_writer = csv.writer(hf)

            summary_writer.writerow(CSV_HEADER)
            changed_writer.writerow(CSV_HEADER)
            history_writer.writerow(HISTORY_HEADER)

            for var_id, history in variant_histories.items():
                latest = variant_latest.get(var_id, {})

                # Sort history by date (should already be chronological,
                # but ensure correctness)
                history.sort(key=lambda x: x[0])

                # ---- Classification change analysis ----
                # Get category for each snapshot where variant had a classification.
                # Skip empty/Other classifications from snapshots where extraction
                # failed, so they don't create false P->Other->P changes.
                categories_seen = []
                for h in history:
                    cat = categorize_classification(h[1])
                    if h[1]:  # Only include if classification string is non-empty
                        categories_seen.append(cat)

                # If no valid classifications found, treat as single "Other"
                if not categories_seen:
                    categories_seen = [categorize_classification('')]

                # Deduplicate consecutive identical categories to get transition path
                # e.g. [P, P, VUS, VUS, P, P] -> [P, VUS, P]
                deduped = [categories_seen[0]]
                for c in categories_seen[1:]:
                    if c != deduped[-1]:
                        deduped.append(c)

                classification_changed = len(deduped) > 1
                change_type = ""
                num_changes = 0
                had_oscillation = False

                if classification_changed:
                    change_type = "->".join(deduped)
                    num_changes = len(deduped) - 1
                    # FIX v5.1: Oscillation = any category reappears in deduped
                    # e.g. [P, VUS, P] has P twice -> oscillation
                    had_oscillation = len(deduped) != len(set(deduped))

                # Build VCV ID
                try:
                    vcv_id = f"VCV{int(var_id):09d}"
                except (ValueError, TypeError):
                    vcv_id = f"VCV{var_id}"

                gene = latest.get('gene_symbol', '')
                current_cls = latest.get('classification', '')
                current_cat = categorize_classification(current_cls)
                first_valid_idx = next((i for i, h in enumerate(history) if h[1]), 0)
                first_cls = history[first_valid_idx][1] if history else ''
                first_date = history[first_valid_idx][0] if history else ''

                row = [
                    var_id, vcv_id, gene,
                    latest.get('variant_type', ''),
                    latest.get('molecular_consequence', ''),
                    current_cls, current_cat,
                    latest.get('review_status', ''),
                    first_cls, first_date,
                    classification_changed, change_type, num_changes, had_oscillation,
                    latest.get('num_submissions', 0),
                    latest.get('num_pmids', 0),
                    latest.get('num_acmg', 0),
                    latest.get('hgvs_c', ''),
                    latest.get('hgvs_p', ''),
                    latest.get('dbsnp_id', ''),
                    latest.get('chromosome', ''),
                    latest.get('start', '') or '',
                    latest.get('stop', '') or '',
                    latest.get('assembly', ''),
                ]

                summary_writer.writerow(row)

                if classification_changed:
                    changed_writer.writerow(row)
                    changed_count += 1
                else:
                    stable_count += 1

                if had_oscillation:
                    oscillation_count += 1

                # Write history rows (skip entries with empty classification)
                for h_date, h_cls, h_review, h_nsub in history:
                    if not h_cls:
                        continue
                    history_writer.writerow([
                        var_id, vcv_id, gene, h_date,
                        h_cls, categorize_classification(h_cls),
                        h_review, h_nsub
                    ])

                # Update stats accumulators
                classification_dist[current_cat] += 1
                if classification_changed:
                    change_type_dist[change_type] += 1
                if gene:
                    gene_dist[gene] += 1
                nsubs = latest.get('num_submissions', 0)
                submission_total += nsubs
                if nsubs > submission_max:
                    submission_max = nsubs
                submission_list_for_median.append(nsubs)

                processed += 1
                if processed % 500000 == 0:
                    self.logger.info(
                        f"  Written {processed:,} / {total_variants:,} variants..."
                    )

        # Compute final stats
        if submission_list_for_median:
            submission_list_for_median.sort()
            median_subs = submission_list_for_median[len(submission_list_for_median) // 2]
            mean_subs = submission_total / len(submission_list_for_median)
        else:
            median_subs = 0
            mean_subs = 0

        stats = {
            "summary": {
                "total_variants": processed,
                "variants_with_changes": changed_count,
                "variants_stable": stable_count,
                "variants_with_oscillation": oscillation_count,
            },
            "current_classification_distribution": dict(classification_dist),
            "change_type_distribution": dict(change_type_dist.most_common(50)),
            "genes": dict(gene_dist.most_common(100)),
            "submissions_per_variant": {
                "mean": round(mean_subs, 2),
                "median": median_subs,
                "max": submission_max,
            },
        }

        with open(stats_path, 'w') as f:
            json.dump(stats, f, indent=2, default=str)

        self.logger.info(f"\nResults saved:")
        self.logger.info(f"  Summary CSV: {summary_path}")
        self.logger.info(f"  Changed variants: {changed_path} ({changed_count:,})")
        self.logger.info(f"  History CSV: {history_path}")
        self.logger.info(f"  Statistics: {stats_path}")
        self.logger.info(f"\nVariants with changes: {changed_count:,}")
        self.logger.info(f"Variants stable: {stable_count:,}")
        self.logger.info(f"Variants with oscillation: {oscillation_count:,}")

        return stats


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="ClinVar FTP-Only Analysis v5.1 - Memory-Efficient with Checkpointing"
    )

    parser.add_argument("--cache-dir", default="./clinvar_cache",
                       help="Directory with cached FTP snapshots")
    parser.add_argument("--output-dir", "-o", default="./clinvar_output",
                       help="Output directory")
    parser.add_argument("--start-year", type=int, default=2012,
                       help="Start year (default: 2012)")
    parser.add_argument("--end-year", type=int, default=2027,
                       help="End year (default: 2027)")
    parser.add_argument("--workers", "-w", type=int, default=8,
                       help="Number of parallel workers (default: 8)")
    parser.add_argument("--batch-size", "-b", type=int, default=500000,
                       help="(Unused in v5.1, kept for CLI compatibility)")

    args = parser.parse_args()

    print("=" * 70)
    print("ClinVar FTP-Only Analysis v5.1")
    print("Memory-Efficient with Checkpointing & Streaming Output")
    print("=" * 70)

    analyzer = ClinVarFTPAnalyzerV5(
        cache_dir=args.cache_dir,
        output_dir=args.output_dir
    )

    try:
        stats = analyzer.run(
            start_year=args.start_year,
            end_year=args.end_year,
            max_workers=args.workers,
        )

        if stats:
            print("\n" + "=" * 70)
            print("FINAL SUMMARY")
            print("=" * 70)
            for key, value in stats["summary"].items():
                print(f"  {key}: {value:,}")
            return 0
        return 1

    except KeyboardInterrupt:
        print("\nInterrupted by user. Checkpoints saved - can resume later.")
        return 1
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
