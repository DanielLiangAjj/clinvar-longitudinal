#!/usr/bin/env python
"""R4: manuscript Figure 4 (XML structural evolution), rendered with the module's house style
(as the original was) after the panel a/d/e layout fixes."""
import json, sys
from pathlib import Path
REV = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses")
EXT = REV / "external_dependencies" / "ClinVar Research Project"
sys.path.insert(0, str(EXT / "scripts"))
import clinvar_visualizations as V
V.set_nature_style()
xml = json.load(open(EXT / "data" / "xml_structure_analysis.json"))
hist = json.load(open("/Users/danielliang/Downloads/ClinVar_Research_Project/figures/.cache/history_agg.json"))
V.figure5(xml, REV / "figures" / "r4", hist)
