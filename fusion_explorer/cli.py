"""CLI entry points for fusion analysis and patient-level clinical integration."""

import argparse
import csv
import json
from pathlib import Path
import sqlite3
import sys

from .core import Filters, analyze, import_callset
from .clinical import ClinicalFilters, combine_clinical, import_clinical, import_crosswalk


def fusion_filter_arguments(command):
    command.add_argument("--cohort", choices=["TCGA", "GTEx", "ALL", "UNKNOWN"], default="TCGA")
    command.add_argument("--tissue", choices=["all", "tumor", "normal", "unknown"], default="all")
    command.add_argument("--cancer", action="append", default=[], help="Repeat for multiple cancer types")
    command.add_argument("--gene", default="", help="Exact partner-gene match, case-insensitive")
    command.add_argument("--fusion", default="", help="Directed GENE1::GENE2 match")
    command.add_argument("--min-junction", type=int, default=0)
    command.add_argument("--min-spanning", type=int, default=0)
    command.add_argument("--min-ffpm", type=float, default=0)


def fusion_filters(args):
    return Filters(args.cohort, args.tissue, tuple(args.cancer), args.gene,
                   args.fusion, args.min_junction, args.min_spanning, args.min_ffpm)


def parser():
    root = argparse.ArgumentParser(description="TCGA Fusion Explorer — basic callset analysis")
    commands = root.add_subparsers(dest="command", required=True)
    ingest = commands.add_parser("import", help="Stream an existing CSV/TSV[.gz] into SQLite")
    ingest.add_argument("--input", required=True, type=Path)
    ingest.add_argument("--db", default=Path("data/fusions.sqlite"), type=Path)
    ingest.add_argument("--provenance", default="", help="Source URL, publication, version or notes")
    run = commands.add_parser("analyze", help="Import optionally, then analyze and export")
    run.add_argument("--input", type=Path, help="Optional callset to import before analysis")
    run.add_argument("--provenance", default="")
    run.add_argument("--db", default=Path("data/fusions.sqlite"), type=Path)
    run.add_argument("--source-id", type=int, help="Default: latest source; cannot be combined with --input")
    run.add_argument("--out", default=Path("results"), type=Path)
    fusion_filter_arguments(run)
    run.add_argument("--top", type=int, default=20)

    clinical = commands.add_parser("clinical-import", help="Import a patient-level clinical CSV/TSV[.gz]")
    clinical.add_argument("--input", type=Path, required=True)
    clinical.add_argument("--db", type=Path, default=Path("data/fusions.sqlite"))
    clinical.add_argument("--provenance", default="")

    crosswalk = commands.add_parser("mapping-import", help="Import a verified sample_id → patient_id crosswalk")
    crosswalk.add_argument("--input", type=Path, required=True)
    crosswalk.add_argument("--db", type=Path, default=Path("data/fusions.sqlite"))
    crosswalk.add_argument("--source-id", type=int, help="Fusion source; default latest")
    crosswalk.add_argument("--provenance", required=True, help="Crosswalk origin and verification notes")

    combine = commands.add_parser("clinical-combine", help="Left join clinical records and export a match audit")
    combine.add_argument("--db", type=Path, default=Path("data/fusions.sqlite"))
    combine.add_argument("--source-id", type=int, help="Fusion source; default latest")
    combine.add_argument("--clinical-source-id", type=int, help="Clinical source; default latest")
    combine.add_argument("--mapping-source-id", type=int, help="Explicit crosswalk selection; default canonical IDs only")
    combine.add_argument("--out", type=Path, default=Path("results/clinical"))
    fusion_filter_arguments(combine)
    combine.add_argument("--min-age", type=float, help="Inclusive age in years")
    combine.add_argument("--max-age", type=float, help="Inclusive age in years")
    combine.add_argument("--stage", action="append", default=[], help="Exact stage; repeat to include multiple stages")
    combine.add_argument("--sex", choices=["male", "female", "other"], default="")
    combine.add_argument("--vital-status", choices=["alive", "dead"], default="")
    combine.add_argument("--matched-only", action="store_true")
    combine.add_argument("--race", action="append", default=[], help="Race label; repeat to include multiple categories, including unknown")
    combine.add_argument("--histology", action="append", default=[], help="Source histology label or ICD-O-3 code, including unknown")
    combine.add_argument("--age-unknown-only", action="store_true", help="Only patients with missing age; mutually exclusive with age ranges")
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "import":
            result = import_callset(args.db, args.input, args.provenance)
        elif args.command == "clinical-import":
            result = import_clinical(args.db, args.input, args.provenance)
        elif args.command == "mapping-import":
            result = import_crosswalk(args.db, args.input, args.source_id, args.provenance)
        elif args.command == "clinical-combine":
            report = combine_clinical(
                args.db, args.out, fusion_filters(args),
                ClinicalFilters(args.min_age, args.max_age, tuple(args.stage),
                                args.sex, args.vital_status, args.matched_only,
                                tuple(args.race), tuple(args.histology), args.age_unknown_only),
                args.source_id, args.clinical_source_id, args.mapping_source_id)
            result = {"match_audit": report["match_audit_before_clinical_filters"],
                      "summary": report["summary_after_clinical_filters"],
                      "output": str(args.out.resolve())}
        else:
            if args.input and args.source_id is not None:
                raise ValueError("--input and --source-id cannot be combined")
            source_id = args.source_id
            if args.input:
                source_id = import_callset(args.db, args.input, args.provenance)["source_id"]
            filters = fusion_filters(args)
            report = analyze(args.db, args.out, filters, source_id, args.top)
            result = {"summary": report["summary"], "top_fusions": report["top_fusions"],
                      "output": str(args.out.resolve())}
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except (OSError, ValueError, sqlite3.Error, UnicodeError, csv.Error) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

