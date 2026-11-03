"""Streaming ingestion, validated normalization, and SQL aggregation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import re
import sqlite3


ALIASES = {
    "sample_id": ("sample_id", "sample", "Tumor_Sample_Barcode", "Sample_ID", "barcode"),
    "fusion": ("fusion_name", "FusionName", "#FusionName", "fusion", "Fusion"),
    "left_gene": ("LeftGene", "left_gene", "gene1", "Gene1"),
    "right_gene": ("RightGene", "right_gene", "gene2", "Gene2"),
    "cancer_type": ("cancer_type", "tissue_type", "Cancer_Type", "project_id"),
    "cohort": ("data_class", "cohort"),
    "tissue": ("tumor_or_normal", "tissue", "sample_class"),
    "junction_reads": ("Jcount", "JunctionReadCount", "junction_reads"),
    "spanning_reads": ("Scount", "SpanningFragCount", "spanning_reads"),
    "ffpm": ("FFPM", "ffpm"),
    "left_breakpoint": ("LeftBreakpoint", "left_breakpoint"),
    "right_breakpoint": ("RightBreakpoint", "right_breakpoint"),
    "fusion_type": ("PROT_FUSION_TYPE", "fusion_type"),
    "structure_type": ("structure_type",),
}
COLUMNS = tuple(ALIASES) + ("patient_id",)
MISSING = {"", "NA", "N/A", "NAN", "NULL", "NONE", "."}
BARCODE = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-|$)", re.I)
SAMPLE_TYPE = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-(\d{2})", re.I)


@contextmanager
def connect(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=60)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("PRAGMA cache_size=-16384")  # 16 MiB SQLite page cache.
    db.execute("PRAGMA temp_store=FILE")
    try:
        yield db
    finally:
        db.close()


def initialize(db):
    db.executescript("""
        CREATE TABLE IF NOT EXISTS sources (
            id INTEGER PRIMARY KEY, sha256 TEXT UNIQUE NOT NULL,
            filename TEXT NOT NULL, provenance TEXT NOT NULL,
            imported_at TEXT NOT NULL, report_json TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS calls (
            source_id INTEGER NOT NULL REFERENCES sources(id), row_sha256 TEXT NOT NULL,
            sample_id TEXT NOT NULL, fusion TEXT NOT NULL,
            left_gene TEXT NOT NULL, right_gene TEXT NOT NULL,
            cancer_type TEXT NOT NULL, cohort TEXT NOT NULL, tissue TEXT NOT NULL,
            junction_reads INTEGER, spanning_reads INTEGER, ffpm REAL,
            left_breakpoint TEXT, right_breakpoint TEXT, fusion_type TEXT,
            structure_type TEXT, patient_id TEXT,
            UNIQUE(source_id, row_sha256)
        );
        CREATE INDEX IF NOT EXISTS calls_source ON calls(source_id);
        CREATE INDEX IF NOT EXISTS calls_left ON calls(source_id, left_gene);
        CREATE INDEX IF NOT EXISTS calls_right ON calls(source_id, right_gene);
        CREATE INDEX IF NOT EXISTS calls_cancer ON calls(source_id, cancer_type);
    """)


def clean(value):
    value = str(value or "").strip()
    return "" if value.upper() in MISSING else value


def gene(value):
    return clean(value).split("^", 1)[0].upper()


def number(value, label, integer=False):
    value = clean(value)
    if not value:
        return None
    try:
        parsed = float(value)
    except ValueError:
        raise ValueError(f"Invalid {label}: {value!r}") from None
    if not math.isfinite(parsed) or parsed < 0 or (integer and not parsed.is_integer()):
        raise ValueError(f"Invalid {label}: {value!r}")
    return int(parsed) if integer else parsed


def normalize(raw, mapping):
    def get(name):
        return clean(raw.get(mapping.get(name, "")))

    sample = get("sample_id")
    if not sample:
        raise ValueError("Missing sample ID")
    left, right = gene(get("left_gene")), gene(get("right_gene"))
    name = get("fusion")
    partners = [gene(x) for x in name.split("::")]
    if name and (len(partners) != 2 or not all(partners)):
        raise ValueError("Fusion name must use directed GENE1::GENE2 format")
    if not left and partners:
        left = partners[0]
    if not right and len(partners) == 2:
        right = partners[1]
    if not left or not right:
        raise ValueError("Missing fusion partners")
    if name and partners != [left, right]:
        raise ValueError("Fusion name disagrees with left/right gene columns")
    patient = sample[:12].upper() if BARCODE.match(sample) else None
    cohort = get("cohort").upper() or ("TCGA" if patient else "UNKNOWN")
    if cohort == "GTEX":
        cohort = "GTEx"
    tissue = get("tissue").lower()
    tissue = {"tumour": "tumor"}.get(tissue, tissue)
    sample_type = SAMPLE_TYPE.match(sample)
    if not tissue and sample_type:
        code = int(sample_type.group(1))
        tissue = "tumor" if 1 <= code <= 9 or code == 40 else "normal" if 10 <= code <= 19 else "unknown"
    if tissue not in {"tumor", "normal"}:
        tissue = "unknown"
    return {
        "sample_id": sample, "patient_id": patient,
        "fusion": f"{left}::{right}", "left_gene": left, "right_gene": right,
        "cancer_type": get("cancer_type").upper().removeprefix("TCGA-") or "UNKNOWN",
        "cohort": cohort, "tissue": tissue,
        "junction_reads": number(get("junction_reads"), "junction reads", True),
        "spanning_reads": number(get("spanning_reads"), "spanning reads", True),
        "ffpm": number(get("ffpm"), "FFPM"),
        "left_breakpoint": get("left_breakpoint"),
        "right_breakpoint": get("right_breakpoint"),
        "fusion_type": get("fusion_type"), "structure_type": get("structure_type"),
    }


def import_callset(db_path, input_path, provenance="", batch_size=5000):
    """Atomic import. Exact raw-row duplicates are removed within one source.

    SHA-256 identifies the original bytes. All cohorts are stored; query defaults
    select TCGA. Unsupported shorthand barcodes are never guessed into patients.
    """
    path = Path(input_path).expanduser().resolve()
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    sha = digest.hexdigest()
    with connect(db_path) as db:
        initialize(db)
        existing = db.execute("SELECT id, report_json FROM sources WHERE sha256=?", (sha,)).fetchone()
        if existing:
            return {**json.loads(existing["report_json"]), "source_id": existing["id"], "reused": True}
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        report = {"sha256": sha, "rows_read": 0, "imported": 0, "duplicates": 0,
                  "rejected": 0, "missing_patient_id": 0, "missing_support": 0,
                  "unknown_cohort": 0, "unknown_cancer_type": 0,
                  "cohorts": {}, "rejection_examples": [], "reused": False}
        with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
            # Skip metadata comments, but preserve STAR-Fusion's #FusionName header.
            header = next((line for line in handle if line.strip() and
                           (not line.startswith("#") or line.startswith("#FusionName"))), "")
            if not header:
                raise ValueError("Empty input file")
            delimiter = "\t" if "\t" in header else ","
            fields = next(csv.reader([header], delimiter=delimiter))
            if len(fields) != len(set(fields)):
                raise ValueError("Duplicate column names")
            mapping = {name: next((a for a in aliases if a in fields), "")
                       for name, aliases in ALIASES.items()}
            if not mapping["sample_id"] or not (mapping["fusion"] or
                                                (mapping["left_gene"] and mapping["right_gene"])):
                raise ValueError("Input needs a sample ID and fusion name or both partner genes; see README")
            reader = csv.DictReader(handle, fieldnames=fields, delimiter=delimiter)
            with db:
                cursor = db.execute(
                    "INSERT INTO sources(sha256,filename,provenance,imported_at,report_json) VALUES(?,?,?,?,?)",
                    (sha, path.name, provenance, datetime.now(timezone.utc).isoformat(), "{}"))
                source_id = cursor.lastrowid
                placeholders = ",".join("?" for _ in range(len(COLUMNS) + 2))
                sql = f"INSERT OR IGNORE INTO calls(source_id,row_sha256,{','.join(COLUMNS)}) VALUES({placeholders})"
                batch = []

                def flush():
                    before = db.total_changes
                    db.executemany(sql, batch)
                    added = db.total_changes - before
                    report["imported"] += added
                    report["duplicates"] += len(batch) - added
                    batch.clear()

                for row_number, raw in enumerate(reader, start=2):
                    report["rows_read"] += 1
                    try:
                        if None in raw or any(v is None for v in raw.values()):
                            raise ValueError("Row column count differs from header")
                        normalized = normalize(raw, mapping)
                    except ValueError as error:
                        report["rejected"] += 1
                        if len(report["rejection_examples"]) < 20:
                            report["rejection_examples"].append({"record": row_number, "reason": str(error)})
                        continue
                    # Quality counters describe accepted input rows before deduplication.
                    report["missing_patient_id"] += normalized["patient_id"] is None
                    report["missing_support"] += any(normalized[x] is None for x in ("junction_reads", "spanning_reads", "ffpm"))
                    report["unknown_cohort"] += normalized["cohort"] == "UNKNOWN"
                    report["unknown_cancer_type"] += normalized["cancer_type"] == "UNKNOWN"
                    cohort = normalized["cohort"]
                    report["cohorts"][cohort] = report["cohorts"].get(cohort, 0) + 1
                    row_sha = hashlib.sha256(json.dumps(raw, sort_keys=True).encode()).hexdigest()
                    batch.append((source_id, row_sha, *(normalized[x] for x in COLUMNS)))
                    if len(batch) >= batch_size:
                        flush()
                flush()
                if not report["imported"]:
                    raise ValueError(f"No valid fusion calls; import rolled back. Rejections: {report['rejection_examples']}")
                db.execute("UPDATE sources SET report_json=? WHERE id=?", (json.dumps(report), source_id))
        return {**report, "source_id": source_id}


@dataclass(frozen=True)
class Filters:
    cohort: str = "TCGA"
    tissue: str = "all"
    cancers: tuple[str, ...] = ()
    gene: str = ""
    fusion: str = ""
    min_junction: int = 0
    min_spanning: int = 0
    min_ffpm: float = 0

    def where(self, source_id):
        clauses, args = ["source_id=?"], [source_id]
        if self.cohort != "ALL":
            clauses.append("cohort=?")
            args.append(self.cohort)
        if self.tissue != "all":
            clauses.append("tissue=?")
            args.append(self.tissue)
        if self.cancers:
            clauses.append(f"cancer_type IN ({','.join('?' for _ in self.cancers)})")
            args.extend(x.upper().removeprefix("TCGA-") for x in self.cancers)
        if self.gene:
            clauses.append("(left_gene=? OR right_gene=?)")
            args.extend([gene(self.gene)] * 2)
        if self.fusion:
            clauses.append("fusion=?")
            args.append(self.fusion.upper())
        for name, value in (("junction_reads", self.min_junction),
                            ("spanning_reads", self.min_spanning), ("ffpm", self.min_ffpm)):
            if not math.isfinite(value) or value < 0:
                raise ValueError("Support thresholds must be finite and nonnegative")
            if value > 0:
                clauses.append(f"{name}>=?")
                args.append(value)
        return " AND ".join(clauses), args


def source_info(db, source_id=None):
    if source_id is None:
        row = db.execute("SELECT * FROM sources ORDER BY id DESC LIMIT 1").fetchone()
    else:
        row = db.execute("SELECT * FROM sources WHERE id=?", (source_id,)).fetchone()
    if row is None:
        raise ValueError("No matching imported source. Import a callset first.")
    return dict(row)


def analyze(db_path, out_dir, filters=Filters(), source_id=None, top=20):
    """SQL aggregation and streaming export; never materializes the callset."""
    if top < 1:
        raise ValueError("top must be positive")
    if not Path(db_path).is_file():
        raise ValueError("Database does not exist. Import a callset first.")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with connect(db_path) as db:
        initialize(db)
        source = source_info(db, source_id)
        where, args = filters.where(source["id"])
        count_sql = """COUNT(*) calls, COUNT(DISTINCT sample_id) samples,
            COUNT(DISTINCT patient_id) patients_with_canonical_id,
            COUNT(DISTINCT fusion) distinct_fusions,
            COUNT(DISTINCT CASE WHEN cancer_type!='UNKNOWN' THEN cancer_type END) cancer_types"""
        summary = dict(db.execute(f"SELECT {count_sql} FROM calls WHERE {where}", args).fetchone())

        def export(name, query, params):
            cursor = db.execute(query, params)
            with (out / name).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(x[0] for x in cursor.description)
                writer.writerows(cursor)

        export("calls.csv", f"SELECT {','.join(COLUMNS)} FROM calls WHERE {where}", args)
        export("cancer_summary.csv", f"SELECT cancer_type, {count_sql} FROM calls WHERE {where} GROUP BY cancer_type ORDER BY calls DESC, cancer_type", args)
        export("fusion_summary.csv", f"SELECT fusion, COUNT(*) calls, COUNT(DISTINCT sample_id) samples FROM calls WHERE {where} GROUP BY fusion ORDER BY samples DESC, calls DESC, fusion", args)
        export("sample_burden.csv", f"SELECT sample_id, cancer_type, cohort, tissue, COUNT(*) calls, COUNT(DISTINCT fusion) distinct_fusions FROM calls WHERE {where} GROUP BY sample_id,cancer_type,cohort,tissue ORDER BY distinct_fusions DESC,sample_id", args)
        # Exclude the second partner of self-fusions; distinct samples are counted below.
        genes = f"""WITH selected AS (SELECT * FROM calls WHERE {where}),
            occurrences AS (
                SELECT left_gene gene, sample_id FROM selected
                UNION ALL SELECT right_gene gene, sample_id FROM selected WHERE right_gene!=left_gene
            ) SELECT gene, COUNT(*) calls, COUNT(DISTINCT sample_id) samples
              FROM occurrences GROUP BY gene ORDER BY samples DESC,calls DESC,gene"""
        export("gene_summary.csv", genes, args)
        export("tissue_summary.csv", f"SELECT cohort,tissue,COUNT(*) calls,COUNT(DISTINCT sample_id) samples FROM calls WHERE {where} GROUP BY cohort,tissue ORDER BY cohort,tissue", args)
        top_fusions = [dict(r) for r in db.execute(f"SELECT fusion,COUNT(DISTINCT sample_id) samples,COUNT(*) calls FROM calls WHERE {where} GROUP BY fusion ORDER BY samples DESC,calls DESC,fusion LIMIT ?", [*args, top])]
        report = {
            "summary": summary, "filters": filters.__dict__,
            "source": {k: v for k, v in source.items() if k != "report_json"},
            "import_quality": json.loads(source["report_json"]),
            "top_fusions": top_fusions,
            "notes": [
                "Counts describe the supplied callset, not population prevalence or confirmed fusion positivity.",
                "Samples with no reported calls are absent; a tested-sample denominator is not available.",
                "Exact raw-row duplicates are removed within a source; breakpoint-distinct calls remain separate.",
                "Fusion partner order is preserved. Missing support is retained unless a positive threshold is applied.",
                "Shorthand IDs are preserved; only canonical TCGA barcodes produce patient IDs.",
                "Latest imported source is used by default; separate sources are never silently pooled.",
                "Quality counters and cohort counts refer to accepted input rows before deduplication.",
            ],
        }
        (out / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return report

