"""Source-scoped pathology metadata storage and read-only Dashboard queries."""

import csv
import json
from pathlib import Path
from urllib.parse import urlparse

from .clinical import ClinicalFilters, file_sha, table_rows, validate_row
from .core import Filters, connect, initialize, source_info
from .dashboard_data import catalog, readonly, rows, scoped
from .pathology_manifest import LINK_COLUMNS


def import_pathology(db_path, links_path, report_path, source_id=None):
    report = json.loads(Path(report_path).read_text())
    sha = file_sha(links_path)
    with connect(db_path) as db:
        initialize(db)
        source = source_info(db, source_id)
        db.executescript("""
            CREATE TABLE IF NOT EXISTS pathology_sources (
                id INTEGER PRIMARY KEY, fusion_source_id INTEGER NOT NULL REFERENCES sources(id),
                sha256 TEXT NOT NULL, report_sha256 TEXT NOT NULL, report_json TEXT NOT NULL,
                UNIQUE(fusion_source_id,sha256,report_sha256)
            );
            CREATE TABLE IF NOT EXISTS pathology_links (
                pathology_source_id INTEGER NOT NULL REFERENCES pathology_sources(id),
                sample_id TEXT NOT NULL, patient_id TEXT NOT NULL, file_id TEXT NOT NULL,
                metadata_json TEXT NOT NULL, PRIMARY KEY(pathology_source_id,sample_id,file_id)
            );
        """)
        report_sha = file_sha(report_path)
        existing = db.execute("SELECT id FROM pathology_sources WHERE fusion_source_id=? AND sha256=? AND report_sha256=?", (source["id"], sha, report_sha)).fetchone()
        if existing:
            return {"pathology_source_id": existing[0], "reused": True}
        count = 0
        with db, table_rows(links_path) as (fields, records):
            if not set(LINK_COLUMNS).issubset(fields):
                raise ValueError("Incomplete pathology manifest columns")
            pid = db.execute("INSERT INTO pathology_sources(fusion_source_id,sha256,report_sha256,report_json) VALUES(?,?,?,?)",
                             (source["id"], sha, report_sha, json.dumps(report))).lastrowid
            for record in records:
                validate_row(record)
                matching = db.execute("SELECT DISTINCT patient_id,cancer_type,cohort FROM calls WHERE source_id=? AND sample_id=?", (source["id"], record["fusion_sample_id"])).fetchall()
                if not matching or any(r["patient_id"] != record["patient_id"] or r["cancer_type"] != "SKCM" or r["cohort"] != "TCGA" for r in matching):
                    raise ValueError("Pathology patient/sample does not match the selected SKCM fusion source")
                if record["match_level"] != "patient_only_verified" or record["same_specimen_verified"] != "no":
                    raise ValueError("Only patient-level, specimen-unverified links are supported")
                for key in ("case_viewer_entry", "file_details_url"):
                    url = urlparse(record[key])
                    if url.scheme != "https" or url.netloc != "portal.gdc.cancer.gov":
                        raise ValueError("Unexpected GDC link destination")
                db.execute("INSERT INTO pathology_links VALUES(?,?,?,?,?)", (pid, record["fusion_sample_id"], record["patient_id"], record["file_id"], json.dumps(record)))
                count += 1
            if not count or count != report["counts"]["patient_slide_links"]:
                raise ValueError("Manifest count does not agree with the linkage report")
        return {"pathology_source_id": pid, "links": count, "reused": False}


def pathology_source(db_path, fusion_source_id):
    with readonly(db_path) as db:
        if not db.execute("SELECT 1 FROM sqlite_master WHERE name='pathology_sources'").fetchone():
            return None
        source = db.execute("SELECT * FROM pathology_sources WHERE fusion_source_id=? ORDER BY id DESC LIMIT 1", (fusion_source_id,)).fetchone()
        return dict(source) if source else None


def pathology_query(db_path, source_id, filters=Filters(), clinical_source_id=None,
                    mapping_source_id=None, clinical_filters=ClinicalFilters(),
                    patient="", consistency="all", page=0, page_size=50, export_path=None):
    if page < 0 or not 1 <= page_size <= 500:
        raise ValueError("Invalid pathology page")
    if consistency not in {"all", "same_type_hint", "different_type_hint", "unknown"}:
        raise ValueError("Invalid tissue-type consistency filter")
    source = pathology_source(db_path, source_id)
    if source is None:
        return {"available": False}
    with scoped(db_path, source_id, filters, clinical_source_id, mapping_source_id, clinical_filters) as db:
        db.execute("""CREATE TEMP TABLE pathology_scope AS
            SELECT l.* FROM pathology_links l JOIN (SELECT DISTINCT sample_id FROM dashboard_scope) s
            ON s.sample_id=l.sample_id WHERE l.pathology_source_id=?""", (source["id"],))
        clauses, args = [], []
        if patient.strip():
            clauses.append("UPPER(patient_id)=?")
            args.append(patient.strip().upper())
        if consistency != "all":
            clauses.append("json_extract(metadata_json,'$.sample_type_consistency')=?")
            args.append(consistency)
        where = " AND ".join(clauses) or "1=1"
        base = f"SELECT * FROM pathology_scope WHERE {where}"
        summary = rows(db, f"""SELECT COUNT(*) links,COUNT(DISTINCT sample_id) samples,
            COUNT(DISTINCT patient_id) patients,COUNT(DISTINCT file_id) slide_files,
            COUNT(CASE WHEN json_extract(metadata_json,'$.sample_type_consistency')='different_type_hint' THEN 1 END) different_type_links
            FROM ({base})""", args)[0]
        pages = (summary["links"] + page_size - 1) // page_size
        page = min(page, max(0, pages-1))
        records = [json.loads(r[0]) for r in db.execute(f"SELECT metadata_json FROM ({base}) ORDER BY patient_id,sample_id,file_id LIMIT ? OFFSET ?", [*args, page_size, page*page_size])]
        if export_path is not None:
            path = Path(export_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", newline="", encoding="utf-8-sig") as handle:
                writer = csv.DictWriter(handle, fieldnames=LINK_COLUMNS)
                writer.writeheader()
                writer.writerows(json.loads(r[0]) for r in db.execute(f"SELECT metadata_json FROM ({base}) ORDER BY patient_id,sample_id,file_id", args))
        return {"available": True, "source": source, "summary": summary,
                "selected_fusion_samples": db.execute("SELECT COUNT(DISTINCT sample_id) FROM dashboard_scope").fetchone()[0],
                "records": records, "page": page, "pages": pages}
