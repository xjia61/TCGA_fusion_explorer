"""Read-only, bounded Dashboard queries. No Streamlit dependency in this layer."""

from contextlib import contextmanager
import csv
import io
import json
from pathlib import Path
import sqlite3

from .core import COLUMNS, Filters
from .clinical import ClinicalFilters, clinical_join_query, category_expression


@contextmanager
def readonly(path):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise ValueError("Database not found. Import data first.")
    db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True, timeout=30)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA cache_size=-16384")
    db.execute("PRAGMA temp_store=FILE")
    try:
        yield db
    finally:
        db.close()


def rows(db, sql, args=()):
    return [dict(r) for r in db.execute(sql, args)]


def catalog(path):
    with readonly(path) as db:
        tables = {r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if not {"sources", "calls"}.issubset(tables):
            raise ValueError("This is not a Fusion Explorer database.")
        return {
            "fusion_sources": rows(db, "SELECT * FROM sources ORDER BY id DESC"),
            "clinical_sources": rows(db, "SELECT * FROM clinical_sources ORDER BY id DESC") if "clinical_sources" in tables else [],
            "mapping_sources": rows(db, "SELECT * FROM mapping_sources ORDER BY id DESC") if "mapping_sources" in tables else [],
        }


def options(path, source_id, clinical_source_id=None):
    with readonly(path) as db:
        result = {
            "cancers": [r[0] for r in db.execute("SELECT DISTINCT cancer_type FROM calls WHERE source_id=? ORDER BY cancer_type", (source_id,))],
            "cohorts": [r[0] for r in db.execute("SELECT DISTINCT cohort FROM calls WHERE source_id=? ORDER BY cohort", (source_id,))],
            "has_ffpm": bool(db.execute("SELECT 1 FROM calls WHERE source_id=? AND ffpm IS NOT NULL LIMIT 1", (source_id,)).fetchone()),
            "stages": [],
            "races": [], "histologies": [],
        }
        if clinical_source_id is not None:
            result["stages"] = [r[0] for r in db.execute("SELECT DISTINCT stage FROM clinical_records WHERE clinical_source_id=? AND stage IS NOT NULL ORDER BY stage", (clinical_source_id,))]
            for name, key in (("race", "races"), ("histology", "histologies")):
                result[key] = [r[0] for r in db.execute(f"SELECT DISTINCT {category_expression(name)} category FROM clinical_records c WHERE clinical_source_id=? ORDER BY category", (clinical_source_id,))]
        return result


@contextmanager
def scoped(path, source_id, filters, clinical_source_id=None, mapping_source_id=None,
           clinical_filters=ClinicalFilters()):
    """Snapshot the selected calls to disk-backed temporary tables, never main DB."""
    info = catalog(path)
    if not any(r["id"] == source_id for r in info["fusion_sources"]):
        raise ValueError("Fusion source not found.")
    if clinical_source_id is not None and not any(r["id"] == clinical_source_id for r in info["clinical_sources"]):
        raise ValueError("Clinical source not found.")
    if mapping_source_id is not None:
        if clinical_source_id is None or not any(r["id"] == mapping_source_id and r["fusion_source_id"] == source_id for r in info["mapping_sources"]):
            raise ValueError("The crosswalk does not belong to the selected fusion source, or no clinical source is selected.")
    clinical_where, clinical_args = clinical_filters.where()
    if clinical_source_id is None and clinical_where != "1=1":
        raise ValueError("Select a clinical source before applying clinical filters.")
    with readonly(path) as db:
        if clinical_source_id is not None:
            query, args = clinical_join_query(filters, source_id, clinical_source_id, mapping_source_id)
        else:
            where, args = filters.where(source_id)
            query = f"""SELECT *,CASE WHEN cohort='TCGA' THEN patient_id END resolved_patient_id,
                CASE WHEN cohort!='TCGA' THEN 'not_tcga' WHEN patient_id IS NOT NULL THEN 'canonical' ELSE 'unresolved' END patient_id_origin,
                'clinical_not_selected' match_status,NULL clinical_patient_id,NULL clinical_cancer_type,
                NULL age_years,NULL stage,NULL sex,NULL vital_status,'unknown' race,'unknown' histology FROM calls WHERE {where}"""
        db.execute("CREATE TEMP TABLE dashboard_base AS " + query, args)
        db.execute("CREATE TEMP TABLE dashboard_scope AS SELECT * FROM dashboard_base WHERE " + clinical_where, clinical_args)
        yield db


def dashboard_query(path, source_id, filters=Filters(), clinical_source_id=None,
                    mapping_source_id=None, clinical_filters=ClinicalFilters(), page=0, page_size=50):
    if page < 0 or page_size < 1 or page_size > 500:
        raise ValueError("Invalid page or page size")
    with scoped(path, source_id, filters, clinical_source_id, mapping_source_id, clinical_filters) as db:
        summary = rows(db, """SELECT COUNT(*) calls,COUNT(DISTINCT sample_id) samples,
            COUNT(DISTINCT resolved_patient_id) resolved_patients,COUNT(DISTINCT fusion) distinct_fusions,
            COUNT(DISTINCT CASE WHEN match_status='matched' THEN resolved_patient_id END) matched_patients,
            COUNT(CASE WHEN ffpm IS NULL THEN 1 END) missing_ffpm,
            COUNT(DISTINCT CASE WHEN cancer_type!='UNKNOWN' THEN cancer_type END) cancer_types FROM dashboard_scope""")[0]
        audit = rows(db, """SELECT match_status,COUNT(*) sample_records FROM
            (SELECT DISTINCT sample_id,cohort,cancer_type,match_status FROM dashboard_base)
            GROUP BY match_status ORDER BY sample_records DESC""")
        tcga_records = sum(r["sample_records"] for r in audit if r["match_status"] != "not_tcga")
        matched_records = sum(r["sample_records"] for r in audit if r["match_status"] == "matched")
        top = rows(db, """SELECT fusion,COUNT(*) calls,COUNT(DISTINCT sample_id) samples FROM dashboard_scope
            GROUP BY fusion ORDER BY samples DESC,calls DESC,fusion LIMIT 20""")
        cancer = rows(db, """SELECT cancer_type,COUNT(*) calls,COUNT(DISTINCT sample_id) samples,
            COUNT(DISTINCT fusion) distinct_fusions FROM dashboard_scope GROUP BY cancer_type ORDER BY samples DESC,cancer_type""")
        burden = rows(db, """SELECT bucket,COUNT(*) samples FROM (
            SELECT CASE WHEN n=1 THEN '01 · 1' WHEN n<=5 THEN '02 · 2–5' WHEN n<=10 THEN '03 · 6–10'
            WHEN n<=25 THEN '04 · 11–25' WHEN n<=50 THEN '05 · 26–50' ELSE '06 · 51+' END bucket
            FROM (SELECT sample_id,COUNT(DISTINCT fusion) n FROM dashboard_scope GROUP BY sample_id))
            GROUP BY bucket ORDER BY bucket""")
        genes = rows(db, """WITH partners AS (
            SELECT left_gene gene,sample_id FROM dashboard_scope UNION ALL
            SELECT right_gene gene,sample_id FROM dashboard_scope WHERE left_gene!=right_gene)
            SELECT gene,COUNT(DISTINCT sample_id) samples FROM partners GROUP BY gene ORDER BY samples DESC,gene LIMIT 20""")
        heatmap = rows(db, """SELECT fusion,cancer_type,COUNT(DISTINCT sample_id) samples FROM dashboard_scope
            WHERE fusion IN (SELECT fusion FROM dashboard_scope GROUP BY fusion ORDER BY COUNT(DISTINCT sample_id) DESC,fusion LIMIT 12)
            GROUP BY fusion,cancer_type""")
        patient_cte = "WITH patients AS (SELECT DISTINCT resolved_patient_id,age_years,stage,sex,vital_status FROM dashboard_scope WHERE match_status='matched') "
        ages = rows(db, patient_cte + "SELECT CAST(age_years/5 AS INTEGER)*5 age_from,COUNT(*) patients FROM patients WHERE age_years IS NOT NULL GROUP BY age_from ORDER BY age_from")
        stages = rows(db, patient_cte + "SELECT COALESCE(stage,'unknown') stage,COUNT(*) patients FROM patients GROUP BY stage ORDER BY stage")
        vital = rows(db, patient_cte + "SELECT COALESCE(vital_status,'unknown') vital_status,COUNT(*) patients FROM patients GROUP BY vital_status ORDER BY vital_status")
        missing_age = db.execute(patient_cte + "SELECT COUNT(*) FROM patients WHERE age_years IS NULL").fetchone()[0]
        max_page = max(0, (summary["calls"] - 1) // page_size)
        actual_page = min(page, max_page)
        selected_columns = (*COLUMNS, "resolved_patient_id", "match_status", "age_years", "stage", "sex", "vital_status", "race", "histology")
        races = rows(db, "SELECT race,COUNT(DISTINCT resolved_patient_id) patients FROM dashboard_scope WHERE match_status='matched' GROUP BY race ORDER BY patients DESC,race")
        histologies = rows(db, "SELECT histology,COUNT(DISTINCT resolved_patient_id) patients FROM dashboard_scope WHERE match_status='matched' GROUP BY histology ORDER BY patients DESC,histology")
        preview = rows(db, f"SELECT {','.join(selected_columns)} FROM dashboard_scope ORDER BY sample_id,fusion,row_sha256 LIMIT ? OFFSET ?", (page_size, actual_page * page_size))
        return {"summary": summary, "audit": audit, "match_fraction": matched_records/tcga_records if clinical_source_id is not None and tcga_records else None,
                "top_fusions": top, "cancers": cancer, "burden": burden, "genes": genes, "heatmap": heatmap,
                "ages": ages, "missing_age_patients": missing_age, "stages": stages, "vital": vital, "races": races, "histologies": histologies, "preview": preview, "page": actual_page,
                "page_count": max_page+1 if summary["calls"] else 0}


def csv_bytes(records):
    if not records:
        return b""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=records[0].keys())
    writer.writeheader()
    writer.writerows(records)
    return buffer.getvalue().encode("utf-8-sig")
