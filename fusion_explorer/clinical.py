"""Patient-level clinical import and explicit, source-scoped sample crosswalks.

No shorthand barcode guessing, remote fetching, or survival endpoint derivation.
"""

from __future__ import annotations

from contextlib import contextmanager
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

from .core import COLUMNS, Filters, clean, connect, initialize, number, source_info


PATIENT = re.compile(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", re.I)
ALIASES = {
    "patient_id": ("patient_id", "PATIENT_ID", "bcr_patient_barcode", "case_submitter_id", "cases.submitter_id", "submitter_id"),
    "cancer_type": ("cancer_type", "project_id", "project.project_id", "cases.project.project_id", "CANCER_TYPE_ACRONYM"),
    "age_years": ("age_years", "age_at_initial_pathologic_diagnosis", "AGE"),
    "age_days": ("age_at_diagnosis", "diagnoses.age_at_diagnosis"),
    "stage": ("stage", "ajcc_pathologic_stage", "diagnoses.ajcc_pathologic_stage", "AJCC_PATHOLOGIC_TUMOR_STAGE"),
    "sex": ("sex", "gender", "demographic.gender", "SEX"),
    "vital_status": ("vital_status", "demographic.vital_status", "VITAL_STATUS"),
}
CLINICAL_FIELDS = ("patient_id", "cancer_type", "age_years", "stage", "sex", "vital_status")
CLINICAL_COLUMNS = ("clinical_patient_id", "clinical_cancer_type", "age_years", "stage", "sex", "vital_status", "race", "histology")
CATEGORY_ALIASES = {
    "race": ("race", "RACE", "demographic.race"),
    "histology": ("histology", "histological_type", "primary_diagnosis", "diagnoses.primary_diagnosis", "ICD_O_3_HISTOLOGY"),
}
EXTRA_MISSING = {"NOT REPORTED", "NOT AVAILABLE", "UNKNOWN", "[NOT AVAILABLE]", "[NOT APPLICABLE]", "[UNKNOWN]", "--"}


def value(text):
    text = clean(text)
    return "" if text.upper() in EXTRA_MISSING else text


def patient_id(text):
    text = value(text).upper()
    if not PATIENT.fullmatch(text):
        raise ValueError("Expected a canonical patient barcode TCGA-XX-XXXX (not a sample barcode or UUID)")
    return text


def stage_name(text):
    return re.sub(r"^STAGE\s+", "", value(text).upper())


def clinical_category(raw, name):
    for alias in CATEGORY_ALIASES[name]:
        found = value(raw.get(alias, ""))
        if found:
            return found
    return "unknown"


def category_expression(name, alias="c"):
    """Read source categories from retained JSON, including older databases."""
    terms = []
    for key in CATEGORY_ALIASES[name]:
        path = '$."' + key + '"'
        expression = f"TRIM(CAST(json_extract({alias}.raw_json,'{path}') AS TEXT))"
        missing = sorted(EXTRA_MISSING | {"", "NA", "N/A", "NAN", "NULL", "NONE", "."})
        quoted = ",".join("'" + text.replace("'", "''") + "'" for text in missing)
        terms.append(f"CASE WHEN UPPER({expression}) NOT IN ({quoted}) THEN {expression} END")
    return "COALESCE(" + ",".join(terms) + ",'unknown')"


def initialize_clinical(db):
    initialize(db)
    db.executescript("""
        CREATE TABLE IF NOT EXISTS clinical_sources (
            id INTEGER PRIMARY KEY, sha256 TEXT NOT NULL UNIQUE,
            filename TEXT NOT NULL, provenance TEXT NOT NULL,
            imported_at TEXT NOT NULL, report_json TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS clinical_records (
            clinical_source_id INTEGER NOT NULL REFERENCES clinical_sources(id),
            patient_id TEXT NOT NULL, cancer_type TEXT, age_years REAL,
            stage TEXT, sex TEXT, vital_status TEXT, raw_json TEXT NOT NULL,
            PRIMARY KEY (clinical_source_id, patient_id)
        );
        CREATE TABLE IF NOT EXISTS mapping_sources (
            id INTEGER PRIMARY KEY, sha256 TEXT NOT NULL,
            fusion_source_id INTEGER NOT NULL REFERENCES sources(id),
            filename TEXT NOT NULL, provenance TEXT NOT NULL,
            imported_at TEXT NOT NULL, report_json TEXT NOT NULL,
            UNIQUE(sha256, fusion_source_id)
        );
        CREATE TABLE IF NOT EXISTS sample_mappings (
            mapping_source_id INTEGER NOT NULL REFERENCES mapping_sources(id),
            sample_id TEXT NOT NULL, patient_id TEXT NOT NULL,
            PRIMARY KEY (mapping_source_id, sample_id)
        );
        CREATE INDEX IF NOT EXISTS calls_sample ON calls(source_id,sample_id);
    """)


@contextmanager
def table_rows(path):
    """Stream a flat table, including cBioPortal-style metadata comment lines."""
    path = Path(path).expanduser().resolve()
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next((line for line in handle if line.strip() and not line.startswith("#")), "")
        if not header:
            raise ValueError("Empty input table")
        delimiter = "\t" if "\t" in header else ","
        fields = next(csv.reader([header], delimiter=delimiter))
        if len(fields) != len(set(fields)):
            raise ValueError("Duplicate column names")
        yield fields, csv.DictReader(handle, fieldnames=fields, delimiter=delimiter)


def file_sha(path):
    digest = hashlib.sha256()
    with Path(path).expanduser().open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_row(row):
    if None in row or any(v is None for v in row.values()):
        raise ValueError("Row column count differs from header")


def normalize_clinical(raw, fields):
    def get(name):
        candidates = [value(raw[a]) for a in ALIASES[name] if a in fields and value(raw[a])]
        return candidates[0] if candidates else ""

    # An extra patient identifier column must not silently disagree with the chosen one.
    ids = {patient_id(raw[a]) for a in ALIASES["patient_id"] if a in fields and value(raw[a])}
    if len(ids) != 1:
        raise ValueError("Missing or conflicting patient IDs")
    years = number(get("age_years"), "age in years")
    days = number(get("age_days"), "age in days")
    if years is not None and days is not None and abs(years - days / 365.25) > 1:
        raise ValueError("Age-in-years and age-in-days columns disagree by more than one year")
    age = years if years is not None else days / 365.25 if days is not None else None
    sex = get("sex").upper()
    vital = get("vital_status").upper()
    if sex and sex not in {"MALE", "FEMALE", "OTHER"}:
        raise ValueError("Sex must be male, female, other or missing")
    if vital and vital not in {"ALIVE", "DEAD"}:
        raise ValueError("Vital status must be alive, dead or missing; OS_STATUS is not used")
    return (ids.pop(), get("cancer_type").upper().removeprefix("TCGA-") or None,
            age, stage_name(get("stage")) or None, sex or None, vital or None)


def survival_signature(raw):
    """Validate optional explicit endpoints; preserve missingness, never infer events."""
    result = []
    for endpoint in ("os", "pfs"):
        duration = number(value(str(raw.get(endpoint + "_time_months", ""))), endpoint + " time in months")
        event = number(value(str(raw.get(endpoint + "_event", ""))), endpoint + " event", True)
        if event not in (None, 0, 1):
            raise ValueError("Survival events must be 0 (censored) or 1 (event)")
        result.extend((duration, event))
    return tuple(result)


def import_clinical(db_path, input_path, provenance=""):
    """One normalized record per patient per source; conflicts abort atomically."""
    path = Path(input_path).expanduser()
    sha = file_sha(path)
    with connect(db_path) as db:
        initialize_clinical(db)
        existing = db.execute("SELECT id,report_json FROM clinical_sources WHERE sha256=?", (sha,)).fetchone()
        if existing:
            return {**json.loads(existing["report_json"]), "clinical_source_id": existing["id"], "reused": True}
        quality = {"sha256": sha, "rows_read": 0, "imported": 0, "duplicates": 0,
                   "rejected": 0, "rejection_examples": [], "reused": False}
        with table_rows(path) as (fields, rows), db:
            if not any(a in fields for a in ALIASES["patient_id"]):
                raise ValueError("Clinical table requires patient_id or another supported canonical patient ID column")
            source_id = db.execute("INSERT INTO clinical_sources(sha256,filename,provenance,imported_at,report_json) VALUES(?,?,?,?,?)",
                                   (sha, path.name, provenance, datetime.now(timezone.utc).isoformat(), "{}")).lastrowid
            for record, raw in enumerate(rows, start=2):
                quality["rows_read"] += 1
                try:
                    validate_row(raw)
                    normalized = normalize_clinical(raw, fields)
                    endpoint_signature = survival_signature(raw)
                except ValueError as error:
                    quality["rejected"] += 1
                    if len(quality["rejection_examples"]) < 20:
                        quality["rejection_examples"].append({"record": record, "reason": str(error)})
                    continue
                old = db.execute(f"SELECT {','.join(CLINICAL_FIELDS)},raw_json FROM clinical_records WHERE clinical_source_id=? AND patient_id=?", (source_id, normalized[0])).fetchone()
                if old is not None:
                    old_raw = json.loads(old["raw_json"])
                    if (tuple(old)[:-1] != normalized or survival_signature(old_raw) != endpoint_signature
                        or any(clinical_category(old_raw, name).casefold() != clinical_category(raw, name).casefold() for name in CATEGORY_ALIASES)):
                        raise ValueError(f"Conflicting clinical records for {normalized[0]}; entire import rolled back. Supply one curated row per patient.")
                    quality["duplicates"] += 1
                    continue
                db.execute("INSERT INTO clinical_records VALUES(?,?,?,?,?,?,?,?)", (source_id, *normalized, json.dumps(raw, ensure_ascii=False)))
                quality["imported"] += 1
            if not quality["imported"]:
                raise ValueError(f"No valid clinical records; import rolled back. Rejections: {quality['rejection_examples']}")
            db.execute("UPDATE clinical_sources SET report_json=? WHERE id=?", (json.dumps(quality), source_id))
        return {**quality, "clinical_source_id": source_id}


def import_crosswalk(db_path, input_path, fusion_source_id=None, provenance=""):
    """Require source-specific mappings and reject ambiguous/contradictory assignments."""
    if not Path(db_path).is_file():
        raise ValueError("Import a fusion callset before importing its crosswalk")
    path = Path(input_path).expanduser()
    sha = file_sha(path)
    with connect(db_path) as db:
        initialize_clinical(db)
        fusion_source_id = source_info(db, fusion_source_id)["id"]
        existing = db.execute("SELECT id,report_json FROM mapping_sources WHERE sha256=? AND fusion_source_id=?", (sha, fusion_source_id)).fetchone()
        if existing:
            return {**json.loads(existing["report_json"]), "mapping_source_id": existing["id"], "reused": True}
        quality = {"sha256": sha, "fusion_source_id": fusion_source_id,
                   "rows_read": 0, "imported": 0, "duplicates": 0, "reused": False}
        with table_rows(path) as (fields, rows), db:
            if not {"sample_id", "patient_id"}.issubset(fields):
                raise ValueError("Crosswalk requires sample_id and patient_id columns")
            source_id = db.execute("INSERT INTO mapping_sources(sha256,fusion_source_id,filename,provenance,imported_at,report_json) VALUES(?,?,?,?,?,?)",
                                   (sha, fusion_source_id, path.name, provenance, datetime.now(timezone.utc).isoformat(), "{}")).lastrowid
            for record, raw in enumerate(rows, start=2):
                quality["rows_read"] += 1
                validate_row(raw)
                sample = value(raw["sample_id"])
                patient = patient_id(raw["patient_id"])
                if not sample:
                    raise ValueError(f"Missing sample_id in crosswalk record {record}")
                matching = db.execute("SELECT DISTINCT cohort,patient_id FROM calls WHERE source_id=? AND sample_id=?", (fusion_source_id, sample)).fetchall()
                if not matching or any(r["cohort"] != "TCGA" for r in matching):
                    raise ValueError(f"Crosswalk sample {sample!r} must belong to TCGA in the selected fusion source")
                if any(r["patient_id"] and r["patient_id"] != patient for r in matching):
                    raise ValueError(f"Crosswalk contradicts canonical patient barcode for {sample}")
                old = db.execute("SELECT patient_id FROM sample_mappings WHERE mapping_source_id=? AND sample_id=?", (source_id, sample)).fetchone()
                if old:
                    if old["patient_id"] != patient:
                        raise ValueError(f"Ambiguous crosswalk for {sample}; entire import rolled back")
                    quality["duplicates"] += 1
                else:
                    db.execute("INSERT INTO sample_mappings VALUES(?,?,?)", (source_id, sample, patient))
                    quality["imported"] += 1
            if not quality["imported"]:
                raise ValueError("Empty crosswalk; import rolled back")
            db.execute("UPDATE mapping_sources SET report_json=? WHERE id=?", (json.dumps(quality), source_id))
        return {**quality, "mapping_source_id": source_id}


@dataclass(frozen=True)
class ClinicalFilters:
    min_age: float | None = None
    max_age: float | None = None
    stages: tuple[str, ...] = ()
    sex: str = ""
    vital_status: str = ""
    matched_only: bool = False
    races: tuple[str, ...] = ()
    histologies: tuple[str, ...] = ()
    age_unknown_only: bool = False

    def where(self):
        clauses, args = [], []
        for age in (self.min_age, self.max_age):
            if age is not None and (not math.isfinite(age) or age < 0):
                raise ValueError("Age thresholds must be finite and nonnegative")
        if self.min_age is not None and self.max_age is not None and self.min_age > self.max_age:
            raise ValueError("min-age must not exceed max-age")
        if self.age_unknown_only and (self.min_age is not None or self.max_age is not None):
            raise ValueError("Unknown age and a numeric age range cannot be selected together.")
        if self.sex and self.sex.upper() not in {"MALE", "FEMALE", "OTHER"}:
            raise ValueError("Unsupported sex filter")
        if self.vital_status and self.vital_status.upper() not in {"ALIVE", "DEAD"}:
            raise ValueError("Unsupported vital-status filter")
        active = self.matched_only or self.min_age is not None or self.max_age is not None or self.stages or self.sex or self.vital_status or self.races or self.histologies or self.age_unknown_only
        if active:
            clauses.append("match_status='matched'")
        if self.min_age is not None:
            clauses.append("age_years>=?")
            args.append(self.min_age)
        if self.max_age is not None:
            clauses.append("age_years<=?")
            args.append(self.max_age)
        if self.age_unknown_only:
            clauses.append("age_years IS NULL")
        if self.stages:
            stages = tuple(stage_name(s) for s in self.stages)
            if not all(stages):
                raise ValueError("Stage filters must be nonempty")
            clauses.append(f"stage IN ({','.join('?' for _ in stages)})")
            args.extend(stages)
        for name, val in (("sex", self.sex), ("vital_status", self.vital_status)):
            if val:
                clauses.append(f"{name}=?")
                args.append(val.upper())
        for name, values in (("race", self.races), ("histology", self.histologies)):
            if values:
                clauses.append(f"{name} COLLATE NOCASE IN ({','.join('?' for _ in values)})")
                args.extend(value(v) or "unknown" for v in values)
        return " AND ".join(clauses) or "1=1", args


def clinical_join_query(filters, source_id, clinical_source_id, mapping_source_id=None):
    """Shared parameterized join for CLI exports and the read-only Dashboard."""
    fusion_where, fusion_args = filters.where(source_id)
    sql = f"""WITH selected AS (SELECT * FROM calls WHERE {fusion_where}),
            resolved AS (
                SELECT f.*,
                  CASE WHEN f.cohort='TCGA' THEN COALESCE(f.patient_id,m.patient_id) END resolved_patient_id,
                  CASE WHEN f.cohort!='TCGA' THEN 'not_tcga'
                       WHEN f.patient_id IS NOT NULL THEN 'canonical'
                       WHEN m.patient_id IS NOT NULL THEN 'crosswalk' ELSE 'unresolved' END patient_id_origin
                FROM selected f LEFT JOIN sample_mappings m ON m.mapping_source_id=? AND m.sample_id=f.sample_id
            )
            SELECT r.*,
                c.patient_id clinical_patient_id, c.cancer_type clinical_cancer_type,
                c.age_years,c.stage,c.sex,c.vital_status,
                {category_expression('race')} race,{category_expression('histology')} histology,
                CASE WHEN r.cohort!='TCGA' THEN 'not_tcga'
                     WHEN r.resolved_patient_id IS NULL THEN 'unresolved_patient_id'
                     WHEN c.patient_id IS NULL THEN 'clinical_not_found'
                     WHEN r.cancer_type!='UNKNOWN' AND c.cancer_type IS NOT NULL
                          AND r.cancer_type!=c.cancer_type THEN 'cancer_mismatch'
                     ELSE 'matched' END match_status
            FROM resolved r LEFT JOIN clinical_records c
              ON c.clinical_source_id=? AND c.patient_id=r.resolved_patient_id
            """
    return sql, [*fusion_args, mapping_source_id, clinical_source_id]


def combine_clinical(db_path, out_dir, filters=Filters(), clinical_filters=ClinicalFilters(),
                     source_id=None, clinical_source_id=None, mapping_source_id=None):
    """Left join one clinical source; export match audit before clinical filtering.

    Disk-backed temporary tables keep memory bounded and preserve a stable join
    snapshot across all report exports. Fusion calls are never updated in place.
    """
    if not Path(db_path).is_file():
        raise ValueError("Database does not exist. Import fusion and clinical records first.")
    clinical_where, clinical_args = clinical_filters.where()
    with connect(db_path) as db:
        initialize_clinical(db)
        fusion_source = source_info(db, source_id)
        if clinical_source_id is None:
            clinical = db.execute("SELECT * FROM clinical_sources ORDER BY id DESC LIMIT 1").fetchone()
        else:
            clinical = db.execute("SELECT * FROM clinical_sources WHERE id=?", (clinical_source_id,)).fetchone()
        if clinical is None:
            raise ValueError("No matching clinical source. Run clinical-import first.")
        mapping = None
        if mapping_source_id is not None:
            mapping = db.execute("SELECT * FROM mapping_sources WHERE id=? AND fusion_source_id=?", (mapping_source_id, fusion_source["id"])).fetchone()
            if mapping is None:
                raise ValueError("Mapping source not found or belongs to another fusion source")
        query, parameters = clinical_join_query(filters, fusion_source["id"], clinical["id"], mapping_source_id)
        db.execute("CREATE TEMP TABLE combined AS " + query, parameters)
        db.execute("CREATE INDEX combined_patient ON combined(resolved_patient_id)")
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)

        def export(filename, sql, args=()):
            cursor = db.execute(sql, args)
            with (out / filename).open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(c[0] for c in cursor.description)
                writer.writerows(cursor)

        sample_query = "SELECT DISTINCT sample_id,cohort,cancer_type,resolved_patient_id,patient_id_origin,match_status FROM combined"
        export("sample_matches.csv", sample_query + " ORDER BY sample_id,cohort,cancer_type")
        audit = {r["match_status"]: r["n"] for r in db.execute(f"SELECT match_status,COUNT(*) n FROM ({sample_query}) GROUP BY match_status")}
        audit_total = sum(audit.values())
        tcga_total = audit_total - audit.get("not_tcga", 0)
        db.execute(f"CREATE TEMP TABLE filtered AS SELECT * FROM combined WHERE {clinical_where}", clinical_args)
        exported = (*COLUMNS, "resolved_patient_id", "patient_id_origin", "match_status", *CLINICAL_COLUMNS)
        exported_sql = ",".join("COALESCE(age_years,'unknown') age_years" if c == "age_years" else c for c in exported)
        export("calls_clinical.csv", f"SELECT {exported_sql} FROM filtered")
        export("patient_summary.csv", """SELECT resolved_patient_id patient_id,clinical_cancer_type,
            COALESCE(age_years,'unknown') age_years,stage,sex,vital_status,race,histology,COUNT(*) calls,COUNT(DISTINCT sample_id) samples,
            COUNT(DISTINCT fusion) distinct_fusions FROM filtered WHERE match_status='matched'
            GROUP BY resolved_patient_id,clinical_cancer_type,age_years,stage,sex,vital_status,race,histology
            ORDER BY patient_id""")
        export("clinical_summary.csv", """SELECT clinical_cancer_type,stage,sex,vital_status,race,histology,
            COUNT(DISTINCT resolved_patient_id) patients, COUNT(DISTINCT sample_id) samples,
            COUNT(*) calls FROM filtered WHERE match_status='matched'
            GROUP BY clinical_cancer_type,stage,sex,vital_status,race,histology
            ORDER BY clinical_cancer_type,stage,sex,vital_status,race,histology""")
        summary = dict(db.execute("""SELECT COUNT(*) calls,COUNT(DISTINCT sample_id) samples,
            COUNT(DISTINCT CASE WHEN match_status='matched' THEN resolved_patient_id END) matched_patients,
            COUNT(DISTINCT CASE WHEN match_status='matched' THEN sample_id END) matched_samples,
            SUM(CASE WHEN match_status='matched' THEN 1 ELSE 0 END) matched_calls FROM filtered""").fetchone())
        summary["matched_calls"] = summary["matched_calls"] or 0

        def metadata(row):
            return {k: row[k] for k in row.keys() if k != "report_json"}

        report = {
            "fusion_source": metadata(fusion_source), "clinical_source": metadata(clinical),
            "mapping_source": metadata(mapping) if mapping else None,
            "clinical_import_quality": json.loads(clinical["report_json"]),
            "mapping_import_quality": json.loads(mapping["report_json"]) if mapping else None,
            "fusion_filters": asdict(filters), "clinical_filters": asdict(clinical_filters),
            "match_audit_before_clinical_filters": {
                "sample_records": audit_total, "tcga_sample_records": tcga_total,
                "by_status": audit,
                "tcga_matched_fraction": audit.get("matched", 0) / tcga_total if tcga_total else None,
            },
            "summary_after_clinical_filters": summary,
            "notes": [
                "Sample audit counts distinct (sample_id, cohort, cancer_type) combinations after fusion filters but before clinical filters.",
                "Matched fraction is among selected TCGA sample records; it is not fusion prevalence.",
                "Clinical duplicates cannot multiply fusion rows: one record per patient per selected clinical source is enforced.",
                "No shorthand IDs are guessed. Crosswalk use is explicit and scoped to the selected fusion source.",
                "Unmatched and cancer-mismatch rows are retained by default; any clinical filter requires a valid match.",
                "Cancer-mismatch rows retain both cancer labels for audit but are excluded from patient/clinical summaries.",
                "Age_at_diagnosis is interpreted as days; age_years and AGE as years. No survival endpoints are derived.",
                "Patient matching does not establish same-sample or same-section tissue matching.",
            ],
        }
        (out / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return report
