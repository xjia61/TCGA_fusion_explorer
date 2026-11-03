"""Normalize and import public cBioPortal SKCM clinical records, offline."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from .clinical import file_sha, import_clinical, patient_id, survival_signature, value
from .prepare_prad import read_table, write_table


SOURCE_URL = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/skcm_tcga_pan_can_atlas_2018/data_clinical_patient.txt"


def prepare(raw_path, out_dir, db_path):
    raw_path, out = Path(raw_path), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    records, seen = [], set()
    required = ("PATIENT_ID", "CANCER_TYPE_ACRONYM", "AGE", "SEX", "RACE", "OS_STATUS", "OS_MONTHS", "PFS_STATUS", "PFS_MONTHS")
    for row in read_table(raw_path, required):
        pid = patient_id(row["PATIENT_ID"])
        if pid in seen or row["CANCER_TYPE_ACRONYM"] != "SKCM":
            raise ValueError("Expected unique canonical SKCM patients")
        seen.add(pid)
        os_status, pfs_status = value(row["OS_STATUS"]), value(row["PFS_STATUS"])
        os_events = {"0:LIVING": "0", "1:DECEASED": "1", "": ""}
        pfs_events = {"0:CENSORED": "0", "1:PROGRESSION": "1", "": ""}
        if os_status not in os_events or pfs_status not in pfs_events:
            raise ValueError("Unrecognized source OS/PFS status; do not guess event coding")
        histology = value(row.get("ICD_O_3_HISTOLOGY", ""))
        record = {**row, "patient_id": pid, "cancer_type": "SKCM", "age_years": value(row["AGE"]),
                  "sex": value(row["SEX"]), "race": value(row["RACE"]) or "unknown",
                  "histology": "ICD-O-3: " + histology if histology else "unknown",
                  "stage": value(row.get("AJCC_PATHOLOGIC_TUMOR_STAGE", "")),
                  "vital_status": {"0:LIVING": "alive", "1:DECEASED": "dead", "": ""}[os_status],
                  "os_time_months": value(row["OS_MONTHS"]), "os_event": os_events[os_status],
                  "pfs_time_months": value(row["PFS_MONTHS"]), "pfs_event": pfs_events[pfs_status]}
        survival_signature(record)
        records.append(record)
    if not records:
        raise ValueError("No clinical patients")
    path = out/"clinical.csv"
    normalized = ["patient_id", "cancer_type", "age_years", "sex", "race", "histology", "stage", "vital_status",
                  "os_time_months", "os_event", "pfs_time_months", "pfs_event"]
    write_table(path, normalized+[f for f in records[0] if f not in normalized], records)
    manifest = {"study": "skcm_tcga_pan_can_atlas_2018", "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
                "source_url": SOURCE_URL, "source_sha256": file_sha(raw_path), "prepared_sha256": file_sha(path),
                "clinical_patients": len(records),
                "notes": ["Original fields are retained alongside normalized columns.",
                          "AGE is diagnosis age, not inferred disease-onset age. Race/histology missing values remain unknown.",
                          "ICD-O-3 source histology codes are retained; stage subgroups are preserved.",
                          "Explicit event mappings: OS 0:LIVING=0,1:DECEASED=1; PFS 0:CENSORED=0,1:PROGRESSION=1. Time unit: months.",
                          "Incomplete endpoint pairs remain incomplete and are excluded with counts by the survival module.",
                          "SKCM specimens include metastases. Diagnosis-based survival in a specimen-selected cohort can have selection/delayed-entry bias; sample timing is not adjusted."]}
    (out/"clinical_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    provenance = f"{SOURCE_URL}; raw SHA256={manifest['source_sha256']}; explicit OS/PFS normalization; see clinical_manifest.json"
    quality = import_clinical(db_path, path, provenance)
    if quality["rejected"]:
        raise ValueError(f"Clinical records were rejected; review before use: {quality}")
    return {"manifest": manifest, "import": quality}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/skcm/raw/data_clinical_patient.txt"))
    parser.add_argument("--out", type=Path, default=Path("data/skcm/prepared"))
    parser.add_argument("--db", type=Path, default=Path("data/skcm/validation.sqlite"))
    args = parser.parse_args()
    print(json.dumps(prepare(args.input,args.out,args.db),indent=2))
