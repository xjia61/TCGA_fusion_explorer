"""Prepare a local cBioPortal PRAD PanCancer Atlas snapshot for validation.

Run: python3 -m fusion_explorer.prepare_prad --raw-dir data/prad/raw --out data/prad/prepared
Downloads are deliberately separate. Original files remain unmodified.
"""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from .clinical import file_sha, patient_id, table_rows, validate_row, value


STUDY = "prad_tcga_pan_can_atlas_2018"
BASE = f"https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/{STUDY}/"
FILES = ("data_clinical_patient.txt", "data_clinical_sample.txt", "data_sv.txt")


def read_table(path, required):
    with table_rows(path) as (fields, rows):
        if not set(required).issubset(fields):
            raise ValueError(f"{path.name}: missing columns {set(required) - set(fields)}")
        for row in rows:
            validate_row(row)
            yield row


def write_table(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def prepare(raw_dir, out_dir):
    raw, out = Path(raw_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    patients = {}
    for row in read_table(raw / FILES[0], ("PATIENT_ID", "CANCER_TYPE_ACRONYM", "AGE", "SEX", "OS_STATUS")):
        pid = patient_id(row["PATIENT_ID"])
        if row["CANCER_TYPE_ACRONYM"] != "PRAD" or pid in patients:
            raise ValueError("Expected unique PRAD patients in clinical source")
        statuses = {"0:LIVING": "alive", "1:DECEASED": "dead", "": ""}
        status = value(row["OS_STATUS"])
        if status not in statuses:
            raise ValueError(f"Unexpected OS_STATUS {status!r}")
        pfs_status = value(row.get("PFS_STATUS", ""))
        pfs_events = {"0:CENSORED": "0", "1:PROGRESSION": "1", "": ""}
        if pfs_status not in pfs_events:
            raise ValueError(f"Unexpected PFS_STATUS {pfs_status!r}")
        patients[pid] = {
            "patient_id": pid, "cancer_type": "PRAD", "age_years": row["AGE"],
            "sex": row["SEX"], "stage": row.get("AJCC_PATHOLOGIC_TUMOR_STAGE", ""),
            "vital_status": statuses[status], "path_t_stage": row.get("PATH_T_STAGE", ""),
            "source_os_status": row["OS_STATUS"], "source_os_months": row.get("OS_MONTHS", ""),
            "os_time_months": row.get("OS_MONTHS", ""),
            "os_event": {"0:LIVING": "0", "1:DECEASED": "1", "": ""}[status],
            "pfs_time_months": row.get("PFS_MONTHS", ""), "pfs_event": pfs_events[pfs_status],
            "source_pfs_status": row.get("PFS_STATUS", ""),
            "race": value(row.get("RACE", "")) or "unknown",
            "histology": ("ICD-O-3: " + value(row.get("ICD_O_3_HISTOLOGY", ""))) if value(row.get("ICD_O_3_HISTOLOGY", "")) else "unknown",
        }
    samples = {}
    for row in read_table(raw / FILES[1], ("PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE")):
        pid, sample = patient_id(row["PATIENT_ID"]), row["SAMPLE_ID"]
        if pid not in patients or sample in samples or sample[:12] != pid:
            raise ValueError("Sample table has an unknown patient, duplicate sample or barcode disagreement")
        samples[sample] = pid

    fields = ["sample_id", "fusion", "cohort", "cancer_type", "tumor_or_normal",
              "left_gene", "right_gene", "junction_reads", "spanning_reads", "ffpm",
              "left_breakpoint", "right_breakpoint", "fusion_type", "reference_build", "source_event_info"]
    call_count, reporting_samples = 0, set()

    def calls():
        nonlocal call_count
        required = ("Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol", "Site1_Chromosome",
                    "Site1_Position", "Site2_Chromosome", "Site2_Position", "Event_Info",
                    "Tumor_Split_Read_Count", "Tumor_Paired_End_Read_Count", "NCBI_Build")
        for row in read_table(raw / FILES[2], required):
            sample = row["Sample_Id"]
            if sample not in samples or not row["Event_Info"].endswith(" fusion"):
                raise ValueError("SV row lacks a verified sample or is not labeled a fusion event")
            left, right = value(row["Site1_Hugo_Symbol"]), value(row["Site2_Hugo_Symbol"])
            if not left or not right:
                raise ValueError("Missing fusion partners")
            call_count += 1
            reporting_samples.add(sample)
            yield {
                "sample_id": sample, "fusion": f"{left}::{right}", "cohort": "TCGA", "cancer_type": "PRAD",
                "tumor_or_normal": "tumor", "left_gene": left, "right_gene": right,
                "junction_reads": row["Tumor_Split_Read_Count"],
                "spanning_reads": row["Tumor_Paired_End_Read_Count"], "ffpm": "",
                "left_breakpoint": f"chr{row['Site1_Chromosome']}:{row['Site1_Position']}",
                "right_breakpoint": f"chr{row['Site2_Chromosome']}:{row['Site2_Position']}",
                "fusion_type": row.get("Site2_Effect_On_Frame", ""),
                "reference_build": row["NCBI_Build"], "source_event_info": row["Event_Info"],
            }

    # Validate the sample types before labeling this study's calls as tumor.
    with table_rows(raw / FILES[1]) as (_, rows):
        if any(r["SAMPLE_TYPE"] != "Primary" for r in rows):
            raise ValueError("This adapter currently expects primary tumor samples only")
    write_table(out / "fusion_calls.csv", fields, calls())
    clinical_fields = ["patient_id", "cancer_type", "age_years", "stage", "sex", "vital_status",
                       "path_t_stage", "source_os_status", "source_os_months",
                       "os_time_months", "os_event", "pfs_time_months", "pfs_event", "source_pfs_status", "race", "histology"]
    write_table(out / "clinical.csv", clinical_fields, patients.values())
    write_table(out / "sample_patient.csv", ["sample_id", "patient_id"],
                ({"sample_id": sample, "patient_id": pid} for sample, pid in sorted(samples.items())))
    manifest = {
        "study": STUDY, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "raw_sources": [{"filename": name, "url": BASE + name, "sha256": file_sha(raw / name)} for name in FILES],
        "prepared_files": {name: file_sha(out / name) for name in ("fusion_calls.csv", "clinical.csv", "sample_patient.csv")},
        "counts": {"clinical_patients": len(patients), "clinical_samples": len(samples),
                   "fusion_rows": call_count, "samples_with_reported_fusions": len(reporting_samples)},
        "rules": [
            "Same cBioPortal PRAD PanCancer Atlas study for clinical, sample relationships and fusion-labeled SV records.",
            "Site1 and Site2 gene order is preserved; 5-prime/3-prime transcript orientation is not independently verified.",
            "Tumor_Split_Read_Count -> junction_reads; Tumor_Paired_End_Read_Count -> spanning_reads. These are source-specific evidence counts.",
            "FFPM is unavailable and stays missing; positive FFPM filters will exclude these calls.",
            "Explicit OS mapping: 0:LIVING -> event=0, 1:DECEASED -> event=1; OS_MONTHS retained in months.",
            "Explicit source PFS mapping: 0:CENSORED -> event=0, 1:PROGRESSION -> event=1; PFS_MONTHS retained in months. Not relabeled PFI.",
            "PATH_T_STAGE is preserved separately, never substituted for overall AJCC stage.",
            "Samples without a reported fusion are not classified as confirmed negative or necessarily fusion-profiled.",
            "master URLs are mutable; SHA-256 records identify the downloaded snapshot bytes.",
        ],
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.raw_dir, args.out), indent=2))


if __name__ == "__main__":
    main()
