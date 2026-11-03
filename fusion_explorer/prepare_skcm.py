"""Prepare the canonical SKCM fusion/pathology dataset from existing local files."""

import argparse
from pathlib import Path

from .clinical import patient_id
from .core import import_callset
from .prepare_prad import read_table, write_table
from .pathology import import_pathology
from .pathology_manifest import build_manifest


def prepare(raw_dir, out_dir, db_path):
    raw, out = Path(raw_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    samples = {}
    for row in read_table(raw / "data_clinical_sample.txt", ("SAMPLE_ID", "PATIENT_ID", "SAMPLE_TYPE")):
        pid = patient_id(row["PATIENT_ID"])
        if not row["SAMPLE_ID"].startswith(pid+"-") or row["SAMPLE_ID"] in samples:
            raise ValueError("Inconsistent SKCM sample mapping")
        samples[row["SAMPLE_ID"]] = row
    fields = ["sample_id", "fusion", "cohort", "cancer_type", "tumor_or_normal", "junction_reads",
              "spanning_reads", "left_breakpoint", "right_breakpoint", "fusion_type", "reference_build", "source_sample_type"]

    def calls():
        for row in read_table(raw / "data_sv.txt", ("Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol", "Event_Info")):
            sample = samples.get(row["Sample_Id"])
            if sample is None or sample["SAMPLE_TYPE"] not in {"Primary", "Metastasis"} or not row["Event_Info"].lower().endswith(" fusion"):
                raise ValueError("Unsupported SKCM fusion/sample type")
            yield {"sample_id": row["Sample_Id"], "fusion": row["Site1_Hugo_Symbol"]+"::"+row["Site2_Hugo_Symbol"],
                   "cohort": "TCGA", "cancer_type": "SKCM", "tumor_or_normal": "tumor",
                   "junction_reads": row.get("Tumor_Split_Read_Count", ""), "spanning_reads": row.get("Tumor_Paired_End_Read_Count", ""),
                   "left_breakpoint": f"chr{row['Site1_Chromosome']}:{row['Site1_Position']}",
                   "right_breakpoint": f"chr{row['Site2_Chromosome']}:{row['Site2_Position']}",
                   "fusion_type": row.get("Site2_Effect_On_Frame", ""), "reference_build": row.get("NCBI_Build", ""),
                   "source_sample_type": sample["SAMPLE_TYPE"]}
    write_table(out / "fusion_calls.csv", fields, calls())
    manifest_dir = out / "pathology"
    build_manifest(raw/"data_sv.txt", raw/"data_clinical_sample.txt", raw/"gdc_slide_images.json", manifest_dir)
    quality = import_callset(db_path, out/"fusion_calls.csv", "Public cBioPortal SKCM PanCancer Atlas 2018; separate from STAR-Fusion; source hashes in pathology linkage report")
    if quality["rejected"]:
        raise ValueError("SKCM normalized calls were rejected; review import quality before linking pathology")
    linked = import_pathology(db_path, manifest_dir/"fusion_slide_links.csv", manifest_dir/"linkage_report.json", quality["source_id"])
    return {"fusion": quality, "pathology": linked}


if __name__ == "__main__":
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/skcm/raw"))
    parser.add_argument("--out", type=Path, default=Path("data/skcm/prepared"))
    parser.add_argument("--db", type=Path, default=Path("data/skcm/validation.sqlite"))
    args = parser.parse_args()
    print(json.dumps(prepare(args.raw_dir, args.out, args.db), indent=2))
