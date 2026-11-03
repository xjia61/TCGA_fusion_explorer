"""Offline, metadata-only SKCM fusion-to-WSI linkage; never downloads slides."""

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from .clinical import file_sha, patient_id
from .prepare_prad import read_table, write_table


SAMPLE_PREFIX = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2})(?:[A-Z]|-|$)", re.I)
LINK_COLUMNS = ["fusion_sample_id", "patient_id", "fusion_sample_type", "fusion_call_rows",
                "file_id", "file_name", "slide_submitter_id", "slide_sample_prefix",
                "slide_sample_type_hint", "sample_type_consistency", "match_level",
                "same_specimen_verified", "slide_kind_hint", "access", "file_size_bytes",
                "case_id", "case_viewer_entry", "file_details_url", "download_url"]


def build_manifest(fusion_path, samples_path, slides_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    samples = {}
    for row in read_table(Path(samples_path), ("PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE")):
        pid = patient_id(row["PATIENT_ID"])
        sid = row["SAMPLE_ID"]
        if sid in samples or not sid.startswith(pid + "-") or not SAMPLE_PREFIX.match(sid):
            raise ValueError("Duplicate or noncanonical sample relationship")
        samples[sid] = {"patient_id": pid, "sample_type": row["SAMPLE_TYPE"] or "unknown"}
    fusion_counts = Counter()
    for row in read_table(Path(fusion_path), ("Sample_Id", "Event_Info")):
        if not row["Event_Info"].strip().lower().endswith(" fusion"):
            raise ValueError("Input contains an event not labeled fusion")
        if row["Sample_Id"] not in samples:
            raise ValueError("Fusion sample missing from authoritative sample table")
        fusion_counts[row["Sample_Id"]] += 1
    if not fusion_counts:
        raise ValueError("No fusion records")
    response = json.loads(Path(slides_path).read_text())
    if response.get("warnings"):
        raise ValueError("GDC response contains warnings; review request before matching")
    hits = response["data"]["hits"]
    if response["data"]["pagination"]["total"] != len(hits):
        raise ValueError("Incomplete GDC pagination; retrieve all metadata pages first")
    inventory, by_patient, seen = [], defaultdict(list), set()
    type_hints = {"01": "Primary", "06": "Metastasis", "07": "Additional Metastasis",
                  "02": "Recurrent", "11": "Solid Tissue Normal"}
    for hit in hits:
        fid = hit["file_id"]
        if fid in seen:
            raise ValueError("Duplicate GDC file ID")
        seen.add(fid)
        cases = hit.get("cases", [])
        if len(cases) != 1:
            raise ValueError("Expected one case per slide file; manual review required")
        case = cases[0]
        pid = patient_id(case["submitter_id"])
        entities = [e for e in hit.get("associated_entities", []) if e.get("entity_type") == "slide"]
        slide_id = entities[0].get("entity_submitter_id", "") if len(entities) == 1 else ""
        prefix = SAMPLE_PREFIX.match(slide_id)
        slide_sample = prefix.group(1) if prefix else "unknown"
        if prefix and not slide_sample.startswith(pid + "-"):
            raise ValueError("Slide entity barcode conflicts with the GDC case")
        sample_type = type_hints.get(slide_sample[-2:], "unknown") if prefix else "unknown"
        kind = "diagnostic" if re.search(r"-DX\d+$", slide_id) else "tissue_section" if re.search(r"-TS[A-Z\d]+$", slide_id) else "unknown"
        record = {"patient_id": pid, "file_id": fid, "file_name": hit["file_name"],
                  "slide_submitter_id": slide_id or "unknown", "slide_sample_prefix": slide_sample,
                  "slide_sample_type_hint": sample_type, "slide_kind_hint": kind,
                  "access": hit["access"], "file_size_bytes": hit.get("file_size", ""),
                  "case_id": case["case_id"], "case_viewer_entry": "https://portal.gdc.cancer.gov/cases/" + case["case_id"],
                  "file_details_url": "https://portal.gdc.cancer.gov/files/" + fid,
                  "download_url": "https://api.gdc.cancer.gov/data/" + fid if hit["access"] == "open" else ""}
        inventory.append(record)
        by_patient[pid].append(record)
    links, unmatched = [], []
    for sid, count in sorted(fusion_counts.items()):
        sample = samples[sid]
        pid = sample["patient_id"]
        if pid not in by_patient:
            unmatched.append({"fusion_sample_id": sid, "patient_id": pid,
                              "fusion_sample_type": sample["sample_type"], "reason": "no_slide_in_GDC_snapshot"})
        for slide in by_patient.get(pid, []):
            hint = slide["slide_sample_type_hint"]
            consistent = "unknown" if hint == "unknown" or sample["sample_type"] == "unknown" else "same_type_hint" if hint == sample["sample_type"] else "different_type_hint"
            links.append({**slide, "fusion_sample_id": sid, "fusion_sample_type": sample["sample_type"],
                          "fusion_call_rows": count, "sample_type_consistency": consistent,
                          "match_level": "patient_only_verified", "same_specimen_verified": "no"})
    write_table(out / "fusion_slide_links.csv", LINK_COLUMNS, links)
    write_table(out / "slide_inventory.csv", list(inventory[0]) if inventory else ["file_id"], inventory)
    write_table(out / "unmatched_samples.csv", ["fusion_sample_id", "patient_id", "fusion_sample_type", "reason"], unmatched)
    # Three case entries for manual preview; select one source file per patient.
    previews, preview_patients = [], set()
    for row in sorted(links, key=lambda r: (r["sample_type_consistency"] != "same_type_hint", r["patient_id"], r["file_id"])):
        if row["patient_id"] not in preview_patients and row["access"] == "open":
            previews.append(row)
            preview_patients.add(row["patient_id"])
        if len(previews) == 3:
            break
    write_table(out / "preview_cases.csv", LINK_COLUMNS, previews)
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "fusion_study": "skcm_tcga_pan_can_atlas_2018", "pathology_project": "TCGA-SKCM",
        "sources": {
            "fusion": "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/skcm_tcga_pan_can_atlas_2018/data_sv.txt",
            "samples": "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/skcm_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
            "slides_endpoint": "https://api.gdc.cancer.gov/files",
            "slides_filters": {"cases.project.project_id": "TCGA-SKCM", "data_type": "Slide Image"},
        },
        "inputs": [{"filename": Path(p).name, "sha256": file_sha(p)} for p in (fusion_path, samples_path, slides_path)],
        "counts": {"fusion_call_rows": sum(fusion_counts.values()), "fusion_samples": len(fusion_counts),
                   "fusion_patients": len({samples[s]["patient_id"] for s in fusion_counts}),
                   "GDC_slide_files": len(hits), "GDC_slide_patients": len(by_patient),
                   "matched_fusion_samples": len(fusion_counts)-len(unmatched),
                   "unmatched_fusion_samples": len(unmatched), "patient_slide_links": len(links),
                   "unique_linked_slide_files": len({r["file_id"] for r in links}),
                   "matched_patients": len({r["patient_id"] for r in links}),
                   "same_specimen_verified_links": 0},
        "link_type_hints": dict(Counter(r["sample_type_consistency"] for r in links)),
        "fusion_sample_types": dict(Counter(samples[s]["sample_type"] for s in fusion_counts)),
        "WSI_bytes_downloaded": 0,
        "notes": [
            "Fusion sample-to-patient relationships come from the same cBioPortal study's clinical sample table.",
            "Slide file-to-case and slide entity IDs come from GDC file metadata. Only exact canonical patient IDs are joined.",
            "All links are patient-level only. Shared sample-type hints never establish identical specimen, tissue block or section.",
            "Slide sample-type hints are derived from TCGA sample codes; unknown values stay unknown.",
            "Case links are entry points to review available images, not guaranteed direct links to a particular slide viewer.",
            "No slide image, thumbnail or pathology-report PDF was downloaded or visually inspected.",
            "This canonical cBioPortal callset differs from the earlier 465-sample abbreviated STAR-Fusion subset; they are not pooled.",
        ],
    }
    (out / "linkage_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fusion", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--slides", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_manifest(args.fusion, args.samples, args.slides, args.out), indent=2))


if __name__ == "__main__":
    main()
