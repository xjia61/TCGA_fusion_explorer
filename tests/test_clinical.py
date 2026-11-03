import csv
import gzip
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

from fusion_explorer.clinical import (
    ClinicalFilters, combine_clinical, import_clinical, import_crosswalk,
)
from fusion_explorer.core import Filters, import_callset


class ClinicalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.db = self.root / "test.sqlite"
        self.fusion = self.root / "fusions.csv"
        self.clinical = self.root / "clinical.csv"
        self.crosswalk = self.root / "mapping.csv"
        self.write(self.fusion, ["sample_id", "fusion", "cohort", "cancer_type"], [
            ["TCGA-ZZ-0001-01A", "A::B", "TCGA", "BRCA"],
            ["TCGA-ZZ-0001-01A", "C::D", "TCGA", "BRCA"],
            ["TCGA-ZZ-0001-11A", "A::B", "TCGA", "BRCA"],
            ["OV-ZZ-0002-TP", "A::B", "TCGA", "OV"],
            ["TCGA-ZZ-0003-01A", "A::B", "TCGA", "LUAD"],
            ["TCGA-ZZ-0004-01A", "A::B", "TCGA", "PRAD"],
            ["GTEX-TEST", "A::B", "GTEx", "Breast"],
        ])
        self.source = import_callset(self.db, self.fusion)["source_id"]
        self.fields = ["patient_id", "cancer_type", "age_years", "stage", "sex", "vital_status"]
        self.clinical_rows = [
            ["TCGA-ZZ-0001", "BRCA", "50", "Stage IIA", "female", "alive"],
            ["TCGA-ZZ-0002", "OV", "60", "III", "female", "dead"],
            ["TCGA-ZZ-0004", "BRCA", "70", "IV", "male", "alive"],
        ]
        self.write(self.clinical, self.fields, self.clinical_rows)
        self.clinical_id = import_clinical(self.db, self.clinical, "synthetic test data")["clinical_source_id"]
        self.write(self.crosswalk, ["sample_id", "patient_id"], [["OV-ZZ-0002-TP", "TCGA-ZZ-0002"]])

    def tearDown(self):
        self.temp.cleanup()

    def write(self, path, fields, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            writer.writerows(rows)

    def read(self, filename, out="out"):
        with (self.root / out / filename).open() as handle:
            return list(csv.DictReader(handle))

    def combine(self, **kwargs):
        return combine_clinical(self.db, self.root / "out", **kwargs)

    def mapping(self):
        return import_crosswalk(self.db, self.crosswalk, self.source, "synthetic verified fixture")["mapping_source_id"]

    def test_left_join_audit_and_no_row_multiplication(self):
        report = self.combine(filters=Filters(cohort="ALL"))
        audit = report["match_audit_before_clinical_filters"]
        self.assertEqual(audit["sample_records"], 6)
        self.assertEqual(audit["by_status"], {"matched": 2, "unresolved_patient_id": 1,
                                           "clinical_not_found": 1, "cancer_mismatch": 1, "not_tcga": 1})
        self.assertEqual(audit["tcga_matched_fraction"], 2 / 5)
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 7)
        self.assertEqual(report["summary_after_clinical_filters"]["matched_patients"], 1)
        patients = self.read("patient_summary.csv")
        self.assertEqual(len(patients), 1)
        self.assertEqual((patients[0]["samples"], patients[0]["calls"], patients[0]["distinct_fusions"]), ("2", "3", "2"))
        self.assertEqual(len(self.read("calls_clinical.csv")), 7)

    def test_explicit_crosswalk_resolves_without_mutating_fusions(self):
        mapping = self.mapping()
        report = self.combine(mapping_source_id=mapping)
        self.assertEqual(report["summary_after_clinical_filters"]["matched_patients"], 2)
        sample = next(r for r in self.read("sample_matches.csv") if r["sample_id"] == "OV-ZZ-0002-TP")
        self.assertEqual((sample["patient_id_origin"], sample["resolved_patient_id"], sample["match_status"]), ("crosswalk", "TCGA-ZZ-0002", "matched"))
        with sqlite3.connect(self.db) as db:
            self.assertIsNone(db.execute("SELECT patient_id FROM calls WHERE sample_id='OV-ZZ-0002-TP'").fetchone()[0])
        # A previous mapping import never opts subsequent analyses in implicitly.
        self.assertEqual(self.combine()["summary_after_clinical_filters"]["matched_patients"], 1)

    def test_clinical_filters_use_years_exact_stages_and_valid_matches(self):
        report = self.combine(mapping_source_id=self.mapping(), clinical_filters=ClinicalFilters(
            min_age=55, max_age=65, stages=("stage III",), sex="female", vital_status="dead"))
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 1)
        self.assertEqual(self.read("patient_summary.csv")[0]["patient_id"], "TCGA-ZZ-0002")
        # Match audit is retained before narrowing by clinical fields.
        self.assertEqual(report["match_audit_before_clinical_filters"]["sample_records"], 5)

    def test_matched_only_excludes_cancer_mismatch(self):
        report = self.combine(clinical_filters=ClinicalFilters(matched_only=True))
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 3)
        self.assertTrue(all(r["match_status"] == "matched" for r in self.read("calls_clinical.csv")))

    def test_missing_age_is_not_zero_and_empty_filter_results_have_headers(self):
        self.write(self.clinical, self.fields, [["TCGA-ZZ-0001", "BRCA", "", "", "", ""]])
        import_clinical(self.db, self.clinical)
        self.assertEqual(self.combine()["summary_after_clinical_filters"]["matched_calls"], 3)
        report = self.combine(clinical_filters=ClinicalFilters(min_age=0))
        self.assertEqual(report["summary_after_clinical_filters"]["matched_calls"], 0)
        self.assertEqual(self.read("calls_clinical.csv"), [])
        self.assertTrue((self.root / "out" / "patient_summary.csv").read_text().startswith("patient_id,"))

    def test_days_to_years_and_gzip_gdc_flattened_aliases(self):
        path = self.root / "gdc.tsv.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("cases.submitter_id\tcases.project.project_id\tdiagnoses.age_at_diagnosis\tdiagnoses.ajcc_pathologic_stage\tdemographic.gender\tdemographic.vital_status\nTCGA-ZZ-0001\tTCGA-BRCA\t18262.5\tStage IIA\tfemale\talive\n")
        imported = import_clinical(self.db, path)
        self.combine(clinical_source_id=imported["clinical_source_id"])
        self.assertEqual(self.read("patient_summary.csv")[0]["age_years"], "50.0")
        self.assertEqual(self.read("patient_summary.csv")[0]["stage"], "IIA")

    def test_clinical_idempotence_and_duplicate_rows(self):
        self.assertTrue(import_clinical(self.db, self.clinical)["reused"])
        self.write(self.clinical, self.fields, self.clinical_rows * 2)
        result = import_clinical(self.db, self.clinical)
        self.assertEqual((result["imported"], result["duplicates"]), (3, 3))
        self.assertEqual(self.combine()["summary_after_clinical_filters"]["calls"], 6)

    def test_conflicting_clinical_patient_aborts_entire_source(self):
        conflicting = self.clinical_rows[0][:]
        conflicting[2] = "99"
        self.write(self.clinical, self.fields, self.clinical_rows + [conflicting])
        with self.assertRaisesRegex(ValueError, "Conflicting clinical"):
            import_clinical(self.db, self.clinical)
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM clinical_sources").fetchone()[0], 1)
            self.assertEqual(db.execute("SELECT COUNT(*) FROM clinical_records").fetchone()[0], 3)

    def test_invalid_clinical_rows_are_reported(self):
        bad = [r[:] for r in self.clinical_rows]
        bad[0][0] = "not-a-patient"
        bad[1][2] = "inf"
        bad[2][5] = "1:DECEASED"
        self.write(self.clinical, self.fields, self.clinical_rows[:1] + bad)
        result = import_clinical(self.db, self.clinical)
        self.assertEqual((result["imported"], result["rejected"]), (1, 3))

    def test_all_invalid_clinical_rolls_back(self):
        self.write(self.clinical, self.fields, [["TCGA-ZZ-0001-01A", "BRCA", "50", "II", "female", "alive"]])
        with self.assertRaisesRegex(ValueError, "No valid clinical"):
            import_clinical(self.db, self.clinical)
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM clinical_sources").fetchone()[0], 1)

    def test_ambiguous_crosswalk_aborts(self):
        self.write(self.crosswalk, ["sample_id", "patient_id"], [
            ["OV-ZZ-0002-TP", "TCGA-ZZ-0002"], ["OV-ZZ-0002-TP", "TCGA-ZZ-9999"]])
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            self.mapping()
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM sample_mappings").fetchone()[0], 0)
            self.assertEqual(db.execute("SELECT COUNT(*) FROM mapping_sources").fetchone()[0], 0)

    def test_canonical_override_gtex_and_unknown_crosswalk_rejected(self):
        for sample, message in (("TCGA-ZZ-0001-01A", "contradicts"),
                                ("GTEX-TEST", "must belong"), ("NONEXISTENT", "must belong")):
            with self.subTest(sample=sample):
                self.write(self.crosswalk, ["sample_id", "patient_id"], [[sample, "TCGA-ZZ-9999"]])
                with self.assertRaisesRegex(ValueError, message):
                    self.mapping()

    def test_crosswalk_idempotence_and_scope(self):
        mapping = self.mapping()
        self.assertTrue(import_crosswalk(self.db, self.crosswalk, self.source)["reused"])
        self.write(self.fusion, ["sample_id", "fusion", "cohort", "cancer_type"], [["OV-ZZ-0002-TP", "X::Y", "TCGA", "OV"]])
        new_source = import_callset(self.db, self.fusion)["source_id"]
        with self.assertRaisesRegex(ValueError, "another fusion source"):
            self.combine(source_id=new_source, mapping_source_id=mapping)

    def test_invalid_age_filter_rejected(self):
        for filters in (ClinicalFilters(min_age=-1), ClinicalFilters(max_age=float("nan")),
                        ClinicalFilters(min_age=80, max_age=20)):
            with self.assertRaises(ValueError):
                self.combine(clinical_filters=filters)

    def test_clinical_sources_can_be_selected_explicitly(self):
        self.write(self.clinical, self.fields, [["TCGA-ZZ-0002", "OV", "60", "III", "female", "dead"]])
        import_clinical(self.db, self.clinical)
        self.assertEqual(self.combine()["summary_after_clinical_filters"]["matched_patients"], 0)
        self.assertEqual(self.combine(clinical_source_id=self.clinical_id)["summary_after_clinical_filters"]["matched_patients"], 1)

    def test_sql_parameters_do_not_interpret_filter_content(self):
        report = self.combine(clinical_filters=ClinicalFilters(stages=("IIA' OR 1=1 --",)))
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 0)

    def test_race_histology_unknown_and_case_insensitive_filters(self):
        self.write(self.clinical, self.fields + ["RACE", "histological_type"], [
            self.clinical_rows[0] + ["White", "Adenocarcinoma"],
            self.clinical_rows[1] + ["", "Not Reported"],
            self.clinical_rows[2] + ["Asian", "Other"],
        ])
        import_clinical(self.db, self.clinical)
        report = self.combine(clinical_filters=ClinicalFilters(races=("white",), histologies=("adenocarcinoma",)))
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 3)
        self.assertEqual(self.read("calls_clinical.csv")[0]["race"], "White")
        report = self.combine(mapping_source_id=self.mapping(), clinical_filters=ClinicalFilters(races=("unknown",), histologies=("unknown",)))
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 1)
        self.assertEqual(self.read("patient_summary.csv")[0]["histology"], "unknown")

    def test_age_unknown_filter_export_and_conflict(self):
        row = self.clinical_rows[0][:]
        row[2] = ""
        self.write(self.clinical, self.fields, [row])
        import_clinical(self.db, self.clinical)
        report = self.combine(clinical_filters=ClinicalFilters(age_unknown_only=True))
        self.assertEqual(report["summary_after_clinical_filters"]["calls"], 3)
        self.assertEqual(self.read("calls_clinical.csv")[0]["age_years"], "unknown")
        with self.assertRaises(ValueError):
            ClinicalFilters(min_age=0, age_unknown_only=True).where()

    def test_conflicting_race_records_abort(self):
        self.write(self.clinical, self.fields + ["race"], [
            self.clinical_rows[0] + ["White"], self.clinical_rows[0] + ["Asian"]])
        with self.assertRaisesRegex(ValueError, "Conflicting clinical"):
            import_clinical(self.db, self.clinical)

    def test_cli_import_mapping_and_combine(self):
        def run(*args):
            result = subprocess.run([sys.executable, "-m", "fusion_explorer", *args, "--db", str(self.db)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            return json.loads(result.stdout)

        imported = run("clinical-import", "--input", str(self.clinical))
        mapping = run("mapping-import", "--input", str(self.crosswalk), "--provenance", "synthetic fixture")
        result = run("clinical-combine", "--clinical-source-id", str(imported["clinical_source_id"]),
                     "--mapping-source-id", str(mapping["mapping_source_id"]), "--matched-only", "--out", str(self.root / "out"))
        self.assertEqual(result["summary"]["matched_patients"], 2)


if __name__ == "__main__":
    unittest.main()
