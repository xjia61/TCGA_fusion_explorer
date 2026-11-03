import csv
from pathlib import Path
import tempfile
import unittest

from fusion_explorer.prepare_prad import prepare
from fusion_explorer.core import import_callset
from fusion_explorer.clinical import import_clinical, combine_clinical


class PradPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.patient_fields = ["PATIENT_ID", "CANCER_TYPE_ACRONYM", "AGE", "SEX", "OS_STATUS", "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_T_STAGE"]
        self.patient_row = ["TCGA-ZZ-0001", "PRAD", "60", "Male", "0:LIVING", "", "T3A"]
        self.write("data_clinical_patient.txt", self.patient_fields, [self.patient_row])
        self.write("data_clinical_sample.txt", ["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE"], [["TCGA-ZZ-0001", "TCGA-ZZ-0001-01", "Primary"]])
        self.write("data_sv.txt", ["Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol", "Site1_Chromosome", "Site1_Position", "Site2_Chromosome", "Site2_Position", "Event_Info", "Tumor_Split_Read_Count", "Tumor_Paired_End_Read_Count", "NCBI_Build"],
                   [["TCGA-ZZ-0001-01", "A", "B", "1", "100", "2", "200", "A-B fusion", "4", "8", "GRCh37"]])

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, fields, rows):
        with (self.raw / name).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields)
            writer.writerows(rows)

    def test_prepared_fields_and_end_to_end_join(self):
        prepared = self.root / "prepared"
        manifest = prepare(self.raw, prepared)
        self.assertEqual(manifest["counts"]["fusion_rows"], 1)
        with (prepared / "clinical.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["vital_status"], "alive")
        self.assertEqual(row["stage"], "")
        self.assertEqual(row["path_t_stage"], "T3A")
        db = self.root / "test.sqlite"
        import_callset(db, prepared / "fusion_calls.csv")
        import_clinical(db, prepared / "clinical.csv")
        report = combine_clinical(db, self.root / "results")
        self.assertEqual(report["summary_after_clinical_filters"]["matched_patients"], 1)

    def test_inconsistent_patient_sample_relationship_rejected(self):
        self.write("data_clinical_sample.txt", ["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE"], [["TCGA-ZZ-0001", "TCGA-ZZ-9999-01", "Primary"]])
        with self.assertRaisesRegex(ValueError, "barcode disagreement"):
            prepare(self.raw, self.root / "prepared")

    def test_unrecognized_os_status_not_guessed(self):
        self.patient_row[4] = "SOMETHING ELSE"
        self.write("data_clinical_patient.txt", self.patient_fields, [self.patient_row])
        with self.assertRaisesRegex(ValueError, "Unexpected OS_STATUS"):
            prepare(self.raw, self.root / "prepared")


if __name__ == "__main__":
    unittest.main()
