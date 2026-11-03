import csv
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from fusion_explorer.core import import_callset
from fusion_explorer.clinical import ClinicalFilters, combine_clinical
from fusion_explorer.prepare_skcm_clinical import prepare


class SkcmClinicalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.raw = self.root/"patient.tsv"
        self.db = self.root/"test.sqlite"
        self.fields = ["PATIENT_ID","CANCER_TYPE_ACRONYM","AGE","SEX","RACE","OS_STATUS","OS_MONTHS","PFS_STATUS","PFS_MONTHS","AJCC_PATHOLOGIC_TUMOR_STAGE","ICD_O_3_HISTOLOGY"]
        self.records = [
            ["TCGA-ZZ-0001","SKCM","50","Male","","1:DECEASED","10","1:PROGRESSION","5","STAGE IIA","8720/3"],
            ["TCGA-ZZ-0002","SKCM","","Female","Asian","0:LIVING","","","","",""],
        ]
        self.write()

    def tearDown(self):
        self.temp.cleanup()

    def write(self):
        with self.raw.open("w",newline="") as handle:
            writer=csv.writer(handle,delimiter="\t")
            writer.writerow(self.fields);writer.writerows(self.records)

    def test_normalization_and_patient_join_with_unknowns(self):
        result=prepare(self.raw,self.root/"prepared",self.db)
        self.assertEqual(result["import"]["imported"],2)
        calls=self.root/"calls.csv"
        calls.write_text("sample_id,fusion,cohort,cancer_type\nTCGA-ZZ-0001-06,A::B,TCGA,SKCM\nTCGA-ZZ-0002-01,C::D,TCGA,SKCM\n")
        import_callset(self.db,calls)
        report=combine_clinical(self.db,self.root/"out",clinical_filters=ClinicalFilters(races=("unknown",)))
        self.assertEqual(report["summary_after_clinical_filters"]["matched_patients"],1)
        with sqlite3.connect(self.db) as db:
            row=db.execute("SELECT stage,raw_json FROM clinical_records WHERE patient_id='TCGA-ZZ-0001'").fetchone()
        raw=json.loads(row[1])
        self.assertEqual(row[0],"IIA")
        self.assertEqual(raw["os_event"],"1")
        self.assertEqual(raw["pfs_time_months"],"5")
        self.assertEqual(raw["histology"],"ICD-O-3: 8720/3")
        self.assertIn("OS_STATUS",raw)

    def test_reimport_is_idempotent(self):
        prepare(self.raw,self.root/"prepared",self.db)
        self.assertTrue(prepare(self.raw,self.root/"prepared",self.db)["import"]["reused"])

    def test_wrong_cancer_and_unknown_event_coding_rejected(self):
        for column,bad in ((1,"PRAD"),(5,"UNRECOGNIZED")):
            old=self.records[0][column]
            self.records[0][column]=bad;self.write()
            with self.assertRaises(ValueError):
                prepare(self.raw,self.root/"prepared",self.db)
            self.records[0][column]=old


if __name__ == "__main__":
    unittest.main()
