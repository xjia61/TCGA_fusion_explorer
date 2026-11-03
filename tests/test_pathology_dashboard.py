import csv
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fusion_explorer.core import Filters, import_callset
from fusion_explorer.pathology import import_pathology, pathology_query
from fusion_explorer.pathology_manifest import build_manifest


class PathologyDashboardTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.db = self.root/"test.sqlite"
        (self.root/"samples.tsv").write_text("PATIENT_ID\tSAMPLE_ID\tSAMPLE_TYPE\nTCGA-ZZ-0001\tTCGA-ZZ-0001-06\tMetastasis\nTCGA-ZZ-0002\tTCGA-ZZ-0002-01\tPrimary\n")
        (self.root/"sv.tsv").write_text("Sample_Id\tEvent_Info\nTCGA-ZZ-0001-06\tA-B Fusion\nTCGA-ZZ-0002-01\tC-D Fusion\n")
        (self.root/"calls.csv").write_text("sample_id,fusion,cohort,cancer_type\nTCGA-ZZ-0001-06,A::B,TCGA,SKCM\nTCGA-ZZ-0002-01,C::D,TCGA,SKCM\n")
        (self.root/"slides.json").write_text(json.dumps({"data":{"pagination":{"total":1},"hits":[{
            "file_id":"file-1","file_name":"test.svs","file_size":100,"access":"open",
            "cases":[{"case_id":"case-1","submitter_id":"TCGA-ZZ-0001"}],
            "associated_entities":[{"entity_type":"slide","entity_submitter_id":"TCGA-ZZ-0001-01A-01-DX1"}]}]},"warnings":{}}))
        build_manifest(self.root/"sv.tsv",self.root/"samples.tsv",self.root/"slides.json",self.root/"manifest")
        self.links = self.root/"manifest/fusion_slide_links.csv"
        self.report = self.root/"manifest/linkage_report.json"
        import_callset(self.db,self.root/"calls.csv")
        import_pathology(self.db,self.links,self.report,1)

    def tearDown(self):
        self.temp.cleanup()

    def test_patient_and_fusion_filters_export_and_tissue_mismatch(self):
        result = pathology_query(self.db,1,Filters(gene="A"),export_path=self.root/"export.csv")
        self.assertEqual(result["summary"]["patients"],1)
        self.assertEqual(result["summary"]["different_type_links"],1)
        self.assertEqual(result["records"][0]["same_specimen_verified"],"no")
        with (self.root/"export.csv").open(encoding="utf-8-sig") as handle:
            self.assertEqual(len(list(csv.DictReader(handle))),1)
        self.assertEqual(pathology_query(self.db,1,Filters(gene="C"))["summary"]["links"],0)
        self.assertEqual(pathology_query(self.db,1,patient="TCGA-ZZ-0002")["summary"]["links"],0)

    def test_idempotence_and_source_isolation(self):
        self.assertTrue(import_pathology(self.db,self.links,self.report,1)["reused"])
        (self.root/"new.csv").write_text("sample_id,fusion,cohort,cancer_type\nTCGA-ZZ-0001-06,A::C,TCGA,SKCM\n")
        new=import_callset(self.db,self.root/"new.csv")["source_id"]
        self.assertFalse(pathology_query(self.db,new)["available"])

    def test_wrong_patient_import_rejected_atomically(self):
        content=self.links.read_text().replace("TCGA-ZZ-0001,","TCGA-ZZ-9999,")
        self.links.write_text(content)
        with self.assertRaisesRegex(ValueError,"does not match"):
            import_pathology(self.db,self.links,self.report,1)
        self.assertEqual(pathology_query(self.db,1)["summary"]["links"],1)

    @unittest.skipUnless(importlib.util.find_spec("streamlit"),"Install dashboard dependencies")
    def test_dashboard_pathology_tab_and_patient_filter(self):
        from streamlit.testing.v1 import AppTest
        with patch.dict(os.environ,{"FUSION_EXPLORER_DB":str(self.db)}):
            app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/"dashboard.py")).run(timeout=30)
            self.assertEqual(len(app.exception),0)
            self.assertIn("Pathology slides",[t.label for t in app.tabs])
            self.assertEqual(next(m for m in app.metric if m.label=="Linked patients").value,"1")
            next(t for t in app.text_input if t.label=="Pathology patient ID").set_value("TCGA-ZZ-0002")
            next(b for b in app.button if b.label=="Filter pathology list").click()
            app.run(timeout=30)
            self.assertEqual(len(app.exception),0)
            self.assertEqual(next(m for m in app.metric if m.label=="Linked patients").value,"0")


if __name__ == "__main__":
    unittest.main()
