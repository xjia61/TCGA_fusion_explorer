import csv
import json
from pathlib import Path
import tempfile
import unittest

from fusion_explorer.pathology_manifest import build_manifest


class PathologyManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fusion = self.root / "sv.tsv"
        self.samples = self.root / "samples.tsv"
        self.slides = self.root / "slides.json"
        self.fusion.write_text("Sample_Id\tEvent_Info\nTCGA-ZZ-0001-06\tA-B Fusion\nTCGA-ZZ-0002-01\tC-D fusion\n")
        self.samples.write_text("SAMPLE_ID\tPATIENT_ID\tSAMPLE_TYPE\nTCGA-ZZ-0001-06\tTCGA-ZZ-0001\tMetastasis\nTCGA-ZZ-0002-01\tTCGA-ZZ-0002\tPrimary\n")
        self.response = {"data": {"pagination": {"total": 1}, "hits": [
            {"file_id": "file-test", "file_name": "test.svs", "access": "open", "file_size": 100,
             "cases": [{"case_id": "case-test", "submitter_id": "TCGA-ZZ-0001"}],
             "associated_entities": [{"entity_type": "slide", "entity_submitter_id": "TCGA-ZZ-0001-01A-01-DX1"}]}]}, "warnings": {}}
        self.write()

    def tearDown(self):
        self.temp.cleanup()

    def write(self):
        self.slides.write_text(json.dumps(self.response))

    def build(self):
        return build_manifest(self.fusion, self.samples, self.slides, self.root / "out")

    def test_patient_only_links_and_unmatched_samples(self):
        report = self.build()
        self.assertEqual(report["counts"]["matched_fusion_samples"], 1)
        self.assertEqual(report["counts"]["unmatched_fusion_samples"], 1)
        self.assertEqual(report["counts"]["same_specimen_verified_links"], 0)
        self.assertEqual(report["WSI_bytes_downloaded"], 0)
        with (self.root / "out/fusion_slide_links.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["match_level"], "patient_only_verified")
        self.assertEqual(row["sample_type_consistency"], "different_type_hint")
        self.assertEqual(row["same_specimen_verified"], "no")
        self.assertTrue(row["case_viewer_entry"].endswith("case-test"))

    def test_incomplete_pagination_rejected(self):
        self.response["data"]["pagination"]["total"] = 2
        self.write()
        with self.assertRaisesRegex(ValueError, "pagination"):
            self.build()

    def test_conflicting_slide_identity_rejected(self):
        self.response["data"]["hits"][0]["associated_entities"][0]["entity_submitter_id"] = "TCGA-ZZ-9999-01A-01-DX1"
        self.write()
        with self.assertRaisesRegex(ValueError, "conflicts"):
            self.build()

    def test_fusion_without_authoritative_sample_mapping_rejected(self):
        self.fusion.write_text("Sample_Id\tEvent_Info\nSKCM-ZZ-0001-TM\tA-B Fusion\n")
        with self.assertRaisesRegex(ValueError, "authoritative"):
            self.build()


if __name__ == "__main__":
    unittest.main()
