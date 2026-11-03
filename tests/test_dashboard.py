import csv
import hashlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fusion_explorer.core import Filters, analyze, import_callset
from fusion_explorer.clinical import ClinicalFilters, combine_clinical, import_clinical
from fusion_explorer.dashboard_data import catalog, options, dashboard_query

try:
    from streamlit.testing.v1 import AppTest
except ImportError:
    AppTest = None


ROOT = Path(__file__).resolve().parents[1]


class DashboardTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.db = self.root / "test.sqlite"
        fusion = self.root / "calls.csv"
        clinical = self.root / "clinical.csv"
        self.write(fusion, ["sample_id", "fusion", "cohort", "cancer_type"], [
            ["TCGA-ZZ-0001-01", "A::B", "TCGA", "BRCA"],
            ["TCGA-ZZ-0001-01", "C::D", "TCGA", "BRCA"],
            ["TCGA-ZZ-0002-01", "A::B", "TCGA", "PRAD"],
            ["GTEX-DEMO", "A::B", "GTEx", "Breast"],
        ])
        self.write(clinical, ["patient_id", "age_years", "stage", "sex", "vital_status"], [
            ["TCGA-ZZ-0001", "50", "II", "female", "alive"],
            ["TCGA-ZZ-0002", "60", "III", "male", "dead"],
        ])
        import_callset(self.db, fusion)
        import_clinical(self.db, clinical)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, path, fields, records):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            writer.writerows(records)

    def test_dashboard_matches_cli_and_deduplicates_patients(self):
        result = dashboard_query(self.db, 1, clinical_source_id=1)
        report = combine_clinical(self.db, self.root / "out")
        self.assertEqual(result["summary"]["calls"], report["summary_after_clinical_filters"]["calls"])
        self.assertEqual(result["summary"]["matched_patients"], 2)
        self.assertEqual(sum(r["patients"] for r in result["ages"]), 2)
        self.assertEqual(result["match_fraction"], 1)
        self.assertEqual(result["top_fusions"][0], {"fusion": "A::B", "calls": 2, "samples": 2})
        self.assertEqual(len(result["cancers"]), 2)

    def test_clinical_filters_and_unfiltered_match_audit(self):
        result = dashboard_query(self.db, 1, clinical_source_id=1,
                                 clinical_filters=ClinicalFilters(min_age=55, sex="male"))
        self.assertEqual(result["summary"]["calls"], 1)
        self.assertEqual(result["preview"][0]["sample_id"], "TCGA-ZZ-0002-01")
        self.assertEqual(result["audit"], [{"match_status": "matched", "sample_records": 2}])

    def test_pagination_empty_results_and_no_main_database_writes(self):
        before = hashlib.sha256(self.db.read_bytes()).hexdigest()
        first = dashboard_query(self.db, 1, page_size=1)
        second = dashboard_query(self.db, 1, page=1, page_size=1)
        last = dashboard_query(self.db, 1, page=999, page_size=1)
        self.assertNotEqual(first["preview"], second["preview"])
        self.assertEqual(last["page"], 2)
        empty = dashboard_query(self.db, 1, Filters(gene="MISSING"))
        self.assertEqual(empty["preview"], [])
        self.assertEqual(empty["page_count"], 0)
        self.assertEqual(hashlib.sha256(self.db.read_bytes()).hexdigest(), before)

    def test_options_missing_fields_and_invalid_sources(self):
        self.assertFalse(options(self.db, 1)["has_ffpm"])
        self.assertEqual(options(self.db, 1, 1)["stages"], ["II", "III"])
        self.assertEqual(len(catalog(self.db)["fusion_sources"]), 1)
        with self.assertRaises(ValueError):
            dashboard_query(self.db, 999)
        with self.assertRaises(ValueError):
            dashboard_query(self.db, 1, clinical_filters=ClinicalFilters(min_age=1))
        with self.assertRaises(ValueError):
            dashboard_query(self.db, 1, clinical_source_id=1, mapping_source_id=999)
        with self.assertRaises(ValueError):
            catalog(self.root / "absent.sqlite")
        self.assertFalse((self.root / "absent.sqlite").exists())

    @unittest.skipIf(AppTest is None, "Install dashboard dependencies to run Streamlit interaction checks")
    def test_streamlit_initial_state_filter_and_empty_results(self):
        with patch.dict(os.environ, {"FUSION_EXPLORER_DB": str(self.db)}):
            app = AppTest.from_file(str(ROOT / "dashboard.py")).run(timeout=30)
            self.assertEqual(len(app.exception), 0)
            self.assertEqual(app.metric[0].value, "3")
            self.assertEqual(app.metric[3].value, "2")
            self.assertEqual(next(w for w in app.multiselect if w.label == "Race").options, ["unknown"])
            self.assertEqual(next(w for w in app.multiselect if w.label == "Histology").options, ["unknown"])
            gene = next(w for w in app.text_input if w.label == "Partner gene")
            gene.set_value("C")
            next(b for b in app.button if b.label == "Apply filters").click()
            app.run(timeout=30)
            self.assertEqual(len(app.exception), 0)
            self.assertEqual(app.metric[0].value, "1")
            next(w for w in app.text_input if w.label == "Partner gene").set_value("ABSENT")
            next(b for b in app.button if b.label == "Apply filters").click()
            app.run(timeout=30)
            self.assertEqual(len(app.exception), 0)
            self.assertEqual(app.metric[0].value, "0")


if __name__ == "__main__":
    unittest.main()
