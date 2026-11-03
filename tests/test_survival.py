import csv
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fusion_explorer.core import Filters, import_callset
from fusion_explorer.clinical import ClinicalFilters, import_clinical
from fusion_explorer.survival import SurvivalSettings, endpoint_record, fit_group, survival_catalog, survival_query


HAS_LIFELINES = importlib.util.find_spec("lifelines") is not None


class EndpointTests(unittest.TestCase):
    def test_explicit_events_missing_and_invalid(self):
        self.assertEqual(endpoint_record({"os_time_months": "0", "os_event": "0"}, "OS"), ((0, 0), None))
        self.assertEqual(endpoint_record({"vital_status": "DEAD"}, "OS"), (None, "missing_endpoint"))
        for time, event in ((-1, 1), ("inf", 0), (1, 2), (1, "bad")):
            self.assertEqual(endpoint_record({"os_time_months": time, "os_event": event}, "OS")[1], "invalid_endpoint")


@unittest.skipUnless(HAS_LIFELINES, "Install dashboard dependencies for survival fitting")
class SurvivalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.db = self.root / "test.sqlite"
        f = self.root / "fusion.csv"
        c = self.root / "clinical.csv"
        self.write(f, ["sample_id", "fusion", "cohort", "cancer_type"], [
            ["TCGA-ZZ-0001-01", "A::B", "TCGA", "PRAD"],
            ["TCGA-ZZ-0001-01", "C::D", "TCGA", "PRAD"],
            ["TCGA-ZZ-0002-01", "A::B", "TCGA", "PRAD"],
            ["TCGA-ZZ-0003-01", "C::D", "TCGA", "PRAD"],
            ["TCGA-ZZ-0004-01", "E::F", "TCGA", "PRAD"],
        ])
        self.clinical_path = c
        self.fields = ["patient_id", "os_time_months", "os_event", "pfs_time_months", "pfs_event"]
        self.records = [["TCGA-ZZ-0001", 1, 1, 1, 1], ["TCGA-ZZ-0002", 2, 0, 1, 1],
                        ["TCGA-ZZ-0003", 2, 1, 2, 0], ["TCGA-ZZ-0004", 4, 0, 4, 0]]
        self.write(c, self.fields, self.records)
        import_callset(self.db, f)
        import_clinical(self.db, c)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, path, fields, records):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            writer.writerows(records)

    def test_km_tied_event_and_censor_hand_calculation(self):
        summary, curve, risks, censors = fit_group("test", [(1, 1), (2, 1), (2, 0)], 2)
        self.assertAlmostEqual(curve[1]["survival"], 2/3)
        self.assertAlmostEqual(curve[-1]["survival"], 1/3)
        self.assertEqual(summary["median_months"], 2)
        self.assertEqual(summary["landmark_at_risk"], 2)
        self.assertTrue(all(0 <= r["lower"] <= r["survival"] <= r["upper"] <= 1 for r in curve))
        self.assertEqual(censors[0]["censored"], 1)

    def test_all_censored_empty_group_and_no_extrapolation(self):
        summary, curve, _, _ = fit_group("test", [(1, 0), (2, 0)], 60)
        self.assertIsNone(summary["median_months"])
        self.assertIsNone(summary["landmark_survival"])
        self.assertTrue(all(r["survival"] == 1 for r in curve))
        self.assertEqual(fit_group("empty", [])[0]["patients"], 0)

    def test_target_other_patient_dedup_and_endpoints(self):
        cat = survival_catalog(self.db, 1, 1)
        self.assertEqual(cat["endpoints"]["OS"], {"patients": 4, "events": 2})
        result = survival_query(self.db, 1, 1, SurvivalSettings(("A::B",)))
        self.assertEqual([r["patients"] for r in result["groups"]], [2, 2])
        self.assertEqual([r["events"] for r in result["groups"]], [1, 1])
        pfs = survival_query(self.db, 1, 1, SurvivalSettings(("A::B",), endpoint="PFS"))
        self.assertEqual([r["events"] for r in pfs["groups"]], [2, 0])
        json.dumps(result, allow_nan=False)

    def test_overlapping_groups_excluded_explicitly(self):
        result = survival_query(self.db, 1, 1, SurvivalSettings(("A::B", "C::D"), mode="fusion_groups"))
        self.assertEqual(result["excluded"]["overlapping_selected_fusions"], 1)
        self.assertEqual(result["excluded"]["outside_selected_fusions"], 1)
        self.assertEqual([r["patients"] for r in result["groups"]], [1, 1])

    def test_no_outcome_conditioning_or_normal_tissue_analysis(self):
        for f, cf in ((Filters(), ClinicalFilters(vital_status="alive")),
                      (Filters(tissue="normal"), ClinicalFilters())):
            with self.assertRaises(ValueError):
                survival_query(self.db, 1, 1, SurvivalSettings(("A::B",)), f, cf)

    def test_conflicting_survival_records_roll_back(self):
        changed = self.records[0][:]
        changed[2] = 0
        self.write(self.clinical_path, self.fields, self.records + [changed])
        with self.assertRaisesRegex(ValueError, "Conflicting clinical"):
            import_clinical(self.db, self.clinical_path)

    def test_missing_endpoint_exclusion_is_reported(self):
        self.records[0][1] = ""
        self.write(self.clinical_path, self.fields, self.records)
        clinical_id = import_clinical(self.db, self.clinical_path)["clinical_source_id"]
        result = survival_query(self.db, 1, clinical_id, SurvivalSettings(("A::B",)))
        self.assertEqual(result["excluded"]["missing_endpoint"], 1)
        self.assertEqual(sum(r["patients"] for r in result["groups"]), 3)

    @unittest.skipUnless(importlib.util.find_spec("streamlit"), "Install Streamlit for the UI check")
    def test_survival_tab_and_endpoint_switch(self):
        from streamlit.testing.v1 import AppTest
        dashboard = Path(__file__).resolve().parents[1] / "dashboard.py"
        with patch.dict(os.environ, {"FUSION_EXPLORER_DB": str(self.db)}):
            app = AppTest.from_file(str(dashboard)).run(timeout=30)
            self.assertEqual(len(app.exception), 0)
            self.assertIn("Fusion survival analysis", [t.label for t in app.tabs])
            next(s for s in app.selectbox if s.label == "Survival endpoint").select("PFS")
            next(b for b in app.button if b.label == "Calculate survival curves").click()
            app.run(timeout=30)
            self.assertEqual(len(app.exception), 0)
            self.assertEqual(app.session_state["survival_settings"].endpoint, "PFS")


if __name__ == "__main__":
    unittest.main()
