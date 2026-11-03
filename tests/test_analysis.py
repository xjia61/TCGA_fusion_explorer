import csv
import gzip
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

from fusion_explorer.core import Filters, analyze, import_callset


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.db = self.root / "test.sqlite"
        self.input = self.root / "calls.tsv"
        self.fields = ["data_class", "tissue_type", "sample", "fusion_name", "LeftGene",
                       "RightGene", "Jcount", "Scount", "FFPM", "tumor_or_normal", "LeftBreakpoint"]
        self.a = ["TCGA", "BRCA", "TCGA-AA-0001-01A", "A::B", "A^ENSG1", "B^ENSG2", "4", "2", "0.5", "tumor", "chr1:10:+"]
        self.rows = [self.a, self.a[:],
                     ["TCGA", "BRCA", "TCGA-AA-0002-01A", "A::B", "A", "B", "1", "0", "0.1", "tumor", "chr1:11:+"],
                     ["TCGA", "BRCA", "TCGA-AA-0001-11A", "B::A", "B", "A", "2", "1", "0.2", "normal", "chr1:12:+"],
                     ["GTEx", "Breast", "GTEX-TEST", "A::C", "A", "C", "3", "1", "0.3", "normal", "chr1:13:+"],
                     ["TCGA", "OV", "OV-24-1428-TP", "A::A", "A", "A", "", "", "", "tumor", ""],
                     ["TCGA", "OV", "TCGA-AA-0003-01A", "A::B", "A", "B", "-2", "1", "1", "tumor", ""]]
        self.write(self.rows)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, rows, fields=None, path=None):
        with (path or self.input).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields or self.fields)
            writer.writerows(rows)

    def read(self, name, out=None):
        with ((out or self.root / "out") / name).open() as handle:
            return list(csv.DictReader(handle))

    def test_import_quality_and_idempotence(self):
        quality = import_callset(self.db, self.input, batch_size=2)
        self.assertEqual((quality["rows_read"], quality["imported"], quality["duplicates"], quality["rejected"]), (7, 5, 1, 1))
        self.assertEqual(quality["missing_patient_id"], 2)
        self.assertEqual(quality["missing_support"], 1)
        self.assertTrue(import_callset(self.db, self.input)["reused"])
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM calls").fetchone()[0], 5)

    def test_counts_filters_and_direction(self):
        import_callset(self.db, self.input)
        report = analyze(self.db, self.root / "out")
        self.assertEqual(report["summary"], {"calls": 4, "samples": 4, "patients_with_canonical_id": 2, "distinct_fusions": 3, "cancer_types": 2})
        self.assertEqual(self.read("gene_summary.csv")[0], {"gene": "A", "calls": "4", "samples": "4"})
        self.assertEqual(len(self.read("fusion_summary.csv")), 3)
        report = analyze(self.db, self.root / "filtered", Filters(tissue="tumor", cancers=("TCGA-BRCA",), gene="a", min_junction=2))
        self.assertEqual(report["summary"]["calls"], 1)
        self.assertEqual(report["top_fusions"][0]["fusion"], "A::B")
        reverse = analyze(self.db, self.root / "reverse", Filters(fusion="b::a"))
        self.assertEqual(reverse["summary"]["calls"], 1)

    def test_missing_support_excluded_only_when_threshold_positive(self):
        import_callset(self.db, self.input)
        report = analyze(self.db, self.root / "out", Filters(min_ffpm=0.01))
        self.assertEqual(report["summary"]["calls"], 3)

    def test_breakpoint_distinct_calls_retained_sample_frequency_deduplicated(self):
        second = self.a[:]
        second[-1] = "chr1:999:+"
        self.write([self.a, second])
        import_callset(self.db, self.input)
        report = analyze(self.db, self.root / "out")
        self.assertEqual(report["summary"]["calls"], 2)
        self.assertEqual(report["top_fusions"], [{"fusion": "A::B", "samples": 1, "calls": 2}])

    def test_empty_filter_writes_headers_and_zero_counts(self):
        import_callset(self.db, self.input)
        report = analyze(self.db, self.root / "out", Filters(gene="NOT_PRESENT"))
        self.assertEqual(report["summary"]["calls"], 0)
        self.assertEqual(self.read("calls.csv"), [])
        self.assertTrue((self.root / "out" / "calls.csv").read_text().startswith("sample_id,"))

    def test_invalid_import_rolls_back(self):
        self.write([["TCGA", "OV", "TCGA-AA-0003-01A", "A::B", "A", "B", "NaN", "1", "inf", "tumor", ""]])
        with self.assertRaisesRegex(ValueError, "No valid"):
            import_callset(self.db, self.input)
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM sources").fetchone()[0], 0)

    def test_gzip_and_star_fusion_header(self):
        path = self.root / "star.tsv.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("# metadata\n#FusionName\tsample_id\tJunctionReadCount\tSpanningFragCount\tFFPM\nA::B\tTCGA-AA-0001-01A\t2\t1\t0.3\n")
        quality = import_callset(self.db, path)
        self.assertEqual(quality["imported"], 1)
        report = analyze(self.db, self.root / "out")
        self.assertEqual(report["summary"]["calls"], 1)
        self.assertEqual(self.read("calls.csv")[0]["tissue"], "tumor")

    def test_sources_are_not_silently_pooled(self):
        first = import_callset(self.db, self.input)
        self.write([self.a])
        second = import_callset(self.db, self.input)
        self.assertNotEqual(first["source_id"], second["source_id"])
        self.assertEqual(analyze(self.db, self.root / "out")["summary"]["calls"], 1)
        self.assertEqual(analyze(self.db, self.root / "old", source_id=first["source_id"])["summary"]["calls"], 4)

    def test_unknown_sample_not_silently_tcga(self):
        self.write([["other", "A::B"]], ["sample_id", "fusion"])
        import_callset(self.db, self.input)
        self.assertEqual(analyze(self.db, self.root / "out")["summary"]["calls"], 0)
        self.assertEqual(analyze(self.db, self.root / "all", Filters(cohort="ALL"))["summary"]["calls"], 1)

    def test_invalid_numeric_thresholds(self):
        for value in (-1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                Filters(min_ffpm=value).where(1)

    def test_conflicting_partners_and_malformed_rows_are_rejected(self):
        conflict = self.a[:]
        conflict[3] = "C::B"
        self.write([self.a, conflict, self.a[:-1]])
        quality = import_callset(self.db, self.input)
        self.assertEqual(quality["imported"], 1)
        self.assertEqual(quality["rejected"], 2)
        self.assertIn("disagrees", quality["rejection_examples"][0]["reason"])

    def test_missing_required_header_fails_without_source(self):
        self.write([["A::B"]], ["fusion"])
        with self.assertRaisesRegex(ValueError, "sample ID"):
            import_callset(self.db, self.input)
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM sources").fetchone()[0], 0)

    def test_cli_success_and_error(self):
        result = subprocess.run([sys.executable, "-m", "fusion_explorer", "analyze", "--input", str(self.input), "--db", str(self.db), "--out", str(self.root / "out")], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["summary"]["calls"], 4)
        error = subprocess.run([sys.executable, "-m", "fusion_explorer", "analyze", "--db", str(self.db), "--min-ffpm", "nan", "--out", str(self.root / "bad")], capture_output=True, text=True)
        self.assertEqual(error.returncode, 2)
        self.assertIn("finite", error.stderr)


if __name__ == "__main__":
    unittest.main()

